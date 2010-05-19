/***************************************************************************
 *   Copyright (C) 2004 by Mario Juric                                     *
 *   mjuric@astro.Princeton.EDU                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "config.h"

#include "../pipeline.h"
#include "Bond2010_gpu.cu.h"
#include "transform.h"

#include "spline.h"
#include "analysis.h"
#include "io.h"
#include <fstream>

#include <astro/system/config.h>
#include <astro/useall.h>

using namespace Bond2010;

// os_Bond2010 -- Generate kinematics based on Bond et al. (in prep)
class os_Bond2010 : public osink, os_Bond2010_data
{	
	interval_list icomp_thin, icomp_thick, icomp_halo;

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
		virtual const std::string &name() const { static std::string s("Bond2010"); return s; }
		virtual double ordering() const { return ord_kinematics; }
		virtual bit_map getAffectedComponents() const
		{
			bit_map ret = icomp_thin;
			ret |= icomp_thick;
			ret |= icomp_halo;
			return ret;
		}

		os_Bond2010() : osink()
		{
			prov.insert("vcyl");
			req.insert("comp");
			req.insert("XYZ");
		}
};
extern "C" opipeline_stage *create_module_bond2010() { return new os_Bond2010(); }	// Factory; called by opipeline_stage::create()

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

extern os_Bond2010_data os_Bond2010_par;

size_t os_Bond2010::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- component tags exist in input
	//	- geocentric XYZ coordinates exist in input

	// fetch prerequisites
	cint_t   &comp  = in.col<int>("comp");
	cint_t   &hidden= in.col<int>("hidden");
	cfloat_t &XYZ   = in.col<float>("XYZ");
	cfloat_t &vcyl   = in.col<float>("vcyl");

	comp_thin = icomp_thin;
	comp_thick = icomp_thick;
	comp_halo = icomp_halo;
	cuxUploadConst("os_Bond2010_par", static_cast<os_Bond2010_data&>(*this));	// for GPU execution
	os_Bond2010_par = static_cast<os_Bond2010_data&>(*this);			// for CPU execution

	CALL_KERNEL(os_Bond2010_kernel, otable_ks(begin, end), rng, comp, hidden, XYZ, vcyl);
	return nextlink->process(in, begin, end, rng);
}

static int split_fvec(std::vector<float>& arr, const std::string &text)
{
	std::vector<float>::value_type tmp;
	std::stringstream ss(text);

	arr.clear();
	while(ss >> tmp) { arr.push_back(tmp); }
	return arr.size();
}

static std::vector<float> split_fvec(const std::string &text)
{
	std::vector<float> ret;
	split_fvec(ret, text);
	return ret;
}

static void fvecToFarray(const std::vector<float>& src, farray5& dst)
{
	for(int i=0; i != src.size(); i++)
	{
		dst[i] = src[i];
	}
	for(int i = src.size(); i < 5; i++)
	{
		dst[i] = 0;
	}
}

OSTREAM(const farray5 &v)
{
	for(int i = 0; i != 4; i++) { out << v[i] << " "; }
	return out << v[4];
}

OSTREAM(const vmoments &v)
{
	out << "\n";
	out << "    mean_0: " << v.m[0] << "\n";
	out << "    mean_1: " << v.m[1] << "\n";
	out << "    mean_2: " << v.m[2] << "\n";
	out << "    s11^.5: " << v.ss[0] << "\n";
	out << "       s12: " << v.ss[1] << "\n";
	out << "       s13: " << v.ss[2] << "\n";
	out << "    s22^.5: " << v.ss[3] << "\n";
	out << "       s23: " << v.ss[4] << "\n";
	out << "    s33^.5: " << v.ss[5] << "\n";
	return out;
}

bool os_Bond2010::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// Component IDs
	read_component_map(icomp_thin,  cfg, "comp_thin");
	read_component_map(icomp_thick, cfg, "comp_thick");
	read_component_map(icomp_halo,  cfg, "comp_halo");

	// Load disk position and orientation
//	Rg = cfg.get("Rg");
	std::string ctr, orient;
	cfg.get(ctr,    "center",      "galplane");
	cfg.get(orient, "orientation", "galplane");
	load_transform(&T.x, M, ctr, orient, cfg);

	typedef std::vector<float> fvec;
	fvec 	vPhi1_fvec, vPhi2_fvec, vR_fvec, vZ_fvec,
		sigmaPhi1_fvec, sigmaPhi2_fvec, sigmaR_fvec, sigmaZ_fvec, sigmaRPhi_fvec, sigmaPhiZ_fvec, sigmaRZ_fvec,
		HvPhi_fvec, HvR_fvec, HvTh_fvec,
		HsigmaPhi_fvec, HsigmaR_fvec, HsigmaTh_fvec, HsigmaRPhi_fvec, HsigmaThPhi_fvec, HsigmaRTh_fvec;

	// Defaults are from Table 1 of Bond et al. 2010 (arXiv:0909.0013v2)
	// These are (a,b,c,d,e) coefficients of a + b*|Z|^c + d*Rcyl^e
	cfg.get(vR_fvec           , "vR"           , split_fvec("  0   0   0    0 0"));
	cfg.get(vPhi1_fvec        , "vPhi"         , split_fvec("-194 19.2 1.25 0 0"));
	cfg.get(vZ_fvec           , "vZ"           , split_fvec("  0   0   0    0 0"));
	cfg.get(sigmaR_fvec       , "sigmaR"       , split_fvec(" 40   5   1.5  0 0"));
	cfg.get(sigmaRPhi_fvec    , "sigmaRPhi"    , split_fvec("  0   0   0    0 0"));
	cfg.get(sigmaRZ_fvec      , "sigmaRZ"      , split_fvec("  0   0   0    0 0"));
	cfg.get(sigmaPhi1_fvec    , "sigmaPhi1"    , split_fvec(" 12   1.8 2    0 0"));
	cfg.get(sigmaPhi2_fvec    , "sigmaPhi2"    , split_fvec(" 34   1.2 2    0 0"));
	cfg.get(sigmaPhiZ_fvec    , "sigmaPhiZ"    , split_fvec("  0   0   0    0 0"));
	cfg.get(sigmaZ_fvec       , "sigmaZ"       , split_fvec(" 25   4   1.5  0 0"));

	// first disk gaussian
	fvecToFarray(vR_fvec,        disk_1.m[0]);
	fvecToFarray(vPhi1_fvec,     disk_1.m[1]);
	fvecToFarray(vZ_fvec,        disk_1.m[2]);
	fvecToFarray(sigmaR_fvec,    disk_1.ss[0]);
	fvecToFarray(sigmaRPhi_fvec, disk_1.ss[1]);
	fvecToFarray(sigmaRZ_fvec,   disk_1.ss[2]);
	fvecToFarray(sigmaPhi1_fvec, disk_1.ss[3]);
	fvecToFarray(sigmaPhiZ_fvec, disk_1.ss[4]);
	fvecToFarray(sigmaZ_fvec,    disk_1.ss[5]);

	// second disk gaussian
	disk_2 = disk_1;
	fvecToFarray(sigmaPhi2_fvec, disk_2.ss[3]);
	float DeltavPhi;
	cfg.get(DeltavPhi    , "DeltavPhi"    , 34.0f);
	disk_2.m[1][0] += DeltavPhi;
	cfg.get(fk           , "fk"           ,  3.0f);
	fk = fk / (1. + fk);	// renormalize to probability of drawing from the first gaussian

	// The meaning of the coefficients is the same as above.
	// Defaults reproduce the Bond et al. 2010. model for the halo
	cfg.get(HvR_fvec          , "HvR"          , split_fvec("0   0 0 0 0"));
	cfg.get(HvPhi_fvec        , "HvPhi"        , split_fvec("0   0 0 0 0"));
	cfg.get(HvTh_fvec         , "HvTh"         , split_fvec("0   0 0 0 0"));
	cfg.get(HsigmaR_fvec      , "HsigmaR"      , split_fvec("141 0 0 0 0"));
	cfg.get(HsigmaRPhi_fvec   , "HsigmaRPhi"   , split_fvec("0   0 0 0 0"));
	cfg.get(HsigmaRTh_fvec    , "HsigmaRTh"    , split_fvec("0   0 0 0 0"));
	cfg.get(HsigmaPhi_fvec    , "HsigmaPhi"    , split_fvec("85  0 0 0 0"));
	cfg.get(HsigmaThPhi_fvec  , "HsigmaThPhi"  , split_fvec("0   0 0 0 0"));
	cfg.get(HsigmaTh_fvec     , "HsigmaTh"     , split_fvec("75  0 0 0 0"));

	fvecToFarray(HvR_fvec,         halo.m[0]);
	fvecToFarray(HvTh_fvec,        halo.m[1]);
	fvecToFarray(HvPhi_fvec,       halo.m[2]);
	fvecToFarray(HsigmaR_fvec,     halo.ss[0]);
	fvecToFarray(HsigmaRTh_fvec,   halo.ss[1]);
	fvecToFarray(HsigmaRPhi_fvec,  halo.ss[2]);
	fvecToFarray(HsigmaTh_fvec,    halo.ss[3]);
	fvecToFarray(HsigmaThPhi_fvec, halo.ss[4]);
	fvecToFarray(HsigmaPhi_fvec,   halo.ss[5]);

	// Output model parameters
	MLOG(verb1) << "Kinematics: Bond+2010 for components " << icomp_thin << " (thin), " << icomp_thick << " (thick), " << icomp_halo << " (halo)" << "   ## " << instanceName();
	MLOG(verb2) << "Component IDs (thn, thk, hl): " << icomp_thin << " " << icomp_thick << " " << icomp_halo;
	MLOG(verb2) << "Disk gaussian normalizations: " << fk << " : " << (1-fk);
	MLOG(verb2) << "Second disk gaussian offset:  " << DeltavPhi;
	MLOG(verb2) << "Disk #1 (a + b*|Z|^c + d*Rcyl^e) coeffs (0=Rcyl, 1=Phi, 2=Z) :   " << disk_1;
	MLOG(verb2) << "Disk #2 (a + b*|Z|^c + d*Rcyl^e) coeffs (0=Rcyl, 1=Phi, 2=Z) :   " << disk_2;
	MLOG(verb2) << "Halo    (a + b*|Z|^c + d*Rcyl^e) coeffs (0=r, 1=Theta, 2=Phi):   " << halo;

	return true;
}
