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
#include "kinTMIII_gpu.cu.h"

#include "spline.h"
#include "analysis.h"
#include "io.h"
#include <fstream>

#include <astro/system/config.h>
#include <astro/useall.h>

// os_kinTMIII -- Generate kinematics based on Bond et al. (in prep)
class os_kinTMIII : public osink, os_kinTMIII_data
{	
	float DeltavPhi;
	interval_list icomp_thin, icomp_thick, icomp_halo;
	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
		virtual const std::string &name() const { static std::string s("kinTMIII"); return s; }
		virtual double ordering() const { return ord_kinematics; }
		virtual bit_map getAffectedComponents() const
		{
			bit_map ret = icomp_thin;
			ret |= icomp_thick;
			ret |= icomp_halo;
			return ret;
		}

		os_kinTMIII() : osink()
		{
			prov.insert("vcyl");
			req.insert("comp");
			req.insert("XYZ");
		}
};
extern "C" opipeline_stage *create_module_kintmiii() { return new os_kinTMIII(); }	// Factory; called by opipeline_stage::create()

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

extern os_kinTMIII_data os_kinTMIII_par;

size_t os_kinTMIII::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input

	// fetch prerequisites
	cint_t   &comp  = in.col<int>("comp");
	cfloat_t &XYZ   = in.col<float>("XYZ");
	cfloat_t &vcyl   = in.col<float>("vcyl");

	comp_thin = icomp_thin;
	comp_thick = icomp_thick;
	comp_halo = icomp_halo;
	cuxUploadConst("os_kinTMIII_par", static_cast<os_kinTMIII_data&>(*this));	// for GPU execution
	os_kinTMIII_par = static_cast<os_kinTMIII_data&>(*this);			// for CPU execution

	CALL_KERNEL(os_kinTMIII_kernel, otable_ks(begin, end), rng, comp, XYZ, vcyl);
	return nextlink->process(in, begin, end, rng);
}

int split_fvec(std::vector<float>& arr, const std::string &text)
{
	std::vector<float>::value_type tmp;
	std::stringstream ss(text);

	arr.clear();
	while(ss >> tmp) { arr.push_back(tmp); }
	return arr.size();
}

std::vector<float> split_fvec(const std::string &text)
{
	std::vector<float> ret;
	split_fvec(ret, text);
	return ret;
}

void fvecToFarray(const std::vector<float>& src, farray5& dst)
{
	for (int i=0;i<5;i++) {
		dst[i]=src[i];
		}
}

void farray_to_iarray(farray5& fa, iarray5& ia)
{
	for (int i=0;i<5;i++)
		ia[i]=int(fa[i]*100);
}

void farray_to_i8array(farray5& fa, i8array5& ia)
{
	for (int i=0;i<5;i++)
		ia[i]=char(fa[i]/10);
}

template<typename T> inline OSTREAM(const std::vector<T> &v) { FOREACH(v) { out << *i << " "; }; return out; }

bool os_kinTMIII::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// Component IDs
	read_component_map(icomp_thin, cfg, "comp_thin");
	read_component_map(icomp_thick, cfg, "comp_thick");
	read_component_map(icomp_halo, cfg, "comp_halo");
// 	cfg.get(comp_thin,  "comp_thin",  0);
// 	cfg.get(comp_thick, "comp_thick", 1);
// 	cfg.get(comp_halo,  "comp_halo",  2);

	// Distance to the Galactic center
	Rg = cfg.get("Rg");

	cfg.get(fk           , "fk"           , 3.0f);
	cfg.get(DeltavPhi    , "DeltavPhi"    , 34.0f);
	fk = fk / (1. + fk);	// renormalize to probability of drawing from the first gaussian

	typedef std::vector<float> fvec;
	fvec 	vPhi1_fvec, vPhi2_fvec, vR_fvec, vZ_fvec,
		sigmaPhiPhi1_fvec, sigmaPhiPhi2_fvec, sigmaRR_fvec, sigmaZZ_fvec, sigmaRPhi_fvec, sigmaZPhi_fvec, sigmaRZ_fvec,
		HvPhi_fvec, HvR_fvec, HvZ_fvec,
		HsigmaPhiPhi_fvec, HsigmaRR_fvec, HsigmaZZ_fvec, HsigmaRPhi_fvec, HsigmaZPhi_fvec, HsigmaRZ_fvec;

	cfg.get(vR_fvec           , "vR"           , split_fvec("0 0 0 0 0"));	
	cfg.get(vPhi1_fvec        , "vPhi"         , split_fvec("-194 19.2 1.25 0 0"));	
	cfg.get(vZ_fvec           , "vZ"           , split_fvec("0 0 0 0 0"));	
	cfg.get(sigmaRR_fvec      , "sigmaRR"      , split_fvec("40 5 1.5 0 0"));	
	cfg.get(sigmaRPhi_fvec    , "sigmaRPhi"    , split_fvec("0 0 0 0 0"));	
	cfg.get(sigmaRZ_fvec      , "sigmaRZ"      , split_fvec("0 0 0 0 0"));	
	cfg.get(sigmaPhiPhi1_fvec , "sigmaPhiPhi1" , split_fvec("12 1.8 2 0 11"));		// dynamically changeable
	cfg.get(sigmaPhiPhi2_fvec , "sigmaPhiPhi2" , split_fvec("34 1.2 2 0 0"));
	cfg.get(sigmaZPhi_fvec    , "sigmaZPhi"    , split_fvec("0 0 0 0 0"));	
	cfg.get(sigmaZZ_fvec      , "sigmaZZ"      , split_fvec("25 4 1.5 0 0"));	

	fvecToFarray(vR_fvec, vR);
	fvecToFarray(vPhi1_fvec, vPhi1);
	fvecToFarray(vZ_fvec, vZ);
	fvecToFarray(sigmaRR_fvec, sigmaRR);
	fvecToFarray(sigmaRPhi_fvec, sigmaRPhi);
	fvecToFarray(sigmaRZ_fvec, sigmaRZ);
	fvecToFarray(sigmaPhiPhi1_fvec, sigmaPhiPhi1);
	fvecToFarray(sigmaPhiPhi2_fvec, sigmaPhiPhi2);
	fvecToFarray(sigmaZPhi_fvec, sigmaZPhi);
	fvecToFarray(sigmaZZ_fvec, sigmaZZ);

	cfg.get(HvR_fvec          , "HvR"          , split_fvec("0 0 0 0 0"));	
	cfg.get(HvPhi_fvec        , "HvPhi"        , split_fvec("0 0 0 0 0"));	
	cfg.get(HvZ_fvec          , "HvZ"          , split_fvec("0 0 0 0 0"));	
	cfg.get(HsigmaRR_fvec     , "HsigmaRR"     , split_fvec("135 0 0 0 0"));	
	cfg.get(HsigmaRPhi_fvec   , "HsigmaRPhi"   , split_fvec("0 0 0 0 0"));	
	cfg.get(HsigmaRZ_fvec     , "HsigmaRZ"     , split_fvec("0 0 0 0 0"));	
	cfg.get(HsigmaPhiPhi_fvec , "HsigmaPhiPhi" , split_fvec("85 0 0 0 0"));	
	cfg.get(HsigmaZPhi_fvec   , "HsigmaZPhi"   , split_fvec("0 0 0 0 0"));	
	cfg.get(HsigmaZZ_fvec     , "HsigmaZZ"     , split_fvec("85 0 0 0 0"));	

	fvecToFarray(HvR_fvec, HvR);
	fvecToFarray(HvPhi_fvec, HvPhi);
	fvecToFarray(HvZ_fvec, HvZ);
	fvecToFarray(HsigmaRR_fvec, HsigmaRR);
	fvecToFarray(HsigmaRPhi_fvec, HsigmaRPhi);
	fvecToFarray(HsigmaRZ_fvec, HsigmaRZ);
	fvecToFarray(HsigmaPhiPhi_fvec, HsigmaPhiPhi);
	fvecToFarray(HsigmaZPhi_fvec, HsigmaZPhi);
	fvecToFarray(HsigmaZZ_fvec, HsigmaZZ);

 	vPhi2 = vPhi1;
 	vPhi2[0] += DeltavPhi;

	// Output model parameters
	MLOG(verb2) << "Component IDs (thn, thk, hl): " << icomp_thin << " " << icomp_thick << " " << icomp_halo;

	MLOG(verb2) << "Disk gaussian normalizations: " << fk << " : " << (1-fk);
	MLOG(verb2) << "Second disk gaussian offset:  " << DeltavPhi;

	MLOG(verb2) << "vR coefficients:              " << vR_fvec;
	MLOG(verb2) << "vZ coefficients:              " << vZ_fvec;
	MLOG(verb2) << "sigmaRR coefficients:         " << sigmaRR_fvec;
	MLOG(verb2) << "sigmaRPhi coefficients:       " << sigmaRPhi_fvec;
	MLOG(verb2) << "sigmaRZ coefficients:         " << sigmaRZ_fvec;
	MLOG(verb2) << "sigmaPhiPhi1 coefficients:    " << sigmaPhiPhi1_fvec;
	MLOG(verb2) << "sigmaPhiPhi2 coefficients:    " << sigmaPhiPhi2_fvec;
	MLOG(verb2) << "sigmaZPhi coefficients:       " << sigmaZPhi_fvec;
	MLOG(verb2) << "sigmaZZ coefficients:         " << sigmaZZ_fvec;

	MLOG(verb2) << "HvR coefficients:             " << HvR_fvec;
	MLOG(verb2) << "HvZ coefficients:             " << HvZ_fvec;
	MLOG(verb2) << "HsigmaRR coefficients:        " << HsigmaRR_fvec;
	MLOG(verb2) << "HsigmaRPhi coefficients:      " << HsigmaRPhi_fvec;
	MLOG(verb2) << "HsigmaRZ coefficients:        " << HsigmaRZ_fvec;
	MLOG(verb2) << "HsigmaPhiPhi coefficients:    " << HsigmaPhiPhi_fvec;
	MLOG(verb2) << "HsigmaZPhi coefficients:      " << HsigmaZPhi_fvec;
	MLOG(verb2) << "HsigmaZZ coefficients:        " << HsigmaZZ_fvec;

	return true;
}
