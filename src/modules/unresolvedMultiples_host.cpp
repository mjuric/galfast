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

#include "galfast_config.h"

#include "../pipeline.h"
#include "unresolvedMultiples_gpu.cu.h"

#include "analysis.h"
#include "spline.h"
#include <fstream>

#include <astro/system/config.h>
#include <astro/useall.h>

// os_unresolvedMultiples -- Generate unresolved muliple systems
class os_unresolvedMultiples : public osink
{
	protected:
		std::string absmagSys;
		multiplesAlgorithms::algo algo;			// algorithm for magnitude assignment to secondaries

		cuxTexture<float> secProb, cumLF, invCumLF;		// probability and LF texture data
	public:
		virtual bool runtime_init(otable &t);
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
		virtual const std::string &name() const { static std::string s("unresolvedMultiples"); return s; }
		//virtual int priority() { return PRIORITY_STAR; } // ensure this is placed near the beginning of the pipeline
		virtual double ordering() const { return ord_multiples; }

		os_unresolvedMultiples() : osink()
		{
			req.insert("absmag");
			req.insert("comp");
		}
};
extern "C" opipeline_stage *create_module_unresolvedmultiples() { return new os_unresolvedMultiples(); }	// Factory; called by opipeline_stage::create()

bool os_unresolvedMultiples::runtime_init(otable &t)
{
	// Not ready until absmag is available
	if(!osink::runtime_init(t)) { return false; }

	// by default, absmagSys1 is aliased to absmag. Drop this alias, as we're going to
	// provide absmagSys1
	t.drop_column("M1");

	// output absolute magnitudes
	otable::columndef &col = t.getColumn("absmag");
	const std::string &absmag = col.getPrimaryName();
	absmagSys = absmag + "Sys";
	std::string band = col.get_property("band");
	std::string absmagSysDef = absmagSys + "[2]{class=magnitude;alias=absmagSys;band=" + band + ";fieldNames=0:M1,1:M2;}";
	t.use_column(absmagSysDef);

	// number of components present
	std::string ncompDef = absmagSys + "Ncomp{type=int;fmt=%1d;}";
	t.use_column(ncompDef);

	return true;
}

DECLARE_TEXTURE(secProb,  float, 1, cudaReadModeElementType);
DECLARE_TEXTURE(cumLF,    float, 1, cudaReadModeElementType);
DECLARE_TEXTURE(invCumLF, float, 1, cudaReadModeElementType);

size_t os_unresolvedMultiples::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input
	//	- all stars are main sequence
	cint_t   &comp  = in.col<int>("comp");
	cint_t   &hidden = in.col<int>("hidden");
	cfloat_t &M     = in.col<float>("absmag");
	cfloat_t &Msys  = in.col<float>(absmagSys);
	cint_t   &ncomp = in.col<int>(absmagSys+"Ncomp");

	{
		cuxTextureBinder
			t1(::secProb,	secProb),
			t2(::cumLF,	cumLF),
			t3(::invCumLF,	invCumLF);

		CALL_KERNEL(os_unresolvedMultiples_kernel, otable_ks(begin, end), applyToComponents, rng, Msys.width(), M, Msys, ncomp, comp, hidden, algo);
	}
	
	return nextlink->process(in, begin, end, rng);
}

bool os_unresolvedMultiples::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	read_component_map(applyToComponents, cfg);

	std::string LFfile, binaryFractionFile, strAlgo;
	cfg.get(LFfile, "lumfunc", "");
	cfg.get(binaryFractionFile, "fraction_file", "");
	strAlgo = cfg.get("algorithm");

	// decide on the secondary assignment algorithm
	using namespace multiplesAlgorithms;
	if(strAlgo == "LF_M2_gt_M1") 		{ algo = LF_M2_GT_M1; }
	else if(strAlgo == "LF")		{ algo = LF; }
	else if(strAlgo == "equal_mass")	{ algo = EQUAL_MASS; }
	else { THROW(EAny, "Unknow secondary mag. assignment algorithm '" + strAlgo + "'"); }

	// Load binary fraction
	if(!binaryFractionFile.empty())
	{
		secProb = load_and_resample_texture_1D(binaryFractionFile.c_str(), 64);
	}
	else
	{
		// 100% binary fraction across all plausible absolute magnitudes
		secProb = load_constant_texture_1D(1, -100, +100);
	}

	// Load luminosity function
	std::vector<double> x, y;
	if(!LFfile.empty())
	{
		text_input_or_die(datain, LFfile);
		::load(datain, x, 0, y, 1);
	}
	else
	{
		// generate uniform LF extending over a plausible range of
		// absolute magnitudes
		x.push_back(-100); y.push_back(1.);
		x.push_back(+100); y.push_back(1.);
	}

	// Construct cumulative distribution (the normalized integral of
	// piecewise linearly interpolated luminosify function)
	const int NPIX = 256;
	spline lf; lf.construct(x, y);
	double dx = (x.back() - x.front()) / (NPIX-1);
	std::vector<double> ycum(NPIX), xcum(NPIX);
	xcum[0] = x.front(); ycum[0] = 0;
	double yprev = lf(x.front());
	FOR(1, NPIX)
	{
		double xx = x.front() + i*dx;
		double yy = lf(xx);

		double dy = yy - yprev;
		double dA = (yprev + 0.5*dy)*dx;	// increase in area from y[i-1] to y[i]

		xcum[i] = xx;
		ycum[i] = ycum[i-1] + dA;

		yprev = yy;
		//std::cerr << xcum[i] << " " << ycum[i] << "\n";
	}
	double norm = ycum.back();
	FOR(0, ycum.size()) { ycum[i] /= norm; }
	//FOR(0, ycum.size()) { std::cerr << xcum[i] << " " << ycum[i] << "\n"; }

	// NOTE: WARNING: because of resampling, invCumLF(cumLF(x)) != x,
	// so DONT EVER DEPEND ON IT!
	cumLF    = construct_texture_by_resampling_1D(&xcum[0], &ycum[0], xcum.size(), NPIX);
	invCumLF = construct_texture_by_resampling_1D(&ycum[0], &xcum[0], xcum.size(), NPIX);
	//FOR(0, xcum.size()) { std::cerr << xcum[i] << " " << ycum[i] << " " << cumLFManager.sample(xcum[i]) << "\n"; }

// 	for(float u=0; u <=1; u += 0.01)
// 	{
// 		std::cerr << u << "\t" << invCumLFManager.sample(u) << "\n";
// 	}

	return true;
}
