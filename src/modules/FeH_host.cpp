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

#include "module_lib.h"
#include "../pipeline.h"
#include "FeH_gpu.cu.h"

#include <astro/system/config.h>
#include <astro/useall.h>

// os_FeH -- Generate Fe/H based on Ivezic et al (2008)
class os_FeH : public osink, os_FeH_data
{
public:
	interval_list icomp_thin, icomp_thick, icomp_halo;

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("FeH"); return s; }
	virtual double ordering() const { return ord_feh; }
	virtual bit_map getAffectedComponents() const
	{
		bit_map ret = icomp_thin;
		ret |= icomp_thick;
		ret |= icomp_halo;
		return ret;
	}

	os_FeH() : osink()
	{
		prov.insert("FeH");
		req.insert("comp");
		req.insert("XYZ");
	}
};
extern "C" opipeline_stage *create_module_feh() { return new os_FeH(); }	// Factory; called by opipeline_stage::create()

DECLARE_KERNEL(os_FeH_kernel(otable_ks ks, os_FeH_data par, gpu_rng_t rng, cint_t::gpu_t comp, cfloat_t::gpu_t XYZ, cfloat_t::gpu_t FeH))
size_t os_FeH::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactic XYZ coordinates exist in input
	//	- all stars are main sequence

	// fetch prerequisites
	cint_t   &comp  = in.col<int>("comp");
	cfloat_t &XYZ   = in.col<float>("XYZ");
	cfloat_t &FeH   = in.col<float>("FeH");

	comp_thin = icomp_thin;
	comp_thick = icomp_thick;
	comp_halo = icomp_halo;
	CALL_KERNEL(os_FeH_kernel, otable_ks(begin, end), *this, rng, comp, XYZ, FeH);
	return nextlink->process(in, begin, end, rng);
}

bool os_FeH::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	cfg.get(A[0],     "A0",     0.63f);
	cfg.get(sigma[0], "sigma0", 0.20f);
	cfg.get(offs[0],  "offs0",  0.00f);
	cfg.get(A[1],     "A1",     0.37f);
	cfg.get(sigma[1], "sigma1", 0.20f);
	cfg.get(offs[1],  "offs1",  0.14f);
	
	cfg.get(Hmu,      "Hmu",     500.f);
	cfg.get(muInf,    "muInf",  -0.82f);
	cfg.get(DeltaMu,  "deltaMu", 0.55f);

	// renormalize disk gaussian amplitudes to sum up to 1
	double sumA = A[0] + A[1];
	A[0] /= sumA; A[1] /= sumA;

	// Halo configuration
	cfg.get(sigma[2], "sigmaH",  0.30f);
	cfg.get(offs[2],  "offsH",  -1.46f);

	// Component IDs
	read_component_map(icomp_thin, cfg, "comp_thin");
	read_component_map(icomp_thick, cfg, "comp_thick");
	read_component_map(icomp_halo, cfg, "comp_halo");
// 	comp_thin  = componentMap.seqIdx(cfg.get("comp_thin"));
// 	comp_thick = componentMap.seqIdx(cfg.get("comp_thick"));
// 	comp_halo  = componentMap.seqIdx(cfg.get("comp_halo"));

	// Output model parameters
	MLOG(verb2) << "Component IDs (thin, thick, halo):        "<< icomp_thin << " " << icomp_thick << " " << icomp_halo;
	MLOG(verb2) << "Normalized disk amplitudes  (A[0], A[1]): "<< A[0] << " " << A[1];
	MLOG(verb2) << "Disk sigma          (sigma[0], sigma[1]): "<< sigma[0] << " " << sigma[1];
	MLOG(verb2) << "Disk offsets          (offs[0], offs[1]): "<< offs[0] << " " << offs[1];
	MLOG(verb2) << "Disk median Z dep. (muInf, deltaMu, Hmu): "<< muInf << " " << DeltaMu << " " << Hmu;
	MLOG(verb2) << "Halo gaussian              (muH, sigmaH): "<< offs[2] << " " << sigma[2];

	return true;
}
