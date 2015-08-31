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

#include "module_lib.h"
#include "../pipeline.h"
#include "GaussianFeH_gpu.cu.h"

#include <numeric>

#include <astro/system/config.h>
#include <astro/useall.h>

// os_GaussianFeH -- Generate Fe/H based on Ivezic et al (2008)
class os_GaussianFeH : public osink
{
protected:
	os_GaussianFeH_data par;

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("GaussianFeH"); return s; }
	virtual double ordering() const { return ord_feh; }

	os_GaussianFeH() : osink()
	{
		prov.insert("FeH");
		prov.insert("FeH_comp{type=int; fmt=%3d;}");
		req.insert("comp");
		req.insert("XYZ");
	}
};
extern "C" opipeline_stage *create_module_gaussianfeh() { return new os_GaussianFeH(); }	// Factory; called by opipeline_stage::create()

DECLARE_KERNEL(os_GaussianFeH_kernel(otable_ks ks, bit_map applyToComponents, os_GaussianFeH_data par, gpu_rng_t rng,
		cint_t::gpu_t comp, cint_t::gpu_t hidden,
		cfloat_t::gpu_t XYZ,
		cint_t::gpu_t   FeH_comp,
		cfloat_t::gpu_t FeH));

size_t os_GaussianFeH::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactic XYZ coordinates exist in input
	//	- all stars are main sequence

	// fetch prerequisites
	cint_t   &comp  = in.col<int>("comp");
	cint_t   &hidden= in.col<int>("hidden");
	cfloat_t &XYZ   = in.col<float>("XYZ");
	cfloat_t &FeH   = in.col<float>("FeH");
	cint_t   &FeH_comp = in.col<int>("FeH_comp");

	CALL_KERNEL(os_GaussianFeH_kernel, otable_ks(begin, end), applyToComponents, par, rng, comp, hidden, XYZ, FeH_comp, FeH);
	return nextlink->process(in, begin, end, rng);
}

bool os_GaussianFeH::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	read_component_map(applyToComponents, cfg);

	std::vector<float> mean  = cfg.get("mean");
	std::vector<float> sigma = cfg.get("sigma");

	ASSERT(mean.size() > 0);				// TODO: change to throw
	ASSERT(mean.size()  < os_GaussianFeH_data::MAXCOMP);	// TODO: change to throw
	ASSERT(sigma.size() == mean.size());			// TODO: change to throw
	par.ncomp = mean.size();
	FOR(0, par.ncomp) { par.mean[i] = mean[i]; }
	FOR(0, par.ncomp) { par.sigma[i] = sigma[i]; }

	// load the weights of each Gaussian, and compute the cumulative probability vector pcum
	std::vector<double> wt = par.ncomp > 1 ? cfg.get("wt") : std::vector<double>(1, 1.);
	ASSERT(wt.size() == par.ncomp);
	std::partial_sum(wt.begin(), wt.end(), par.pcum);
	ASSERT(par.pcum[par.ncomp-1] > 0);			// TODO: change to throw
	FOR(0, par.ncomp) { par.pcum[i] /= par.pcum[par.ncomp-1]; }

	// Output model parameters
	MLOG(verb1) << "Metallicity: (multi)Gaussian for components " << applyToComponents << "   ## " << instanceName();
	FOR(0, par.ncomp)
	{
		MLOG(verb2) << "           : (mu, sigma, pcum) = " << par.mean[i] << ", " << par.sigma[i] << ", " << par.pcum[i];
	}

	return true;
}
