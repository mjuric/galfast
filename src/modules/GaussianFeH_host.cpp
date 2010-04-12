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
#include "GaussianFeH_gpu.cu.h"

#include <astro/system/config.h>
#include <astro/useall.h>

// os_GaussianFeH -- Generate Fe/H based on Ivezic et al (2008)
class os_GaussianFeH : public osink
{
protected:
	float mean, sigma;

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("GaussianFeH"); return s; }
	virtual double ordering() const { return ord_feh; }

	os_GaussianFeH() : osink()
	{
		prov.insert("FeH");
		req.insert("comp");
		req.insert("XYZ");
	}
};
extern "C" opipeline_stage *create_module_gaussianfeh() { return new os_GaussianFeH(); }	// Factory; called by opipeline_stage::create()

DECLARE_KERNEL(os_GaussianFeH_kernel(otable_ks ks, bit_map applyToComponents, float mean, float sigma, gpu_rng_t rng,
		cint_t::gpu_t comp,
		cfloat_t::gpu_t XYZ,
		cfloat_t::gpu_t FeH));

size_t os_GaussianFeH::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactic XYZ coordinates exist in input
	//	- all stars are main sequence

	// fetch prerequisites
	cint_t   &comp  = in.col<int>("comp");
	cfloat_t &XYZ   = in.col<float>("XYZ");
	cfloat_t &FeH   = in.col<float>("FeH");

	CALL_KERNEL(os_GaussianFeH_kernel, otable_ks(begin, end), applyToComponents, mean, sigma, rng, comp, XYZ, FeH);
	return nextlink->process(in, begin, end, rng);
}

bool os_GaussianFeH::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	read_component_map(applyToComponents, cfg);

	mean  = cfg.get("mean");
	sigma = cfg.get("sigma");

	// Output model parameters
	MLOG(verb1) << "Metallicity: Gaussian for components " << applyToComponents << "   ## " << instanceName();
	MLOG(verb2) << "           : (mu, sigma) = " << mean      << ", " << sigma;

	return true;
}
