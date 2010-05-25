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
#include "fixedFeH_gpu.cu.h"

#include <astro/system/config.h>
#include <astro/useall.h>

// os_fixedFeH -- Generate a fixed Fe/H
class os_fixedFeH : public osink
{
	protected:
		float fixedFeH;

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
		virtual const std::string &name() const { static std::string s("fixedFeH"); return s; }
		virtual double ordering() const { return ord_feh; }

		os_fixedFeH() : osink(), fixedFeH(0)
		{
			prov.insert("FeH");
		}
};
extern "C" opipeline_stage *create_module_fixedfeh() { return new os_fixedFeH(); }	// Factory; called by opipeline_stage::create()

size_t os_fixedFeH::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input
	//	- all stars are main sequence
	cfloat_t &FeH   = in.col<float>("FeH");
	cint_t  &comp   = in.col<int>("comp");
	cint_t  &hidden = in.col<int>("hidden");

	CALL_KERNEL(os_fixedFeH_kernel, otable_ks(begin, end), applyToComponents, fixedFeH, comp, hidden, FeH);
	return nextlink->process(in, begin, end, rng);
}

bool os_fixedFeH::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	read_component_map(applyToComponents, cfg);

	if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	cfg.get(fixedFeH, "FeH", 0.f);

	return true;
}
