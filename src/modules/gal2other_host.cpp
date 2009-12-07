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
#include "gal2other_gpu.cu.h"

#include "spline.h"
#include "analysis.h"
#include "io.h"
#include <fstream>

#include <astro/system/config.h>
#include <astro/useall.h>

// convert between coordinate systems
class os_gal2other : public osink
{
public:
	int coordsys;

public:
	size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("gal2other"); return s; }
	virtual double ordering() const { return ord_database; }

	os_gal2other() : osink(), coordsys(GAL)
	{
		req.insert("lb");
	}
};
extern "C" opipeline_stage *create_module_gal2other() { return new os_gal2other(); }	// Factory; called by opipeline_stage::create()

size_t os_gal2other::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	cdouble_t &lb   = in.col<double>("lb");

	if(coordsys == EQU)
	{
		cdouble_t &out = in.col<double>("radec");
		CALL_KERNEL(os_gal2other_kernel, otable_ks(begin, end), coordsys, lb, out);
	}

	return nextlink->process(in, begin, end, rng);
}

bool os_gal2other::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	//if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	std::string cs;
	cfg.get(cs, "coordsys", "gal");
	     if(cs == "gal") { coordsys = GAL; /* -- noop -- */ }
	else if(cs == "equ") { coordsys = EQU; prov.insert("radec[2]"); }
	else { THROW(EAny, "Unknown coordinate system (" + cs + ") requested."); }

	return true;
}
