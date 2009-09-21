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
#include "vel2pm_gpu.cu.h"

#include "spline.h"
#include "analysis.h"
#include "io.h"
#include <fstream>

#include <astro/system/config.h>
#include <astro/useall.h>

// os_vel2pm -- Convert velocities to proper motions
class os_vel2pm : public osink, public os_vel2pm_data
{
protected:
	std::string output_col_name;

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("vel2pm"); return s; }

	os_vel2pm() : osink()
	{
		coordsys=GAL;
		req.insert("lb");
		req.insert("XYZ");
		req.insert("vcyl");
	}
};
	
extern "C" opipeline_stage *create_module_vel2pm() { return new os_vel2pm(); }	// Factory; called by opipeline_stage::create()

size_t os_vel2pm::process(otable &in, size_t begin, size_t end, rng_t &rng)
{ 
	// ASSUMPTIONS:
	//	vcyl() velocities are in km/s, XYZ() distances in parsecs
	//
	// OUTPUT:
	//	Proper motions in mas/yr for l,b directions in pm[0], pm[1]
	//	Radial velocity in km/s in pm[2]
	cdouble_t &lb0 = in.col<double>("lb");
	cfloat_t  &XYZ  = in.col<float>("XYZ");
	cfloat_t  &vcyl = in.col<float>("vcyl");
	cfloat_t  &pmlb = in.col<float>(output_col_name);

	CALL_KERNEL(os_vel2pm_kernel, otable_ks(begin, end), *this, rng, lb0, XYZ, vcyl,pmlb);
	return nextlink->process(in, begin, end, rng);
}

bool os_vel2pm::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	//if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	std::string cs;
	cfg.get(cs, "coordsys", "gal");
	     if(cs == "gal") { coordsys = GAL; output_col_name = "pmlb"; }
	else if(cs == "equ") { coordsys = EQU; output_col_name = "pmradec"; }
	else { THROW(EAny, "Unknown coordinate system (" + cs + ") requested."); }
	prov.insert(output_col_name);

	// LSR and peculiar Solar motion
	cfg.get(vLSR, "vLSR", -220.0f);
	cfg.get(u0,   "u0",    -10.0f);
	cfg.get(v0,   "v0",     -5.3f);
	cfg.get(w0,   "w0",      7.2f);

	return true;
}
