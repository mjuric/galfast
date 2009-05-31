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

#include "gpc_cpp.h"

#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <sstream>

#include "simulate.h"
#include "projections.h"
#include "model.h"
#include "paralax.h"
#include "analysis.h"
//#include "dm.h"
#include "io.h"
#include "gpu.h"
#include "simulate_base.h"

#include <vector>
#include <map>
#include <string>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <astro/io/format.h>
#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

using namespace boost::lambda;
namespace ct = column_types;

#include "observe.h"


//size_t os_vel2pm::process(otable &in, size_t begin, size_t end, rng_t &rng)
DECLARE_KERNEL(
	os_vel2pm_kernel(
		otable_ks ks, os_vel2pm_data par, gpu_rng_t rng, 
		ct::cdouble::gpu_t lb0, 
		ct::cfloat::gpu_t XYZ,
		ct::cfloat::gpu_t vcyl,  
		ct::cfloat::gpu_t pmlb))

size_t os_vel2pm::process(otable &in, size_t begin, size_t end, rng_t &rng)
{ 
	// ASSUMPTIONS:
	//	vcyl() velocities are in km/s, XYZ() distances in parsecs
	//
	// OUTPUT:
	//	Proper motions in mas/yr for l,b directions in pm[0], pm[1]
	//	Radial velocity in km/s in pm[2]
	using namespace column_types;
	ct::cdouble &lb0 = in.col<double>("lb");
	ct::cfloat  &XYZ  = in.col<float>("XYZ");
	ct::cfloat  &vcyl = in.col<float>("vcyl");
	ct::cfloat  &pmlb = in.col<float>(output_col_name);

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
