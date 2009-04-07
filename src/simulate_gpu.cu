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

#ifndef __CONFIG_H
#define __CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <math.h>

#include <astro/constants.h>

#include "simulate_base.h"
#include "column.h"
#include "gpu.h"

namespace ct = column_types;
KERNEL(
	ks,
	os_FeH_kernel(otable_ks ks, os_FeH_data par, gpu_rng_t rng, ct::cint::gpu_t comp, ct::cfloat::gpu_t XYZ, ct::cfloat::gpu_t FeH),
	os_FeH_kernel,
	(ks, par, rng, comp, XYZ, FeH)
)
{
	uint32_t row = ks.row();
	if(row == (uint32_t)(-1)) { return; }
	rng.load(ks);

	float feh;
	switch(comp[row])
	{
		case 0: // BahcallSoneira_model::THIN:
		case 1: // BahcallSoneira_model::THICK:
		{
			// choose the gaussian to draw from
			float p = rng.uniform()*(par.A[0]+par.A[1]);
			int i = p < par.A[0] ? 0 : 1;

			// calculate mean
			float muD = par.muInf + par.DeltaMu*exp(-fabs(XYZ(row, 2))/par.Hmu);		// Bond et al. A2
			float aZ = muD - 0.067f;

			// draw
			feh = rng.gaussian(par.sigma[i]) + aZ + par.offs[i];
		} break;
		case 2: //BahcallSoneira_model::HALO:
			feh = par.offs[2] + rng.gaussian(par.sigma[2]);
			break;
		default:
			//THROW(ENotImplemented, "We should have never gotten here");
			feh = -9999.f;
			break;
	}
	FeH[row] = feh;
}


// equgal - Equatorial to Galactic coordinates
using namespace peyton;
typedef double Radians;
static const double angp = ctn::d2r * 192.859508333; //  12h 51m 26.282s (J2000)
static const double dngp = ctn::d2r * 27.128336111;  // +27d 07' 42.01" (J2000)
static const double l0 = ctn::d2r * 32.932;
static const double ce = 0.88998740217659689; // cos(dngp)
static const double se = 0.45598511375586859; // sin(dngp)

inline __device__ double2 galequ(const double2 lb)
{
	const double cb = cos(lb.y);
	const double sb = sin(lb.y);
	const double cl = cos(lb.x-l0);
	const double sl = sin(lb.x-l0);

//	// TODO: These should be precomputed constants
//	const double ce = cos(dngp);
//	const double se = sin(dngp);

	double2 r;
	r.x = atan2(
			cb*cl,
			sb*ce-cb*se*sl
		) + angp;
	r.y = asin(cb*ce*sl + sb*se);

	while(r.x < 0.) { r.x += ctn::pi2; }
	return r;
}

KERNEL(
	ks,
	os_gal2other_kernel(otable_ks ks, int coordsys, ct::cdouble::gpu_t lb0, ct::cdouble::gpu_t out),
	os_gal2other_kernel,
	(ks, coordsys, lb0, out)
)
{
	uint32_t row = ks.row();
	if(row == (uint32_t)(-1)) { return; }

	double2 lb, ret;

	// convert to radians
	lb.x = lb0(row, 0) * ctn::d2r;
	lb.y = lb0(row, 1) * ctn::d2r;

	// rotate to output coordinate system
	switch(coordsys)
	{
	case EQU:
		ret = galequ(lb);
		break;
	default:
		ret.x = ret.y = -9999.;
		break;
	}

	// convert to degrees
	out(row, 0) = ret.x / ctn::d2r;
	out(row, 1) = ret.y / ctn::d2r;
}
