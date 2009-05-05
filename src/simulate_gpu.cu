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
	ks, 3*4,
	os_FeH_kernel(
		otable_ks ks, os_FeH_data par, gpu_rng_t rng, 
		ct::cint::gpu_t comp, 
		ct::cfloat::gpu_t XYZ, 
		ct::cfloat::gpu_t FeH),
	os_FeH_kernel,
	(ks, par, rng, comp, XYZ, FeH)
)
{
	rng.load(ks);
	uint32_t tid = threadID();
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
//		FeH[row] = -((threadID() % 300) / 100.f);
// 		float temp = -rng.uniform();
// 		FeH[row] = temp;
// 		continue;

//		feh = -comp[row]; FeH[row] = feh; continue;
		float feh;
		int component = comp[row];
#if 1
		if(component < 2)
		{
			// choose the gaussian to draw from
			float p = rng.uniform()*(par.A[0]+par.A[1]);
			int i = p < par.A[0] ? 0 : 1;

			// calculate mean
			float muD = par.muInf + par.DeltaMu*exp(-fabs(XYZ(row, 2))/par.Hmu);		// Bond et al. A2
			float aZ = muD - 0.067f;

			// draw
			feh = rng.gaussian(par.sigma[i]) + aZ + par.offs[i];
		}
		else if(component == 2)
		{
			feh = par.offs[2] + rng.gaussian(par.sigma[2]);
		}
		else
		{
			feh = -9999.f;
		}
#else
		/*
			This chunk of code appears to compile incorrectly with CUDA 2.1 nvcc (V0.2.1221).
			Expected result of running with below: feh=1 when comp=1, 0 otherwise.
			Actual result: feh=comp when comp=1,2
		*/
		feh = 0.f;
		switch(component)
		{
/*			case 2: //BahcallSoneira_model::HALO:
				feh = -2.f;
				break;*/
			case 3: // BahcallSoneira_model::THIN:
				feh = -3.0f;
				break;
			case 1: // BahcallSoneira_model::THICK:
				feh = component;
				break;
			case 4: // BahcallSoneira_model::THIN:
				feh = -4.0f;
				break;
/*			default:
				//THROW(ENotImplemented, "We should have never gotten here");
				feh = -9.f;
				break;*/
		}
#endif
		FeH[row] = feh;
	}
	rng.store(ks);
}

inline __device__ double sqrDouble(double x) { return x*x;}

static const double degToRadCoef=0.0174532925;//PI/180
inline __device__ double radDouble(double degrees) { return degrees*degToRadCoef; }

inline __device__ void vel_cyl2xyz(
		float &vx, float &vy, float &vz, 
		const float vr, const float vphi, const float vz0, 
		const float X, const float Y)
{
	// convert galactocentric cylindrical to cartesian velocities
	float rho = sqrt(X*X + Y*Y);
	float cphi = X / rho;
	float sphi = Y / rho;

	vx = -sphi * vr + cphi * vphi;
	vy =  cphi * vr + sphi * vphi;
	vz = vz0;
}

// convert cartesian galactocentric velocities to vl,vr,vb wrt. the observer
inline __device__ void vel_xyz2lbr(
		float &vl, float &vr, float &vb, 
		const float vx, const float vy, const float vz, 
		const double l, const double b)
{
	double cl, sl, cb, sb;
	cl = cos(l);
	sl = sin(l);
	cb = cos(b);
	sb = sin(b);

	float tmp;
	vl  =  vx*sl  - vy*cl;
	tmp =  vx*cl  + vy*sl;
	vr  = -cb*tmp + vz*sb;
	vb  =  sb*tmp + vz*cb;
}

/// !!!!! 3*4 FIX-.. calculate correct value !!!!!!!!!
KERNEL(
	ks, 3*4, 
	os_vel2pm_kernel(
		otable_ks ks, os_vel2pm_data par, gpu_rng_t rng, 
		ct::cdouble::gpu_t lb0, 
		ct::cfloat::gpu_t XYZ,
		ct::cfloat::gpu_t vcyl,  
		ct::cfloat::gpu_t pmlb),
	os_vel2pm_kernel,
	(ks, par, rng, lb0, XYZ, vcyl,pmlb)
)
{
	rng.load(ks);
	uint32_t tid = threadID();

	for(uint32_t row=ks.row_begin(); row < ks.row_end(); row++)
	{
		// fetch prerequisites
		double l = radDouble(lb0(row, 0));
		double b = radDouble(lb0(row, 1));
		float X = XYZ(row, 0);
		float Y = XYZ(row, 1);
		float Z = XYZ(row, 2);
		float vx = vcyl(row, 0);
		float vy = vcyl(row, 1);
		float vz = vcyl(row, 2);

		// convert the velocities from cylindrical to galactocentric cartesian system
		float pm[3];
		vel_cyl2xyz(pm[0], pm[1], pm[2],   vx, vy, vz,   X, Y);

		// switch to Solar coordinate frame
		pm[0] -= par.u0;
		pm[1] -= par.v0 + par.vLSR;
		pm[2] -= par.w0;

		// convert to velocities wrt. the observer
		vel_xyz2lbr(pm[0], pm[1], pm[2],   pm[0], pm[1], pm[2],  l, b);

		// convert to proper motions
		float D = sqrt(sqrDouble(X) + sqrDouble(Y) + sqrDouble(Z));
		pm[0] /= 4.74 * D*1e-3;	// proper motion in mas/yr (4.74 km/s @ 1kpc is 1mas/yr)
		pm[1] /= 4.74 * D*1e-3;

		// rotate to output coordinate system
		switch(par.coordsys)
		{
		case GAL:
//			array_copy(s.pmlb(), pm, 3);
			pmlb(row, 0) = pm[0];
			pmlb(row, 1) = pm[1];
			pmlb(row, 2) = pm[2];
			break;
		case EQU:
			//THROW(EAny, "Output in equatorial system not implemented yet.");
			//array_copy(s.pmradec(), pm, 3);
			break;
		default:
			//THROW(EAny, "Unknown coordinate system [id=" + str(coordsys) + "] requested");
			break;
		}
	}

	rng.store(ks);
}


// equgal - Equatorial to Galactic coordinates
using namespace peyton;
typedef double Radians;
static const double angp = ctn::d2r * 192.859508333; //  12h 51m 26.282s (J2000)
//static const double dngp = ctn::d2r * 27.128336111;  // +27d 07' 42.01" (J2000)
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
	ks, 0,
	os_gal2other_kernel(otable_ks ks, int coordsys, ct::cdouble::gpu_t lb0, ct::cdouble::gpu_t out),
	os_gal2other_kernel,
	(ks, coordsys, lb0, out)
)
{
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
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
}

KERNEL(
	ks, 0,
	os_fixedFeH_kernel(otable_ks ks, float fixedFeH, ct::cfloat::gpu_t FeH),
	os_fixedFeH_kernel,
	(ks, fixedFeH, FeH)
)
{
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		FeH[row] = fixedFeH;
	}
}

#if __CUDACC__
__global__ void clockInstruction(int *dest, int *clk, int a, int b)
{
	clock_t t0 = clock();
	*dest = a + b;
	clock_t t1 = clock();
	*clk = t1 - t0;
}
#endif
