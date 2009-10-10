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

#ifndef vel2pm_gpu_cu_h__
#define vel2pm_gpu_cu_h__

//
// Module data passed to device kernel
//
struct os_vel2pm_data
{
	int coordsys;		// output coordinate system (GAL, EQU)

	float vLSR,		// Local standard of rest velocity
 	      u0, v0, w0;	// Solar peculiar motion

 	float Rg;		// Distance to the Galactic center
};

//
// Device kernel implementation
//
#if !__CUDACC__ && !BUILD_FOR_CPU

	DECLARE_KERNEL(
		os_vel2pm_kernel(
			otable_ks ks, os_vel2pm_data par, gpu_rng_t rng, 
			cdouble_t::gpu_t lb0,
			cfloat_t::gpu_t XYZ,
			cfloat_t::gpu_t vcyl,
			cfloat_t::gpu_t pmlb
		)
	);

#else // #if !__CUDACC__ && !BUILD_FOR_CPU

	#include <astro/constants.h>
	#include <astro/types.h>

	__device__ inline double rad(double deg) { return 0.017453292519943295769  * deg; } // convert deg to rad
	__device__ inline float  radf(float deg) { return 0.017453292519943295769f * deg; } // convert deg to rad

	/* convert cartesian to celestial coordinates */
	__device__ void xyz2lbr(peyton::Radians &l, peyton::Radians &b, float &r, float x, float y, float z)
	{
		using namespace peyton;	// for mathematical constants

		r = sqrt(x*x + y*y + z*z);
		b = asin(z / r);
		l = atan2(y, x);
		if(l < 0) l += ctn::twopi;
	}

	/* convert celestial to cartesian coordinates */
	__device__ void lbr2xyz(float &x, float &y, float &z, peyton::Radians l, peyton::Radians b, float r)
	{
		x = r*cos(l)*cos(b);
		y = r*sin(l)*cos(b);
		z = r*sin(b);
	}

	// convert cartesian _Galactic_ (not galactocentric!!) velocities to vl,vr,vb wrt. the observer (velocity units are unimportant)
	__device__ void vel_xyz2lbr(float &vl, float &vb, float &vr, const float vx, const float vy, const float vz, const float l, const float b)
	{
		float cl, sl, cb, sb;
		cl = cosf(l);
		sl = sinf(l);
		cb = cosf(b);
		sb = sinf(b);

		float tmp;
		vl  = -vx*sl  + vy*cl;
		tmp =  vx*cl  + vy*sl;
		vr  =  cb*tmp + vz*sb;
		vb  = -sb*tmp + vz*cb;
	}

	// convert vl,vr,vb velocities wrt. the observer to cartesian _Galactic_ (not galactocentric!!) velocities (velocity units are unimportant)
	__device__ void vel_lbr2xyz(float &vx, float &vy, float &vz, const float vl, const float vb, const float vr, const float l, const float b)
	{
		float cl, sl, cb, sb;
		cl = cosf(l);
		sl = sinf(l);
		cb = cosf(b);
		sb = sinf(b);

		float tmp = sb*vb - cb*vr;
		vx = -sl*vl - cl*tmp;
		vy =  cl*vl - sl*tmp;
		vz =  cb*vb + sb*vr;
	}

	// convert galactocentric cylindrical to galactocentric cartesian velocities
	__device__ inline void vel_cyl2xyz(float &vx, float &vy, float &vz, const float vr, const float vphi, const float vz0, const float X, const float Y)
	{
		float rho = sqrtf(X*X + Y*Y);
		float cphi = X / rho;
		float sphi = Y / rho;

		vx = -sphi * vphi + cphi * vr;
		vy =  cphi * vphi + sphi * vr;
		vz = vz0;
	}

	/*
		Convert proper motions in one coordinate system to de. The coordinate
		system orientations are given by specifying the pole of the source coordinate
		system in destination coordinates, and the longitude of ascending node
		of the source coordinate system, in source coordinates.

		transform_pmf: uses single-precision floating point arithmetic.

		Input:		l,b in radians; vl,vb (units not important)
		Output:		vlout, vbout in destination coordinate system
	*/
	__device__ void transform_pm(float &vlout, float &vbout, float l, float b, float vl, float vb,
		const double ce, const double se, const double l0)
	{
		double cb = cos(b);
		double sb = sin(b);
		double cl = cos(l - l0);
		double sl = sin(l - l0);

		double tmp1 = cb * se - ce * sb * sl; // ??
		double tmp2 = ce * sb - cb * se * sl; // \propto cos(alpha)
		double tmp4 = cb * cl;                // \propto sin(alpha)
		double tmp3 = sb * se + cb * ce * sl; // sin(delta)

		if(fabs(tmp3) > 0.9999999)
		{
			// we're practically at the north/south pole.
			// better signal that computation here will be incorrect,
			// than give an incorrect result.
			vlout = vbout = 999.99;
		}
		else
		{
			double denom1 = 1. / ( tmp4*tmp4 + tmp2*tmp2 );
			double cd = sqrt(1. - tmp3*tmp3); // cos(delta);

			vl /= cb;
			vlout = (-ce * cl * vb + cb * tmp1 * vl) * denom1;
			vbout = (tmp1 * vb  + cb * ce * cl * vl) / cd;
			vlout *= cd;
		}
	}

	__device__ void transform_pmf(float &vlout, float &vbout, float l, float b, float vl, float vb,
		const float ce, const float se, const float l0)
	{
		float cb = cosf(b);
		float sb = sinf(b);
		float cl = cosf(l - l0);
		float sl = sinf(l - l0);

		float tmp1 = cb * se - ce * sb * sl; // ??
		float tmp2 = ce * sb - cb * se * sl; // \propto cos(alpha)
		float tmp3 = sb * se + cb * ce * sl; // sin(delta)
		float cd = sqrtf(1 - sqrf(tmp3)); // cos(delta);

		vl /= cb;
		vlout = (-ce * cl * vb + cb * tmp1 * vl) / ( sqrf(cb * cl) + sqrf(tmp2) );
		vbout = (tmp1 * vb  + cb * ce * cl * vl) / cd;
		vlout *= cd;
	}

	/*
		Convert proper motions in galactic coordinate system to proper motions in equatorial
		coordinate system. Uses single-precision floating point arithmetic.

		Input:		l,b in radians; vl,vb (units not important)
		Output:		vra,vdec (in units of vl,vb)
	*/
	__device__ void pm_galequf(float &vra, float &vdec, float l, float b, float vl, float vb)
	{
		using namespace float_galequ_constants;
		transform_pm(vra, vdec, l, b, vl, vb, ce, se, l0);
	}

	__device__ void pm_galequ(float &vra, float &vdec, float l, float b, float vl, float vb)
	{
		using namespace galequ_constants;
		transform_pm(vra, vdec, l, b, vl, vb, ce, se, l0);
	}

	__device__ void pm_equgalf(float &vl, float &vb, float ra, float dec, float vra, float vdec)
	{
		using namespace float_galequ_constants;

		// transform_pm(vl, vb, ra, dec, vra, vdec, ce, se, angp - 0.5*ctn::pi);
		transform_pm(vl, vb, ra, dec, vra, vdec, ce, se, angp - halfpi);
	}

	__device__ void pm_equgal(float &vl, float &vb, float ra, float dec, float vra, float vdec)
	{
		using namespace galequ_constants;
		using namespace peyton;	// for mathematical constants

		// transform_pm(vl, vb, ra, dec, vra, vdec, ce, se, angp - 0.5*ctn::pi);
		transform_pm(vl, vb, ra, dec, vra, vdec, ce, se, angp - halfpi);
	}

	__device__ float3 vcyl2pm(float l, float b, float vx, float vy, float vz, float X, float Y, float Z, os_vel2pm_data &par)
	{
	#if 0
		// These must produce (mu_l, mu_b, v_radial) = (12.72 mas/yr, -14.75 mas/yr, -27.34 km/s)
		X = 8100; Y = 100; Z = 3000;
		vx = 10; vy = 50; vz = -30;
		l = atan2(-Y, 8000-X);
		b = asin(Z / sqrt(sqrf(8000-X) + sqrf(Y) + sqrf(Z)));
		#if __DEVICE_EMULATION__
		printf("%f %f\n", deg(l), deg(b));
		#endif
	#endif
		// convert the velocities from cylindrical to galactocentric cartesian system
		float3 pm;	// pm.x == mu_l, pm.y == mu_b, pm.z == radial velocity (km/s)
		vel_cyl2xyz(pm.x, pm.y, pm.z,   vx, vy, vz,   X, Y);

		// switch to Solar coordinate frame
		pm.x -= par.u0;
		pm.y -= par.v0 + par.vLSR;
		pm.z -= par.w0;

		// convert to velocities wrt. the observer
		vel_xyz2lbr(pm.x, pm.y, pm.z,   -pm.x, -pm.y, pm.z,  l, b);

		// convert to proper motions
		float D = sqrtf(sqrf(8000.f-X) + sqrf(Y) + sqrf(Z));
		pm.x /= 4.74 * D*1e-3;	// proper motion in mas/yr (4.74 km/s @ 1kpc is 1mas/yr)
		pm.y /= 4.74 * D*1e-3;

		return pm;
	}

	KERNEL(
		ks, 0, 
		os_vel2pm_kernel(
			otable_ks ks, os_vel2pm_data par, gpu_rng_t rng, 
			cdouble_t::gpu_t lb0,
			cfloat_t::gpu_t XYZ,
			cfloat_t::gpu_t vcyl,
			cfloat_t::gpu_t pmout),
		os_vel2pm_kernel,
		(ks, par, rng, lb0, XYZ, vcyl, pmout)
	)
	{

		for(uint32_t row=ks.row_begin(); row < ks.row_end(); row++)
		{
			// fetch prerequisites
			// NOTE: Single precision should be sufficient here for (l,b)
			float l = radf(lb0(row, 0));
			float b = radf(lb0(row, 1));
			float X = par.Rg - XYZ(row, 0);
			float Y = -XYZ(row, 1);
			float Z = XYZ(row, 2);
			float vx = vcyl(row, 0);
			float vy = vcyl(row, 1);
			float vz = vcyl(row, 2);

			float3 pm = vcyl2pm(l, b, vx, vy, vz, X, Y, Z, par);

			// rotate to output coordinate system
			switch(par.coordsys)
			{
			case GAL:
				break;
			case EQU:
				pm_galequ(pm.x, pm.y, l, b, pm.x, pm.y);
				break;
			default:
				pm.x = pm.y = pm.z = -9.99;
				break;
			}

			pmout(row, 0) = pm.x;
			pmout(row, 1) = pm.y;
			pmout(row, 2) = pm.z;
		}
	}

#endif // #else (!__CUDACC__ && !BUILD_FOR_CPU)

#endif // vel2pm_gpu_cu_h__
