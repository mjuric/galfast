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

/**
	Constants and common device functions shared by all modules
*/

#ifndef module_lib_h__
#define module_lib_h__

	#include "gpu.h"
	#include <astro/constants.h>

// 	#ifdef __CUDACC__
// 		// distance to the Galactic center
// 		__device__ __constant__ float Rg_gpu;
// 		__device__ inline float Rg() { return Rg_gpu; }
// 	#else
// 		// galactic center distance, CPU version
// 		extern double Rg_impl;
// 		inline double Rg() { return Rg_impl; }
// 	#endif

	// useful math constants
	static const double dbl_pi  = 3.14159265358979323846264338;
	static const double dbl_d2r = 0.01745329251994329576923691; // (pi/180.0)
	static const double dbl_r2d = 57.2957795130823208767981548; // (180.0/pi)
	static const float flt_pi  = 3.14159265358979323846264338f;
	static const float flt_d2r = 0.01745329251994329576923691f; // (pi/180.0)
	static const float flt_r2d = 57.2957795130823208767981548f; // (180.0/pi)

	namespace cudacc
	{
		inline __device__ float deg(float rd) { return flt_r2d * rd; }
		inline __device__ float rad(float dg) { return flt_d2r * dg; }
		template<typename T>
			__host__ __device__ inline __device__ T sqr(const T x) { return x*x; }
	}

	static const float ABSMAG_NOT_PRESENT = 99.999f;

	static const int GAL = 0;
	static const int EQU = 1;

	namespace galequ_constants
	{
		static const double angp = peyton::ctn::d2r * 192.859508333; //  12h 51m 26.282s (J2000)
		static const double dngp = peyton::ctn::d2r * 27.128336111;  // +27d 07' 42.01" (J2000)
		static const double l0   = peyton::ctn::d2r * 32.932;	// galactic longitude of ascending node of galactic coordinate system (where b=0, dec=0)
		static const double ce   = 0.88998740217659689; // cos(dngp)
		static const double se   = 0.45598511375586859; // sin(dngp)

		static const double halfpi = peyton::ctn::halfpi;
	};

	namespace float_galequ_constants
	{
		static const float angp = (float)galequ_constants::angp; //  12h 51m 26.282s (J2000)
		static const float dngp = (float)galequ_constants::dngp;  // +27d 07' 42.01" (J2000)
		static const float ce   = (float)galequ_constants::ce;
		static const float se   = (float)galequ_constants::se;
		static const float l0   = (float)galequ_constants::l0;

		static const float halfpi = (float)galequ_constants::halfpi;
	};

	//
	// Store a unit vector in terms of its sines/cosines
	//
	struct direction
	{
		float cl, cb, sl, sb;

		// Compute cartesian xyz position given direction and distance
		// Note: x=y=z=0 is the position of the Sun. The galactic center
		//       lies on the positive x axis. This *OPPOSITE* from the
		//       Juric et al. (2008) heliocentric convention, but more
		//       consistent with tradition (and math).
		__device__ inline float3 xyz(const float d) const
		{
			float3 ret;

			ret.x = d*cl*cb;
			ret.y = d*sl*cb;
			ret.z =    d*sb;

			return ret;
		};

		direction() {}
		direction(double l_, double b_)				// l,b are in radians
		: cl(cos(l_)), cb(cos(b_)), sl(sin(l_)), sb(sin(b_))
		{ }
	};

	float inline __device__ sqrf(float x) { return x*x; }

#endif
