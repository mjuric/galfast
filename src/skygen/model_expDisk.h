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

#ifndef expDisk_h__
#define expDisk_h__

#include "skygen.h"
#include <cstdlib>

// luminosity function texture reference
DEFINE_TEXTURE(expDiskLF, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

// ccw rotation of vector (x,y) by angle alpha
// to rotate in cw direction, just flip the sign of sina
__device__ inline void rotate2D(float &xout, float &yout, float x, float y, float cosa, float sina)
{
	xout = cosa * x - sina * y;
	yout = sina * x + cosa * y;
}

// exponential disk model
struct ALIGN(16) expDisk : public modelConcept
{
public:
	struct ALIGN(16) host_state_t
	{
		cuxTexture<float> lf;
	};

public:
	float f;		// overal factor computed by the setup routines that ensures rho(x_lf,y_lf,z_lf)=f where f is the value from the config file
	float l, h;		// scale length, heights
	float r_cut2;		// cylindrical radius of sharp density cutoff

#if 0
	float Rg;		// distance from the Sun to Galactic center
	float rsun, zsun;	// cylindric galactocentric r and z coordinates of the Sun
	float asin, acos;	// sin/cos of the angle between the Sun and the plane of the Galaxy as viewed from the Galactic center (i.e., sin(a) = zsun / rsun)
#else
	float M[3][3];		// Galactic->disk coordinate system rotation matrix. Use with transform().
	float3 T;		// Negative disk center offset (in Galactic coordinates). Use with transform().
#endif
public:
	struct state
	{
		float rho;
	};
	void load(host_state_t &hstate, const peyton::system::Config &cfg);
	void prerun(host_state_t &hstate, bool draw);
	void postrun(host_state_t &hstate, bool draw);
	bool hint_absmag(host_state_t &hstate, float &M0, float &M1) const;

public:
	__device__ float rho(float x, float y, float z) const
	{
#if 0
		printf("Coordinates before translation, rotation: %.9f %.9f %.9f\n", xx, yy, zz);
#endif

#if 0
		// transform to galactocentric coordinate system
		x = Rg-x;
		rotate2D(x, z, x, z, acos, asin);
#else
		float3 v = {x, y, z};
		v = transform(v, T, M);
#endif

#if 0
		printf("Coordinates after translation, rotation: %.9f %.9f %.9f\n", v.x, v.y, v.z);
		printf("% 11.9f % 11.9f % 11.9f\n", M[0][0], M[0][1], M[0][2]);
		printf("% 11.9f % 11.9f % 11.9f\n", M[1][0], M[1][1], M[1][2]);
		printf("% 11.9f % 11.9f % 11.9f\n", M[2][0], M[2][1], M[2][2]);
		abort();
#endif
#if 0
		printf("Coordinates after translation, rotation: %.9f %.9f %.9f\n", x, y, z);
		printf("% 11.9f % 11.9f\n", acos, -asin);
		printf("% 11.9f % 11.9f\n", asin,  acos);
		abort();
#endif

		// check for hard cutoffs
		float r2 = v.x*v.x + v.y*v.y;
		if(r_cut2 && r2 + v.z*v.z > r_cut2) { return 0.f; }

		// compute the profile
		float r = sqrtf(r2);
		float rho = f * expf(-r/l - fabsf(v.z)/h);

		return rho;
	}

public:
	__device__ void setpos(state &s, float x, float y, float z) const
	{
		s.rho = rho(x, y, z);
	}

	__device__ float rho(state &s, float M) const
	{
		float phi = TEX1D(expDiskLF, M);
		return phi * s.rho;
	}
};

MODEL_IMPLEMENTATION(expDisk);

#endif // #ifndef expDisk_h__
