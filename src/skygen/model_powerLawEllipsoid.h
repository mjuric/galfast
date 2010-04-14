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

#ifndef powerLawEllipsoid_h__
#define powerLawEllipsoid_h__

#include "skygen.h"

// luminosity function texture reference
DEFINE_TEXTURE(powerLawEllipsoidLF, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

// exponential disk model
struct ALIGN(16) powerLawEllipsoid : public modelConcept
{
public:
	struct ALIGN(16) host_state_t
	{
		cuxTexture<float> lf;
	};

public:
	float3 T;		// Translation parameter for transform() (negative distance to ellipsoid center, in Galactic coordinates)
	float f, n, ba, ca;	// Density normalization, power law index, b/a and c/a axes ratios
	float rminSq, rmaxSq;	// Minimum/maximum radius from the center of the ellipsoid where the density is nonzero (the square of)
	float rot[3][3];	// Rotation matrix of the ellipsoid (transforms points from Galactic to ellipsoid coordinate system)

	int dummy[1];		// for padding (nvcc vs. gcc structure alignment issues)

protected:
	void load_geometry(host_state_t &hstate, const peyton::system::Config &cfg);

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
	__host__ __device__ float rho(float x, float y, float z) const
	{
		using namespace cudacc;

#if 0
		printf("Coordinates before translation, rotation: %.9f %.9f %.9f\n", x, y, z);
#endif

		// translate and rotate into ellipsoid-centric frame
#if 0
		float3 v = { c.x - x, c.y - y, z - c.z };
		v = matmul3d(rot, v);
#else
		float3 v = {x, y, z};
		v = transform(v, T, rot);
#endif

#if 0
		printf("Coordinates after translation, rotation: %.9f %.9f %.9f\n", v.x, v.y, v.z);
		printf("% 11.9f % 11.9f % 11.9f\n", rot[0][0], rot[0][1], rot[0][2]);
		printf("% 11.9f % 11.9f % 11.9f\n", rot[1][0], rot[1][1], rot[1][2]);
		printf("% 11.9f % 11.9f % 11.9f\n", rot[2][0], rot[2][1], rot[2][2]);
		abort();
#endif

		// test for whether we're within the ellipsoid
		float DSq = sqr(v.x) + sqr(v.y) + sqr(v.z);
		if(rmaxSq && (DSq < rminSq || DSq >= rmaxSq)) { return 0.f; }

		// compute the density
		float rSq = sqr(v.x) + sqr(v.y/ba) + sqr(v.z/ca);
		float rho = f*powf(rSq, 0.5f*n);
		return rho;
	}

public:
	__device__ void setpos(state &s, float x, float y, float z) const
	{
		s.rho = rho(x, y, z);
	}

	__device__ float rho(state &s, float M) const
	{
		float phi = TEX1D(powerLawEllipsoidLF, M);
		return phi * s.rho;
	}
};

MODEL_IMPLEMENTATION(powerLawEllipsoid);

#endif // #ifndef powerLawEllipsoid_h__
