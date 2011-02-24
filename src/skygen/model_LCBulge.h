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

#ifndef LCBulge_h__
#define LCBulge_h__

#include "skygen.h"

#include "model_powerLawEllipsoid.h"

// luminosity function texture reference
DEFINE_TEXTURE(LCBulgeLF, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

// exponential disk model
struct ALIGN(16) LCBulge : public modelConcept
{
public:
	struct ALIGN(16) host_state_t
	{
		cuxTexture<float> lf;
	};

	struct state
	{
		float rho;
	};

public:
	// Exponential ellipsoid parameters for
	//
	//	rho = LF(M) * f * exp(-t/l)
	//	l = ( x^n + (y/ba)^n + (y/ca)^n )^(1/n)
	//
	// model. Note that LF(M) * f is precomputed on load.

	float ba, ca;		// ellipsoid axes ratios
	float n;		// power law index of the exponent
	float l;		// exponential scale

	// Positioning
	float3 T;		// Translation parameter for transform() (negative distance to ellipsoid center, in Galactic coordinates)
	float rot[3][3];	// Rotation matrix of the ellipsoid (transforms points from Galactic to ellipsoid coordinate system)

public:
	void load(host_state_t &hstate, const peyton::system::Config &cfg);
	void prerun(host_state_t &hstate, bool draw);
	void postrun(host_state_t &hstate, bool draw);
	bool hint_absmag(host_state_t &hstate, float &M0, float &M1) const;

public:
	__host__ __device__ float rho(float x, float y, float z) const
	{
		using namespace cudacc;

		// translate and rotate into ellipsoid-centric frame
		float3 v = {x, y, z};
		v = transform(v, T, rot);

		// compute the exponent and density
		float t = powf(powf(v.x, n) + powf(v.y / ba, n) + powf(v.z / ca, n), 1./n);
		float rho = expf(-t / l);

		return rho;
	}

public:
	__device__ void setpos(state &s, float x, float y, float z) const
	{
		s.rho = rho(x, y, z);
	}

	__device__ float rho(state &s, float M) const
	{
		float phi = TEX1D(LCBulgeLF, M);
		return phi * s.rho;
	}
};

MODEL_IMPLEMENTATION(LCBulge);

#endif // #ifndef LCBulge_h__
