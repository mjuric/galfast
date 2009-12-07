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

#ifndef brokenPowerLaw_h__
#define brokenPowerLaw_h__

#include "skygen.h"

#include "model_powerLawEllipsoid.h"

// luminosity function texture reference
//DEFINE_TEXTURE(powerLawEllipsoidLF, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

// exponential disk model
struct ALIGN(16) brokenPowerLaw : public powerLawEllipsoid
{
public:
	struct state : public powerLawEllipsoid::state {};

public:
	// Note: almost all parameters are inherited from powerLawEllipsoid

	static const int MAXPOWERLAWS = 10;
	int nbreaks;				// number of power law breaks
	float   nArr[MAXPOWERLAWS];		// power law indices
	float	fArr[MAXPOWERLAWS];		// piece-wise normalizations (computed in load())
	float rbreakSq[MAXPOWERLAWS-1];		// ellipsoid radii below which they kick in

public:
	void load(host_state_t &hstate, const peyton::system::Config &cfg);

public:
	__host__ __device__ float rho(float x, float y, float z) const
	{
		using namespace cudacc;

		// translate and rotate into ellipsoid-centric frame
		float3 v = { c.x - x, c.y - y, z - c.z };
		v = matmul3d(rot, v);

		// test for whether we're within the ellipsoid
		float DSq = sqr(v.x) + sqr(v.y) + sqr(v.z);
		if(rmaxSq && (DSq < rminSq || DSq >= rmaxSq)) { return 0.f; }

		// find the power law that applies in this region of space
		int k = 0;
		float rSq = sqr(v.x) + sqr(v.y/ba) + sqr(v.z/ca);
		while(k != nbreaks && rSq > rbreakSq[k])
		{
			k++;
		}

		float rho = f * fArr[k] * powf(rSq, 0.5f*nArr[k]);
		return rho;
	}

public:
	__device__ void setpos(state &s, float x, float y, float z) const
	{
		s.rho = rho(x, y, z);
	}

	__device__ float rho(state &s, float M) const
	{
		return powerLawEllipsoid::rho(s, M);
	}
};

MODEL_IMPLEMENTATION(brokenPowerLaw);

#endif // #ifndef brokenPowerLaw_h__
