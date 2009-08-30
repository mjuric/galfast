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

// luminosity function texture reference
DEFINE_TEXTURE(expDiskLF, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

// exponential disk model
struct ALIGN(16) expDisk : public modelConcept
{
public:
	struct ALIGN(16) host_state_t
	{
		texptr<float> lf;
	};

protected:
	float f, l, h, z0;
	float r_cut2;

	int comp;

public:
	struct state
	{
		float rho;
	};
	void load(host_state_t &hstate, const peyton::system::Config &cfg);
	void prerun(host_state_t &hstate, bool draw);
	void postrun(host_state_t &hstate, bool draw);

protected:
	__device__ float rho(float x, float y, float z, float M) const
	{
		float r2 = x*x + y*y;
		if(r2 + z*z > r_cut2) { return 0.f; }

		float r = sqrtf(r2);
		float rho = f * expf((Rg()-r)/l  + (fabsf(z0) - fabsf(z + z0))/h);

		return rho;
	}

public:
	__device__ void setpos(state &s, float x, float y, float z) const
	{
		s.rho = rho(x, y, z, 0.f);
	}

	__device__ float rho(state &s, float M) const
	{
		float phi = TEX1D(expDiskLF, M);
		return phi * s.rho;
	}

	__device__ int component(float x, float y, float z, float M, gpuRng::constant &rng) const
	{
		return comp;
	}
};

MODEL_IMPLEMENTATION(expDisk);

#endif // #ifndef expDisk_h__
