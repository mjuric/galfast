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

#ifndef model_densityCube_h__
#define model_densityCube_h__

#include "skygen.h"

// luminosity function texture
DEFINE_TEXTURE(densityCubeLF, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

// 3D density texture
DEFINE_TEXTURE(densityCubeTex, float, 3, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

// exponential disk model
struct ALIGN(16) densityCube : public modelConcept
{
public:
	// this remains on the host
	struct ALIGN(16) host_state_t
	{
		cuxTexture<float, 1> lf;	// LF
		cuxTexture<float, 3> den;	// density cube
	};

protected:
	// uploaded to a GPU __constant__
	
	float f;	// scale factor by which to multiply sampled densities
//	int comp;	// component ID for this model

	float Rg;	// distance to the Galactic center

public:
	struct state
	{
		float rho;	// sampled density, before multiplying by LF
	};
	void load(host_state_t &hstate, const peyton::system::Config &cfg);
	void prerun(host_state_t &hstate, bool draw);
	void postrun(host_state_t &hstate, bool draw);
	bool hint_absmag(host_state_t &hstate, float &M0, float &M1) const;

public:
	__device__ void setpos(state &s, float x, float y, float z) const
	{
		x = Rg - x;	// convert to galactocentric
		y *= -1;	// ...

		s.rho = f * TEX3D(densityCubeTex, x, y, z);
	}

	__device__ float rho(state &s, float M) const
	{
		float phi = TEX1D(densityCubeLF, M);
		return phi * s.rho;
	}
};

MODEL_IMPLEMENTATION(densityCube);

#endif // #ifndef model_densityCube_h__
