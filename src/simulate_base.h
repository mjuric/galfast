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

#ifndef __simulate_base_h
#define __simulate_base_h

#include "gpu.h"

/**
	Things "C"-ish enough for CUDA to swallow (and to upload to GPU) go into this header.
*/

// returns numeric_limits::epsilon() for the type of x
#define EPSILON_OF(x) std::numeric_limits<typeof(x)>::epsilon()

static const float ABSMAG_NOT_PRESENT = 99.999f;

struct os_FeH_data
{
	float A[2], sigma[3], offs[3];
	float Hmu, muInf, DeltaMu;
	int comp_thin, comp_thick, comp_halo;
};

struct ALIGN(16) os_photometry_data
{
	int ncolors, bidx;	// number of colors, bootstrap band index

	float FeH0, dFeH;	// FeH ...
	float Mr0, dMr;		// ... and Mr texture coordinates (TODO: fully convert to cuxTexture)

	uint32_t comp0, comp1;	// component ID range [comp0, comp1) to which this photometry module will be asigning magnitudes

	static const int N_REDDENING = 17;
	float reddening[N_REDDENING];	// reddening coefficients for the loaded bands (NOTE: hardcoded maximum of 17 bands (16 colors))
};

struct os_vel2pm_data
{
	int coordsys;

	float vLSR,		// Local standard of rest velocity
 	      u0, v0, w0;	// Solar peculiar motion
};

static const int GAL = 0;
static const int EQU = 1;

struct farray5
{
	float data[5];
	
	__device__ float& operator [] (int i) { return data[i]; }
	__device__ const float& operator [] (int i) const { return data[i]; }
};

struct iarray5
{
	short int data[5];
	
	__device__ short int& operator [] (int i) { return data[i]; } 
	__device__ const short int& operator [] (int i) const { return data[i]; } 
};

struct i8array5
{
	char data[5];
	
	__device__ char& operator [] (int i) { return data[i]; } 
	__device__ const char& operator [] (int i) const { return data[i]; }
};

struct os_kinTMIII_data
{
	int comp_thin, comp_thick, comp_halo;
	float fk, DeltavPhi;
	farray5 	vPhi1, vPhi2, vR, vZ,
		sigmaPhiPhi1, sigmaPhiPhi2, sigmaRR, sigmaZZ, sigmaRPhi, sigmaZPhi, sigmaRZ,
		HvPhi, HvR, HvZ,
		HsigmaPhiPhi, HsigmaRR, HsigmaZZ, HsigmaRPhi, HsigmaZPhi, HsigmaRZ;
};

namespace peyton { namespace system { class Config; }};
class otable;
class osink;
class opipeline;

struct skygenConfig;
struct skypixel;
struct skyConfigInterface
{
	virtual bool init(
		otable &t,
		const peyton::system::Config &cfg,	// model cfg file
		const skygenConfig &sc,
		const skypixel *pixels) = 0;
	virtual void initRNG(rng_t &rng) = 0;		// initialize the random number generator from CPU RNG
	virtual double integrateCounts() = 0;		// return the expected starcounts contributed by this model
	virtual void setDensityNorm(float norm) = 0;
	virtual size_t run(otable &in, osink *nextlink) = 0;
	virtual ~skyConfigInterface() {};
};

namespace multiplesAlgorithms
{
	enum algo {
		LF_M2_GT_M1	= 1,
		LF		= 2,
		EQUAL_MASS	= 3
	};
}

#endif
