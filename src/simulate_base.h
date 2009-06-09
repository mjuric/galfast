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

static const float ABSMAG_NOT_PRESENT = 99.999f;

struct os_FeH_data
{
	float A[2], sigma[3], offs[3];
	float Hmu, muInf, DeltaMu;
};

struct os_photometry_data
{
	int ncolors, bidx;
	float FeH0, dFeH;
	float Mr0, dMr;
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
	
	__device__ float& operator [] (int i) { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
	__device__ const float& operator [] (int i) const { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
};

static const int BahcallSoneira_model_THIN = 0, BahcallSoneira_model_THICK = 1, BahcallSoneira_model_HALO = 2;

struct iarray5
{
	short int data[5];
	
	__device__ short int& operator [] (int i) { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
	__device__ const short int& operator [] (int i) const { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
};

struct i8array5
{
	char data[5];
	
	__device__ char& operator [] (int i) { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
	__device__ const char& operator [] (int i) const { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
};

struct os_kinTMIII_data
{	
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

struct skyConfigInterface
{
	virtual bool init(
		const peyton::system::Config &cfg,
		const peyton::system::Config &pdf_cfg,
		const peyton::system::Config &foot_cfg,
		const peyton::system::Config &model_cfg,
		otable &t,
		opipeline &pipe) = 0;
	virtual size_t run(otable &in, osink *nextlink, rng_t &rng) = 0;
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
