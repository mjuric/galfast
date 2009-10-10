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

#ifndef __CONFIG_H
#define __CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <cmath>
#include <limits>
#include <cassert>
#include "skygen.h"

#include <gsl/gsl_statistics_float.h>

using namespace cudacc;

__constant__ gpuRng::constant rng;	// GPU RNG
__constant__ lambert proj[2];		// Projection definitions for the two hemispheres

/********************************************************/

DEFINE_TEXTURE(ext_north, float, 3, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
DEFINE_TEXTURE(ext_south, float, 3, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D texture resampler, used mainly to debug extinction maps. Should be moved to a separate file.
//////////////////////////////////////////////////////////////////////////////////////////////////////////

__constant__ lambert texProj;

template<bool deproject>
__global__ void resample_extinction_kernel(gptr<float4, 1> out,
	float2 xrange, float2 yrange, float2 DMrange,
	int nx, int ny, int nDM)
{
	float dx  = ( xrange.y -  xrange.x) / (nx-1);
	float dy  = ( yrange.y -  yrange.x) / (ny-1);
	float dDM = (DMrange.y - DMrange.x) / (nDM-1);

	int nthreads = gridDim.x * blockDim.x;
	int at = blockDim.x * blockIdx.x + threadIdx.x;
	int end = nx * ny * nDM;

	while(at < end)
	{
		int k = at / (nx * ny);
		int at2 = at - k * nx * ny;
		int j = at2 / nx;
		int i = at2 - j * nx;

		float x  =  xrange.x + dx  * (i + 0.);
		float y  =  yrange.x + dy  * (j + 0.);
		float DM = DMrange.x + dDM * (k + 0.);

		float v;
		if(deproject)
		{
			float xim, yim;
			texProj.project(xim, yim, direction(x, y));
			if(!(-2.f < xim && xim < 2.f) || !(-2.f < yim && yim < 2.f))
			{
				// back away from the pole
				xim = 2.0;
				yim = 0.f;
			}
			v = TEX3D(ext_north, xim, yim, DM);
			x = deg(x);
			y = deg(y);
		}
		else
		{
			v = TEX3D(ext_north, x, y, DM);
		}

		out(at) = make_float4(x, y, DM, v);
		at += nthreads;
	}
}

cuxSmartPtr<float4> resample_extinction_texture(cuxTexture<float, 3> &tex, float2 crange[3], int npix[3], lambert *proj)
{
	cuxTextureBinder tb(ext_northManager, tex);

	cuxSmartPtr<float4> out(npix[0] * npix[1] * npix[2]);
	if(proj != NULL)
	{
		cuxUploadConst(texProj, *proj);
	}

	// TODO: There's really no good reason for this to be hardcoded...
	int nblocks = 30;
	int nthreads = 128;
	
	if(proj != NULL)
		resample_extinction_kernel<true><<<nblocks, nthreads>>>(out, crange[0], crange[1], crange[2], npix[0], npix[1], npix[2]);
	else
		resample_extinction_kernel<false><<<nblocks, nthreads>>>(out, crange[0], crange[1], crange[2], npix[0], npix[1], npix[2]);

	return out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// Get the extinction at projected coordinates (X,Y), for hemisphere specified
// by projIdx
__device__ float sampleExtinction(int projIdx, float X, float Y, float DM)
{
	float Am;

	if(projIdx == 0)
		Am = TEX3D(ext_north, X, Y, DM);
	else
		Am = TEX3D(ext_south, X, Y, DM);

	return Am;
}

/**
	Compute the quantities that are change only if a distance bin
	is changed (and will be cached otherwise).

	This includes the distance D, extinction, and 3D coordinates XYZ.
*/
template<typename T>
__device__ float3 skygenGPU<T>::compute_pos(float &D, float &Am, float M, const int im, const pencilBeam &pix) const
{
	float m = m0 + im*dm;
	float DM = m - M;

	D = powf(10.f, 0.2f*DM + 1.f);
	Am = sampleExtinction(pix.projIdx, pix.X, pix.Y, DM);

	return pix.xyz(D);
}

__device__ int ijToDiagIndex(int i, int j, const int x, const int y)
{
	int idx;

	int l = x < y ? x : y-1; // "less" (smaller dimension)
	int m = x < y ? y-1 : x; // "more" (bigger dimension)

	int d = i+j;
	idx = d*(d+1)/2 + i;
	if(l < d)
	{
		int dp = d - l;
		idx -= dp*(dp+1)/2;
		if(d > m)
		{
			int dm = d - m;
			idx -= dm*(dm+1)/2;
		}
	}
	return idx;
} 

/**
	Given a 'diagonal' linear index in (X,Y,M,m) space, decompose it into ilb, iM, im
	and k indices that can be used to compute physical coordinates.
*/
__device__ bool diagIndexToIJ(int &ilb, int &i, int &j, int &k, const int x, const int y)
{
	int kmax = x*y;
	if(k >= kmax)
	{
		ilb += k / kmax;
		k %= kmax;
	}

	int l = x < y ? x : y; // smaller dimension
	int m = x < y ? y : x; // bigger dimension
	int d;

	if(2*k < l*(l-1))
	{
		d = int(0.5*(sqrtf(1+8*k)-1));
		i = k - d*(d+1)/2;
	}
	else if(2*k < l*(2*m+1-l))
	{
		if(x >= y)
		{
			int ka = k - y*(y-1)/2;
			d = int(ka / y) + (y-1);
			i = (ka % y) + d - (y-1);
		}
		else
		{
			int ka = k - x*(x-1)/2;
			d = int(ka / x) + (x-1);
			i = ka % x;
		}
	}
	else
	{
		int ka = x*y - k - 1;
		int dd = int(0.5*(sqrtf(1+8*ka)-1));
		d = (x+y-2) - dd;
		i = (x-1) - (ka-dd*(dd+1)/2);
	}

	j = d - i;
	return true;
}

/**
	Diagonally advance the index in (XY,M,m) space. Since most models can be decomposed
	into LF(M)*den(X,Y,DM), this advancement usually leaves the thread in the same
	distance bin, allowing the (usually expensive) computation of den() to be cached.
*/
template<typename T>
__device__ bool skygenGPU<T>::advance(int &ilb, int &i, int &j, pencilBeam &dir, const int x, const int y) const
{
	i++; j--;

	if(j < 0) // slipped out through the top
	{
		if(i < y)
		{
			// This is the upper-left triangle
			j += i+1; i = 0;
		}
		else
		{
			// Mid range (hit only if x > y)
			i -= (y-1); j += y;
		}
	}
	else if(i >= x) // slipped out through the right edge
	{
		if(x > y)
		{
			// Bottom triangle (hit only if x > y)
			i -= y - j - 2; j = y-1;
		} else {
			if(j + 1 < y - x)
			{
				// Mid range (only if y > x)
				i -= x; j += x + 1;
			}
			else
			{
				// Bottom triangle (only if y > x)
				i = x - y + j + 2; j = y-1;
			}
		}
		if(i == x)
		{
			ilb++;
			if(ilb != npixels)
			{
				dir = pixels(ilb);
			}
		}
	}
	else
	{
		return false;	// no need to recompute distance
	}

	return true; // recompute the distance
}

/**
	Draw up to ndraw stars in magnitude bin (M,im) in pencil beam pix.

	Take care there's enough space in the output table for the generated
	stars. If there isn't, ndraw will be != 0 upon return.
*/
template<typename T>
__device__ void skygenGPU<T>::draw_stars(int &ndraw, const float &M, const int &im, const pencilBeam &pix) const
{
	if(!ndraw) { return; }

	int idx = atomicAdd(nstars.ptr, ndraw);

	for(; ndraw && idx < stopstars; idx++)
	{
		// Draw the distance and absolute magnitude
		float Mtmp, mtmp, DM;
		stars.M(idx) = Mtmp = M + dM*(rng.uniform() - 0.5f);
		mtmp = m0 + dm*(im + rng.uniform() - 0.5f);
		DM = mtmp - Mtmp;
		stars.DM(idx) = DM;
		float D = powf(10, 0.2f*DM + 1.f);

		// Draw the position within the pixel
		float x = pix.X, y = pix.Y;
		x += pix.dx*(rng.uniform() - 0.5f);
		y += pix.dx*(rng.uniform() - 0.5f);
		stars.projIdx(idx) = pix.projIdx;
		stars.projXY(idx, 0) = x;
		stars.projXY(idx, 1) = y;

		// Draw extinction
		stars.Am(idx) = sampleExtinction(pix.projIdx, x, y, DM);

		// Transform projected coordinates to (l,b), in degrees
		double l, b;
		proj[pix.projIdx].deproject(l, b, x, y);
		direction dir2(l, b);		// Note: we do this _before_ converting l,b to degrees
		l *= dbl_r2d;
		if(l < 0.) l += 360.;
		if(l > 360.) l -= 360.;
		b *= dbl_r2d;
		stars.lb(idx, 0) = l;
		stars.lb(idx, 1) = b;

		// Compute and store the 3D position (XYZ)
		float3 pos = dir2.xyz(D);
		stars.XYZ(idx, 0) = pos.x;
		stars.XYZ(idx, 1) = pos.y;
		stars.XYZ(idx, 2) = pos.z;

		// Store the component ID
		//stars.comp(idx) = model.component(pos.x, pos.y, pos.z, Mtmp, rng);
		stars.comp(idx) = model.component();

		ndraw--;
	}
}

/**
	The main skygen kernel (effectively, the entry point).

	Depending on the flag 'draw', it either draws the stars, or computes the number
	of stars in the given footprint. Launched from skygenHost<T>::run() and
	skygenHost<T>::integrateCounts().
*/
static const float POGSON = 0.4605170185988091f;
static const int block = 10;

template<typename T>
template<int draw>
__device__ void skygenGPU<T>::kernel() const
{
	int ilb, im, iM, k, ndraw = 0;
	float3 pos;
	pencilBeam pix;
	float Am;

	double count = 0.f, countCovered = 0.f;
	float maxCount1 = 0.;
	int bc = 0;					// block counter (decreses from block..0)
	float D;
	typename T::state ms;				// internal model class' state

	// Initialize (or load previously stored) execution state
	int tid = threadID();
	ilb = 0;
	k = block*(tid-nthreads);
	if(draw)
	{
		if(ks.continuing(tid))
		{
			ks.load(tid,   ilb, im, iM, k, bc, pos, D, pix, Am, ms, ndraw);
		}
		if(ilb >= npixels) { return; } // this thread has already finished
		rng.load(tid);

		// finish previous draw that didn't complete before the space ran out
		if(ndraw)
		{
			float M = M1 - iM*dM;
			draw_stars(ndraw, M, im, pix);
		}
	}

	//
	// Crawl through (X,Y,M,m) space, sample the densities and
	// either sum them up or draw the stars.
	//
	// We crawl through this space by incrementing a linear 'diagonal index' k
	// (TODO: explain this better).
	//
	// To evenly distribute work, while still maintaining some locality (i.e.,
	// not moving between distance bins often), we crawl in blocks of size 'block'
	// and then jump block*nthreads ahead.
	//
	while(ndraw == 0)
	{
		// advance the index in (X,Y,M,m) space (indexed by (ilb,iM,im), or linear index k)
		bool moved;	// whether we moved to a different distance bin
		if(bc == 0)	// jump over block*nthreads, or advance by 1?
		{
			// jump block*nthreads
		 	bc = block;
			k += block*nthreads;
			diagIndexToIJ(ilb, im, iM, k, nm, nM);
			if(ilb >= npixels) { break; }

			pix = pixels(ilb);
			moved = true;
		}
		else
		{
			// advance by 1
			moved = advance(ilb, im, iM, pix, nm, nM);
		}
		bc--;

		// compute the absolute magnitude and position for this pixel
		float M = M1 - iM*dM;
		if(moved)
		{
			if(ilb >= npixels) { break; }

			// We moved to a new distance bin. Recompute and reset.
			pos = compute_pos(D, Am,     M, im, pix);
			model.setpos(ms, pos.x, pos.y, pos.z);
		}

		// Test if this location has been extincted away.
		float m = m0 + dm*im + Am;
		if(m > m1) { continue; }

		// compute the density in this pixel
		float rho = model.rho(ms, M);
		rho *= norm;
		rho *= D*D*D;			  // multiply by volume (part one)
		rho *= pix.dA * POGSON * dm * dM; // multiply by volume (part two). TODO: Optimize this by precomputing

		if(draw) // draw stars
		{
			if(ndraw == 0)
			{
				ndraw = rng.poisson(rho);
			}

			draw_stars(ndraw, M, im, pix);
		}
		else
		{
			// sum up the stars in the volume
			count += rho;
			countCovered += rho * pix.coveredFraction;
			if(maxCount1 < rho) { maxCount1 = rho; }
#if 1
			// add this sample to the correct histogram bin
			// NOTE: This is basically debugging info, and can be thrown
			// out at a later date.
			int rhoBin = round((log10f(rho) - lrho0) / dlrho);
			if(rhoBin < 0) { rhoBin = 0; }
			if(rhoBin >= nhistbins) { rhoBin = (nhistbins-1); }
			rhoHistograms(nthreads*rhoBin + tid)++;
#endif
		}
	};

	tid = threadID();
	if(!draw)
	{
		counts(tid) = count;
		countsCovered(tid) = countCovered;
		maxCount(tid) = maxCount1;
	}
	else
	{
		// store execution state
		ks.store(tid,   ilb, im, iM, k, bc, pos, D, pix, Am, ms, ndraw);
		rng.store(tid);
	}
}

// default kernels (do nothing)
// NOTE: CUDA compatibility -- in principle, we could only _declare_, but 
// not define these functions to ensure they can never be instantiated without
// specialization. However, nvcc then fails to compile this code.
template<typename T> __global__ void integrate_counts_kernel() {  }
template<typename T> __global__ void draw_sources_kernel() {  }

// Launch the apropriate kernel. The launched kernels are specialized
// by the particular models.
template<typename T>
void skygenHost<T>::compute(bool draw)
{
	cuxErrCheck( cudaThreadSynchronize() );

	kernelRunSwatch.start();
	if(draw)
	{
		draw_sources_kernel<T><<<gridDim, blockDim, shb>>>();
	}
	else
	{
		integrate_counts_kernel<T><<<gridDim, blockDim, shb>>>();
	}
	cuxErrCheck( cudaThreadSynchronize() );
	kernelRunSwatch.stop();
}

#include "model_J08.h"
#include "model_expDisk.h"
#include "model_densityCube.h"
#include "model_powerLawEllipsoid.h"
