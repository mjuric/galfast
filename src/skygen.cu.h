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

//
// GPU implementation
//

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <cmath>
#include <limits>
#include <cassert>
#include "skygen.h"

#include <gsl/gsl_statistics_float.h>

__device__ __constant__ gpuRng::constant rng;
__device__ __constant__ skyConfigGPU<expModel> expModelSky;
__device__ __constant__ skyConfigGPU<expDisk> expDiskSky;
__device__ __constant__ lambert proj[2];

#if 0
lfTextureManager texLFMgr(texLF);
lfParams lfTextureManager::set(float *cpu_lf, int lfLen, float M0, float M1, float dM)
{
	free();

	// Upload luminosity function as texture
	par = make_lfParams(M0, M1, dM);
	cuxErrCheck( cudaMallocArray( &lfArray, &texref.channelDesc, lfLen, 1));
	cuxErrCheck( cudaMemcpyToArray( lfArray, 0, 0, cpu_lf, lfLen*sizeof(float), cudaMemcpyHostToDevice));
	bind();

	return par;
}

void lfTextureManager::bind()
{
	//fprintf(stderr, "Binding luminosity function texture.\n");
	assert(lfArray);
	cuxErrCheck( cudaBindTextureToArray(&texref, lfArray, &texref.channelDesc) );
}

void lfTextureManager::free()
{
	if(lfArray)
	{
		cudaUnbindTexture(&texref);
		cudaFreeArray(lfArray);
		lfArray = NULL;
	}
}
#endif

void expModel::prerun(host_state_t &hstate, bool draw)
{
	// bind the luminosity function texture to texture reference
	hstate.lf.bind_texture(expModelLF);
}

void expDisk::prerun(host_state_t &hstate, bool draw)
{
	// bind the luminosity function texture to texture reference
	hstate.lf.bind_texture(expDiskLF);
}

void expModel::postrun(host_state_t &hstate, bool draw)
{
	if(draw)
	{
		hstate.lf.unbind_texture(expModelLF);
	}
}

void expDisk::postrun(host_state_t &hstate, bool draw)
{
	if(draw)
	{
		hstate.lf.unbind_texture(expDiskLF);
	}
}

template<typename T>
void skyConfig<T>::download(bool draw)
{
	if(draw)
	{
		this->nstars.download(&stars_generated, 1);

		// this is for debugging purposes mostly
		int *ilb = new int[this->nthreads];
		int *im = new int[this->nthreads];
		int *iM = new int[this->nthreads];
		delete [] cpu_state;
		cpu_state = new int3[this->nthreads];
		this->ks.ilb.download(ilb, this->nthreads);
		this->ks.im.download(im,   this->nthreads);
		this->ks.iM.download(iM,   this->nthreads);
		for(int i=0; i != this->nthreads; i++)
		{
			cpu_state[i] = make_int3(ilb[i], im[i], iM[i]);
		}
		delete [] ilb; delete [] im, delete [] iM;
	}
	else
	{
		float *cpu_counts = new float[this->nthreads];
		float *cpu_countsCovered = new float[this->nthreads];
		this->counts.download(cpu_counts, this->nthreads);
		this->countsCovered.download(cpu_countsCovered, this->nthreads);

		// sum up the total expected number of stars
		nstarsExpectedToGenerate = 0;
		nstarsExpected = 0;
		for(int i=0; i != this->nthreads; i++)
		{
			nstarsExpectedToGenerate += cpu_counts[i];
			nstarsExpected += cpu_countsCovered[i];
		}
		delete [] cpu_counts;
		delete [] cpu_countsCovered;

		int *cpu_rhoHistograms = new int[this->nthreads*this->nhistbins];
		this->rhoHistograms.download(cpu_rhoHistograms, this->nthreads*this->nhistbins);

		// sum up the total
		cpu_hist = new int[this->nhistbins];
		memset(cpu_hist, 0, sizeof(float)*this->nhistbins);
		for(int i=0; i != this->nthreads; i++)
		{
			for(int j=0; j != this->nhistbins; j++)
			{
				cpu_hist[j] += cpu_rhoHistograms[this->nthreads*j + i];
			}
		}
		delete [] cpu_rhoHistograms;

		// Download list of maximum densities encountered by each thread
		delete [] cpu_maxCount;
		cpu_maxCount = new float[this->nthreads];
		this->maxCount.download(cpu_maxCount, this->nthreads);
	}
}

#if 0
static const float inf  =  std::numeric_limits<float>::infinity();
static const float snan =  std::numeric_limits<float>::signaling_NaN();
static const int smallest_float_int = 0x00000001; // hex representation of the smallest representable positive 32-bit float (IEEE 754)
static const float feps = *(float *)&smallest_float_int;
#endif

/********************************************************/

static const float POGSON = 0.4605170185988091f;
static const int block = 10;

#if _EMU_DEBUG
float mmax = 22;
float mmin = 15;
#endif

template<typename T>
__device__ float3 skyConfigGPU<T>::compute_pos(float &D, float M, const int im, const direction &dir) const
{
	float m = m0 + im*dm;
#if 0 && _EMU_DEBUG
	if(m > mmax) {
		mmax = m;
		printf("mmax = %f, im=%d, M=%f\n", mmax, im, M);
	}
	if(m < mmin) {
		mmin = m;
		printf("mmin = %f, im=%d, M=%f\n", mmin, im, M);
	}
#endif
	D = powf(10, 0.2f*(m - M) + 1.);
	return position(dir, D);
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

template<typename T>
__device__ bool skyConfigGPU<T>::advance(int &ilb, int &i, int &j, skypixel &dir, const int x, const int y) const
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

#if _EMU_DEBUG
long long voxelsVisited = 0;
#endif

template<typename T>
template<int draw>
__device__ void skyConfigGPU<T>::kernel() const
{
	int ilb, im, iM;
	float3 pos;
	skypixel pix;
	int k;

	double count = 0.f, countCovered = 0.f;
	float maxCount1 = 0.;
	int bc = 0;
	float D;
	typename T::state ms;

#if _EMU_DEBUG
	long long voxelsVisited1 = 0;
#endif

	// Initialize (or load previously stored) execution state
	int tid = threadID();
	ilb = 0;
	k = block*(tid-nthreads);
	if(draw)
	{
		if(ks.continuing(tid))
		{
			ks.load(tid,   ilb, im, iM, k, bc, pos, D, pix, ms);
		}
		if(ilb >= npixels) { return; } // this thread has already finished
		rng.load(tid);
	}

	while(true)
	{
		// advance the index
		bool moved;
		if(bc == 0)
		{
		 	bc = block;
			k += block*nthreads;
			diagIndexToIJ(ilb, im, iM, k, nm, nM);
			if(ilb >= npixels) { break; }

			pix = pixels(ilb);
			moved = true;
		}
		else
		{
			moved = advance(ilb, im, iM, pix, nm, nM);
		}

#if _EMU_DEBUG
		{
			voxelsVisited1++;

			// verify diagonal indices were computed correctly
			int idx = ijToDiagIndex(im, iM, nm, nM);
			int kk = k + (block-bc);
			if(idx != kk)
			{
				printf("ERROR: idx != k: %d != %d (i=%d, j=%d, ilb=%d)", idx, kk, im, iM, ilb);
				abort();
			}

			// verify this pixel was supposed to be processed by this thread
			uint64_t gidx = ilb*nm*nM + idx;	// global index
			int lidx = (gidx % (nthreads*block)) - tid*block;	// local index
			if(0 > lidx || lidx >= block)
			{
				printf("ERROR: This thread (tid=%d) should not be processing pixel %ld!\n", tid, gidx);
				abort();
			}

			// check there was no change in distance unless 'moved' is true
			if(!moved)
			{
				float D2;
				float M = M1 - iM*dM;
				float3 pos2 = compute_pos(D2, M, im, pix);
				if(fabsf(D/D2-1) > 1e-5)
				{
					printf("ERROR: Unexpected distance while not moving! Old D = %f, new D = %f\n", D, D2);
					abort();
				}
			}
		}
#endif
		bc--;

		// compute absolute magnitude and position for this pixel
		float M = M1 - iM*dM;
		if(moved)
		{
			if(ilb >= npixels) { break; }

			// we moved to a new distance bin. Recompute and reset.
			pos = compute_pos(D, M, im, pix);
			model.setpos(ms, pos.x, pos.y, pos.z);
		}

		// TODO: test for extinction here

		// compute the density in this pixel
		float rho = model.rho(ms, M);
#if __DEVICE_EMULATION__
//		printf("dN=%f dA=%f dm=%f dM=%f D=%f\n", rho, dA, dm, dM, D);
#endif
		rho *= D*D*D; // multiply by volume (part one)
		rho *= pix.dA * POGSON * dm * dM; // multiply by volume (part two)
#if __DEVICE_EMULATION__
//		printf("dN=%g dA=%f dm=%f dM=%f\n", rho, dA, dm, dM); abort();
#endif
		if(draw) // draw stars
		{
			int k = rng.poisson(rho);

			if(k)
			{
				int idx = atomicAdd(nstars.ptr, k);

				do
				{
					float Mtmp, mtmp, DMtmp;
					stars.M(idx) = Mtmp = M + dM*(rng.uniform() - 0.5f);
//					stars.m[idx] = mtmp = m0 + dm*(im + rng.uniform() - 0.5f);
					mtmp = m0 + dm*(im + rng.uniform() - 0.5f);
					DMtmp = mtmp - Mtmp;
					stars.DM(idx) = DMtmp;
					float D = powf(10, 0.2f*DMtmp + 1.f);

					float x, y;
					stars.projIdx(idx) = pix.projIdx;
					proj[pix.projIdx].convert(pix, x, y);
					x += pix.dx*(rng.uniform() - 0.5f);
					y += pix.dx*(rng.uniform() - 0.5f);

					double l, b;
					proj[pix.projIdx].inverse(x, y, l, b);
					direction dir2(l, b);		// Note: we do this _before_ converting l,b to degrees

					l *= dbl_r2d;
					if(l < 0.) l += 360.;
					if(l > 360.) l -= 360.;
					b *= dbl_r2d;
					stars.lb(idx, 0) = l;
					stars.lb(idx, 1) = b;

					float3 pos = position(dir2, D);
					stars.XYZ(idx, 0) = pos.x;
					stars.XYZ(idx, 1) = pos.y;
					stars.XYZ(idx, 2) = pos.z;
					stars.comp(idx) = model.component(pos.x, pos.y, pos.z, Mtmp, rng);

// 					printf("%3d %5d %13.8f %13.8f %6.3f %6.3f %10.2f %10.2f %10.2f %3d\n",
// 						sizeof(s), idx, deg(s.l), deg(s.b), s.M, s.m, s.pos.x, s.pos.y, s.pos.z, s.comp);

					idx++;
					k--;
				} while(k);

				if(idx >= stopstars)
				{
					// memory buffer filled. Exit.
					break;
				}
			}
		}
		else
		{
			count += rho;
			countCovered += rho * pix.coveredFraction;
			if(maxCount1 < rho) { maxCount1 = rho; }
#if 1
			// add this sample to the correct histogram bin
			int rhoBin = round((log10f(rho) - lrho0) / dlrho);
			if(rhoBin < 0) { rhoBin = 0; }
			if(rhoBin >= nhistbins) { rhoBin = (nhistbins-1); }
			rhoHistograms(nthreads*rhoBin + tid)++;
			//if(rho > 1e-2f) { printf("%g (D=%f) -> bin %d\n", rho, D, rhoBin); }
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
		ks.store(tid,   ilb, im, iM, k, bc, pos, D, pix, ms);
		rng.store(tid);
	}

#if _EMU_DEBUG
	int locked;
	do {
		locked = atomicCAS(lock.ptr, 0, 1);
	} while(locked != 0);
	voxelsVisited += voxelsVisited1;
	lock[0] = 0;

	#if 0
		if(draw) {
			printf("%d stars after thread %d/%d done (in %d voxels).\n", locked, tid, nthreads, voxelsVisited1);
		} else {
			printf("Thread %d/%d done (sampled %d voxels).\n", tid, nthreads, voxelsVisited1);
		}
	#endif
#endif
}

// default kernels (do nothing)
template<typename T> __global__ void compute_sky() { }
template<typename T> __global__ void draw_sky() { }

// expModel overriden kernels
template<> __global__ void compute_sky<expModel>() { expModelSky.kernel<0>(); }
template<> __global__ void draw_sky<expModel>() { expModelSky.kernel<1>(); }

// expModel overriden kernels
template<> __global__ void compute_sky<expDisk>() { expDiskSky.kernel<0>(); }
template<> __global__ void draw_sky<expDisk>() { expDiskSky.kernel<1>(); }

template<typename T>
void skyConfig<T>::compute(bool draw)
{
	cuxErrCheck( cudaThreadSynchronize() );

	swatch.reset();
	swatch.start();
	if(draw)
	{
		draw_sky<T><<<gridDim, blockDim, shb>>>();
	}
	else
	{
		compute_sky<T><<<gridDim, blockDim, shb>>>();
	}
	cuxErrCheck( cudaThreadSynchronize() );
	swatch.stop();
}

// Explicit instantiations
template class skyConfig<expModel>;
template class skyConfig<expDisk>;
