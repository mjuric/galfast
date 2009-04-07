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

#include "config.h"

#include <iostream>
#include <fstream>

#include <astro/exceptions.h>
#include <astro/util.h>
#include <astro/math.h>
#include <astro/system/log.h>
#include <astro/io/format.h>
#include <astro/useall.h>

#include <vector>
#include "gpu.h"

stopwatch kernelRunSwatch;

// CUDA emulation for the CPU
// Used by CPU versions of CUDA kernels
__TLS char shmem[16384];
namespace gpuemu	// prevent collision with nvcc's symbols
{
	__TLS uint3 blockIdx;
	__TLS uint3 threadIdx;
	__TLS uint3 blockDim;	// Note: uint3 instead of dim3, because __TLS variables have to be PODs
	__TLS uint3 gridDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
}

__TLS int  active_compute_device;
GPUMM gpuMMU;


#if HAVE_CUDA || !ALIAS_GPU_RNG
gpu_rng_t::gpu_rng_t(rng_t &rng)
{
	float v = rng.uniform();
	seed = *(uint32_t*)&v;
}
#endif

#if HAVE_CUDA

#include <cuda_runtime.h>


///////////////////////////////////////////////////////////
// CUDA helpers

//
// Find the dimensions (bx,by) of a 2D grid of blocks that 
// has as close to nblocks blocks as possible
//
void find_best_factorization(unsigned int &bx, unsigned int &by, int nblocks)
{
	bx = -1;
	int best_r = 100000;
	for(int bytmp = 1; bytmp != 65536; bytmp++)
	{
		int r  = nblocks % bytmp;
		if(r < best_r && nblocks / bytmp < 65535)
		{
			by = bytmp;
			bx = nblocks / bytmp;
			best_r = r;
			
			if(r == 0) { break; }
			bx++;
		}
	}
	if(bx == -1) { std::cerr << "Unfactorizable?!\n"; exit(-1); }
}

//
// Given a total number of threads, their memory requirements, and the
// number of threadsPerBlock, compute the optimal allowable grid dimensions.
// Returns false if the requested number of threads are impossible to fit to
// shared memory.
//
bool calculate_grid_parameters(dim3 &gridDim, int threadsPerBlock, int neededthreads, int dynShmemPerThread, int staticShmemPerBlock)
{
	const int shmemPerMP =  16384;

	int dyn_shared_mem_required = dynShmemPerThread*threadsPerBlock;
	int shared_mem_required = staticShmemPerBlock + dyn_shared_mem_required;
	if(shared_mem_required > shmemPerMP) { return false; }

	// calculate the total number of threads
	int nthreads = neededthreads;
	int over = neededthreads % threadsPerBlock;
	if(over) { nthreads += threadsPerBlock - over; } // round up to multiple of threadsPerBlock

	// calculate the number of blocks
	int nblocks = nthreads / threadsPerBlock;
	if(nthreads % threadsPerBlock) { nblocks++; }

	// calculate block dimensions so that there are as close to nblocks blocks as possible
	find_best_factorization(gridDim.x, gridDim.y, nblocks);
	gridDim.z = 1;

//	MLOG(verb2) << "Grid parameters: tpb(" << threadsPerBlock << "), nthreads(" << neededthreads << "), shmemPerTh(" << dynShmemPerThread <<
//			"), shmemPerBlock(" << staticShmemPerBlock << ") --> grid(" << gridDim[0] << ", " << gridDim[1] << ", " << gridDim[2] << ")";

	return true;
}

static int cuda_initialized = 0;
static int cuda_enabled = 0;
bool gpuExecutionEnabled(const char *kernel)
{
	return cuda_enabled;
}

bool cuda_init()
{
	if(cuda_initialized) { return true; }

	// get requested device from environment
	const char *devStr = getenv("CUDA_DEVICE");
	int dev = devStr == NULL ? 0 : atoi(devStr);
	cudaError err;

	// disable GPU acceleration
	if(dev == -1)
	{
		cuda_initialized = 1;
		cuda_enabled = 0;
		return true;
	}
#if 0
	// get device properties
	cudaDeviceProp deviceProp;
	err = cudaGetDeviceProperties(&deviceProp, dev);
	if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); return false; }

	// use the device
	MLOG(verb1) << io::format("Using GPU Device %d: \"%s\"") << dev << deviceProp.name;
#endif
	err = cudaSetDevice(dev);
	if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); return false; }

	cuda_initialized = 1;
	cuda_enabled = 1;
	return true;
}


//////////////////////////////////////////////

size_t GPUMM::allocated() const
{
	size_t total = 0;
	FOREACH(gpuPtrs)
	{
		total += i->second.ptr.memsize();
	}
	return total;
}

// garbage collection: free GPU memory that has been synced to the host or released
void GPUMM::gc()
{
	std::vector<void*> toDelete;
	FOREACH(gpuPtrs)
	{
		gpu_ptr &g = i->second;
		if(g.lastop != SYNCED_TO_HOST || g.lastop != RELEASED_TO_HOST) { continue; }

		cudaFree(g.ptr.base);
		toDelete.push_back(i->first);
	}
	FOREACH(toDelete)
	{
		gpuPtrs.erase(*i);
	}
}

// copies the data to GPU device, allocating GPU memory if necessary
xptr GPUMM::syncToDevice(const xptr &hptr)
{
	// get GPU buffer
	gpu_ptr &g = gpuPtrs[hptr.base];
	if(!g.ptr.elementSize())
	{
		// check if any GC is needed
		size_t total = allocated();
		size_t newmem = hptr.memsize();
		if(total + newmem > gc_treshold) { gc(); }

		// allocate GPU memory
		g.ptr = hptr;
#if 0
		cudaError err = cudaMallocPitch(&g.ptr.base, &g.ptr.pitch(), hptr.ncols() * hptr.elementSize(), hptr.nrows());
#else
		cudaError err = cudaMalloc((void**)&g.ptr.base, hptr.memsize());
#endif
		if(err != cudaSuccess)
		{
			MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err) << "\n";
		}
	}

	// sync GPU with host, if needed
	if(g.lastop != SYNCED_TO_DEVICE)
	{
//		MLOG(verb1) << "Syncing to device (" << hptr.memsize() << " bytes)";
#if 0
		cudaError err = cudaMemcpy2D(g.ptr.base, g.ptr.pitch(), hptr.base, hptr.pitch(), hptr.width()*hptr.elementSize(), hptr.height(), cudaMemcpyHostToDevice);
#else
		cudaError err = cudaMemcpy(g.ptr.base, hptr.base, hptr.memsize(), cudaMemcpyHostToDevice);
#endif
		if(err != cudaSuccess)
		{
			MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err) << "\n";
		}
		g.lastop = SYNCED_TO_DEVICE;
	}

	return g.ptr;
}

void GPUMM::syncToHost(xptr &hptr)
{
	// get GPU buffer
	if(!gpuPtrs.count(hptr.base)) { return; }
	gpu_ptr &g = gpuPtrs[hptr.base];

	if(g.lastop != SYNCED_TO_HOST)
	{
#if 0
		cudaError err = cudaMemcpy2D(hptr.base, hptr.pitch(), g.ptr.base, g.ptr.pitch(), hptr.width()*hptr.elementSize(), hptr.height(), cudaMemcpyDeviceToHost);
#else
		cudaError err = cudaMemcpy(hptr.base, g.ptr.base, hptr.memsize(), cudaMemcpyDeviceToHost);
#endif
		if(err != cudaSuccess)
		{
			MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err) << "\n";
		}
		g.lastop = SYNCED_TO_HOST;
	}
}

#endif // HAVE_CUDA

#ifdef HAVE_CUDA
cudaError_t cudaConfigureCall(dim3 gridDim, dim3 blockDim, size_t sharedMem, cudaStream_t stream);
cudaError_t cudaSetupArgument(const void *arg, size_t size, size_t offset);
cudaError_t cudaLaunch(const char *entry);
struct cudaFuncAttributes
{
	size_t constSizeBytes;
	size_t localSizeBytes;
	int maxThreadsPerBlock;
	int numRegs;
	size_t sharedSizeBytes;
};
typedef int cuda_stream_t;
cudaError_t cudaFuncGetAttributes(cudaFuncAttributes *attr, const char *func)
{
	return cudaSuccess;
}
struct emptyArg {};

#define CUDA_RETURN_ON_FAIL(x) \
	{ cudaError_t ret_54e843 = (x); if(ret_54e843 != cudaSuccess) { return ret_54e843; } }

namespace nv
{
	struct kernel
	{
		dim3 gridDim, blockDim;
		size_t sharedMem;
		cuda_stream_t stream;
		size_t nthreads;
		const char *name;
		cudaError_t err;

		kernel(const char *kernel_name_, size_t nthreads_, size_t sharedMem_ = 0, cuda_stream_t stream_ = -1)
		: name(kernel_name_), nthreads(nthreads_), sharedMem(sharedMem_), stream(stream_)
		{
			int dev;
			cudaGetDevice(&dev);
			cudaDeviceProp dprop;
			cudaGetDeviceProperties(&dprop, dev);

			cudaFuncAttributes attr;
			err = cudaFuncGetAttributes(&attr, name);
			if(err != cudaSuccess) { return; }

			// compute the number of threads per block
			unsigned int nbreg = dprop.regsPerBlock / attr.numRegs; // Threads per block limit due to number of registers
			unsigned int nbmem = (dprop.sharedMemPerBlock - attr.sharedSizeBytes) / sharedMem; // Threads per block limit due to shared mem. size
			blockDim.x = attr.maxThreadsPerBlock;
			blockDim.x = std::min(blockDim.x, nbreg);
			blockDim.x = std::min(blockDim.x, nbmem);

			// compute grid dimensions
			calculate_grid_parameters(gridDim, blockDim.x, nthreads, sharedMem, attr.sharedSizeBytes);
		}
	};
};

template<typename T>
inline cudaError_t bindKernelParam(const T &p, size_t &offs)
{
	cudaError_t ret = cudaSetupArgument(&p, sizeof(T), offs);
	offs += sizeof(T);
	return ret;
}

template<>
inline cudaError_t  bindKernelParam(const emptyArg &p, size_t &offs)
{
	return cudaSuccess;
}

template<typename T1, typename T2, typename T3>
cudaError_t callKernel3(const nv::kernel &kc, const T1 &v1, const T2 &v2, const T3 &v3)
{
	// setup launch configuration
	if(kc.err) { return kc.err; }
	CUDA_RETURN_ON_FAIL( cudaConfigureCall(kc.gridDim, kc.blockDim, kc.sharedMem, kc.stream) );

	// push parameters to the stack
	size_t offs = 0;
	bindKernelParam(v1, offs);
	bindKernelParam(v2, offs);
	bindKernelParam(v3, offs);

	// launch the kernel
	return cudaLaunch(kc.name);
}

void test_cuda_caller()
{
	int var1;
	double var2;
	dim3 var3;

	int nthreads = 10000;
	callKernel3(nv::kernel("test_kernel", nthreads), var1, var2, var3);
}
#endif
