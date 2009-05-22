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
#include "model.h"
#include "analysis.h"

#if HAVE_CUDA
#include <cuda_runtime.h>
#endif // HAVE_CUDA

xptrng::ptr_desc *xptrng::ptr_desc::null = NULL;

xptrng::ptr_desc::ptr_desc(size_t es, size_t width, size_t height, size_t pitch)
	: m_elementSize(es), m_dim(width, height, 1), m_pitch(pitch), m_data(NULL)
{
	refcnt = 1;
	masterDevice = -1;
	deviceDataPointers[-1] = m_data = new char[memsize()];
}
xptrng::ptr_desc::~ptr_desc()
{
	delete [] ((char*)m_data);
	FOREACH(deviceDataPointers)
	{
		if(i->first < 0) { continue; }

		cuxErrCheck( cudaFree(i->second) );
	}
	FOREACH(cudaArrayPointers)
	{
		cuxErrCheck( cudaFreeArray(i->second) );
	}
}
void *xptrng::ptr_desc::syncToDevice(int dev)
{
	if(dev == -2)
	{
		dev = gpuGetActiveDevice();
	}

	if(masterDevice != dev)
	{
		// check if this device-to-device copy. If so, do the copy via host
		if(masterDevice != -1 && dev != -1)
		{
			syncToDevice(-1); // this will change masterDevice to -1
		}

		// allocate/copy to device
		cudaError err;

		// determine destination and copy direction
		cudaMemcpyKind dir = cudaMemcpyDeviceToHost;
		if(dev != -1)
		{
			dir = cudaMemcpyHostToDevice;

			// allocate device space (if unallocated)
			if(!deviceDataPointers.count(dev))
			{
				cuxErrCheck( cudaMalloc(&deviceDataPointers[dev], memsize()) );
			}
		}
		void *dest = deviceDataPointers[dev];
		void *src = deviceDataPointers[masterDevice];

		// do the copy
		cuxErrCheck( cudaMemcpy(dest, src, memsize(), dir) );

		// record new master device
		masterDevice = dev;
	}
	return deviceDataPointers[masterDevice];
}
cudaArray *xptrng::ptr_desc::getCUDAArray(cudaChannelFormatDesc &channelDesc, int dev, bool forceUpload)
{
	syncToDevice(-1);	// ensure the data is on the host

	if(dev == -2)
	{
		dev = gpuGetActiveDevice();
		assert(dev >= 0);
	}

	cudaArray* cu_array;
	cudaError err;

	if(cudaArrayPointers.count(dev))
	{
		cu_array = cudaArrayPointers[dev];
	}
	else
	{
		// autocreate
		cuxErrCheck( cudaMallocArray(&cu_array, &channelDesc, m_dim.x, m_dim.y) );

		cudaArrayPointers[dev] = cu_array;
		forceUpload = true;
	}

	if(forceUpload)
	{
		cuxErrCheck( cudaMemcpy2DToArray(cu_array, 0, 0, m_data, m_pitch, m_dim.x*m_elementSize, m_dim.y, cudaMemcpyHostToDevice) );
	}

	return cu_array;
}


stopwatch kernelRunSwatch;
gpu_rng_t::persistent_rng gpu_rng_t::gpuRNG;

gpu_prng::mwc &gpu_rng_t::persistent_rng::get(rng_t &seeder)
{
	if(!gpuRNG)
	{
		std::string file = datadir() + "/safeprimes32.txt";

		gpuRNG = new gpu_prng::mwc;
		gpuRNG->srand((uint32_t)(seeder.uniform()*(1<<24)), 1<<16, file.c_str());
	}
	return *gpuRNG;
}

// CUDA emulation for the CPU
// Used by CPU versions of CUDA kernels
__TLS char impl_shmem[16384];
namespace gpuemu	// prevent collision with nvcc's symbols
{
	__TLS uint3 blockIdx;
	__TLS uint3 threadIdx;
	__TLS uint3 blockDim;	// Note: uint3 instead of dim3, because __TLS variables have to be PODs
	__TLS uint3 gridDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
}

__TLS int  active_compute_device;

#if 0
#if HAVE_CUDA || !ALIAS_GPU_RNG
struct rng_mwc
{
	static uint32_t nstreams;	// number of allocated/initialized streams
	static uint32_t *cpu_streams;	// Pointer to streams state on the CPU (used in CPU mode)
#if HAVE_CUDA
	static uint32_t *gpu_streams;	// Pointer to streams state on the device (used in GPU mode)
	static bool onDevice;		// Whether the master copy is on the device (GPU)
#endif

	static void init(rng_t &rng)
	{
		static bool initialized = false;
		if(initialized) { return; }
		initialized = true;

		nstreams = 1<<16;
		cpu_streams = new uint32_t[3*nstreams];

		// initialize CPU streams
		text_input_or_die(in, datadir() + "/safeprimes32.txt");
		std::vector<int> primes;
		load(in, primes, 0);
		if(primes.size() < nstreams)
		{
			THROW(EAny, "Insufficient number of safe primes in " + datadir() + "/safeprimes32.txt");
		}
		for(int i = 0; i != nstreams; i++)
		{
			cpu_streams[i] = primes[i];	// multiplier
			cpu_streams[  nstreams + i] = (int)(rng.uniform() * cpu_streams[i]);	// initial carry (nas to be < multiplier)
			float v = rng.uniform();
			cpu_streams[2*nstreams + i] = *(uint32_t*)&v;	// initial x
		}

		DLOG(verb1) << "Initialized " << nstreams << " multiply-with-carry RNG streams";
	}

	static void checkInit()
	{
		if(cpu_streams) return;

		MLOG(verb1) << "ERROR: Must call rng_mwc::init before using GPU random number generator";
		abort();
	}
	static uint32_t statebytes()
	{
		return sizeof(uint32_t)*3*nstreams;
	}

	static uint32_t *gpuStreams()
	{
		checkInit();
#if HAVE_CUDA
		// sync with device, if on CPU
		if(!onDevice)
		{
			cudaError err;

			// allocate device space (if unallocated)
			if(gpu_streams == NULL)
			{
				err = cudaMalloc((void**)&gpu_streams, statebytes());
				if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); abort(); }
			}

			// copy to device
			err = cudaMemcpy(gpu_streams, cpu_streams, statebytes(), cudaMemcpyHostToDevice);
			if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); abort(); }

			onDevice = true;
		}

		return gpu_streams;
#else
		MLOG(verb1) << "ERROR: We should have never gotten here with CUDA support disabled!";
		abort();
#endif
	}

	static uint32_t *cpuStreams()
	{
		checkInit();
#if HAVE_CUDA
		if(onDevice)
		{
			// wait for currently executing kernels to finish
			cudaError err = cudaThreadSynchronize();
			if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); abort(); }

			// copy to host
			err = cudaMemcpy(cpu_streams, gpu_streams, statebytes(), cudaMemcpyDeviceToHost);
			if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); abort(); }

			onDevice = false;
		}
#endif
		return cpu_streams;
	}
};
uint32_t rng_mwc::nstreams = 0;
uint32_t *rng_mwc::cpu_streams = NULL;
#if HAVE_CUDA
uint32_t *rng_mwc::gpu_streams = NULL;
bool rng_mwc::onDevice = false;
#endif

gpu_rng_t::gpu_rng_t(rng_t &rng)
{
	rng_mwc::init(rng);

	nstreams = rng_mwc::nstreams;
	streams = gpuGetActiveDevice() < 0 ? rng_mwc::cpuStreams() : rng_mwc::gpuStreams();
}
#endif
#endif

#if 0
//////////////////////////////////////////////
#define ENABLE_PAGELOCKED 0
void xptr::alloc(size_t eSize, size_t ncol, size_t nrow, size_t ptch)
{
	if(eSize == (size_t)-1) { eSize = elementSize(); } else { m_elementSize = eSize; }
	if(ncol == (size_t)-1) { ncol = ncols(); } else { dim[0] = ncol; }
	if(nrow == (size_t)-1) { nrow = nrows(); } else { dim[1] = nrow; }
	if(ptch == (size_t)-1) { ptch = pitch(); } else { m_pitch[0] = ptch; }

	free();
#if ENABLE_PAGELOCKED
	std::cerr << "*************** PAGELOCKED ALLOC OF " << memsize() << " bytes.\n";
	cudaMallocHost((void **)&base, memsize());
#else
	base = new char[memsize()];
#endif
}
	
void xptr::free()
{
#if ENABLE_PAGELOCKED
	if(base != NULL)
	{
		std::cerr << "*************** PAGELOCKED FREE\n";
		cudaFreeHost(base);
	}
#else
	delete [] base;
#endif
	base = NULL;
}
//////////////////////////////////////////////
#endif

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

	MLOG(verb2) << "Grid parameters: tpb(" << threadsPerBlock << "), nthreads(" << neededthreads << "), shmemPerTh(" << dynShmemPerThread <<
			"), shmemPerBlock(" << staticShmemPerBlock <<
			") --> grid(" << gridDim.x << ", " << gridDim.y << ", " << gridDim.z <<
			") block(" << threadsPerBlock << ", 1, 1)" <<
			" " << shared_mem_required << " shmem/block (" << (float)shared_mem_required / shmemPerMP << ")" << 
			" == total of " << nthreads << " threads.";

	return true;
}

const char *cpuinfo()
{
	static char buf[1000];
	FILE *f = popen("cat /proc/cpuinfo | grep 'model name' | head -n 1 | awk -F': ' '{ print $2}'", "r");
	fgets(buf, 1000, f);
	pclose(f);

	int len = strlen(buf);
	if(len && buf[len-1] == '\n') buf[len-1] = 0;
	return buf;
}


void abort_on_cuda_error(cudaError err)
{
	if(err == cudaSuccess) { return; }

	MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err);
	abort();
}

void cuxErrCheck(cudaError err)
{
	if(err != cudaSuccess)
	{
		abort_on_cuda_error(err);
//		throw cuxException(err);
	}
}

#if HAVE_CUDA
static int cuda_initialized = 0;
static int cuda_enabled = 0;
bool gpuExecutionEnabled(const char *kernel)
{
/*	if(strcmp(kernel, "os_kinTMIII_kernel") == 0)
	{
		return false;
	}*/
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

		MLOG(verb1) << "Using CPU: \"" << cpuinfo() << "\"";
		return true;
	}
#if !CUDA_DEVEMU
	// get device properties
	cudaDeviceProp deviceProp;
	err = cudaGetDeviceProperties(&deviceProp, dev);
	if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); return false; }

	// use the device
	MLOG(verb1) << io::format("Using CUDA Device %d: \"%s\"") << dev << deviceProp.name;
#else
	MLOG(verb1) << "Using CUDA Device Emulation";
#endif
	err = cudaSetDevice(dev);
	if(err != cudaSuccess) { MLOG(verb1) << "CUDA Error: " << cudaGetErrorString(err); return false; }

	cuda_initialized = 1;
	cuda_enabled = 1;
	return true;
}

#endif // HAVE_CUDA

#if HAVE_CUDA
#if CUDART_VERSION < 2020
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
#endif

struct emptyArg {};

#define CUDA_RETURN_ON_FAIL(x) \
	{ cudaError_t ret_54e843 = (x); if(ret_54e843 != cudaSuccess) { return ret_54e843; } }

namespace nv
{
	struct kernel
	{
		dim3 gridDim, blockDim;
		size_t sharedMem;
		cudaStream_t stream;
		size_t nthreads;
		const char *name;
		cudaError_t err;

		kernel(const char *kernel_name_, size_t nthreads_, size_t sharedMem_ = 0, cudaStream_t stream_ = -1)
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
