#ifndef __gpu_h
#define __gpu_h

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

// Possible states
// #define HAVE_CUDA 1 -- CUDA support is on
//	#define BUILD_FOR_CPU 1	-- build CPU versions of kernels
//
// Note: __CUDACC__ will never be defined when BUILD_FOR_CPU is defined
//

// Thread local storage -- use for shared memory emulation in CPU mode
//#define __TLS __thread
#define __TLS

#if __CUDACC__
extern __shared__ char shmem[];
#else
extern __TLS char shmem[16384];
#endif

#if HAVE_CUDA
	#include <cuda_runtime.h>

	bool cuda_init();
#else
	// Emulate CUDA types and keywords

	#define __device__
	#define __host__
	#define __constant__

	struct float4 { float x, y, z, w; };
	struct double2 { double x, y; };
	struct uint3 { unsigned int x, y, z; };
	struct uint4 { unsigned int x, y, z, w; };

	struct dim3
	{
		unsigned int x, y, z;
		#if defined(__cplusplus)
		dim3(unsigned int x = 1, unsigned int y = 1, unsigned int z = 1) : x(x), y(y), z(z) {}
		dim3(uint3 v) : x(v.x), y(v.y), z(v.z) {}
		operator uint3(void) { uint3 t; t.x = x; t.y = y; t.z = z; return t; }
		#endif /* __cplusplus */
	};

	/*
	* Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
	*
	* NOTICE TO USER:
	*
	* This source code is subject to NVIDIA ownership rights under U.S. and
	* international Copyright laws.  Users and possessors of this source code
	* are hereby granted a nonexclusive, royalty-free license to use this code
	* in individual and commercial software.
	*
	* NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
	* CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
	* IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
	* REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
	* MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
	* IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
	* OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
	* OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
	* OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
	* OR PERFORMANCE OF THIS SOURCE CODE.
	*
	* U.S. Government End Users.   This source code is a "commercial item" as
	* that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
	* "commercial computer  software"  and "commercial computer software
	* documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
	* and is provided to the U.S. Government only as a commercial end item.
	* Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
	* 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
	* source code with only those rights set forth herein.
	*
	* Any use of this source code in individual and commercial software must
	* include, in the user documentation and internal comments to the code,
	* the above Disclaimer and U.S. Government End Users Notice.
	*/

	enum cudaError
	{
	cudaSuccess                           =      0,   ///< No errors
	cudaErrorMissingConfiguration         =      1,   ///< Missing configuration error
	cudaErrorMemoryAllocation             =      2,   ///< Memory allocation error
	cudaErrorInitializationError          =      3,   ///< Initialization error
	cudaErrorLaunchFailure                =      4,   ///< Launch failure
	cudaErrorPriorLaunchFailure           =      5,   ///< Prior launch failure
	cudaErrorLaunchTimeout                =      6,   ///< Launch timeout error
	cudaErrorLaunchOutOfResources         =      7,   ///< Launch out of resources error
	cudaErrorInvalidDeviceFunction        =      8,   ///< Invalid device function
	cudaErrorInvalidConfiguration         =      9,   ///< Invalid configuration
	cudaErrorInvalidDevice                =     10,   ///< Invalid device
	cudaErrorInvalidValue                 =     11,   ///< Invalid value
	cudaErrorInvalidPitchValue            =     12,   ///< Invalid pitch value
	cudaErrorInvalidSymbol                =     13,   ///< Invalid symbol
	cudaErrorMapBufferObjectFailed        =     14,   ///< Map buffer object failed
	cudaErrorUnmapBufferObjectFailed      =     15,   ///< Unmap buffer object failed
	cudaErrorInvalidHostPointer           =     16,   ///< Invalid host pointer
	cudaErrorInvalidDevicePointer         =     17,   ///< Invalid device pointer
	cudaErrorInvalidTexture               =     18,   ///< Invalid texture
	cudaErrorInvalidTextureBinding        =     19,   ///< Invalid texture binding
	cudaErrorInvalidChannelDescriptor     =     20,   ///< Invalid channel descriptor
	cudaErrorInvalidMemcpyDirection       =     21,   ///< Invalid memcpy direction
	cudaErrorAddressOfConstant            =     22,   ///< Address of constant error
	cudaErrorTextureFetchFailed           =     23,   ///< Texture fetch failed
	cudaErrorTextureNotBound              =     24,   ///< Texture not bound error
	cudaErrorSynchronizationError         =     25,   ///< Synchronization error
	cudaErrorInvalidFilterSetting         =     26,   ///< Invalid filter setting
	cudaErrorInvalidNormSetting           =     27,   ///< Invalid norm setting
	cudaErrorMixedDeviceExecution         =     28,   ///< Mixed device execution
	cudaErrorCudartUnloading              =     29,   ///< CUDA runtime unloading
	cudaErrorUnknown                      =     30,   ///< Unknown error condition
	cudaErrorNotYetImplemented            =     31,   ///< Function not yet implemented
	cudaErrorMemoryValueTooLarge          =     32,   ///< Memory value too large
	cudaErrorInvalidResourceHandle        =     33,   ///< Invalid resource handle
	cudaErrorNotReady                     =     34,   ///< Not ready error
	cudaErrorInsufficientDriver           =     35,   ///< CUDA runtime is newer than driver
	cudaErrorSetOnActiveProcess           =     36,   ///< Set on active process error
	cudaErrorNoDevice                     =     38,   ///< No available CUDA device
	cudaErrorStartupFailure               =   0x7f,   ///< Startup failure
	cudaErrorApiFailureBase               =  10000    ///< API failure base
	};

	enum cudaMemcpyKind
	{
	cudaMemcpyHostToHost          =   0,      ///< Host   -> Host
	cudaMemcpyHostToDevice        =   1,      ///< Host   -> Device
	cudaMemcpyDeviceToHost        =   2,      ///< Device -> Host
	cudaMemcpyDeviceToDevice      =   3       ///< Device -> Device
	};

	/**
	* Channel format kind
	*/
	/*DEVICE_BUILTIN*/
	enum cudaChannelFormatKind
	{
	cudaChannelFormatKindSigned           =   0,      ///< Signed channel format
	cudaChannelFormatKindUnsigned         =   1,      ///< Unsigned channel format
	cudaChannelFormatKindFloat            =   2,      ///< Float channel format
	cudaChannelFormatKindNone             =   3,      ///< No channel format
	};

	/**
	* CUDA Channel format descriptor
	*/
	/*DEVICE_BUILTIN*/
	struct cudaChannelFormatDesc
	{
	int                        x; ///< x
	int                        y; ///< y
	int                        z; ///< z
	int                        w; ///< w
	enum cudaChannelFormatKind f; ///< Channel format kind
	};
	/*        END NVIDIA CODE          */

	inline const char *cudaGetErrorString(cudaError err) { return "CUDA Error: CUDA Error when no CUDA used (?!)"; }

	inline cudaError cudaMalloc(void** devPtr, size_t count)
	{
		*devPtr = malloc(count);
		return cudaSuccess;
	}
	inline cudaError cudaFree(void* devPtr)
	{
		free(devPtr);
		return cudaSuccess;
	}
	inline cudaError cudaMemcpy(void* dst, const void* src, size_t count, cudaMemcpyKind kind)
	{
		memcpy(dst, src, count);
		return cudaSuccess;
	}

	/**
	* CUDA array
	*/
	/*DEVICE_BUILTIN*/
	struct cudaArray;

	inline cudaError cudaFreeArray(cudaArray* array ) { assert(0); }
	inline cudaError cudaMallocArray(cudaArray** array, const struct cudaChannelFormatDesc* desc, size_t width, size_t height ) { assert(0); }
	inline cudaError cudaMemcpy2DToArray(cudaArray* dstArray, size_t dstX, size_t dstY, const void* src, size_t spitch, size_t width, size_t height, enum cudaMemcpyKind kind) { assert(0); }
#endif

void abort_on_cuda_error(cudaError err);
#define CUDA_ASSERT(err) \
	if((err) != cudaSuccess) { abort_on_cuda_error(err); }

// Based on NVIDIA's LinuxStopWatch class
// TODO: Should be moved to libpeyton
class stopwatch
{
protected:
	struct timeval  start_time;	// Start of measurement
	float  diff_time;		// Time difference between the last start and stop
	float  total_time;		// TOTAL time difference between starts and stops
	bool running;			// flag if the stop watch is running
	int clock_sessions;		// Number of times clock has been started and stopped (for averaging)

public:
	stopwatch() :
		start_time(),
		diff_time(0.0),
		total_time(0.0),
		running(false),
		clock_sessions(0)
	{ }

	// Start time measurement
	void start()
	{
		gettimeofday( &start_time, 0);
		running = true;
	}

	// Stop time measurement and increment add to the current diff_time summation
	// variable. Also increment the number of times this clock has been run.
	void stop()
	{
		diff_time = getDiffTime();
		total_time += diff_time;
		running = false;
		clock_sessions++;
	}

	// Reset the timer to 0. Does not change the timer running state but does
	// recapture this point in time as the current start time if it is running.
	void reset()
	{
		diff_time = 0;
		total_time = 0;
		clock_sessions = 0;
		if( running )
		{
			gettimeofday( &start_time, 0);
		}
	}

	// Time in sec. after start. If the stop watch is still running (i.e. there
	// was no call to stop()) then the elapsed time is returned added to the
	// current diff_time sum, otherwise the current summed time difference alone
	// is returned.
	float getTime() const
	{
		// Return the TOTAL time to date
		float retval = total_time;
		if(running)
		{
			retval += getDiffTime();
		}
		return 0.001 * retval;
	}

	// Time in msec. for a single run based on the total number of COMPLETED runs
	// and the total time.
	float getAverageTime() const
	{
		return 0.001 * total_time/clock_sessions;
	}

	int nSessions() const
	{
		return clock_sessions;
	}

private:

	// helpers functions

	float getDiffTime() const
	{
		struct timeval t_time;
		gettimeofday( &t_time, 0);

		// time difference in milli-seconds
		return  (float) (1000.0 * ( t_time.tv_sec - start_time.tv_sec)
				+ (0.001 * (t_time.tv_usec - start_time.tv_usec)) );
	}
};

// Rounds up v to nearest integer divisable by mod
inline int roundUpModulo(int v, int mod)
{
	int r = v % mod;
	int pitch = r ? v + (mod-r) : v;
	return pitch;
}


// For CPU versions of GPU algorithms
#if !__CUDACC__
namespace gpuemu // prevent collision with nvcc's symbols
{
	extern __TLS uint3 blockIdx;
	extern __TLS uint3 threadIdx;
	extern __TLS uint3 blockDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
	extern __TLS uint3 gridDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
}
using namespace gpuemu;
#endif

inline __device__ uint32_t threadID()
{
	// this supports 3D grids with 1D blocks of threads
	// NOTE: This could/should be optimized to use __mul24 (but be careful not to overflow a 24-bit number!)
	// Number of cycles (I think...): 4*4 + 16 + 3*4
#if 0 && __CUDACC__
	// This below is untested...
	const uint32_t id =
		  threadIdx.x
		+ __umul24(blockDim.x, blockIdx.x)
		+ __umul24(blockDim.x, blockIdx.y) * gridDim.x
		+ __umul24(blockDim.x, blockIdx.z) * __umul24(gridDim.x, gridDim.y);
#else
	// 16 + 16 + 16 cycles (assuming FMAD)
	const uint32_t id = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
#endif
	return id;
}

/* Support structures */
struct kernel_state
{
	__host__ __device__ uint32_t nthreads() const
	{
		uint32_t nthreads;
		nthreads = m_end - m_begin;
		nthreads = nthreads / m_step + (nthreads % m_step ? 1 : 0);
		return nthreads;
	}

	template<typename T>
	__device__ T* sharedMemory() const
	{
		return (T*)(shmem + sharedBytesPerThread*threadIdx.x);
	}
	uint32_t shmemBytes() const { return sharedBytesPerThread; }

	/////////////////////////
	uint32_t m_begin, m_step, m_end;
	uint32_t sharedBytesPerThread;

	kernel_state(uint32_t b, uint32_t e, uint32_t s, uint32_t shbytes = 0)
		: m_begin(b), m_end(e), m_step(s), sharedBytesPerThread(shbytes)
	{
	}

//	__device__ uint32_t row() const { uint32_t row = threadID(); return row < nthreads() ? beg + row : (uint32_t)-1; }
	__device__ uint32_t row_begin() const { return m_begin + m_step * threadID(); }
	__device__ uint32_t row_end()   const { uint32_t tmp = m_begin + m_step * (threadID()+1); return tmp <= m_end ? tmp : m_end; }
};

typedef kernel_state otable_ks;

/*  Support macros  */

extern stopwatch kernelRunSwatch;

bool calculate_grid_parameters(dim3 &gridDim, int threadsPerBlock, int neededthreads, int dynShmemPerThread, int staticShmemPerBlock);

//
// Support for run-time selection of execution on GPU or CPU
//
#if HAVE_CUDA
	bool gpuExecutionEnabled(const char *kernel);

	extern __TLS int  active_compute_device;		// helpers
	inline int gpuGetActiveDevice() { return active_compute_device; }

	struct activeDevice
	{
		int prev_active_device;

		activeDevice(int dev)
		{
			prev_active_device = active_compute_device;
			active_compute_device = dev;
		}

		~activeDevice()
		{
			active_compute_device = prev_active_device;
		}
	};
#else
	inline bool gpuExecutionEnabled(const char *kernel) { return false; }
	inline int gpuGetActiveDevice() { return -1; }
	struct activeDevice
	{
		activeDevice(int dev) {}
	};
#endif

#if HAVE_CUDA
	#define DECLARE_KERNEL(kDecl) \
		void cpulaunch_##kDecl; \
		void gpulaunch_##kDecl;

	#define CALL_KERNEL(kName, ...) \
		{ \
			swatch.start(); \
			activeDevice dev(gpuExecutionEnabled(#kName)? 0 : -1); \
			if(gpuGetActiveDevice() < 0) \
			{ \
				cpulaunch_##kName(__VA_ARGS__); \
			} \
			else \
			{ \
				gpulaunch_##kName(__VA_ARGS__); \
			} \
			swatch.stop(); \
			static bool firstTime = true; if(firstTime) { swatch.reset(); kernelRunSwatch.reset(); firstTime = false; } \
		}
#else
	// No CUDA
	#define DECLARE_KERNEL(kDecl) \
		void cpulaunch_##kDecl;

	#define CALL_KERNEL(kName, ...) \
	{ \
		swatch.start(); \
		cpulaunch_##kName(__VA_ARGS__); \
		swatch.stop(); \
		static bool firstTime = true; if(firstTime) { swatch.reset(); kernelRunSwatch.reset(); firstTime = false; } \
	}
#endif

#if HAVE_CUDA && !BUILD_FOR_CPU
	// Building GPU kernels
	#define KERNEL(ks, shmemPerThread, kDecl, kName, kArgs) \
		__global__ void gpu_##kDecl; \
		void gpulaunch_##kDecl \
		{ \
			int dynShmemPerThread = ks.shmemBytes() ? ks.shmemBytes() : shmemPerThread;  /* built in the algorithm */ \
		        int staticShmemPerBlock = 96;   /* read from .cubin file */ \
		        int threadsPerBlock = 192;      /* TODO: This should be computed as well */ \
			dim3 gridDim; \
			calculate_grid_parameters(gridDim, threadsPerBlock, ks.nthreads(), dynShmemPerThread, staticShmemPerBlock); \
			\
			kernelRunSwatch.start(); \
			gpu_##kName<<<gridDim, threadsPerBlock, threadsPerBlock*dynShmemPerThread>>>kArgs; \
			cudaError err = cudaThreadSynchronize();\
			if(err != cudaSuccess) { abort(); } \
			kernelRunSwatch.stop(); \
		} \
		__global__ void gpu_##kDecl
#endif

#if !HAVE_CUDA || BUILD_FOR_CPU
	// Building CPU kernels

/*	#if BUILD_FOR_CPU
		// Building CPU kernels in CUDA-enabled binary*/
		#define KERNEL_NAME(kDecl) cpulaunch_##kDecl
// 	#else
// 		// Building CPU kernels only
// 		#define KERNEL_NAME(kDecl) kDecl
// 	#endif

	#define KERNEL(ks, shmemPerThread, kDecl, kName, kArgs) \
		void cpu_##kDecl; \
		void KERNEL_NAME(kDecl) \
		{ \
			int dynShmemPerThread = shmemPerThread;      /* built in the algorithm */ \
		        int staticShmemPerBlock = 96;   /* read from .cubin file */ \
		        int threadsPerBlock = 192;      /* TODO: This should be computed as well */ \
			calculate_grid_parameters((dim3 &)gridDim, threadsPerBlock, ks.nthreads(), dynShmemPerThread, staticShmemPerBlock); \
			\
			kernelRunSwatch.start(); \
			threadIdx.x = threadIdx.y = threadIdx.z = 0; \
			blockIdx.x = blockIdx.y = blockIdx.z = 0; \
			blockDim.x = threadsPerBlock; blockDim.y = blockDim.z = 1; \
			for(uint32_t __i=0; __i != ks.nthreads(); __i++) \
			{ \
				if(0) { MLOG(verb1) << "t(" << threadIdx.x << "), b(" << blockIdx.x << "," << blockIdx.y << "," << blockIdx.z << ")"; } \
				if(0) { MLOG(verb1) << "  db(" << blockDim.x << "," << blockDim.y << "," << blockDim.z << ")"; } \
				if(0) { MLOG(verb1) << "  dg(" << gridDim.x << "," << gridDim.y << "," << gridDim.z << ")"; } \
				cpu_##kName kArgs; \
				threadIdx.x++; \
				if(threadIdx.x == blockDim.x) { threadIdx.x = 0; blockIdx.x++; } \
				if(blockIdx.x  == gridDim.x) { blockIdx.x = 0;  blockIdx.y++; } \
				if(blockIdx.y  == gridDim.y) { blockIdx.y = 0;  blockIdx.z++; } \
			} \
			kernelRunSwatch.stop(); \
		} \
		void cpu_##kDecl
#endif

// thin random number generator abstraction
struct rng_t
{
	virtual float uniform() = 0;
	virtual float gaussian(const float sigma) = 0;
	virtual ~rng_t() {}
	// interface compatibility with gpu_rng_t
	void load(const otable_ks &o) {}
};
#define ALIAS_GPU_RNG 0

// GPU random number generator abstraction
#if !HAVE_CUDA && ALIAS_GPU_RNG
typedef rng_t &gpu_rng_t;
#endif

#include <gsl/gsl_rng.h>
#if !__CUDACC__
#include <iostream>
#endif

inline uint32_t rng_mwc(uint32_t *xc)
{
	#define c (xc[0])
	#define x (xc[1])
	#define a (xc[2])

	uint64_t xnew = (uint64_t)a*x + c;
//	printf("%016llx\n", xnew);
	c = xnew >> 32;
	x = (xnew << 32) >> 32;
	return x;

	#undef c
	#undef x
	#undef a
}

#if HAVE_CUDA || !ALIAS_GPU_RNG
struct gpu_rng_t
{
	uint32_t *streams;	// pointer to RNG stream states vector
	uint32_t nstreams;	// number of initialized streams

//	gpu_rng_t(uint32_t s) : seed(s) { }
	gpu_rng_t(rng_t &rng);		// initialization constructor from existing rng

	__device__ float uniform() const
	{
		/*
			Marsaglia's Multiply-With-Carry RNG. For theory and details see:
			
				http://www.stat.fsu.edu/pub/diehard/cdrom/pscript/mwc1.ps
				http://www.ms.uky.edu/~mai/RandomNumber
				http://www.ast.cam.ac.uk/~stg20/cuda/random/index.html
		*/
		#define a  (((uint32_t*)shmem)[               threadIdx.x])
		#define c  (((uint32_t*)shmem)[  blockDim.x + threadIdx.x])
		#define xn (((uint32_t*)shmem)[2*blockDim.x + threadIdx.x])

		uint64_t xnew = (uint64_t)a*xn + c;
		c = xnew >> 32;
		xn = (xnew << 32) >> 32;
		return 2.32830643708e-10f * xn;

		#undef a
		#undef c
		#undef xn
	}

	__device__ float uniform_pos() const
	{
		float x;
		do { 
			x = uniform(); 
			} 
		while (x == 0.f);
		return x;
	}

	__device__ float gaussian(const float sigma)
	{
		float x, y, r2;

		do
		{
			/* choose x,y in uniform square (-1,-1) to (+1,+1) */
			x = -1.f + 2.f * uniform_pos();
			y = -1.f + 2.f * uniform_pos();

			/* see if it is in the unit circle */
			r2 = x * x + y * y;
		}
		while (r2 > 1.0f || r2 == 0.f);

		/* Box-Muller transform */
		return sigma * y * sqrt (-2.0f * logf (r2) / r2);
	}

	__device__ bool load(kernel_state &ks)
	{
#if 0
		((int32_t*)shmem)[threadIdx.x] = seed + threadID();
		if(!rng) rng = gsl_rng_alloc(gsl_rng_default);
#else
#if 0
		gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
		gsl_rng_set(rng, seed + threadID());
		((int32_t*)shmem)[threadIdx.x] = gsl_rng_get(rng);
		gsl_rng_free(rng);
#else
		// load the RNG state
		uint32_t tid = threadID();
		if(tid >= nstreams)
		{
			// we should somehow abort the entire kernel here
#if !__CUDACC__
			ASSERT(tid >= nstreams) {
				std::cerr << "threadID= " << tid << " >= nstreams=" << nstreams << "\n";
			}
#endif
			return false;
		};
		((uint32_t*)shmem)[               threadIdx.x] = streams[             tid];
		((uint32_t*)shmem)[  blockDim.x + threadIdx.x] = streams[  nstreams + tid];
		((uint32_t*)shmem)[2*blockDim.x + threadIdx.x] = streams[2*nstreams + tid];
		return true;
#endif
		//std::cerr << seed << " " << threadID() << " " << ((int32_t*)shmem)[threadIdx.x] << "\n";
#endif
	}

	__device__ void store(kernel_state &ks)
	{
#if 0
		if(threadIdx.x == 0)
		{
			// This "stores" the RNG state by storing the current
			// state of threadID=0 RNG, and disregarding the rest.
			seed = ((uint32_t *)shmem)[0];
		}
#else
		// store the RNG state
		uint32_t tid = threadID();
		streams[             tid] = ((uint32_t*)shmem)[               threadIdx.x];
		streams[  nstreams + tid] = ((uint32_t*)shmem)[  blockDim.x + threadIdx.x];
		streams[2*nstreams + tid] = ((uint32_t*)shmem)[2*blockDim.x + threadIdx.x];
#endif
	}

/*	void srand(uint32_t s)
	{
		seed = s;
	}*/
};
#endif

#endif
