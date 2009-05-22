#ifndef __gpu_h
#define __gpu_h

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

//////////////////////////////////////////////////////////////////////////
// Possible states of important defines:
//
// #define HAVE_CUDA 1 -- CUDA support is on
//	#define BUILD_FOR_CPU 1	-- build CPU versions of kernels
//
// Note: __CUDACC__ will never be defined when BUILD_FOR_CPU is defined
//
//////////////////////////////////////////////////////////////////////////


// Thread local storage -- use for shared memory emulation in CPU mode
//#define __TLS __thread
#define __TLS

//////////////////////////////////////////////////////////////////////////
// CUDA (or emulation) APIs
//////////////////////////////////////////////////////////////////////////
#if HAVE_CUDA
	#include <cuda_runtime.h>
#else
	#include "cuda_emulation.h"
#endif

//////////////////////////////////////////////////////////////////////////
// Shared memory access and CPU emulation
//////////////////////////////////////////////////////////////////////////
#if __CUDACC__
	extern __shared__ char impl_shmem[];
#else
	// For CPU versions of GPU algorithms
	extern __TLS char impl_shmem[16384];
	namespace gpuemu // prevent collision with nvcc's symbols
	{
		extern __TLS uint3 blockIdx;
		extern __TLS uint3 threadIdx;
		extern __TLS uint3 blockDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
		extern __TLS uint3 gridDim;		// Note: uint3 instead of dim3, because __TLS variables have to be PODs
	}
	using namespace gpuemu;
#endif
// template<typename T> inline __device__ T* shmem()              { return (T*)impl_shmem; }
// template<typename T> inline __device__ T& shmem(const int idx) { return ((T*)impl_shmem)[idx]; }
// 
// // explicit instantiations are necessary to work around a bug in nvcc, that appears to
// // forget __device__ when automatically instantiating the templates above (?)
// template<> inline __device__ uint32_t& shmem(const int idx) { return ((uint32_t*)impl_shmem)[idx]; }
// template<> inline __device__      int& shmem(const int idx) { return      ((int*)impl_shmem)[idx]; }

#define shmem(type) ((type*)impl_shmem)

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
extern stopwatch kernelRunSwatch;

//////////////////////////////////////////////////////////////////////////
// Useful host and device functions
//////////////////////////////////////////////////////////////////////////

// Rounds up v to nearest integer divisable by mod
inline int roundUpModulo(int v, int mod)
{
	int r = v % mod;
	int pitch = r ? v + (mod-r) : v;
	return pitch;
}

// Computes the linear ID of the thread
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

bool cuda_init();
bool calculate_grid_parameters(dim3 &gridDim, int threadsPerBlock, int neededthreads, int dynShmemPerThread, int staticShmemPerBlock);

//////////////////////////////////////////////////////////////////////////
// cux/tptr/gptr/hptr abstraction layer
//////////////////////////////////////////////////////////////////////////
#include "cuda_cux.h"

//////////////////////////////////////////////////////////////////////////
// Support structures and macros
//////////////////////////////////////////////////////////////////////////
struct kernel_state
{
	__host__ __device__ uint32_t nthreads() const
	{
		uint32_t nthreads;
		nthreads = m_end - m_begin;
		nthreads = nthreads / m_step + (nthreads % m_step ? 1 : 0);
		return nthreads;
	}

	// Return the pointer to this thread's piece of shared memory
	template<typename T>
	__device__ T* sharedMemory() const
	{
		return (T*)(shmem(char) + sharedBytesPerThread*threadIdx.x);
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

//////////////////////////////////////////////////////////////////////////
// Run-time selection of GPU/CPU execution
//////////////////////////////////////////////////////////////////////////
#if HAVE_CUDA
	// Return true if <kernel> should execute on GPU
	bool gpuExecutionEnabled(const char *kernel);

	extern __TLS int  active_compute_device;
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
	//
	// Not compiling with CUDA -- always run on CPU
	//
	inline bool gpuExecutionEnabled(const char *kernel) { return false; }
	inline int gpuGetActiveDevice() { return -1; }
	struct activeDevice
	{
		activeDevice(int dev) {}
	};
#endif

//
// Building with CUDA support. Declarations and runtime stubs for
// calling the GPU or CPU kernel versions.
//
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
#endif

//
// No CUDA at all. Only declare, build and call host-based kernels
//
#if !HAVE_CUDA
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

//
// Have CUDA and building GPU kernel.
//
#if HAVE_CUDA && !BUILD_FOR_CPU
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

//
// No CUDA, or building the CPU version of the kernel.
//
#if !HAVE_CUDA || BUILD_FOR_CPU
	#define KERNEL(ks, shmemPerThread, kDecl, kName, kArgs) \
		void cpu_##kDecl; \
		void cpulaunch_##kDecl \
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

//////////////////////////////////////////////////////////////////////////
// GPU random number generators
//////////////////////////////////////////////////////////////////////////
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

#if HAVE_CUDA || !ALIAS_GPU_RNG
	#include "cuda_rng.h"
	//typedef gpu_prng::mwc gpu_rng_t;
	typedef prngs::gpu::mwc gpu_prng_impl;
	typedef prngs::cpu::mwc cpu_prng_impl;
	struct gpu_rng_t : public gpu_prng_impl
	{
		struct persistent_rng
		{
			gpu_prng_impl gpuRNG;
			cpu_prng_impl cpuRNG;
			enum { EMPTY, GPU, CPU } state;

			persistent_rng() : state(EMPTY) { }
			~persistent_rng()
			{
				gpuRNG.free();
				cpuRNG.free();
			}

			gpu_prng_impl &get(rng_t &seeder);
		};
		static persistent_rng gpuRNG;

		gpu_rng_t(rng_t &seeder)
		{
			(gpu_prng_impl&)*this = gpuRNG.get(seeder);
		}
	};
#endif

#if 0
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

#endif
