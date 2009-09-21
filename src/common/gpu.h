#ifndef __gpu_h
#define __gpu_h

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdint.h>

//////////////////////////////////////////////////////////////////////////
// Possible states of important defines:
//
// #define HAVE_CUDA 1 -- CUDA support is on
//	#define BUILD_FOR_CPU 1	-- build CPU versions of kernels
//
// Note: __CUDACC__ will never be defined when BUILD_FOR_CPU is defined
//
//////////////////////////////////////////////////////////////////////////

#include "cux.h"

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

	kernel_state(uint32_t b, uint32_t e, uint32_t s = 0xFFFFFFFF, uint32_t shbytes = 0)
		: m_begin(b), m_end(e), m_step(s == 0xFFFFFFFF ? 256 : s), sharedBytesPerThread(shbytes)
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
			cuxErrCheck( cudaThreadSynchronize() );\
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
	typedef prngs::gpu::mwc gpu_prng_impl;
	typedef prngs::cpu::mwc cpu_prng_impl;

	//
	// Maintains a multi-threaded random number generator instance
	// that persists through the application and is downloaded/uploaded
	// to GPU as needed. It is auto-initialized on first call, using the
	// seed provided by seeder rng, and freed on exit from the application.
	//
	// Usage:
	//	gpu_rng_t bla(seeder)
	//	bla.uniform(); ....
	//
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
		gpu_rng_t() {}
	};
#endif

/**
	rng_gsl_t -- GSL RNG implementation
*/
#include <gsl/gsl_randist.h>
struct rng_gsl_t : public rng_t
{
	bool own;
	gsl_rng *rng;

	rng_gsl_t(gsl_rng *r, bool own_ = false) : rng(r), own(own_) {}
	rng_gsl_t(unsigned long int seed)
		: own(true), rng(NULL)
	{
		rng = gsl_rng_alloc(gsl_rng_default);
		gsl_rng_set(rng, seed);
	}

	virtual ~rng_gsl_t() { if(own && rng) gsl_rng_free(rng); }

	virtual float uniform()
	{
		return (float)gsl_rng_uniform(rng);
	}
	virtual float gaussian(const float sigma) { return (float)gsl_ran_gaussian(rng, sigma); }
};

#endif
