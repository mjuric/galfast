#ifndef __gpu_h
#define __gpu_h

#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

// Possible states
// #define HAVE_CUDA 1 -- CUDA support is on
// 	#define BUILD_FOR_CPU 1	-- build CPU versions of kernels
//
// Note: __CUDACC__ will never be defined when BUILD_FOR_CPU is defined
//

#if HAVE_CUDA
	#include <cuda_runtime.h>

	bool cuda_init();
#else
	// Emulate CUDA types and keywords

	#define __device__
	#define __host__

	struct double2 { double x, y; };
	struct uint3 { unsigned int x, y, z; };

	struct dim3
	{
		unsigned int x, y, z;
		#if defined(__cplusplus)
		dim3(unsigned int x = 1, unsigned int y = 1, unsigned int z = 1) : x(x), y(y), z(z) {}
		dim3(uint3 v) : x(v.x), y(v.y), z(v.z) {}
		operator uint3(void) { uint3 t; t.x = x; t.y = y; t.z = z; return t; }
		#endif /* __cplusplus */
	};
#endif

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

// A pointer type that keeps the information of the type
// and size of the array it's pointing to
struct xptr
{
protected:
	friend struct GPUMM;

	char *base;		// data pointer
	uint32_t m_elementSize;	// size of array element (bytes)
	uint32_t dim[2];	// array dimensions (in elements). dim[0] == ncolumns == width, dim[1] == nrows == height
	uint32_t m_pitch[1];	// array pitch (in bytes). pitch[1] == width of the padded row (in bytes)

public:
	uint32_t elementSize() const { return m_elementSize; }
	uint32_t ncols()   const { return dim[0]; }
	uint32_t nrows()   const { return dim[1]; }
	uint32_t width()   const { return dim[0]; }
	uint32_t height()  const { return dim[1]; }
	void set_height(uint32_t h) { dim[1] = h; }
	uint32_t &pitch()  { return m_pitch[0]; }
	uint32_t pitch()  const { return m_pitch[0]; }
	uint32_t memsize() const { return nrows() * pitch(); }

	template<typename T> const T *get() const { return (const T*)base; }
	template<typename T> T *get() { return (T*)base; }

	operator bool() const { return base; }
	
	xptr(size_t es = 0, size_t ncol = 0, size_t nrow = 1, size_t p = 0) { init(es, ncol, nrow, p); }

	void alloc(size_t eSize = (size_t)-1, size_t ncol = (size_t)-1, size_t nrow = (size_t)-1, size_t ptch = (size_t)-1)
	{
		if(eSize == (size_t)-1) { eSize = elementSize(); } else { m_elementSize = eSize; }
		if(ncol == (size_t)-1) { ncol = ncols(); } else { dim[0] = ncol; }
		if(nrow == (size_t)-1) { nrow = nrows(); } else { dim[1] = nrow; }
		if(ptch == (size_t)-1) { ptch = pitch(); } else { m_pitch[0] = ptch; }

		delete [] base;
		base = new char[memsize()];
	}
	
	void free()
	{
		delete [] base;
		base = NULL;
	}
	~xptr() { }

	void init(size_t es = 0, size_t ncol = 0, size_t nrow = 1, size_t p = 0)
	{
		m_elementSize = es;
		dim[0] = ncol;
		dim[1] = nrow;
		m_pitch[0] = p;
		base = NULL;

		if(memsize()) { alloc(); }
	}
};

#include <map>
struct GPUMM 
{
	static const int gc_treshold = 512*1024*1024;

	static const int NOT_EXIST = -1;
	static const int NEWPTR = 0;
	static const int SYNCED_TO_DEVICE = 1;
	static const int SYNCED_TO_HOST = 2;
	static const int RELEASED_TO_HOST = 3;

protected:
	struct gpu_ptr
	{
		xptr ptr;
		int lastop;

		gpu_ptr() : lastop(NEWPTR) {}
	};
	std::map<void *, gpu_ptr> gpuPtrs;

	size_t allocated() const;
	void gc();	// do garbage collection
public:
#if HAVE_CUDA
	xptr syncToDevice(const xptr &hptr);
	void syncToHost(xptr &hptr);
	int lastOp(const xptr &hptr)
	{
		if(!gpuPtrs.count(hptr.base)) { return NOT_EXIST; }
		return gpuPtrs[hptr.base].lastop;
	}
#else
	xptr syncToDevice(const xptr &hptr) const { return hptr; }
	void syncToHost(xptr &hptr) { }
#endif
};
extern GPUMM gpuMMU;

//#define __TLS __thread
#define __TLS

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
/*	const uint32_t id =
		  threadIdx.x
		+ __umul24(blockDim.x, blockIdx.x)
		+ __umul24(blockDim.x, blockIdx.y) * gridDim.x
		+ __umul24(blockDim.x, blockIdx.z) * __umul24(gridDim.x, gridDim.y)*/
	// 16 + 16 + 16 cycles (assuming MADD)
	const uint32_t id = ((blockIdx.z * gridDim.y + blockIdx.y) * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	return id;
}

/* Support structures */
struct kernel_state
{
	uint32_t m_nthreads;
	__host__ __device__ uint32_t nthreads() const
	{
		return m_nthreads;
	}

	/////////////////////////
	uint32_t begin;
	__device__ uint32_t row() const { uint32_t row = threadID(); return row < nthreads() ? begin + row : (uint32_t)-1; }

	kernel_state(uint32_t b, uint32_t e) : begin(b), m_nthreads(e - b) { }
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
			activeDevice dev(gpuExecutionEnabled(#kName)? 0 : -1); \
			if(gpuGetActiveDevice() < 0) \
			{ \
				cpulaunch_##kName(__VA_ARGS__); \
			} \
			else \
			{ \
				gpulaunch_##kName(__VA_ARGS__); \
			} \
		}
#else
	// No CUDA
	#define DECLARE_KERNEL(kDecl) \
		void cpulaunch_##kDecl;

	#define CALL_KERNEL(kName, ...) \
		cpulaunch_##kName(__VA_ARGS__);
#endif

#if HAVE_CUDA && !BUILD_FOR_CPU
	// Building GPU kernels
	#define KERNEL(ks, kDecl, kName, kArgs) \
		__global__ void gpu_##kDecl; \
		void gpulaunch_##kDecl \
		{ \
			int dynShmemPerThread = 4;      /* built in the algorithm */ \
		        int staticShmemPerBlock = 96;   /* read from .cubin file */ \
		        int threadsPerBlock = 128;      /* TODO: This should be computed as well */ \
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

	#define KERNEL(ks, kDecl, kName, kArgs) \
		void cpu_##kDecl; \
		void KERNEL_NAME(kDecl) \
		{ \
			const int threadsPerBlock = 32; \
			threadIdx.x = 0; \
			blockIdx.x = blockIdx.y = blockIdx.z = 0; \
			blockDim.x = threadsPerBlock; blockDim.y = blockDim.z = 1; \
			gridDim.x = 512; gridDim.y = 512; gridDim.z = 1; \
			kernelRunSwatch.start(); \
			for(uint32_t __i=0; __i != ks.nthreads(); __i++) \
			{ \
				if(0) { MLOG(verb1) << "t(" << threadIdx.x << "), b(" << blockIdx.x << "," << blockIdx.y << "," << blockIdx.z << ")"; } \
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

#if __CUDACC__
extern __shared__ char shmem[];
#else
extern __TLS char shmem[16384];
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

#if HAVE_CUDA || !ALIAS_GPU_RNG
struct gpu_rng_t
{
	uint32_t seed;

	gpu_rng_t(uint32_t s) : seed(s) { }
	gpu_rng_t(rng_t &rng);		// initialization constructor from existing rng

	__device__ float uniform() const
	{
		#define IA 16807
		#define IM 2147483647
		#define AM (1.0f/IM)
		#define IQ 127773
		#define IR 2836
		#define MASK 123459876
		#define idum ((int32_t*)shmem)[threadIdx.x]

		int32_t k;
		float ans;

		idum ^= MASK;			// XORing with MASK allows use of zero and other
		k=idum/IQ;			// simple bit patterns for idum.
		idum=IA*(idum-k*IQ)-IR*k;	// Compute idum=(IA*idum) % IM without over-
		if (idum < 0) idum += IM;	// flows by Schrage's method.
		ans=AM*idum; 			// Convert idum to a floating result.
		idum ^= MASK; 			// Unmask before return.

		#undef idum
		return ans;
	}

	__device__ float uniform_pos() const
	{
		float x;
		do { x = uniform(); } while (x == 0.f);
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

	__device__ void load(kernel_state &ks)
	{
		((int32_t*)shmem)[threadIdx.x] = seed + threadID();
	}

	__device__ void store(kernel_state &ks)
	{
		if(threadIdx.x == 0)
		{
			// This "stores" the RNG state by storing the current
			// state of threadID=0 RNG, and disregarding the rest.
			seed = ((uint32_t *)shmem)[0];
		}
	}

	void srand(uint32_t s)
	{
		seed = s;
	}
};
#endif

#endif
