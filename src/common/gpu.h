#ifndef __gpu_h
#define __gpu_h

#include <time.h>
#include <sys/time.h>

#if HAVE_CUDA
	#include <cuda_runtime.h>
#else
	#define __device__
	#define __host__
	struct double2 { double x, y; };
#endif

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

#if HAVE_CUDA
#include <map>
struct GPUMM 
{
	static const int gc_treshold = 512*1024*1024;

	static const int NOT_EXIST = -1;
	static const int NEWPTR = 0;
	static const int SYNCED_TO_DEVICE = 1;
	static const int SYNCED_TO_HOST = 2;
	static const int RELEASED_TO_HOST = 3;
	struct gpu_ptr
	{
		xptr ptr;
		int lastop;

		gpu_ptr() : lastop(NEWPTR) {}
	};
	std::map<void *, gpu_ptr> gpuPtrs;

	size_t allocated() const;
	void gc();	// do garbage collection
	xptr syncToDevice_aux(const xptr &hptr);
	void syncToHost_aux(xptr &hptr);

	xptr syncToDevice(const xptr &ptr)
	{
		return syncToDevice_aux(ptr);
	}

	void syncToHost(xptr &ptr)
	{
		syncToHost_aux(ptr);
	}

	int lastOp(const xptr &hptr)
	{
		if(!gpuPtrs.count(hptr.base)) { return NOT_EXIST; }
		return gpuPtrs[hptr.base].lastop;
	}
};
extern GPUMM gpuMMU;

#endif

/* Support structures */
struct kernel_state
{
#if HAVE_CUDA
	__device__ uint32_t threadIndex() const
	{
		// this supports 2D grids with 1D blocks of threads
		#if __CUDACC__
		return (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
		#else
		assert(0);
		#endif
	}
#else
	uint32_t m_threadIndex;
	uint32_t threadIndex() const { return m_threadIndex; }
	void set_threadIndex(uint32_t idx) { m_threadIndex = idx; }
#endif

	uint32_t m_nthreads;
	__host__ __device__ uint32_t nthreads() const
	{
		return m_nthreads;
	}

	/////////////////////////
	uint32_t begin;
	__device__ uint32_t row() const { uint32_t row = threadIndex(); return row < nthreads() ? begin + row : (uint32_t)-1; }

	kernel_state(uint32_t b, uint32_t e) : begin(b), m_nthreads(e - b) { }
};

typedef kernel_state otable_ks;

/*  Support macros  */

#ifdef HAVE_CUDA
	#define KERNEL(ks, kDecl, kName, kArgs) \
		__global__ void aux_##kDecl; \
		void kDecl \
		{ \
			int dynShmemPerThread = 4;      /* built in the algorithm */ \
		        int staticShmemPerBlock = 96;   /* read from .cubin file */ \
		        int threadsPerBlock = 64;       /* TODO: This should be computed as well */ \
			int gridDim[3]; \
			calculate_grid_parameters(gridDim, threadsPerBlock, ks.nthreads(), dynShmemPerThread, staticShmemPerBlock); \
			\
			dim3 grid; \
			grid.x = gridDim[0]; grid.y = gridDim[1]; \
			aux_##kName<<<grid, threadsPerBlock, threadsPerBlock*dynShmemPerThread>>>kArgs; \
		} \
		__global__ void aux_##kDecl

	bool calculate_grid_parameters(int gridDim[3], int threadsPerBlock, int neededthreads, int dynShmemPerThread, int staticShmemPerBlock);

#else
	#define DECLARE_KERNEL(kernelName, ...) \
		void kernelName(__VA_ARGS__); \

	#define KERNEL(ks, kDecl, kName, kArgs) \
		void aux_##kDecl; \
		void kDecl \
		{ \
			for(uint32_t __i=0; __i != ks.nthreads(); __i++) \
			{ \
				ks.set_threadIndex(__i); \
				aux_##kName kArgs; \
			} \
		} \
		void aux_##kDecl
#endif

// GPU random number generator abstraction
#if __CUDACC__
extern __shared__ char memory[];
extern __shared__ int32_t mem_int32[];
#endif
#if HAVE_CUDA
struct gpu_rng_t
{
	uint32_t seed;
	gpu_rng_t(uint32_t s) : seed(s) {}

	__device__ float uniform() const
	{
	#if __CUDACC__
		#define IA 16807
		#define IM 2147483647
		#define AM (1.0f/IM)
		#define IQ 127773
		#define IR 2836
		#define MASK 123459876
		#define idum mem_int32[threadIdx.x]

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
	#endif
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
	#if __CUDACC__
		uint32_t *seeds = (uint32_t *)memory;
		seeds[threadIdx.x] = seed + ks.threadIndex();
	#endif
	}

	__device__ void store(kernel_state &ks)
	{
	#if __CUDACC__
		if(threadIdx.x == 0)
		{
			uint32_t *seeds = (uint32_t *)memory;
			seed = seeds[0];
		}
	#endif
	}

	void srand(uint32_t s)
	{
		seed = s;
	}
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

#endif
