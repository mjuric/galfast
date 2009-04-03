#ifndef __gpu_h
#define __gpu_h

#ifndef __CUDACC__
	#define __device__
	#define __host__
#endif

/* Support structures */
struct kernel_state
{
#if __CUDACC__
	__device__ uint32_t threadIndex() const
	{
		// this supports 2D grids with 1D blocks of threads
		return (blockIdx.y * gridDim.x + blockIdx.x) * blockDim.x + threadIdx.x;
	}
#else
	uint32_t m_threadIndex;
	uint32_t threadIndex() const { return m_threadIndex; }
	void set_threadIndex(uint32_t idx) { m_threadIndex = idx; }
#endif

	uint32_t m_nthreads;
	__host__ __device__ uint32_t nthreads() const { return m_nthreads; }
	kernel_state(uint32_t nthr) : m_nthreads(nthr) {}
};

struct otable_ks : public kernel_state
{
	uint32_t begin;
	__device__ uint32_t row() const { uint32_t row = threadIndex(); return row < nthreads() ? begin + row : (uint32_t)-1; }

	otable_ks(uint32_t b, uint32_t e) : begin(b), kernel_state(e - b) { }
};


/*  Support macros  */

#ifdef __CUDACC__
	#define KERNEL(ks, kDecl, kName, kArgs) \
		__global__ void aux_##kDecl; \
		void kDecl \
		{ \
			int dynShmemPerThread = 4;     /* built in the algorithm */ \
		        int staticShmemPerBlock = 32;   /* read from .cubin file */ \
		        int threadsPerBlock = 480; /* TODO: This should be computed as well */ \
			int gridDim[3]; \
			calculate_grid_parameters(gridDim, threadsPerBlock, ks.nthreads(), dynShmemPerThread, staticShmemPerBlock); \
			\
			dim3 grid; \
			grid.x = gridDim[0]; grid.y = gridDim[1]; \
			aux_##kName<<<grid, threadsPerBlock, threadsPerBlock*dynShmemPerThread>>>kArgs; \
		} \
		__global__ void aux_##kDecl
#else
	#define DECLARE_KERNEL(kernelName, ...) \
		void kernelName(__VA_ARGS__); \

	#define KERNEL(ks, kDecl, kCall) \
		void aux_##kDecl; \
		void kDecl \
		{ \
			for(uint32_t __i=0; __i != ks.nthreads(); __i++) \
			{ \
				ks.set_threadIndex(__i); \
				aux_##kCall; \
			} \
		} \
		void aux_##kDecl
#endif

bool calculate_grid_parameters(int gridDim[3], int threadsPerBlock, int neededthreads, int dynShmemPerThread, int staticShmemPerBlock);

// GPU random number generator abstraction
#if __CUDACC__
extern __shared__ char memory[];
extern __shared__ int32_t mem_int32[];
#endif
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
		do { x = uniform(); } while (x == 0);
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
		return sigma * y * sqrt (-2.0 * log (r2) / r2);
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
