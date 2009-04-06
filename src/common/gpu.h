#ifndef __gpu_h
#define __gpu_h

#include <cuda_runtime.h>

/*#ifndef __CUDACC__
	#define __device__
	#define __host__
#endif*/
// A pointer type that keeps the information of the type
// and size of the array it's pointing to
template<typename T>
struct xptr
{
	T *base;		// data pointer
	uint32_t m_elementSize;	// size of array element (bytes)
	uint32_t dim[2];	// array dimensions (in elements). dim[0] == ncolumns == width, dim[1] == nrows == height
	uint32_t m_pitch[1];	// array pitch (in bytes). pitch[1] == width of the padded row (in bytes)

	uint32_t elementSize() const { return m_elementSize; }
	uint32_t ncols()   const { return dim[0]; }
	uint32_t nrows()   const { return dim[1]; }
	uint32_t width()   const { return dim[0]; }
	uint32_t height()  const { return dim[1]; }
	uint32_t &pitch()  { return m_pitch[0]; }
	uint32_t pitch()  const { return m_pitch[0]; }
	uint32_t memsize() const { return nrows() * pitch(); }
	
	const T *get() const { return base; }
	T *get() { return base; }

/*	operator T*() { return base; }
	operator char*() { return (char *)base; }*/
//	operator const xptr<void>&() const { return *(xptr<void> *)this; }
//	template<typename TT> xptr<TT>& cvt() { return *(xptr<TT> *)this; }
// 	T& operator [](uint32_t idx) { return base[idx]; }
// 	T& operator [](uint32_t idx) const { return base[idx]; }

	xptr(size_t es = 0, size_t ncol = 0, size_t nrow = 0, size_t p = 0, T* d = NULL)
	{
		m_elementSize = es;
		dim[0] = ncol;
		dim[1] = nrow;
		m_pitch[0] = p;
		base = d;
	}
};

#ifndef __CXUDACC__
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
		xptr<void> ptr;
		int lastop;

		gpu_ptr() : lastop(NEWPTR) {}
	};
	std::map<void *, gpu_ptr> gpuPtrs;

	size_t allocated() const;
	void gc();	// do garbage collection
	xptr<void> syncToDevice_aux(const xptr<void> &hptr);
	void syncToHost_aux(xptr<void> &hptr);

	template<typename T>
	xptr<T> syncToDevice(const xptr<T> &ptr)
	{
		xptr<void> ret = syncToDevice_aux(*(const xptr<void> *)&ptr);
		return *(xptr<T> *)&ret;
	}

	template<typename T>
	void syncToHost(xptr<T> &ptr)
	{
		syncToHost_aux(*(xptr<void> *)&ptr);
	}

	template<typename T>
	int lastOp(const xptr<T> &hptr)
	{
		return lastOp(*(const xptr<void> *)&hptr);
	}
	
	int lastOp(const xptr<void> &hptr)
	{
		if(!gpuPtrs.count(hptr.base)) { return NOT_EXIST; }
		return gpuPtrs[hptr.base].lastop;
	}
};
extern GPUMM gpuMMU;

#endif

#if 0
template<typename T>
device_ptr
{
	T *data;
	int devId;
};

template typename<T>
struct xptr : public xptr_base
{
	device_ptr<T> get(int devId)
	{
		return mm.get(this, devId);
	}
};

struct gpu_mm
{
	struct {
		void *ptr;
		bool locked;
		bool onGpu;

		gpu_ptr(void ptr_, locked_ = false, ongpu = false) : ptr(ptr_), locked(locked_), onGpu(ongpu) { }
	} gpu_ptr;
	std::map<void *, gpu_ptr> host2gpu;

	template<typename T>
	struct {
		T *gpu_ptr;
		T *cpu_ptr;
		gpu_mm &mm;
		size_t refcnt;

		locked_ptr(gpu_mm &mm_, void *c) : mm(mm), cpu_ptr((T*)c), gpu_ptr((T*)g)
		{
			gpu_mm.host2gpu[c]
			gpu_mm.lock_aux(cpu_ptr);
			refcnt = 1;
		}
		~locked_ptr() { mm.release_aux(cpu_ptr); }
	} locked_ptr;

	template<typename T> locked_ptr<T> lock(T *c)
	{
		return locked_ptr<T>(this, c);
	}
};
#endif

/* Support structures */
struct kernel_state
{
#if 1
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

#ifdef __CUDACC__
	#define KERNEL(ks, kDecl, kName, kArgs) \
		__global__ void aux_##kDecl; \
		void kDecl \
		{ \
			int dynShmemPerThread = 4;     /* built in the algorithm */ \
		        int staticShmemPerBlock = 32;   /* read from .cubin file */ \
		        int threadsPerBlock = 1; /* TODO: This should be computed as well */ \
			int gridDim[3]; \
			calculate_grid_parameters(gridDim, threadsPerBlock, ks.nthreads(), dynShmemPerThread, staticShmemPerBlock); \
			\
			dim3 grid; \
			grid.x = gridDim[0]; grid.y = gridDim[1]; \
			aux_##kName<<<1, 1, threadsPerBlock*dynShmemPerThread>>>kArgs; \
		} \
		__global__ void aux_##kDecl

//			aux_##kName<<<grid, threadsPerBlock, threadsPerBlock*dynShmemPerThread>>>kArgs;

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
