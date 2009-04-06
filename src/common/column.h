#ifndef __column_h
#define __column_h

#include "gpu.h"

template<typename T>
struct column_gpu
{
	T *base;
	uint32_t pitch;

	__device__ T &operator()(const size_t row, const size_t elem)	// 2D table column accessor
	{
		return *((T*)((char*)base + elem * pitch) + row);
	}
	__device__ T &operator[](const size_t i)	// 1D table column accessor (i == the table row)
	{
		return base[i];
	}

	column_gpu(xptr &ptr) : base(ptr.get<T>()), pitch(ptr.pitch()) {}
};

template<typename T>
struct column
{
	typedef column_gpu<T> gpu_t;

	xptr base;
	bool onGPU;

	size_t nrows() const { return base.width(); }
	size_t width() const { return base.height(); }

	column()
	{
		onGPU = false;
	}

	void toHost()
	{
	#if HAVE_CUDA
		if(!onGPU) return;
		gpuMMU.syncToHost(base);
		onGPU = false;
	#endif
	}
	#ifndef __CUDACC__
	T *get() { toHost(); return base.get<T>(); }

	T& operator()(const size_t row, const size_t elem)	// 2D column accessor
	{
		toHost();
		return *((T*)(base.get<char>() + elem * base.pitch()) + row);
	}

	T &operator[](const size_t i)	// 1D column accessor (i == the row)
	{
		toHost();
		return base.get<T>()[i];
	}
	#endif

	// transfer data to GPU
	operator gpu_t()
	{
	#if HAVE_CUDA
		xptr gptr = gpuMMU.syncToDevice(base);
		onGPU = true;
		return gpu_t(gptr);
	#else
		return gpu_t(base);
	#endif
	}
private:
	// prevent copying
	column(const column<T> &);
	column<T>& operator=(const column<T> &a);
};

// convenience typedefs
namespace column_types
{
	typedef column<double> cdouble;
	typedef column<int> cint;
	typedef column<float> cfloat;
};

#endif
