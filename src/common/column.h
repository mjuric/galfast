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

	column_gpu(xptr<T> &ptr) : base(ptr.get()), pitch(ptr.pitch()) {}
};

template<typename T>
struct column
{
	typedef column_gpu<T> gpu_t;

	xptr<T> base;
	bool onGPU;

	column(const xptr<T> &b) : base(b)
	{
		// check if there's a GPU version of this pointer
		onGPU = gpuMMU.lastOp(base) == GPUMM::SYNCED_TO_DEVICE;
	}

	__host__ T& operator()(const size_t row, const size_t elem)	// 2D column accessor
	{
		if(onGPU)
		{
			gpuMMU.syncToHost(base);
			onGPU = false;
		}
		return *((T*)((char*)base.get() + elem * base.pitch()) + row);
	}

	__host__ T &operator[](const size_t i)	// 1D column accessor (i == the row)
	{
		if(onGPU)
		{
			gpuMMU.syncToHost(base);
			onGPU = false;
		}
		return base.get()[i];
	}

	operator gpu_t()
	{
		xptr<T> gptr = gpuMMU.syncToDevice(base);
		onGPU = true;
		return gpu_t(gptr);
	}
};

// convenience typedefs
namespace column_types
{
	typedef column<double> cdouble;
	typedef column<int> cint;
	typedef column<float> cfloat;
};

#endif
