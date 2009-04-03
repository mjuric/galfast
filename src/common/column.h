#ifndef __column_h
#define __column_h

#include "gpu.h"

template<typename T>
struct column
{
	T *base;
	size_t pitch;

	column(void *b, size_t p = 0) : base((T*)b), pitch(p) {}

	__device__ T &operator()(const size_t row, const size_t elem)	// 2D column accessor
	{
		return *((T*)((char*)base + elem * pitch) + row);
	}

	__device__ T &operator[](const size_t i)	// 1D column accessor (i == the row)
	{
		return base[i];
	}
	
	__device__ T &val(const size_t i)	// 1D column accessor (i == the row)
	{
		return base[i];
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
