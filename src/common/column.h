#ifndef __column_h
#define __column_h

#include "gpu2.h"

#if 0
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

	#ifndef __CUDACC__ // because nvcc barfs on explicit template member-function calls
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

	// GPU interface
	operator gpu_t()
	{
		if(gpuGetActiveDevice() >= 0)
		{
			xptr gptr = gpuMMU.syncToDevice(base);
			onGPU = true;
			return gpu_t(gptr);
		}
		else
		{
			toHost();
			return gpu_t(base);
		}
	}
private:
	// prevent copying
	column(const column<T> &);
	column<T>& operator=(const column<T> &a);
};
#endif

#if 1

template<typename T>
struct column : public xptrng::tptr<T>
{
	typedef xptrng::gptr<T> gpu_t;

	// NOTE: WARNING: The meaning of width/height is being _inverted_ here
	// to bring it in line with the typical table metaphore. In memory, however
	// the array is stored such that nrows() is its width. This allows the GPU
	// to coalesce memory accesses to adjacent rows of the table.
	size_t nrows() const { return xptrng::tptr<T>::width(); }
	size_t width() const { return xptrng::tptr<T>::height(); }

	column() : xptrng::tptr<T>(0, 1) { }
	void reshape(const column<T> &t)
	{
		(xptrng::tptr<T> &)(*this) = t.clone();
	}
	void resize(size_t nrows, size_t width, int es = 0)
	{
		if(!es) { es = this->elementSize(); }
		(xptrng::tptr<T>&)(*this) = xptrng::tptr<T>(nrows, width, es);
	}
	T *get() const { return this->syncToHost(); }
	column<T>& operator=(void *ptr)
	{
		assert(ptr == NULL); // can only be used to set it to 0 (for now)
		(xptrng::tptr<T>&)(*this) = xptrng::tptr<T>();
		return *this;
	}
// // 	column(size_t nrows, size_t width)
// // 		: xptrng::tptr<T>(nrows, width)
// // 	{}
// // 	column<T>& operator=(const column<T> &a)
// // 	{
// // 		return (xptrng::tptr<T>&)(*this) = a;
// // 	}

// 	#ifndef __CUDACC__ // because nvcc barfs on explicit template member-function calls
// 	T *get() { syncToHost(); return base.get<T>(); }
// 
// 	T& operator()(const size_t row, const size_t elem)	// 2D column accessor
// 	{
// 		syncToHost();
// 		return *((T*)(base.get<char>() + elem * base.pitch()) + row);
// 	}
// 
// 	T &operator[](const size_t i)	// 1D column accessor (i == the row)
// 	{
// 		syncToHost();
// 		return base.get<T>()[i];
// 	}
// 	#endif

	// GPU interface
	operator gpu_t()
	{
/*		if(gpuGetActiveDevice() >= 0)
		{
			xptr gptr = gpuMMU.syncToDevice(base);
			onGPU = true;
			return xptrng::make_gptr<T>(gptr.get<void>(), gptr.pitch());
		}
		else
		{
			toHost();
			return xptrng::make_gptr<T>(base.get<void>(), base.pitch());
		}*/
		return (xptrng::tptr<T>&)(*this);
	}
private:
	// prevent copying
	column(const column<T> &);
	//column<T>& operator=(const column<T> &a);
};
#endif

// convenience typedefs
namespace column_types
{
	typedef column<double> cdouble;
	typedef column<int> cint;
	typedef column<float> cfloat;
};

#endif
