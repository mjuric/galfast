#ifndef __column_h
#define __column_h

#include "gpu.h"

template<typename T>
struct column : public xptrng::tptr<T>
{
	typedef xptrng::gptr<T> gpu_t;
	typedef xptrng::hptr<T> host_t;

	// NOTE: WARNING: The meaning of width/height is being _inverted_ here
	// to bring it in line with the typical table metaphore. In memory, however
	// the array is stored such that nrows() is its width. This allows the GPU
	// to coalesce memory accesses to adjacent rows of the table.
	// NOTE #2: The std::max() in width() is there to represent the 1D arrays
	// as 2D arrays with 2nd dimension equal to 1
	uint32_t nrows() const { return xptrng::tptr<T>::width(); }
	uint32_t width() const { return std::max(xptrng::tptr<T>::height(), 1U); }

//	column() : xptrng::tptr<T>(0, 1) { }
	column() { }
	void reshape(const column<T> &t)
	{
		(xptrng::tptr<T> &)(*this) = t.clone();
	}
	void resize(size_t nrows, size_t width, int es = 0)
	{
		if(!es) { es = this->elementSize(); }
		(xptrng::tptr<T>&)(*this) = xptrng::tptr<T>(nrows, width, 0., es);
	}
	T *get() const { return this->syncToHost(); }
	column<T>& operator=(void *ptr)
	{
		assert(ptr == NULL); // can only be used to set it to 0 (for now)
		(xptrng::tptr<T>&)(*this) = NULL; //xptrng::tptr<T>();
		return *this;
	}

	// GPU interface
	operator gpu_t() { return (xptrng::tptr<T>&)(*this); }
	operator host_t() { return (xptrng::tptr<T>&)(*this); }
private:
	// prevent copying
	column(const column<T> &);
	//column<T>& operator=(const column<T> &a);
};

// convenience typedefs
namespace column_types
{
	typedef column<double> cdouble;
	typedef column<int> cint;
	typedef column<float> cfloat;
};

#endif
