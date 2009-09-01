#ifndef __column_h
#define __column_h

#include "gpu.h"

template<typename T>
struct column : public xptrng::xptr<T>
{
public:
	typedef xptrng::xptr<T> xptr_t;

	typedef xptrng::gptr<T, 2> gpu_t;
	typedef xptrng::hptr<T, 2> host_t;

public:
	// NOTE: WARNING: The meaning of width/height is being _inverted_ here
	// to bring it in line with the typical table metaphore. In memory, however
	// the array is stored such that nrows() is its width. This allows the GPU
	// to coalesce memory accesses to adjacent rows of the table.
	uint32_t nrows() const { return xptr_t::width(); }
	uint32_t width() const { return xptr_t::height(); }

	// Change the width/height/element size of this column to match that
	// of the argument column
	void reshape(const column<T> &t)
	{
		(xptr_t &)(*this) = xptr_t(t.nrows(), t.width(), 1, t.elementSize());
	}

	// Resize the column to have the given number of rows, width and
	// element size
	void resize(uint32_t nrows, uint32_t width = 0, uint32_t es = 0)
	{
		if(!es) { es = this->elementSize(); }
		if(!width) { width = this->width(); }
		(xptr_t&)(*this) = xptr_t(nrows, width, 1, es);
	}

	// Return the base of the host memory containing the column data
	T *get() const { return this->syncToHost(); }

public:
	column() { }

private:
	// prevent copying
	column(const column<T> &);
};

// convenience typedefs
typedef column<double> cdouble_t;
typedef column<int> cint_t;
typedef column<float> cfloat_t;

#endif
