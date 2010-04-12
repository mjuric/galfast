#ifndef __column_h
#define __column_h

#include "gpu.h"
#include <vector>

/**
	bit_map -- A simple fixed-width bit-map abstraction.

	Can/should be thought of as an array of bools with some extra
	syntactic suggar added.

	NOTE: No constructors are declared on the device!
*/
typedef std::vector<std::pair<uint32_t, uint32_t> > interval_list;

#ifndef __CUDACC__
#include <iostream>
#include <astro/macros.h>

inline std::ostream &operator <<(std::ostream &out, const interval_list &il)
{
	bool first = true;
	FOREACH(il)
	{
		if(first)
		{
			out << "{ ";
			first = false;
		}
		else
		{
			out << ", ";
		}

		if(i->first != i->second)
		{
			out << "[" << i->first << "-" << i->second << "]";
		}
		else
		{
			out << i->first;
		}
	}
	out << " }";
	return out;
}

#endif

struct bit_map
{
protected:
	static const int MAXBITS = 64;	// Note: if you change MAXBITS, make sure to change the bounds of bits, and isset/set
	uint32_t bits[2];
public:
	__device__ __host__ int isset(int bit) const
	{
		if(bit < 32) { return (1U << bit) & bits[0]; }
		bit -= 32;
		{ return (1U << bit) & bits[1]; }
	}
	__device__ __host__ void set(int bit, int value = true)
	{
		if(value) {
			// set bit
			if(bit < 32) { bits[0] |= (1 << bit); return; }
			bit -= 32;
			{ bits[1] |= (1 << bit); }
		} else {
			// clear bit
			if(bit < 32) { bits[0] &= ~(1 << bit); return; }
			bit -= 32;
			{ bits[1] &= ~(1 << bit); return; }
		}
	}
public:
#ifndef __CUDACC__
	bit_map(const interval_list &il);
	bit_map() { }
#endif

public:
	__device__ __host__ void set_all(int value = true) { for(int i=0; i != MAXBITS; i++) { set(i, value); } }

	static int maxbits() { return MAXBITS; }

	bool isanyset() { for(int i=0; i != MAXBITS; i++) { if(isset(i)) { return true; } }; return false; }
	bool operator ==(const bit_map &b)
	{
		for(int i = 0; i != MAXBITS; i++)
		{
			if(isset(i) != b.isset(i)) { return false; }
		}
		return true;
	}
	bit_map &operator=(const bit_map &a)
	{
		for(int i = 0; i != MAXBITS; i++)
		{
			set(i, a.isset(i));
		}
		return *this;
	}
	bit_map &operator|=(const bit_map &a)
	{
		for(int i = 0; i != MAXBITS; i++)
		{
			set(i, a.isset(i) || isset(i));
		}
		return *this;
	}
};

/**
	column -- Adaptation of cuxSmartPointer to table column semantics

	Used internally for storage of columns in otable. One shouldn't need
	to use this class directly.
*/
template<typename T>
struct column : public cuxSmartPtr<T>
{
public:
	typedef cuxSmartPtr<T> cuxSmartPtr_t;

	typedef gptr<T, 2> gpu_t;
	typedef hptr<T, 2> host_t;

public:
	// NOTE: WARNING: The meaning of width/height is being _inverted_ here
	// to bring it in line with the typical table metaphore. In memory, however
	// the array is stored such that nrows() is its width. This allows the GPU
	// to coalesce memory accesses to adjacent rows of the table.
	uint32_t nrows() const { return cuxSmartPtr_t::width(); }
	uint32_t width() const { return cuxSmartPtr_t::height(); }

	// Change the width/height/element size of this column to match that
	// of the argument column
	void reshape(const column<T> &t)
	{
		(cuxSmartPtr_t &)(*this) = cuxSmartPtr_t(t.nrows(), t.width(), 1, t.elementSize());
	}

	// Resize the column to have the given number of rows, width and
	// element size
	void resize(uint32_t nrows, uint32_t width = 0, uint32_t es = 0)
	{
		if(!es) { es = this->elementSize(); }
		if(!width) { width = this->width(); }
		(cuxSmartPtr_t&)(*this) = cuxSmartPtr_t(nrows, width, 1, es);
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
typedef column<double>	cdouble_t;
typedef column<int>	cint_t;
typedef column<float>	cfloat_t;

typedef column<double>::gpu_t	gcdouble_t;
typedef column<int>::gpu_t	gcint_t;
typedef column<float>::gpu_t	gcfloat_t;

#endif
