/***************************************************************************
 *   Copyright (C) 2004 by Mario Juric                                     *
 *   mjuric@astro.Princeton.EDU                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef _gpu2_h__
#define _gpu2_h__

//#include "gpu.h"
#include <assert.h>
#include <map>
#include <set>
#include <algorithm>

// CUDA API wrappers
struct cuxException
{
	cudaError err;
	cuxException(cudaError err_) : err(err_) {}
	const char *msg() const { return cudaGetErrorString(err); }
};

void cuxErrCheck_impl(cudaError err, const char *fun, const char *file, const int line);
#define cuxErrCheck(expr) \
	cuxErrCheck_impl(expr, __PRETTY_FUNCTION__, __FILE__, __LINE__)

template<typename T>
	T *cuxNew(uint32_t size = 1)
	{
		T *v;
		cuxErrCheck( cudaMalloc((void**)&v, sizeof(T)*size) );
		return v;
	}

template<typename T>
	void cuxFree(T *v)
	{
		cuxErrCheck( cudaFree(v) );
	}

template<typename T, int dim = 1>
	struct array_ptr
	{
		T *ptr;
		uint32_t extent[dim-1];	// extent[0] == width of a row of data, in bytes
					// extent[1] == number of rows of data in a slice of a 3D data cube
					// extent[2] == number of slices in a 3D data sub-cube of a 4D data cube (etc...)

		// Access
		__host__ __device__ T &operator()(const uint32_t x, const uint32_t y, const uint32_t z) const	// 3D accessor (valid only if dim = 3)
		{
			return *((T*)((char*)ptr + y * (extent[0] + z * extent[1])) + x);
		}
		__host__ __device__ T &operator()(const uint32_t x, const uint32_t y) const	// 2D accessor (valid only if dim >= 2)
		{
			return *((T*)((char*)ptr + y * extent[0]) + x);
		}
		__host__ __device__ T &operator()(const uint32_t i) const			// 1D accessor (valid only if dim >= 2)
		{
			return ptr[i];
		}
		__host__ __device__ operator bool() const
		{
			return ptr != NULL;
		}
		__host__ __device__ array_ptr<T, dim> &operator=(void *)			// allow the setting of pointer to NULL
		{
			ptr = NULL;
			return *this;
		}
	};

template<typename T>
	inline array_ptr<T, 2> make_array_ptr(T *data, uint32_t pitch)
	{
		array_ptr<T, 2> ptr;
		ptr.ptr = data;
		ptr.extent[0] = pitch;
		return ptr;
	}

template<typename T>
	inline array_ptr<T, 3> make_array_ptr(T *data, uint32_t pitch, uint32_t ydim)
	{
		array_ptr<T, 3> ptr;
		ptr.ptr = data;
		ptr.extent[0] = pitch;
		ptr.extent[1] = ydim;
		return ptr;
	}

template<typename T, int dim = 1>
	struct cux_ptr : public array_ptr<T, dim>
	{
		void   upload(const T* src, int count = 1)  { if(!this->ptr) { alloc(count); } cuxErrCheck( cudaMemcpy(this->ptr,  src, count*sizeof(T), cudaMemcpyHostToDevice) ); }
		void download(T* dest, int count = 1) const {                            cuxErrCheck( cudaMemcpy(dest, this->ptr, count*sizeof(T), cudaMemcpyDeviceToHost) ); }
		void alloc(int count) { this->ptr = cuxNew<T>(count); }
		void free() { cuxFree(this->ptr); }

		cux_ptr<T> operator =(T *p) { this->ptr = p; return *this; }
		cux_ptr<T> operator =(const cux_ptr<T> &p) { this->ptr = p.ptr; return *this; }
	};

template<typename T>
	void cuxUploadConst(T &dest, const T &source)
	{
		unsigned size = sizeof(source);
		cuxErrCheck( cudaMemcpyToSymbol(dest, &source, size) );
	}

template<typename T>
	void cuxUploadConst(const char *symbol, const T &source)
	{
		unsigned size = sizeof(source);
		//fprintf(stderr, "cuxUploadConst: Uploading %u bytes to symbol %s.\n", size, symbol);
		cuxErrCheck( cudaMemcpyToSymbol(symbol, &source, size) );
	}

////////////////////////////////////////////////////////////////////////////////////////////////
// Simplified texturing support
////////////////////////////////////////////////////////////////////////////////////////////////

struct spline;
struct cuxTextureManager
{
public:
	struct textureParameters
	{
		float x0, inv_dx;
	};

protected:
	textureParameters par;
	const char *parSymbol;
	textureReference &texref;
	cudaArray *texdata;

	spline *cputex;

	textureParameters make_textureParameters(float x0, float inv_dx)
	{
		textureParameters ret = { x0, inv_dx };
		return ret;
	}
public:
	cuxTextureManager(textureReference &tr, const char *tp)
		: texref(tr), parSymbol(tp), texdata(0)
	{
		par.x0 = 0.f; par.inv_dx = 1.f;
	}

	float sample(float x) const;

	void load(const char *fn, int nsamples);
	void construct(double *x, double *y, int ndata, int nsamples);
	void bind();

	void free();
	~cuxTextureManager() { free(); }
protected:
	void set(float *cpu_lf, int lfLen, float M0, float dM);
};

#if !HAVE_CUDA || BUILD_FOR_CPU
	template<typename T>
	typename T::value_type tex1D(T &t, float x)
	{
		return t.sample(x);
	}

	#define DEFINE_TEXTURE(name) \
		extern cuxTextureManager name##Manager; \
		cuxTextureManager &name = name##Manager;
#else
	#define DEFINE_TEXTURE_REFERENCE(name) \
		texture<float, 1, cudaReadModeElementType>       tex_##name(false, cudaFilterModeLinear, cudaAddressModeClamp); \
		struct textureSampler_##name : public cuxTextureManager::textureParameters \
		{ \
			__device__ inline float sample(float x) const \
			{ \
				float val, ix; \
				ix  = (x  -  x0) * inv_dx + 0.5; \
				val = tex1D(tex_##name, ix); \
				return val; \
			} \
		}; \
		__device__ __constant__ textureSampler_##name name

	#define DEFINE_TEXTURE_MANAGER(name) \
		cuxTextureManager name##Manager(tex_##name, #name)

	#define DEFINE_TEXTURE(name) \
		DEFINE_TEXTURE_REFERENCE(name); \
		DEFINE_TEXTURE_MANAGER(name);

	// support for textureManagers as class members
	#define DECLARE_MEMBER_TEXTURE_MANAGER(name) \
		cuxTextureManager name##Manager(tex_##name, #name)

	#define INIT_MEMBER_TEXTURE_MANAGER(name) \
		name##Manager(tex_##name, #name)
#endif

#define DECLARE_TEXTURE(name) \
	extern cuxTextureManager name##Manager;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

namespace xptrng
{
	// Pointer to n-D array in global GPU memory
	//	- Used to obtain a pointer to and access device memory (in kernels)
	//	- Thin veneer over array_ptr, to facilitate obtaining by cast from xptr<>
	//	- Must be obtained via cast from xptr<>
	//	- Remains valid until a hptr<> or a texture_bind is called for the parent xptr<>
	template<typename T, int dim = 2>
	struct gptr : public array_ptr<T, dim>
	{
		__device__ __host__ gptr<T, dim> &operator=(void *)			// allow the setting of pointer to NULL
		{
			this->ptr = NULL;
			return *this;
		}
	};

	// Pointer to n-D array in host memory
	//	- Used to obtain a pointer to and access host memory
	//	- Thin veneer over array_ptr, to facilitate obtaining by cast from xptr<>
	//	- Must be obtained via cast from xptr<>
	//	- Remains valid until a gptr<> or a texture_bind is called for the parent xptr<>
	template<typename T, int dim = 2>
	struct hptr : public array_ptr<T, dim>
	{
		hptr<T, dim> &operator=(void *)			// allow the setting of pointer to NULL
		{
			this->ptr = NULL;
			return *this;
		}
	};

	// Inner part, low level implementation of xptr<T>. See the documentation of xptr<T>
	// for more details.
	//     - only meant to be used from xptr
	//     - points to a well-formed block 3D block of elements of size m_elementSize,
	//	 with byte dimensions in m_data.extent[0-2], and logical width given in
	//	 m_width.
	//     - is reference counted, auto de-allocates on all devices upon final release
	struct xptr_impl_t
	{
	public:
		// data members
		array_ptr<char, 4> m_data;	// master copy of the data (can be on host or device, depending on onDevice)
		char *slave;			// slave copy of the data  (can be on host or device, depending on onDevice)
		cudaArray* cuArray;		// CUDA array copy of the data

		bool onDevice;			// true if the "master copy" of the data is on the device
		bool cleanCudaArray;		// true if the last access operation was obtaining a reference to cudaArray

		uint32_t m_elementSize;		// size of the array element (bytes)
		uint32_t m_width;		// logical width of the array (in elements)

		int refcnt;			// number of xptrs pointing to this m_impl

		std::set<textureReference *> boundTextures;	// list of textures bound to the cuArray

	public:
		// constructors
		xptr_impl_t(size_t es, size_t pitch, size_t width, size_t height = 1, size_t depth = 1);
		~xptr_impl_t();

		// aux methods
		operator bool() const					// true if the pointer is considered non-null (effectively defines what non-null means)
		{
			return m_data.extent[0] != 0;
		}
		uint32_t memsize() const				// number of bytes allocated
		{
			uint32_t size = m_data.extent[0]*m_data.extent[1]*m_data.extent[2];
			return size;
		}

		// reference management
		xptr_impl_t *addref() { ++refcnt; return this; }
		int release()
		{
			--refcnt;
			if(refcnt == 0) { delete this; return 0; }
			return refcnt;
		}

		// upload/download to/from the device
		void *syncTo(bool device);
		void *syncToDevice() { return syncTo(true); }
		void *syncToHost() { return syncTo(false); }

		// texture access
		void bind_texture(textureReference &texref);
		void unbind_texture(textureReference &texref);

	private:
		// garbage collection facilities
		struct allocated_pointers : public std::set<xptr_impl_t *>
		{
			~allocated_pointers();
		};
		static allocated_pointers all_xptrs;
		static void global_gc();

	private:
		cudaArray *getCUDAArray(cudaChannelFormatDesc &channelDesc);

		// garbage collection -- release all unused copies of this pointer
		// on devices other than master
		void gc();

		// ensure no shenanigans (no copy constructor & operator)
		xptr_impl_t(const xptr_impl_t &);
		xptr_impl_t& operator=(const xptr_impl_t &);
	};

	// Smart GPU/CPU memory pointer with on-demand garbage collection
	//	- Points to a sized (up to three-dimensional) block of memory
	//	- At any given time, the block is either on GPU, CPU, or bound to texture
	//	- To access the block on host/device, obtain a hptr<> or gptr<> via cast.
	//	  Obtaining either moves the block to the apropriate device, and invalidates
	//	  any previously obtained pointers of the other type. The pointer is then
	//	  "locked" to that particular device.
	//	- Accessing the data with operator() is equivalent to obtaining a hptr<>
	//	- Blocks remain allocated on device and host (for speed) until an OOM
	//	  condition is reached on the device. Garbage collection will then remove
	//	  all unlocked device pointers, and retry the allocation.
	//	- The block can be bound to texture using bind_texture, and must be unbound
	//	  before the pointer is deallocated
	//	- The pointer is reference-counted, and automatically releases all memory
	//	  when the reference count falls to zero.
	template<typename T>
	struct xptr
	{
	protected:
		xptr_impl_t *m_impl;

	public:
		uint32_t elementSize() const { return m_impl->m_elementSize; }	// size of the stored element (typically, sizeof(T))
		uint32_t width()   const { return m_impl->m_width; }		// logical width of the data array
		uint32_t height()  const { return m_impl->m_data.extent[1]; }	// logical height of the data array
		uint32_t depth()   const { return m_impl->m_data.extent[2]; }	// logical depth of the data array
		uint32_t pitch()   const { return m_impl->m_data.extent[0]; }	// byte width of the data array
		uint32_t size()    const { return width()*height()*depth(); }	// logical size (number of elements) in the data array

		xptr(const xptr& t)				// construct a copy of existing pointer
		{
			m_impl = t.m_impl->addref();
		}
		xptr &operator =(const xptr &t)			// construct a copy of existing pointer
		{
			if(t.m_impl != m_impl)
			{
				m_impl->release();
				m_impl = t.m_impl->addref();
			}
			return *this;
		}
		~xptr()
 		{
 			m_impl->release();
 		}

		xptr<T> clone(bool copyData = false) const	// make the exact copy of this pointer, optionally copying the data as well
		{
			xptr<T> ret(new xptr_impl_t(elementSize(), pitch(), width(), height(), depth()));
			if(copyData)
			{
				syncToHost();
				memcpy(ret.m_impl->m_data.ptr, m_impl->m_data.ptr, m_impl->memsize());
			}
			return ret;
		}

		// comparisons/tests
 		operator bool() const			// return true if this is not a null pointer
 		{
 			return (bool)(*m_impl);
 		}

		// host data accessors (note: use hptr<> interface if possible, for speed)
		T& operator()(const uint32_t x, const uint32_t y, const uint32_t z) const	// 3D accessor (valid only if dim = 3)
		{
			if(m_impl->onDevice || !m_impl->m_data.ptr) { syncToHost(); }
			return (*(array_ptr<T, 3> *)(&m_impl->m_data))(x, y, z);
		}
		T& operator()(const uint32_t x, const uint32_t y) const	// 2D accessor (valid only if dim >= 2)
		{
			if(m_impl->onDevice || !m_impl->m_data.ptr) { syncToHost(); }
			return (*(array_ptr<T, 2> *)(&m_impl->m_data))(x, y);
		}
		T& operator()(const uint32_t x) const			// 1D accessor (valid only if dim >= 2)
		{
			if(m_impl->onDevice || !m_impl->m_data.ptr) { syncToHost(); }
			return (*(array_ptr<T, 1> *)(&m_impl->m_data))(x);
		}

		// texture access
		void bind_texture(textureReference &texref)
		{
			m_impl->bind_texture(texref);
		}
		void unbind_texture(textureReference &texref)
		{
			m_impl->unbind_texture(texref);
		}
	public:
		template<int dim>
			operator hptr<T, dim>()			// request access to data on the host
			{
				syncToHost();
				return *(hptr<T, dim> *)(&m_impl->m_data);
			}

		template<int dim>
			operator gptr<T, dim>()			// request access to data on the GPU
			{
				syncToDevice();
				return *(gptr<T, dim> *)(&m_impl->m_data);
			}

		xptr(uint32_t width = 0, uint32_t height = 1, uint32_t depth = 1, int elemSize = sizeof(T), uint32_t align = 128)
		{
			m_impl = new xptr_impl_t(elemSize, roundUpModulo(elemSize*width, align), width, height, depth);
		}


	protected:
		// multi-device support. Don't call these directly; use gptr instead.
		T* syncToHost() const { return (T*)const_cast<xptr<T> *>(this)->m_impl->syncToHost(); }
		T* syncToDevice() { return (T*)m_impl->syncToDevice(); }

		// protected constructors
		xptr(xptr_impl_t *d)
		{
			m_impl = d;
		}
	};

	template<typename T>
	void copy(xptr<T> &p, const T* data)	// copy data from a C pointer to an xptr<>. Assumes the argument has the required number of elements.
	{
		if(data == NULL) { return; }

		hptr<T, 3> v = p;
		for(int i = 0; i != p.width(); i++)
		{
			for(int j = 0; j != p.height(); j++)
			{
				for(int k = 0; k != p.depth(); k++)
				{
					v(i, j, k) = *data;
					data++;
				}
			}
		}
	}
};
using namespace xptrng;

#endif // _gpu2_h__
