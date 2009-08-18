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
			return *((T*)((char*)ptr + y * extent[0] + z * extent[0] * extent[1]) + x);
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
#ifdef __CUDACC__
		cuxErrCheck( cudaMemcpyToSymbol(dest, &source, size) );
#elif BUILD_FOR_CPU
		memcpy(&dest, &source, size);
#else
		assert(0);
		//#error cuxUploadConst can be used only in .cu files, or when BUILD_FOR_CPU is defined
#endif
	}

template<typename T>
	void cuxUploadConst(const char *symbol, const T &source)
	{
		unsigned size = sizeof(source);
		//fprintf(stderr, "cuxUploadConst: Uploading %u bytes to symbol %s.\n", size, symbol);
		cuxErrCheck( cudaMemcpyToSymbol(symbol, &source, size) );
	}

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
		// iterator-ish interface (host)
		struct iterator
		{
			xptr<T> *d;
			uint32_t i, j, k;

			iterator(xptr<T> *d_ = NULL, uint32_t i_ = 0, uint32_t j_ = 0, uint32_t k_ = 0)
			: d(d_), i(i_), j(j_), k(k_)
			{}

			T& operator *()  const { return (*d)(i, j, k); }
			T& operator ->() const { return (*d)(i, j, k); }
			iterator &operator ++()	// prefix
			{
				if(++i >= d->width())
				{
					i = 0;
					if(++j >= d->height())
					{
						j = 0;
						++k;
					}
				}
				return *this;
			}
			bool operator==(const iterator &a) const
			{
				return a.i == i && a.j == j && a.k == k;
			}
			bool operator!=(const iterator &a) const
			{
				return !(a == *this);
			}
		};
		iterator begin() { return iterator(this, 0, 0, 0); }
		iterator end() { return iterator(this, 0, 0, depth()); }
	public:
		uint32_t elementSize() const { return m_impl->m_elementSize; }		// size of the stored element (typically, sizeof(T))
		uint32_t width()       const { return m_impl->m_width; }		// logical width of the data array
		uint32_t height()      const { return m_impl->m_data.extent[1]; }	// logical height of the data array
		uint32_t depth()       const { return m_impl->m_data.extent[2]; }	// logical depth of the data array
		uint32_t pitch()       const { return m_impl->m_data.extent[0]; }	// byte width of the data array
		uint32_t size()        const { return width()*height()*depth(); }	// logical size (number of elements) in the data array
		uint32_t extent(int n) const						// logical size (number of elements) of the requested dimension (0=first, 1=second, ...)
			{ return n == 0 ? width() : m_impl->m_data.extent[n]; }

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
				ret.syncToHost();
				memcpy(ret.m_impl->m_data.ptr, m_impl->m_data.ptr, m_impl->memsize());
			}
			return ret;
		}

		// comparisons/tests
 		operator bool() const			// return true if this is not a null pointer
 		{
 			return (bool)(*m_impl);
 		}

		// host data accessors (note: use hptr<> interface if maximum speed is needed)
		T& operator()(const uint32_t x, const uint32_t y, const uint32_t z) const	// 3D accessor (valid only if dim = 3)
		{
			assert(x < width());
			assert(y < height());
			assert(z < depth());
			if(m_impl->onDevice || !m_impl->m_data.ptr) { syncToHost(); }
			return (*(array_ptr<T, 3> *)(&m_impl->m_data))(x, y, z);
		}
		T& operator()(const uint32_t x, const uint32_t y) const	// 2D accessor (valid only if dim >= 2)
		{
			assert(x < width());
			assert(y < height());
			if(m_impl->onDevice || !m_impl->m_data.ptr) { syncToHost(); }
			return (*(array_ptr<T, 2> *)(&m_impl->m_data))(x, y);
		}
		T& operator()(const uint32_t x) const			// 1D accessor (valid only if dim >= 2)
		{
			assert(x < width());
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

		xptr(uint32_t width = 0, uint32_t height = 1, uint32_t depth = 1, uint32_t elemSize = 0xffffffff, uint32_t align = 128)
		{
			if(elemSize == 0xffffffff) { elemSize = sizeof(T); }
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

////////////////////////////////////////////////////////////////////////////////////////////////
// Simplified texturing support
////////////////////////////////////////////////////////////////////////////////////////////////

/*
	Work around CUDA defficiency with some built-in struct alignments.

	CUDA header files declare some structs (float2 being an example) with
	__builtin_align() attribute that resolves to __align__ only when using
	CUDACC. This makes those structure's alignments different in nvcc compiled
	object code, compared to GCC's. Example: with nvcc, float2 is 8-byte
	aligned; on gcc, it's 4-byte aligned (given its members are all 4-byte
	aligned floats). Therefore, a struct that has a float2 member may be
	packed differently on gcc and nvcc. Example: struct { float a; float2 b; };
	On nvcc, &b = &a + 8 (in bytes). On gcc, &b = &a + 4 (bytes).

	This cludge works around the problem by deriving an aligned type from
	the problematic CUDA type. It should be used instead of the CUDA type
	in structures where this problem may occur.
*/
struct ALIGN(8) afloat2 : public float2
{
	afloat2& operator =(const float2 &a) { (float2&)*this = a; return *this; }
};

// CPU emulation & management
template<int dim>
	struct cuxTexCoords
	{
		afloat2 tc[dim];
		__host__ const float2 &operator[](int i) const { return tc[i]; }
		__host__       float2 &operator[](int i)       { return tc[i]; }
	};

struct cuxTextureInterface
{
	virtual void bind(const void *data_, const float2 *texcoord) = 0;
	virtual void unbind() = 0;
};

struct cuxTextureBinder
{
	cuxTextureInterface &tex;

	template<typename T>
	cuxTextureBinder(cuxTextureInterface &tex_, const xptr<T> &data, const float2 *texcoord)
		: tex(tex_)
	{
		tex.bind(&data, texcoord);
	}

	~cuxTextureBinder()
	{
		tex.unbind();
	}
};


template<typename T, int dim, enum cudaTextureReadMode mode>
	struct cuxTexture : public cuxTextureInterface
	{
	public:
		xptr<T>	data;
		textureReference &texref;
		cuxTexCoords<dim> tc;
		const char *tcSymbolName;
	public:
		cuxTexture(textureReference &texref_, const char *tcSymbolName)
			: texref(texref_), tcSymbolName(tcSymbolName)
		{
		}

		virtual void bind(const void *xptr_data_, const float2 *texcoord)
		{
			bind(*(const xptr<T> *)xptr_data_, texcoord);
		}

		void bind(const xptr<T> &data_, const float2 *texcoord)
		{
			unbind();

			data = data_;
			for(int i=0; i != dim; i++) { tc[i] = texcoord[i]; }

			cuxUploadConst(tcSymbolName, tc);
			data.bind_texture(texref);
		}

		virtual void unbind()
		{
			if(!data) { return; }

			data.unbind_texture(texref);
			data = 0U;
		}

		T tex1D(float x) const		// sample a 1D texture
		{
			// FIXME: implement interpolation, clamp modes, normalized coordinates
			uint32_t i = (uint32_t)x;
			if(i < 0) { i = 0; }
			if(i >= data.width()) { i = data.width()-1; }
			return data(i);
		}
	};

template<typename T, int dim, enum cudaTextureReadMode mode>
	T tex1D(const cuxTexture<T, dim, mode> &texref, float x)
	{
		return texref.tex1D(x);
	}

// Sampler routines
template<typename T, typename Texref>
	inline __device__ T sample_impl(Texref r, float x, float2 tc)
	{
		float xi = (x - tc.x) * tc.y + 0.5f;
		T v = tex1D(r, xi);
#if __DEVICE_EMULATION__
//		printf("phi=%f\n", v);
#endif
		return v;
	}

template<typename T, typename Texref>
	inline __device__ T sample_impl(Texref r, float x, float y, float2 tcx, float2 tcy)
	{
		float xi = (x - tcx.x) * tcx.y + 0.5f;
		float yi = (y - tcy.x) * tcy.y + 0.5f;
		T v = tex2D(r, xi, yi);
		return v;
	}

template<typename T, typename Texref>
	inline __device__ T sample_impl(Texref r, float x, float y, float z, float2 tcx, float2 tcy, float2 tcz)
	{
		float xi = (x - tcx.x) * tcx.y + 0.5f;
		float yi = (y - tcy.x) * tcy.y + 0.5f;
		float zi = (z - tcz.x) * tcz.y + 0.5f;
		T v = tex3D(r, xi, yi, zi);
		return v;
	}

// texture class (GPU and CPU)
#if !__CUDACC__
	#define DEFINE_TEXTURE(name, T, dim, mode, norm, fMode, aMode) \
		extern cuxTexture<T, dim, mode> name##Manager; \
		static cuxTexture<T, dim, mode> &name = name##Manager

	#define DECLARE_TEXTURE(name, T, dim, mode) \
		DEFINE_TEXTURE(name, T, dim, mode, dummy, dummy, dummy)

	#define TEX1D(name, x)       sample(name, x)
	#define TEX2D(name, x, y)    sample(name, x, y)
	#define TEX3D(name, x, y, z) sample(name, x, y, z)

	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(cuxTexture<T, 1, mode> r, float x)
		{
			return sample_impl<T, cuxTexture<T, 1, mode> >(r, x, r.tc[0]);
		}

	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(cuxTexture<T, 2, mode> r, float x, float y)
		{
			return sample_impl<T, cuxTexture<T, 2, mode> >(r, x, y, r.tc[0]);
		}

	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(cuxTexture<T, 3, mode> r, float x, float y, float z)
		{
			return sample_impl<T, cuxTexture<T, 3, mode> >(r, x, y, z, r.tc[0]);
		}
#else
	// real thing
	#define DEFINE_TEXTURE(name, T, dim, mode, norm, fMode, aMode) \
		texture<T, dim, mode> name(norm, fMode, aMode); \
		__constant__ cuxTexCoords<dim> name##TC; \
		cuxTexture<T, dim, mode> name##Manager(name, #name "TC")

	#define TEX1D(name, x)       sample(name, x, name##TC)
	#define TEX2D(name, x, y)    sample(name, x, y, name##TC)
	#define TEX3D(name, x, y, z) sample(name, x, y, z, name##TC)

	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 1, mode> r, float x, float2 tc)
		{
			return sample_impl<T, texture<T, 1, mode> >(r, x, tc);
		}
	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 2, mode> r, float x, float y, float2 tcx, float2 tcy)
		{
			return sample_impl<T, texture<T, 2, mode> >(r, x, y, tcx, tcy);
		}
	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 3, mode> r, float x, float y, float z, float2 tcx, float2 tcy, float2 tcz)
		{
			return sample_impl<T, texture<T, 3, mode> >(r, x, y, z, tcx, tcy, tcz);
		}


	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 1, mode> r, float x, cuxTexCoords<1> tc)
		{
			return sample(r, x, tc[0]);
		}
	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 2, mode> r, float x, float y, float z, cuxTexCoords<2> tc)
		{
			return sample(r, x, y, tc[0], tc[1]);
		}
	template<typename T, enum cudaTextureReadMode mode>
		inline __device__ T sample(texture<T, 3, mode> r, float x, float y, float z, cuxTexCoords<3> tc)
		{
			return sample(r, x, y, z, tc[0], tc[1], tc[2]);
		}
#endif

#endif // _gpu2_h__
