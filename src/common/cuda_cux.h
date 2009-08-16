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
	inline cux_ptr<T> make_cux_ptr(void *data)
	{
		cux_ptr<T> ptr;
		ptr.ptr = data;
		return ptr;
	}

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

// 	/*TEXTURE_TYPE*/
// 	struct cpuTextureReference : public textureReference
// 	{
// 		void *data;
// 	};
// 
// 	template<class T, int dim = 1, enum cudaTextureReadMode mode = cudaReadModeElementType>
// 	struct cpu_texture : public cpuTextureReference
// 	{
// 		typedef T value_type;
// 
// 		__host__ cpu_texture(int                         norm  = 0,
// 				enum cudaTextureFilterMode  fMode = cudaFilterModePoint,
// 				enum cudaTextureAddressMode aMode = cudaAddressModeClamp)
// 		{
// 			normalized     = norm;
// 			filterMode     = fMode;
// 			addressMode[0] = aMode;
// 			addressMode[1] = aMode;
// 			addressMode[2] = aMode;
// 			channelDesc    = cudaCreateChannelDesc<T>();
// 		}
// 
// 		__host__ cpu_texture(int                          norm,
// 				enum cudaTextureFilterMode   fMode,
// 				enum cudaTextureAddressMode  aMode,
// 				struct cudaChannelFormatDesc desc)
// 		{
// 			normalized     = norm;
// 			filterMode     = fMode;
// 			addressMode[0] = aMode;
// 			addressMode[1] = aMode;
// 			addressMode[2] = aMode;
// 			channelDesc    = desc;
// 		}
// 
// 		value_type sample(float x)
// 		{
// 			assert(0);
// 		}
// 	};

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
	// Smart pointer, inner part, low level implementation.
	//     - only meant to be used from tptr and derivatives
	//     - points to a well-formed block of memory with 2D dimensions m_dim, of elements of size m_elementSize,
	//	 with subsequent rows separated by m_pitch bytes.
	//     - is reference counted, auto de-allocates on all devices upon final release
	//     - can upload the data to global device memory, or CUDA array
	struct ptr_desc
	{
	public:
		array_ptr<char, 4> m_data;	// master copy of the data (can be on host or device, depending on onDevice)
		char *slave;			// slave copy of the data  (can be on host or device, depending on onDevice)
		cudaArray* cuArray;		// CUDA array copy of the data

		bool onDevice;			// true if the "master copy" of the data is on the device
		bool cleanCudaArray;		// true if the last access operation was obtaining a reference to cudaArray

/*		static ptr_desc *null;
*/
		uint32_t m_width;		// logical width of the array (in elements)
		uint32_t m_elementSize;		// size of the array element (bytes)

		int refcnt;			// number of xptrs pointing to this desc

		std::set<textureReference *> boundTextures;	// list of textures bound to the cudaArray copy

		uint32_t memsize() const				// number of bytes allocated
		{
			uint32_t size = m_data.extent[0]*m_data.extent[1]*m_data.extent[2];
			return size;
		}

		ptr_desc(size_t es, size_t pitch, size_t width, size_t height = 1, size_t depth = 1);
		~ptr_desc();

// 		static ptr_desc *getnullptr()
// 		{
// 			if(!null) { null = new ptr_desc(0, 0, 0); }
// 			return null;
// 		}

//		operator bool() const { return m_data.extent[0] == 0 && m_data.extent[1] == 0 && m_data.extent[2] == 0; }
		operator bool() const { return m_data.extent[0] == 0; }

		ptr_desc *addref() { ++refcnt; return this; }
		int release()
		{
			--refcnt;
			if(refcnt == 0) { delete this; return 0; }
			return refcnt;
		}

		// upload/download the data to/from the device
		void *syncTo(bool device);
		void *syncToDevice() { return syncTo(true); }
		void *syncToHost() { return syncTo(false); }

		cudaArray *getCUDAArray(cudaChannelFormatDesc &channelDesc);
		// texture access
		void bind_texture(textureReference &texref);
		void unbind_texture(textureReference &texref);
	private:
		// garbage collection -- release all unused copies of this pointer
		// on devices other than master
		void gc();

		// ensure no shenanigans (no copy constructor & operator)
		ptr_desc(const ptr_desc &);
		ptr_desc& operator=(const ptr_desc &);
	};

	// Pointer to n-D array, GPU interface.
	//     - Used to point to device memory, in kernels
	//     - Must be obtained from tptr<>
	template<typename T, int dim = 2>
	struct gptr : public array_ptr<T, dim>
	{
		__device__ __host__ gptr<T, dim> &operator=(void *)			// allow the setting of pointer to NULL
		{
			this->ptr = NULL;
			return *this;
		}
	};

	// Pointer to 2D array, CPU interface.
	//     - Used to access host memory, on the host
	//     - May be obtained from tptr<>
	template<typename T, int dim = 2>
	struct hptr : public array_ptr<T, dim>
	{
		hptr<T, dim> &operator=(void *)			// allow the setting of pointer to NULL
		{
			this->ptr = NULL;
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

	// Smart pointer, public interface.
	//     - attempts to mimic shared_ptr<> semantics
	//     - works in conjunction with gptr<> to provide copies of data on compute devices
	template<typename T>
	struct tptr
	{
		ptr_desc *desc;

		uint32_t elementSize() const { return desc->m_elementSize; }
		uint32_t width()   const { return desc->m_width; }
		uint32_t height()  const { return desc->m_data.extent[1]; }
		uint32_t depth()   const { return desc->m_data.extent[2]; }
		uint32_t ncols()   const { return width(); }
		uint32_t nrows()   const { return height(); }
		uint32_t pitch()   const { return desc->m_data.extent[0]; }
		size_t size() const
		{
			size_t s = width()*height()*depth();
			return s;
		}

/*		tptr()
		{
			desc = ptr_desc::getnullptr()->addref();
		}*/
		tptr(uint32_t width = 0, uint32_t height = 1, uint32_t depth = 1, int elemSize = sizeof(T), uint32_t align = 128)
		{
			desc = new ptr_desc(elemSize, roundUpModulo(elemSize*width, align), width, height, depth);
		}
		tptr(const tptr& t)
		{
			desc = t.desc->addref();
		}
		tptr &operator =(const tptr &t)
		{
			if(t.desc != desc)
			{
				desc->release();
				desc = t.desc->addref();
			}
			return *this;
		}
		~tptr()
 		{
 			desc->release();
 		}

		// memory allocation
		void realloc(const dim3 &dim, int elemSize = sizeof(T), uint32_t align = 128)
		{
			*this = tptr<T>(dim.x, dim.y, dim.z, elemSize, align);
		}
		tptr<T> clone(bool copyData = false) const
		{
			tptr<T> ret(new ptr_desc(elementSize(), pitch(), width(), height(), depth()));
			if(copyData)
			{
				syncToHost();
				memcpy(ret.desc->m_data.ptr, desc->m_data.ptr, desc->memsize());
			}
			return ret;
		}

		void copyFrom(const T* data)
		{
			if(data == NULL) { return; }

			hptr<T, 3> v = *this;
			for(int i = 0; i != width(); i++)
			{
				for(int j = 0; j != height(); j++)
				{
					for(int k = 0; k != depth(); k++)
					{
						v(i, j, k) = *data;
						data++;
					}
				}
			}
		}

		// comparisons/tests
 		operator bool() const
 		{
 			return *desc;
// 			return desc == ptr_desc::getnullptr();
 		}

		// host data accessors (note: use hptr<> interface if possible)
		T& operator()(const uint32_t x, const uint32_t y, const uint32_t z) const	// 3D accessor (valid only if dim = 3)
		{
			if(desc->onDevice || !desc->m_data.ptr) { syncToHost(); }
			return (*(array_ptr<T, 3> *)(&desc->m_data))(x, y, z);
		}
		T& operator()(const uint32_t x, const uint32_t y) const	// 2D accessor (valid only if dim >= 2)
		{
			if(desc->onDevice || !desc->m_data.ptr) { syncToHost(); }
			return (*(array_ptr<T, 2> *)(&desc->m_data))(x, y);
		}
		T& operator()(const uint32_t x) const			// 1D accessor (valid only if dim >= 2)
		{
			if(desc->onDevice || !desc->m_data.ptr) { syncToHost(); }
			return (*(array_ptr<T, 1> *)(&desc->m_data))(x);
		}

		// texture access
		void bind_texture(textureReference &texref)
		{
			desc->bind_texture(texref);
		}
		void unbind_texture(textureReference &texref)
		{
			desc->unbind_texture(texref);
		}
	public:
		template<int dim>
			operator hptr<T, dim>()
			{
				syncToHost();
				return *(hptr<T, dim> *)(&desc->m_data);
			}

		template<int dim>
			operator gptr<T, dim>()
			{
				syncToDevice();
				return *(gptr<T, dim> *)(&desc->m_data);
			}

	protected:
		// multi-device support. Don't call these directly; use gptr instead.
		T* syncToHost() const { return (T*)const_cast<tptr<T> *>(this)->desc->syncToHost(); }
		T* syncToDevice() { return (T*)desc->syncToDevice(); }

		// protected constructors
		tptr(ptr_desc *d)
		{
			desc = d;
		}
	};

};

#endif // _gpu2_h__
