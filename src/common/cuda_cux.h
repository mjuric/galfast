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

	template<typename T>
	struct cux_ptr
	{
		T *ptr;

		void   upload(const T* src, int count = 1)  { if(!ptr) { alloc(count); } cuxErrCheck( cudaMemcpy(ptr,  src, count*sizeof(T), cudaMemcpyHostToDevice) ); }
		void download(T* dest, int count = 1) const {                            cuxErrCheck( cudaMemcpy(dest, ptr, count*sizeof(T), cudaMemcpyDeviceToHost) ); }
		void alloc(int count) { ptr = cuxNew<T>(count); }
		void free() { cuxFree(ptr); }
		cux_ptr<T> operator =(T *p) { ptr = p; return *this; }
		cux_ptr<T> operator =(const cux_ptr<T> &p) { ptr = p.ptr; return *this; }
		
		#if __CUDACC__
		T &operator[](uint32_t idx) const { return ptr[idx]; }
		#endif
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
	#define DEFINE_TEXTURE(name) \
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
		__device__ __constant__ textureSampler_##name name; \
		cuxTextureManager name##Manager(tex_##name, #name);
#endif

#define DECLARE_TEXTURE(name) \
	extern cuxTextureManager name##Manager;

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

namespace xptrng
{
	// Pointer to 2D array, GPU interface.
	//     - Used to point to device memory, in kernels
	//     - May be obtained from tptr<>
	template<typename T>
	struct gptr
	{
		typedef T value_type;

		union {
			char *data;
			T *tdata;
		};
		uint32_t pitch;

		// Access
		__device__ T &operator()(const size_t x, const size_t y) const	// 2D accessor
		{
			return *((T*)(data + y * pitch) + x);
		}
		__device__ T &operator[](const size_t i) const			// 1D accessor
		{
			return ((T*)data)[i];
		}
	};
	template<typename T>
	inline gptr<T> make_gptr(void *data, uint32_t pitch)
	{
		gptr<T> ptr;
		ptr.data = (char*)data;
		ptr.pitch = pitch;
		return ptr;
	}

	// Pointer to 2D array, CPU interface.
	//     - Used to access host memory, on the host
	//     - May be obtained from tptr<>
	template<typename T>
	struct hptr : public gptr<T>
	{
	};
	template<typename T>
	inline hptr<T> make_hptr(void *data, uint32_t pitch)
	{
		hptr<T> ptr;
		ptr.data = (char*)data;
		ptr.pitch = pitch;
		return ptr;
	}

	// Smart pointer, inner part, low level implementation.
	//     - only meant to be used from tptr and derivatives
	//     - points to a well-formed block of memory with 2D dimensions m_dim, of elements of size m_elementSize,
	//	 with subsequent rows separated by m_pitch bytes.
	//     - is reference counted, auto de-allocates on all devices upon final release
	//     - can upload the data to global device memory, or CUDA array
	struct ptr_desc
	{
		friend struct GPUMM;
		static ptr_desc *null;

		const uint32_t m_elementSize;	// size of array element (bytes)
		const dim3     m_dim;		// array dimensions (in elements). dim[0] == ncolumns == width, dim[1] == nrows == height
		const uint32_t m_pitch;		// array pitch (in bytes). pitch[1] == width of the (padded) row (in bytes)

		int refcnt;		// number of xptrs pointing to this desc
		int masterDevice;	// the device holding the master copy of the data (-1 is the host)

		void *m_data;					// pointer to actual host-allocated data
		std::map<int, void*> deviceDataPointers;	// map of device<->device pointer
		std::map<int, cudaArray*> cudaArrayPointers;	// map of device<->device pointer

		uint32_t memsize() const { return m_dim.y * m_pitch; }	// number of bytes allocated

		ptr_desc(size_t es, size_t width, size_t height, size_t pitch);
		~ptr_desc();

		static ptr_desc *getnullptr()
		{
			if(!null) { null = new ptr_desc(0, 0, 0, 0); }
			return null;
		}

		ptr_desc *addref() { ++refcnt; return this; }
		int release()
		{
			--refcnt;
			if(refcnt == 0) { delete this; return 0; }
			return refcnt;
		}

		// multi-device support. Returns the pointer to the copy of the data
		// on the device dev, which becomes the master device for the data.
		void *syncToDevice(int dev = -2);

		cudaArray *getCUDAArray(cudaChannelFormatDesc &channelDesc, int dev = -2, bool forceUpload = false);

	private:
		// ensure no shenanigans (no copy constructor & operator)
		ptr_desc(const ptr_desc &);
		ptr_desc& operator=(const ptr_desc &);
	};

	// Smart pointer, public interface.
	//     - attempts to mimic shared_ptr<> semantics
	//     - works in conjunction with gptr<> to provide copies of data on compute devices
	template<typename T>
	struct tptr
	{
		ptr_desc *desc;

		uint32_t elementSize() const { return desc->m_elementSize; }
		uint32_t width()   const { return desc->m_dim.x; }
		uint32_t height()  const { return desc->m_dim.y; }
		uint32_t ncols()   const { return width(); }
		uint32_t nrows()   const { return height(); }
		uint32_t pitch()   const { return desc->m_pitch; }
		size_t size() const { return width()*height(); }

		tptr()
		{
			desc = ptr_desc::getnullptr()->addref();
		}
		tptr(uint32_t width, uint32_t height, int elemSize = sizeof(T), uint32_t align = 128)
		{
			desc = new ptr_desc(elemSize, width, height, roundUpModulo(elemSize*width, align));
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
		void realloc(uint32_t width, uint32_t height, int elemSize = sizeof(T), uint32_t align = 128)
		{
			*this = tptr<T>(width, height, elemSize, align);
		}
		tptr<T> clone(bool copyData = false) const
		{
			tptr<T> ret(new ptr_desc(elementSize(), width(), height(), pitch()));
			if(copyData)
			{
				syncToHost();
				memcpy(ret.desc->m_data, desc->m_data, desc->memsize());
			}
			return ret;
		}

		// comparisons/tests
 		bool isNull() const
 		{
 			return desc == ptr_desc::getnullptr();
 		}

		// debugging
		void assert_synced(int whichDev)
		{
			assert(desc->masterDevice == whichDev);
		}

		// host data accessors -- deprecated in favor of hptr<> interface
		T &elem(const size_t x, const size_t y)	// 2D accessor
		{
			assert_synced(-1);
			size_t offs = y * pitch() + elementSize()*x;
			assert(offs < desc->memsize());
			return *(T*)((char*)desc->m_data + offs);
		}
		T &elem(const size_t i)			// 1D accessor
		{
			assert_synced(-1);
			assert(i > 0 && i < size());
			return ((T*)desc->m_data)[i];
		}

		cudaArray *getCUDAArray(cudaChannelFormatDesc &channelDesc, int dev = -2, bool forceUpload = false)
		{
			return desc->getCUDAArray(channelDesc, dev, forceUpload);
		}

	public:
		operator hptr<T>()
		{
			return make_hptr<T>(syncToHost(), pitch());
		}
		operator gptr<T>()
		{
			return make_gptr<T>(syncToDevice(), pitch());
		}
		operator cux_ptr<T>()
		{
			return make_cux_ptr<T>(syncToDevice());
		}

	protected:
		// multi-device support. Don't call these directly; use gptr instead.
		T *syncToHost() const { return (T*)const_cast<tptr<T> *>(this)->desc->syncToDevice(-1); }
		T* syncToDevice(int dev = -2) { return (T*)desc->syncToDevice(dev); }

		// protected constructors
		tptr(ptr_desc *d)
		{
			desc = d;
		}
	};

};

#endif // _gpu2_h__
