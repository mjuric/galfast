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

#include "gpu.h"
#include <astro/macros.h>
#include <assert.h>

namespace xptrng
{
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

		ptr_desc(size_t es, size_t width, size_t height, size_t pitch)
			: m_elementSize(es), m_dim(width, height, 1), m_pitch(pitch), m_data(NULL)
		{
			refcnt = 1;
			masterDevice = -1;
			deviceDataPointers[-1] = m_data = new char[memsize()];
		}
		~ptr_desc()
		{
			delete [] ((char*)m_data);
			FOREACH(deviceDataPointers)
			{
				if(i->first < 0) { continue; }

				cudaError err = cudaFree(i->second);
				abort_on_cuda_error(err);
			}
			FOREACH(cudaArrayPointers)
			{
				cudaError err = cudaFreeArray(i->second);
				abort_on_cuda_error(err);
			}
		}

		static ptr_desc *getnullptr()
		{
			if(!null) { null = new ptr_desc(0, 0, 0, 0); }
			return null;
		}

		ptr_desc *addref() { ++refcnt; return this; }
		int release()
		{
			if(--refcnt == 0) { delete this; return 0; }
			return refcnt;
		}

		// multi-device support. Returns the pointer to the copy of the data
		// on the device dev, which becomes the master device for the data.
		void *syncToDevice(int dev = -2)
		{
			if(dev == -2)
			{
				dev = gpuGetActiveDevice();
			}

			if(masterDevice != dev)
			{
				// check if this device-to-device copy. If so, do the copy via host
				if(masterDevice != -1 && dev != -1)
				{
					syncToDevice(-1); // this will change masterDevice to -1
				}

				// allocate/copy to device
				cudaError err;

				// determine destination and copy direction
				cudaMemcpyKind dir = cudaMemcpyDeviceToHost;
				if(dev != -1)
				{
					dir = cudaMemcpyHostToDevice;

					// allocate device space (if unallocated)
					if(!deviceDataPointers.count(dev))
					{
						err = cudaMalloc(&deviceDataPointers[dev], memsize());
						abort_on_cuda_error(err);
					}
				}
				void *dest = deviceDataPointers[dev];
				void *src = deviceDataPointers[masterDevice];

				// do the copy
				err = cudaMemcpy(dest, src, memsize(), dir);
				abort_on_cuda_error(err);

				// record new master device
				masterDevice = dev;
			}
			return deviceDataPointers[masterDevice];
		}

		cudaArray *getCUDAArray(cudaChannelFormatDesc &channelDesc, int dev = -2, bool forceUpload = false)
		{
			syncToDevice(-1);	// ensure the data is on the host

			if(dev == -2)
			{
				dev = gpuGetActiveDevice();
				assert(dev >= 0);
			}
			
			cudaArray* cu_array;
			cudaError err;

			if(cudaArrayPointers.count(dev))
			{
				cu_array = cudaArrayPointers[dev];
			}
			else
			{
				// autocreate
				err = cudaMallocArray(&cu_array, &channelDesc, m_dim.x, m_dim.y);
				CUDA_ASSERT(err);

				cudaArrayPointers[dev] = cu_array;
				forceUpload = true;
			}

			if(forceUpload)
			{
				err = cudaMemcpy2DToArray(cu_array, 0, 0, m_data, m_pitch, m_dim.x*m_elementSize, m_dim.y, cudaMemcpyHostToDevice);
				CUDA_ASSERT(err);
			}

			return cu_array;
		}

		private:
		// ensure no shenanigans (no copy constructor & operator)
		ptr_desc(const ptr_desc &);
		ptr_desc& operator=(const ptr_desc &);
	};

	template<typename T>
	struct tptr
	{
		ptr_desc *desc;

		uint32_t elementSize() const { return sizeof(T); }
		uint32_t width()   const { return desc->m_dim.x; }
		uint32_t height()  const { return desc->m_dim.y; }
		uint32_t ncols()   const { return width(); }
		uint32_t nrows()   const { return height(); }
		uint32_t pitch()   const { return desc->m_pitch; }

		tptr()
		{
			desc = ptr_desc::getnullptr()->addref();
		}
		tptr(uint32_t width, uint32_t height, uint32_t align = 128)
		{
			desc = new ptr_desc(sizeof(T), width, height, roundUpModulo(sizeof(T)*width, align));
		}
		tptr(const tptr& t)
		{
			desc = t.desc->addref();
		}
		tptr &operator =(const tptr &t)
		{
			desc->release();
			desc = t.desc->addref();
			return *this;
		}
		~tptr()
 		{
 			desc->release();
 		}

		// comparisons/tests
		operator bool() const
		{
			return desc != ptr_desc::getnullptr();
		}

		// debugging
		void assert_synced(int whichDev)
		{
			assert(desc->masterDevice == whichDev);
		}

		// accessors
		T &operator()(const size_t x, const size_t y)	// 2D accessor
		{
			assert_synced(-1);
			return *((T*)((char*)desc->m_data + y * pitch()) + x);
		}
		T &operator[](const size_t i)			// 1D accessor
		{
			assert_synced(-1);
			return ((T*)desc->m_data)[i];
		}

		cudaArray *getCUDAArray(cudaChannelFormatDesc &channelDesc, int dev = -2, bool forceUpload = false)
		{
			return desc->getCUDAArray(channelDesc, dev, forceUpload);
		}

		// multi-device support. Don't call these directly; use gptr instead.
		void syncToHost() { desc->syncToDevice(-1); }
		void* syncToDevice(int dev) { return desc->syncToDevice(dev); }
	};

	template<typename T>
	struct gptr
	{
		char *data;
		uint32_t pitch;

		gptr(tptr<T> &ptr)
		{
			pitch = ptr.pitch();
			data = ptr.syncToDevice();
		}

		// Access
		__device__ T &operator()(const size_t x, const size_t y)	// 2D accessor
		{
			return *((T*)(data + y * pitch) + x);
		}
		__device__ T &operator[](const size_t i)			// 1D accessor
		{
			return ((T*)data)[i];
		}
	};

};
