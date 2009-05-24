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

#ifndef _cuda_emulation_h__
#define _cuda_emulation_h__

#include <string.h>
#include <assert.h>

// Emulate CUDA types and keywords

#define __device__
#define __host__
#define __constant__

struct float4 { float x, y, z, w; };
struct double2 { double x, y; };
struct uint3 { unsigned int x, y, z; };
struct uint4 { unsigned int x, y, z, w; };

struct dim3
{
	unsigned int x, y, z;
	#if defined(__cplusplus)
	dim3(unsigned int x = 1, unsigned int y = 1, unsigned int z = 1) : x(x), y(y), z(z) {}
	dim3(uint3 v) : x(v.x), y(v.y), z(v.z) {}
	operator uint3(void) { uint3 t; t.x = x; t.y = y; t.z = z; return t; }
	#endif /* __cplusplus */
};

/*
* Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
*
* NOTICE TO USER:
*
* This source code is subject to NVIDIA ownership rights under U.S. and
* international Copyright laws.  Users and possessors of this source code
* are hereby granted a nonexclusive, royalty-free license to use this code
* in individual and commercial software.
*
* NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
* CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
* IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
* REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
* MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
* IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
* OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
* OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
* OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
* OR PERFORMANCE OF THIS SOURCE CODE.
*
* U.S. Government End Users.   This source code is a "commercial item" as
* that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
* "commercial computer  software"  and "commercial computer software
* documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
* and is provided to the U.S. Government only as a commercial end item.
* Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
* 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
* source code with only those rights set forth herein.
*
* Any use of this source code in individual and commercial software must
* include, in the user documentation and internal comments to the code,
* the above Disclaimer and U.S. Government End Users Notice.
*/

enum cudaError
{
cudaSuccess                           =      0,   ///< No errors
cudaErrorMissingConfiguration         =      1,   ///< Missing configuration error
cudaErrorMemoryAllocation             =      2,   ///< Memory allocation error
cudaErrorInitializationError          =      3,   ///< Initialization error
cudaErrorLaunchFailure                =      4,   ///< Launch failure
cudaErrorPriorLaunchFailure           =      5,   ///< Prior launch failure
cudaErrorLaunchTimeout                =      6,   ///< Launch timeout error
cudaErrorLaunchOutOfResources         =      7,   ///< Launch out of resources error
cudaErrorInvalidDeviceFunction        =      8,   ///< Invalid device function
cudaErrorInvalidConfiguration         =      9,   ///< Invalid configuration
cudaErrorInvalidDevice                =     10,   ///< Invalid device
cudaErrorInvalidValue                 =     11,   ///< Invalid value
cudaErrorInvalidPitchValue            =     12,   ///< Invalid pitch value
cudaErrorInvalidSymbol                =     13,   ///< Invalid symbol
cudaErrorMapBufferObjectFailed        =     14,   ///< Map buffer object failed
cudaErrorUnmapBufferObjectFailed      =     15,   ///< Unmap buffer object failed
cudaErrorInvalidHostPointer           =     16,   ///< Invalid host pointer
cudaErrorInvalidDevicePointer         =     17,   ///< Invalid device pointer
cudaErrorInvalidTexture               =     18,   ///< Invalid texture
cudaErrorInvalidTextureBinding        =     19,   ///< Invalid texture binding
cudaErrorInvalidChannelDescriptor     =     20,   ///< Invalid channel descriptor
cudaErrorInvalidMemcpyDirection       =     21,   ///< Invalid memcpy direction
cudaErrorAddressOfConstant            =     22,   ///< Address of constant error
cudaErrorTextureFetchFailed           =     23,   ///< Texture fetch failed
cudaErrorTextureNotBound              =     24,   ///< Texture not bound error
cudaErrorSynchronizationError         =     25,   ///< Synchronization error
cudaErrorInvalidFilterSetting         =     26,   ///< Invalid filter setting
cudaErrorInvalidNormSetting           =     27,   ///< Invalid norm setting
cudaErrorMixedDeviceExecution         =     28,   ///< Mixed device execution
cudaErrorCudartUnloading              =     29,   ///< CUDA runtime unloading
cudaErrorUnknown                      =     30,   ///< Unknown error condition
cudaErrorNotYetImplemented            =     31,   ///< Function not yet implemented
cudaErrorMemoryValueTooLarge          =     32,   ///< Memory value too large
cudaErrorInvalidResourceHandle        =     33,   ///< Invalid resource handle
cudaErrorNotReady                     =     34,   ///< Not ready error
cudaErrorInsufficientDriver           =     35,   ///< CUDA runtime is newer than driver
cudaErrorSetOnActiveProcess           =     36,   ///< Set on active process error
cudaErrorNoDevice                     =     38,   ///< No available CUDA device
cudaErrorStartupFailure               =   0x7f,   ///< Startup failure
cudaErrorApiFailureBase               =  10000    ///< API failure base
};

enum cudaMemcpyKind
{
cudaMemcpyHostToHost          =   0,      ///< Host   -> Host
cudaMemcpyHostToDevice        =   1,      ///< Host   -> Device
cudaMemcpyDeviceToHost        =   2,      ///< Device -> Host
cudaMemcpyDeviceToDevice      =   3       ///< Device -> Device
};

/**
* Channel format kind
*/
/*DEVICE_BUILTIN*/
enum cudaChannelFormatKind
{
cudaChannelFormatKindSigned           =   0,      ///< Signed channel format
cudaChannelFormatKindUnsigned         =   1,      ///< Unsigned channel format
cudaChannelFormatKindFloat            =   2,      ///< Float channel format
cudaChannelFormatKindNone             =   3,      ///< No channel format
};

/**
* CUDA Channel format descriptor
*/
/*DEVICE_BUILTIN*/
struct cudaChannelFormatDesc
{
int                        x; ///< x
int                        y; ///< y
int                        z; ///< z
int                        w; ///< w
enum cudaChannelFormatKind f; ///< Channel format kind
};
/*******************************************************************/
/*                       END NVIDIA CODE                           */
/*******************************************************************/


inline const char *cudaGetErrorString(cudaError err) { return "CUDA Error: CUDA Error when no CUDA used (?!)"; }

inline cudaError cudaMalloc(void** devPtr, size_t count)
{
	*devPtr = malloc(count);
	return cudaSuccess;
}
inline cudaError cudaFree(void* devPtr)
{
	free(devPtr);
	return cudaSuccess;
}
inline cudaError cudaMemcpy(void* dst, const void* src, size_t count, cudaMemcpyKind kind)
{
	memcpy(dst, src, count);
	return cudaSuccess;
}

/**
* CUDA array
*/
/*DEVICE_BUILTIN*/
struct cudaArray;

inline cudaError cudaFreeArray(cudaArray* array ) { assert(0); }
inline cudaError cudaMallocArray(cudaArray** array, const struct cudaChannelFormatDesc* desc, size_t width, size_t height ) { assert(0); }
inline cudaError cudaMemcpy2DToArray(cudaArray* dstArray, size_t dstX, size_t dstY, const void* src, size_t spitch, size_t width, size_t height, enum cudaMemcpyKind kind) { assert(0); }

inline bool cuda_init() { return true; }

#endif // _cuda_emulation_h__
