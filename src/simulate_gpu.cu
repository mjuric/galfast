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

#ifndef __CONFIG_H
#define __CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <math.h>

#include <astro/constants.h>

#include "simulate_base.h"
#include "column.h"
#include "gpu2.h"

namespace ct = column_types;
KERNEL(
	ks, 3*4,
	os_FeH_kernel(otable_ks ks, os_FeH_data par, gpu_rng_t rng, ct::cint::gpu_t comp, ct::cfloat::gpu_t XYZ, ct::cfloat::gpu_t FeH),
	os_FeH_kernel,
	(ks, par, rng, comp, XYZ, FeH)
)
{
	rng.load(ks);
	uint32_t tid = threadID();
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		float feh;
		switch(comp[row])
		{
			case 0: // BahcallSoneira_model::THIN:
			case 1: // BahcallSoneira_model::THICK:
			{
				// choose the gaussian to draw from
				float p = rng.uniform()*(par.A[0]+par.A[1]);
				int i = p < par.A[0] ? 0 : 1;

				// calculate mean
				float muD = par.muInf + par.DeltaMu*exp(-fabs(XYZ(row, 2))/par.Hmu);		// Bond et al. A2
				float aZ = muD - 0.067f;

				// draw
				feh = rng.gaussian(par.sigma[i]) + aZ + par.offs[i];
			} break;
			case 2: //BahcallSoneira_model::HALO:
				feh = par.offs[2] + rng.gaussian(par.sigma[2]);
				break;
			default:
				//THROW(ENotImplemented, "We should have never gotten here");
				feh = -9999.f;
				break;
		}
		FeH[row] = feh;
	}
	rng.store(ks);
}


// equgal - Equatorial to Galactic coordinates
using namespace peyton;
typedef double Radians;
static const double angp = ctn::d2r * 192.859508333; //  12h 51m 26.282s (J2000)
//static const double dngp = ctn::d2r * 27.128336111;  // +27d 07' 42.01" (J2000)
static const double l0 = ctn::d2r * 32.932;
static const double ce = 0.88998740217659689; // cos(dngp)
static const double se = 0.45598511375586859; // sin(dngp)

inline __device__ double2 galequ(const double2 lb)
{
	const double cb = cos(lb.y);
	const double sb = sin(lb.y);
	const double cl = cos(lb.x-l0);
	const double sl = sin(lb.x-l0);

//	// TODO: These should be precomputed constants
//	const double ce = cos(dngp);
//	const double se = sin(dngp);

	double2 r;
	r.x = atan2(
			cb*cl,
			sb*ce-cb*se*sl
		) + angp;
	r.y = asin(cb*ce*sl + sb*se);

	while(r.x < 0.) { r.x += ctn::pi2; }
	return r;
}

KERNEL(
	ks, 0,
	os_gal2other_kernel(otable_ks ks, int coordsys, ct::cdouble::gpu_t lb0, ct::cdouble::gpu_t out),
	os_gal2other_kernel,
	(ks, coordsys, lb0, out)
)
{
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		double2 lb, ret;

		// convert to radians
		lb.x = lb0(row, 0) * ctn::d2r;
		lb.y = lb0(row, 1) * ctn::d2r;

		// rotate to output coordinate system
		switch(coordsys)
		{
		case EQU:
			ret = galequ(lb);
			break;
		default:
			ret.x = ret.y = -9999.;
			break;
		}

		// convert to degrees
		out(row, 0) = ret.x / ctn::d2r;
		out(row, 1) = ret.y / ctn::d2r;
	}
}

KERNEL(
	ks, 0,
	os_fixedFeH_kernel(otable_ks ks, float fixedFeH, ct::cfloat::gpu_t FeH),
	os_fixedFeH_kernel,
	(ks, fixedFeH, FeH)
)
{
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		FeH[row] = fixedFeH;
	}
}

#include <vector>

#if BUILD_FOR_CPU
extern __TLS std::vector<tptr<float> > *locuses;
extern __TLS std::vector<tptr<uint> >   *flags;
#else
__TLS std::vector<tptr<float> > *locuses;
__TLS std::vector<tptr<uint> >   *flags;
#endif

#if HAVE_CUDA && !BUILD_FOR_CPU
texture<float4, 2, cudaReadModeElementType> color0;
texture<float4, 2, cudaReadModeElementType> color1;
texture<float4, 2, cudaReadModeElementType> color2;
texture<float4, 2, cudaReadModeElementType> color3;
texture<uint4, 2, cudaReadModeElementType> cflags0;
texture<uint4, 2, cudaReadModeElementType> cflags1;
texture<uint4, 2, cudaReadModeElementType> cflags2;
texture<uint4, 2, cudaReadModeElementType> cflags3;

texture<float4, 2, cudaReadModeElementType> *colorTextures[] = { &color0, &color1, &color2, &color3 };
texture<uint4, 2, cudaReadModeElementType> *cflagsTextures[] = { &cflags0, &cflags1, &cflags2, &cflags3 };

__device__ uint fill(float *&colors, uint &flags, const float4 clr, const uint4 f, int &ncolors)
{
	*colors = clr.x; colors++; flags |= f.x; if(--ncolors == 0) return flags;
	*colors = clr.y; colors++; flags |= f.y; if(--ncolors == 0) return flags;
	*colors = clr.z; colors++; flags |= f.z; if(--ncolors == 0) return flags;
	*colors = clr.w; colors++; flags |= f.w; --ncolors; return flags;
}

__device__ uint sampleColors(float *colors, float FeH, float Mr, int ncolors)
{
	float4 clr;
	uint4 f;
	uint flags;

	clr = tex2D(color0, FeH, Mr); f = tex2D(cflags0, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) return flags;
	clr = tex2D(color1, FeH, Mr); f = tex2D(cflags1, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) return flags;
	clr = tex2D(color2, FeH, Mr); f = tex2D(cflags2, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) return flags;
	clr = tex2D(color3, FeH, Mr); f = tex2D(cflags3, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) return flags;
	return 0xFFFFFFFF;
}

#else
uint sampleColors(float *colors, float FeH, float Mr, int ncolors)
{
	int f = (int)FeH;
	int m = (int)Mr;
//	std::cerr << "fm = " << f << " " << m << "   " << FeH << " " << Mr << "\n";

	uint fl = 0;
	for(int ic=0; ic != ncolors; ic++)
	{
		colors[ic] = (*locuses)[ic](f, m);
		fl |= (*flags)[ic](f, m);
	}
	return fl;
}
#endif

#if HAVE_CUDA && BUILD_FOR_CPU
#include <map>
std::map<std::string, xptrng::tptr<float4> > os_photometry_tex_c;
std::map<std::string, xptrng::tptr<uint4> >  os_photometry_tex_f;
void os_photometry_tex_get(const char *id, xptrng::tptr<float4> &c, xptrng::tptr<uint4> &f)
{
	c = os_photometry_tex_c[id];
	f = os_photometry_tex_f[id];
}
void os_photometry_tex_set(const char *id, xptrng::tptr<float4> &c, xptrng::tptr<uint4> &f)
{
	os_photometry_tex_c[id] = c;
	os_photometry_tex_f[id] = f;
}
#endif

#if !BUILD_FOR_CPU || !HAVE_CUDA

void os_photometry_tex_get(const char *id, xptrng::tptr<float4> &c, xptrng::tptr<uint4> &f);
void os_photometry_tex_set(const char *id, xptrng::tptr<float4> &c, xptrng::tptr<uint4> &f);

void os_photometry_set_isochrones(const char *id, std::vector<tptr<float> > *loc, std::vector<tptr<uint> > *flgs)
{
	locuses = loc;
	flags = flgs;

#if HAVE_CUDA
	if(gpuGetActiveDevice() < 0) { return; }

	size_t width  = (*loc)[0].width();
	size_t height = (*loc)[0].height();

	xptrng::tptr<float4> texc;
	xptrng::tptr<uint4>  texf;
	cudaError err;
	int texid = 0;
	for(int i=0; i != loc->size(); i += 4)
	{
		// Get the pre-built arrays cached across kernel calls
		char idx[50];
		sprintf(idx, "%s%d", id, i);
		os_photometry_tex_get(idx, texc, texf);
		if(!texc || !texf)
		{
			texc = xptrng::tptr<float4>(width, height);
			texf =  xptrng::tptr<uint4>(width, height);

			// Pack the lookups to float4
			for(int y=0; y != height; y++)
			{
				for(int x=0; x != width; x++)
				{
					texc(x, y).x =                     (*loc)[i+0](x, y)    ;
					texc(x, y).y = i+1 < loc->size() ? (*loc)[i+1](x, y) : 0;
					texc(x, y).z = i+2 < loc->size() ? (*loc)[i+2](x, y) : 0;
					texc(x, y).w = i+3 < loc->size() ? (*loc)[i+3](x, y) : 0;

					texf(x, y).x =                      (*flgs)[i+0](x, y)    ;
					texf(x, y).y = i+1 < flgs->size() ? (*flgs)[i+1](x, y) : 0;
					texf(x, y).z = i+2 < flgs->size() ? (*flgs)[i+2](x, y) : 0;
					texf(x, y).w = i+3 < flgs->size() ? (*flgs)[i+3](x, y) : 0;
				}
			}
			os_photometry_tex_set(idx, texc, texf);
		}
		printf("%f %f %f %f\n%f %f %f %f\n",
			texc(317, 28).x,
			texc(317, 28).y,
			texc(317, 28).z,
			texc(317, 28).w,
			(*loc)[i+0](317, 28),
			(*loc)[i+1](317, 28),
			(*loc)[i+2](317, 28),
			(*loc)[i+3](317, 28));

		// Bind isochrone array to texture reference
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);
		//cudaArray* cu_array = gpuMMU.mapToCUDAArray(texc, channelDesc);
		cudaArray* cu_array = texc.getCUDAArray(channelDesc);
		// set texture parameters & bind the array to the texture
		colorTextures[texid]->addressMode[0] = cudaAddressModeClamp;
		colorTextures[texid]->addressMode[1] = cudaAddressModeClamp;
		colorTextures[texid]->filterMode = cudaFilterModeLinear;
		colorTextures[texid]->normalized = false;    // access with normalized texture coordinates
		err = cudaBindTextureToArray( *colorTextures[texid], cu_array, channelDesc);
		CUDA_ASSERT(err);

		// Bind flags array to texture reference
		channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindUnsigned);
		//cu_array = gpuMMU.mapToCUDAArray(texf, channelDesc);
		cu_array = texf.getCUDAArray(channelDesc);
		// set texture parameters & bind the array to the texture
		cflagsTextures[texid]->addressMode[0] = cudaAddressModeClamp;
		cflagsTextures[texid]->addressMode[1] = cudaAddressModeClamp;
		cflagsTextures[texid]->filterMode = cudaFilterModePoint;
		cflagsTextures[texid]->normalized = false;    // access with normalized texture coordinates
		err = cudaBindTextureToArray( *cflagsTextures[texid], cu_array, channelDesc);
		CUDA_ASSERT(err);
	}
#endif
}
#endif

typedef ct::cfloat::gpu_t gcfloat;
typedef ct::cint::gpu_t gcint;

#if !__CUDACC__
#include <iostream>
#endif

KERNEL(
	ks, 0,
	os_photometry_kernel(otable_ks ks, os_photometry_data lt, gcint flags, gcfloat bmag, gcfloat Mr, gcfloat mags, gcfloat FeH),
	os_photometry_kernel,
	(ks, lt, flags, bmag, Mr, mags, FeH)
)
{
	float *c = ks.sharedMemory<float>();
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		// construct colors given the absolute magnitude and metallicity
		float fFeH = FeH[row];
		float fMr = Mr[row];
		float iFeH = (fFeH - lt.FeH0) / lt.dFeH;
		float iMr  = (fMr  -  lt.Mr0) / lt.dMr;
		flags[row] = sampleColors(c, iFeH, iMr, lt.ncolors);

		// convert colors to apparent magnitudes
		float mag0 = bmag[row];
		for(int b = 0; b <= lt.ncolors; b++)
		{
			float mag = mag0;
			if(b < lt.bidx) { for(int i=b;       i != lt.bidx; i++) { mag += c[i]; } }
			if(b > lt.bidx) { for(int i=lt.bidx; i != b;    i++) { mag -= c[i]; } }
			mags(row, b) = mag;
		}

/*		if(row == 30)
		{
			for(int i=0; i != lt.ncolors; i++) { std::cerr << "c[" << i << "]=" << c[i] << "\n"; }
			for(int i=0; i != lt.ncolors+1; i++) { std::cerr << "mags[" << i << "]=" << mags(row, i) << "\n"; }
			exit(0);
		}*/
	}
}
