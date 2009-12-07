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

#ifndef photometry_gpu_cu_h__
#define photometry_gpu_cu_h__

#include "photometry.h"

//
// Device kernel implementation
//

DEFINE_TEXTURE(color0, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
DEFINE_TEXTURE(color1, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
DEFINE_TEXTURE(color2, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
DEFINE_TEXTURE(color3, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

DEFINE_TEXTURE(cflags0, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
DEFINE_TEXTURE(cflags1, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
DEFINE_TEXTURE(cflags2, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
DEFINE_TEXTURE(cflags3, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

__device__ uint setFlag(float f, uint cidx)
{
	uint flag = f != 0.f;
	flag <<= cidx;

	return flag;
}

__device__ void fill(float *colors, uint &cidx, uint &flags, const float4 clr, const float4 f, const uint ncolors)
{
	colors[cidx] = clr.x; flags |= setFlag(f.x, cidx); cidx++; if(cidx == ncolors) { return; }
	colors[cidx] = clr.y; flags |= setFlag(f.y, cidx); cidx++; if(cidx == ncolors) { return; }
	colors[cidx] = clr.z; flags |= setFlag(f.z, cidx); cidx++; if(cidx == ncolors) { return; }
	colors[cidx] = clr.w; flags |= setFlag(f.w, cidx); cidx++; 
}

__device__ uint sampleColors(float *colors, float FeH, float Mr, uint ncolors)
{
	float4 clr;
	float4 f;

	uint flags = 0;
	uint cidx = 0;

	clr = TEX2D(color0, FeH, Mr); f = TEX2D(cflags0, FeH, Mr); fill(colors, cidx, flags, clr, f, ncolors); if(cidx == ncolors) { return flags; }
	clr = TEX2D(color1, FeH, Mr); f = TEX2D(cflags1, FeH, Mr); fill(colors, cidx, flags, clr, f, ncolors); if(cidx == ncolors) { return flags; }
	clr = TEX2D(color2, FeH, Mr); f = TEX2D(cflags2, FeH, Mr); fill(colors, cidx, flags, clr, f, ncolors); if(cidx == ncolors) { return flags; }
	clr = TEX2D(color3, FeH, Mr); f = TEX2D(cflags3, FeH, Mr); fill(colors, cidx, flags, clr, f, ncolors); if(cidx == ncolors) { return flags; }

	// We should never reach this point
	return 0xFFFFFFFF;
}

__constant__ os_photometry_data os_photometry_params;

#if BUILD_FOR_CPU
unsigned int __brev(unsigned int a)
{
  a = ((a >>  1) & 0x55555555) + ((a & 0x55555555) <<  1);
  a = ((a >>  2) & 0x33333333) + ((a & 0x33333333) <<  2);
  a = ((a >>  4) & 0x0F0F0F0F) + ((a & 0x0F0F0F0F) <<  4);
  a = ((a >>  8) & 0x00FF00FF) + ((a & 0x00FF00FF) <<  8);
  a = ( a >> 16              ) + ( a               << 16);
  return a;
}
#endif

KERNEL(
	ks, 0,
	os_photometry_kernel(otable_ks ks, bit_map applyToComponents, gcfloat_t Am, gcint_t flags, gcfloat_t DM, gcfloat_t Mr, int nabsmag, gcfloat_t mags, gcfloat_t FeH, gcint_t comp),
	os_photometry_kernel,
	(ks, applyToComponents, Am, flags, DM, Mr, nabsmag, mags, FeH, comp)
)
{
	os_photometry_data &lt = os_photometry_params;
	float *c = ks.sharedMemory<float>();

	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		// should we process this star (based on the galactic model component it belongs to)?
		int cmp = comp(row);
		if(!applyToComponents.isset(cmp)) { continue; }

		// construct colors given the absolute magnitude and metallicity
		float fFeH = FeH(row);

		// If this is a multiple system, the following loop will loop through
		// each system component, compute its luminosity in each band, and
		// sum up the luminosities to obtain the total luminosity of the system
		//
		// For a single star it reduces to computing its luminosity in the
		// photometric system's bands.
		uint flag = 0;
		for(int syscomp = 0; syscomp < nabsmag; syscomp++)
		{
			float fMr = Mr(row, syscomp);
			if(fMr >= ABSMAG_NOT_PRESENT) { break; }

			// get colors of a star with (fFeH, fMr)
			flag |= sampleColors(c, fFeH, fMr, lt.ncolors);

			// compute the luminosities in all bands
			for(int b = 0; b <= lt.ncolors; b++)
			{
				float M = fMr;
				if(b < lt.bidx) { for(int i=b;       i != lt.bidx; i++) { M += c[i]; } }
				if(b > lt.bidx) { for(int i=lt.bidx; i != b;    i++)    { M -= c[i]; } }

				float L = exp10f(-0.4f*M);
				if(syscomp) { L += mags(row, b); }
				mags(row, b) = L;
			}
		}

		// convert the per-color flags to per-band. After this
		// bit ncolors-n is set if mag[n] was obtained by extrapolation
		// (where the least significant bit is n=0)
		bool fl = false;
		for(int i = lt.bidx-1; i >= 0; i--)
		{
			uint mask = 1 << i;
			fl = fl || (flag & mask);	// set fl=true if i-th flag is set
			flag |= fl << i;		// set i-th flag to fl
		}
		fl = false;
		for(int i = lt.bidx; i < lt.ncolors; i++)
		{
			uint mask = 1 << i;
			fl = fl || (flag & mask);	// set fl=true if i-th flag is set
			flag |= fl << i << 1;		// set (i+1)-st flag to fl
		}
		flag &= ~(1 << lt.bidx);			// clear bootstrap band's flag (bootstrapped band is, by definition, not extrapolated)
		flag = __brev(flag);				// reverse the order of the flags
		flag >>= (8*sizeof(flag) - lt.ncolors - 1);
		flags(row) = flag;


		float dm = DM(row);	// distance modulus to the star
		float Am0 = Am(row);	// extinction to the star (in bootstrap band)

		// convert luminosity to apparent magnitude of the system
		// taking extinction and reddening into account
		for(int b = 0; b <= lt.ncolors; b++)
		{
			float Msys = -2.5f * log10f(mags(row, b));
			float msys = dm + Msys;

			float Am = Am0 * lt.reddening[b];

			mags(row, b) = msys + Am;
		}
	}
}

#endif // photometry_gpu_cu_h__
