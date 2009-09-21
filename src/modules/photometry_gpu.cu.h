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

//
// Module data passed to device kernel
//
struct ALIGN(16) os_photometry_data
{
	int ncolors, bidx;	// number of colors, bootstrap band index

	uint32_t comp0, comp1;	// component ID range [comp0, comp1) to which this photometry module will be asigning magnitudes

	static const int N_REDDENING = 17; // This number is set by the number of textures out of which the colors are sampled (4x4=16, currently)
	float reddening[N_REDDENING];	// reddening coefficients for the loaded bands (NOTE: hardcoded maximum of 17 bands (16 colors))
};

//
// Device kernel implementation
//
#if !__CUDACC__ && !BUILD_FOR_CPU

	DECLARE_KERNEL(os_photometry_kernel(otable_ks ks, gcfloat_t Am, gcint_t flags, gcfloat_t DM, gcfloat_t Mr, int nabsmag, gcfloat_t mags, gcfloat_t FeH, gcint_t comp));

	// Textures with "isochrones" (note that we're assuming a
	// single-age population here; no actual dependence on age)
	DECLARE_TEXTURE(color0, float4, 2, cudaReadModeElementType);
	DECLARE_TEXTURE(color1, float4, 2, cudaReadModeElementType);
	DECLARE_TEXTURE(color2, float4, 2, cudaReadModeElementType);
	DECLARE_TEXTURE(color3, float4, 2, cudaReadModeElementType);

	// Textures with extrapolation flags for isochrones -- return
	// nonzero if the isochrone at that point is an extrapolation
	DECLARE_TEXTURE(cflags0, float4, 2, cudaReadModeElementType);
	DECLARE_TEXTURE(cflags1, float4, 2, cudaReadModeElementType);
	DECLARE_TEXTURE(cflags2, float4, 2, cudaReadModeElementType);
	DECLARE_TEXTURE(cflags3, float4, 2, cudaReadModeElementType);

#else // #if !__CUDACC__ && !BUILD_FOR_CPU
	
	DEFINE_TEXTURE(color0, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
	DEFINE_TEXTURE(color1, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
	DEFINE_TEXTURE(color2, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
	DEFINE_TEXTURE(color3, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

	DEFINE_TEXTURE(cflags0, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
	DEFINE_TEXTURE(cflags1, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
	DEFINE_TEXTURE(cflags2, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
	DEFINE_TEXTURE(cflags3, float4, 2, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

	__device__ uint shiftFlag(float f, int ncolors)
	{
		uint flag =   (f == 0.f) ? (0.f) : (1U << ncolors >> 1);

		return flag;
	}

	__device__ uint fill(float *&colors, uint &flags, const float4 clr, const float4 f, int &ncolors)
	{
		*colors = clr.x; colors++; flags |= shiftFlag(f.x, ncolors); if(--ncolors == 0) { return flags; }
		*colors = clr.y; colors++; flags |= shiftFlag(f.y, ncolors); if(--ncolors == 0) { return flags; }
		*colors = clr.z; colors++; flags |= shiftFlag(f.z, ncolors); if(--ncolors == 0) { return flags; }
		*colors = clr.w; colors++; flags |= shiftFlag(f.w, ncolors);    --ncolors;        return flags;
	}

	__device__ uint sampleColors(float *colors, float FeH, float Mr, int ncolors)
	{
		float4 clr;
		float4 f;
		uint flags = 0;

		clr = TEX2D(color0, FeH, Mr); f = TEX2D(cflags0, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) { return flags; }
		clr = TEX2D(color1, FeH, Mr); f = TEX2D(cflags1, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) { return flags; }
		clr = TEX2D(color2, FeH, Mr); f = TEX2D(cflags2, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) { return flags; }
		clr = TEX2D(color3, FeH, Mr); f = TEX2D(cflags3, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) { return flags; }

		// We should never reach this point
		return 0xFFFFFFFF;
	}

	__constant__ os_photometry_data os_photometry_params;

	KERNEL(
		ks, 0,
		os_photometry_kernel(otable_ks ks, gcfloat_t Am, gcint_t flags, gcfloat_t DM, gcfloat_t Mr, int nabsmag, gcfloat_t mags, gcfloat_t FeH, gcint_t comp),
		os_photometry_kernel,
		(ks, Am, flags, DM, Mr, nabsmag, mags, FeH, comp)
	)
	{
		os_photometry_data &lt = os_photometry_params;
		float *c = ks.sharedMemory<float>();

		for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
		{
			// should we process this star (based on the galactic model component it belongs to)?
			int cmp = comp(row);
			if(lt.comp0 > cmp || cmp >= lt.comp1) { continue; }

			// construct colors given the absolute magnitude and metallicity
			float fFeH = FeH(row);

			// If this is a multiple system, the following loop will loop through
			// each system component, compute its luminosity in each band, and
			// sum up the luminosities to obtain the total luminosity of the system
			//
			// For a single star it reduces to computing its luminosity in the
			// photometric system's bands.
			for(int syscomp = 0; syscomp < nabsmag; syscomp++)
			{
				float fMr = Mr(row, syscomp);
				if(fMr >= ABSMAG_NOT_PRESENT) { break; }

				// get colors of a star with (fFeH, fMr)
				int flag = sampleColors(c, fFeH, fMr, lt.ncolors);
				if(syscomp) { flag &= flags(row); }
				flags(row) = flag;

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

#endif // #else (!__CUDACC__ && !BUILD_FOR_CPU)

#endif // photometry_gpu_cu_h__
