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

#ifndef photometry_host_h__
#define photometry_host_h__ 1

//
// Module data passed to device kernel
//
struct ALIGN(16) os_photometry_data
{
	uint32_t ncolors, bidx;	// number of colors, bootstrap band index

	static const int N_REDDENING = 17; // This number is set by the number of textures out of which the colors are sampled (4x4=16, currently)
	float reddening[N_REDDENING];	// reddening coefficients for the loaded bands (NOTE: hardcoded maximum of 17 bands (16 colors))
};

#if !__CUDACC__ && !BUILD_FOR_CPU

	DECLARE_KERNEL(os_photometry_kernel(otable_ks ks, bit_map applyToComponents, gcfloat_t Am, gcint_t flags, gcfloat_t DM, gcfloat_t Mr, int nabsmag, gcfloat_t mags, gcfloat_t FeH, gcint_t comp, gcint_t hidden));

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

#endif // !__CUDACC__ && !BUILD_FOR_CPU

#endif // photometry_host_h__
