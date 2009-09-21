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

#ifndef unresolvedMultiples_gpu_cu_h__
#define unresolvedMultiples_gpu_cu_h__

namespace multiplesAlgorithms
{
	enum algo
	{
		LF_M2_GT_M1	= 1,
		LF		= 2,
		EQUAL_MASS	= 3
	};
}

//
// Device kernel implementation
//

#if !__CUDACC__ && !BUILD_FOR_CPU

	DECLARE_KERNEL(os_unresolvedMultiples_kernel(otable_ks ks, gpu_rng_t rng, int nabsmag, cfloat_t::gpu_t M, cfloat_t::gpu_t Msys, cint_t::gpu_t ncomp, cint_t::gpu_t comp, uint32_t comp0, uint32_t comp1, multiplesAlgorithms::algo algo));

#else

	DEFINE_TEXTURE( secProb, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
	DEFINE_TEXTURE(   cumLF, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);
	DEFINE_TEXTURE(invCumLF, float, 1, cudaReadModeElementType, false, cudaFilterModeLinear, cudaAddressModeClamp);

	__device__ bool draw_companion(float &M2, float M1, multiplesAlgorithms::algo algo, gpu_rng_t &rng)
	{
		// draw the probability that this star has a secondary
		float psec, u;

		psec = TEX1D(secProb, M1);
		u = rng.uniform();
		if(u > psec) { return false; }

		// draw the absolute magnitude of the secondary, subject to requested
		// algorithm
		using namespace multiplesAlgorithms;
		if(algo == EQUAL_MASS) { M2 = M1; return true; }

		float pprim = TEX1D(cumLF, M1);
		u = rng.uniform();
		if(algo == LF_M2_GT_M1)
		{
			// draw subject to the requirement that it is fainter than the primary
			u = pprim + u * (1. - pprim);
		}
		M2 = TEX1D(invCumLF, u);// + rng.gaussian(1.f);
		if(algo == LF_M2_GT_M1 && M2 < M1)
		{
			// This can happen due to resampling of cumLF and invCumLF
			// (see the note in os_unresolvedMultiples::construct)
			M2 = M1;
		}

		return true;
	}

	KERNEL(
		ks, 3*4,
		os_unresolvedMultiples_kernel(otable_ks ks, gpu_rng_t rng, int nabsmag, cfloat_t::gpu_t M, cfloat_t::gpu_t Msys, cint_t::gpu_t ncomp, cint_t::gpu_t comp, uint32_t comp0, uint32_t comp1, multiplesAlgorithms::algo algo),
		os_unresolvedMultiples_kernel,
		(ks, rng, nabsmag, M, Msys, ncomp, comp, comp0, comp1, algo)
	)
	{
		/*
			Input:	M -- absolute magnitude of the primary.
			Output:	Msys[] -- absolute magnitudes of system components.
				M      -- total absolute magnitude of the system
		*/
		rng.load(threadID());
		for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
		{
			int cmp = comp(row);
			if(comp0 > cmp || cmp >= comp1) { continue; }

			float M1 = M(row);
			Msys(row, 0) = M1;
			float Ltot = exp10f(-0.4f*M1);
			int ncomps = 1;			/* Number of components of the multiple system */
			for(int i = 1; i < nabsmag; i++)
			{
				float M2;
				if(draw_companion(M2, M1, algo, rng))
				{
					Msys(row, i) = M2;
					Ltot += exp10f(-0.4f*M2);
					ncomps++;
				}
				else
				{
					Msys(row, i) = ABSMAG_NOT_PRESENT;
				}
			}
			float Mtot = -2.5*log10f(Ltot);
			M(row) = Mtot;
			ncomp(row) = ncomps;
		}
		rng.store(threadID());
	}


#endif // (__CUDACC__ || BUILD_FOR_CPU)

#endif // unresolvedMultiples_gpu_cu_h__
