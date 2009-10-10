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

#ifndef fixedFeH_gpu_cu_h__
#define fixedFeH_gpu_cu_h__

//
// Device kernel implementation
//

#if !__CUDACC__ && !BUILD_FOR_CPU

	DECLARE_KERNEL(os_fixedFeH_kernel(otable_ks ks, float fixedFeH, uint32_t compFirst, uint32_t compLast, cint_t::gpu_t comp, cfloat_t::gpu_t FeH));

#else

	KERNEL(
		ks, 0,
		os_fixedFeH_kernel(otable_ks ks, float fixedFeH, uint32_t compFirst, uint32_t compLast, cint_t::gpu_t comp, cfloat_t::gpu_t FeH),
		os_fixedFeH_kernel,
		(ks, fixedFeH, compFirst, compLast, comp, FeH)
	)
	{
		for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
		{
			int cmp = comp(row);
			if(compFirst > cmp || cmp > compLast) { continue; }

			FeH(row) = fixedFeH;
		}
	}

#endif // (__CUDACC__ || BUILD_FOR_CPU)

#endif // fixedFeH_gpu_cu_h__
