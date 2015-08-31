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

#ifndef GaussianFeH_gpu_cu_h__
#define GaussianFeH_gpu_cu_h__

//
// Module data passed to device kernel
//
struct os_GaussianFeH_data
{
	static const int MAXCOMP = 10;

	int ncomp;
	float mean[MAXCOMP];
	float sigma[MAXCOMP];
	float pcum[MAXCOMP];
};

//
// Device kernel implementation
//
#if (__CUDACC__ || BUILD_FOR_CPU)

KERNEL(
	ks, 3*4,
	os_GaussianFeH_kernel(
		otable_ks ks, bit_map applyToComponents, os_GaussianFeH_data par, gpu_rng_t rng,
		cint_t::gpu_t comp, cint_t::gpu_t hidden,
		cfloat_t::gpu_t XYZ,
		cint_t::gpu_t FeH_comp,
		cfloat_t::gpu_t FeH),
	os_GaussianFeH_kernel,
	(ks, applyToComponents, par, rng, comp, hidden, XYZ, FeH_comp, FeH)
)
{
	uint32_t tid = threadID();
	rng.load(tid);
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		if(hidden(row)) { continue; }

		int cmp = comp(row);
		if(!applyToComponents.isset(cmp)) { continue; }

		int at = 0;
		if(par.ncomp != 1)
		{
			float u = rng.uniform();
			while(par.pcum[at] < u) ++at;
		}
		FeH(row) = par.mean[at] + rng.gaussian(par.sigma[at]);
		FeH_comp(row) = at;
	}
	rng.store(tid);
}

#endif // (__CUDACC__ || BUILD_FOR_CPU)

#endif // GaussianFeH_gpu_cu_h__
