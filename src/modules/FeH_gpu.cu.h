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

#ifndef FEH_gpu_cu_h__
#define FEH_gpu_cu_h__

//
// Module data passed to device kernel
//
struct os_FeH_data
{
	float A[2], sigma[3], offs[3];
	float Hmu, muInf, DeltaMu;
	int comp_thin, comp_thick, comp_halo;
};

//
// Device kernel implementation
//
#if (__CUDACC__ || BUILD_FOR_CPU)

KERNEL(
	ks, 3*4,
	os_FeH_kernel(
		otable_ks ks, os_FeH_data par, gpu_rng_t rng, 
		cint_t::gpu_t comp,
		cfloat_t::gpu_t XYZ,
		cfloat_t::gpu_t FeH),
	os_FeH_kernel,
	(ks, par, rng, comp, XYZ, FeH)
)
{
	uint32_t tid = threadID();
	rng.load(tid);
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		float feh;
		int component = comp(row);
		if(component == par.comp_thin || component == par.comp_thick)
		{
			// choose the gaussian to draw from
			float p = rng.uniform()*(par.A[0]+par.A[1]);
			int i = p < par.A[0] ? 0 : 1;

			// calculate mean
			float muD = par.muInf + par.DeltaMu*exp(-fabs(XYZ(row, 2))/par.Hmu);		// Bond et al. A2
			float aZ = muD - 0.067f;

			// draw
			feh = rng.gaussian(par.sigma[i]) + aZ + par.offs[i];			
		}
		else if(component == par.comp_halo)
		{
			feh = par.offs[2] + rng.gaussian(par.sigma[2]);
		}
		else
		{
			// do nothing (assuming some other module will set Fe/H for this component)
			feh = -100.f;
		}

		if(feh != -100.f)
		{
			FeH(row) = feh;
		}
	}
	rng.store(threadID());
}

#endif // (__CUDACC__ || BUILD_FOR_CPU)

#endif // FEH_gpu_cu_h__
