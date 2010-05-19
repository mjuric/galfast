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

#ifndef gal2other_gpu_cu_h__
#define gal2other_gpu_cu_h__

//
// Module data passed to device kernel, encapsulated in a struct
//

// -- none for this module

//
// Device kernel implementation
//
#if !__CUDACC__ && !BUILD_FOR_CPU

	DECLARE_KERNEL(os_gal2other_kernel(otable_ks ks, int coordsys, cint_t::gpu_t hidden, cdouble_t::gpu_t lb0, cdouble_t::gpu_t out));

#else // #if !__CUDACC__ && !BUILD_FOR_CPU

	inline __device__ double2 galequ(const double2 lb)
	{
		using namespace galequ_constants;	// for ce and se
		using namespace peyton;

		const double cb = cos(lb.y);
		const double sb = sin(lb.y);
		const double cl = cos(lb.x-l0);
		const double sl = sin(lb.x-l0);

		double2 r;
		r.x = atan2(
				cb*cl,
				sb*ce-cb*se*sl
			) + angp;
		r.y = asin(cb*ce*sl + sb*se);

		while(r.x < 0.)        { r.x += ctn::pi2; }
		while(r.x >= ctn::pi2) { r.x -= ctn::pi2; }

		return r;
	}

	KERNEL(
		ks, 0,
		os_gal2other_kernel(otable_ks ks, int coordsys, cint_t::gpu_t hidden, cdouble_t::gpu_t lb0, cdouble_t::gpu_t out),
		os_gal2other_kernel,
		(ks, coordsys, hidden, lb0, out)
	)
	{
		using namespace peyton;	// for mathematical constants

		for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
		{
			if(hidden(row)) { continue; }

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

#endif // #else (!__CUDACC__ && !BUILD_FOR_CPU)

#endif // gal2other_gpu_cu_h__
