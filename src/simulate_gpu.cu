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

#include <stdint.h>
#include "simulate_base.h"
#include "column.h"
#include "gpu.h"

namespace ct = column_types;
KERNEL(
	ks,
	os_FeH_kernel(otable_ks ks, os_FeH_data par, gpu_rng_t rng, ct::cint comp, ct::cfloat XYZ, ct::cfloat FeH),
	os_FeH_kernel,
	(ks, par, rng, comp, XYZ, FeH)
)
{
	uint32_t row = ks.row();
	if(row == (uint32_t)(-1)) { return; }
	rng.load(ks);

	switch(comp.val(row))
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
			FeH[row] = rng.gaussian(par.sigma[i]) + aZ + par.offs[i];
		} break;
		case 2: //BahcallSoneira_model::HALO:
			FeH[row] = par.offs[2] + rng.gaussian(par.sigma[2]);
			break;
		default:
			//THROW(ENotImplemented, "We should have never gotten here");
			FeH[row] = -9999.f;
			break;
	}
}

