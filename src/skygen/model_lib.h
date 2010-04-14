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

#ifndef model_lib_h__
#define model_lib_h__

#include "../modules/module_lib.h"

#include <astro/system/config.h>
#include <astro/math.h>

//
// Loads the luminosity function and computes and sets the normalization based
// on the coordinates where the LF was sampled.
//
// LFSeparableModel concept requirements:
//	*) \rho(x,y,z,M) density must be separable and computed as
//	   \rho(x,y,z) * LF(M)
//      *) \rho(x,y,z) must exist
//      *) A non-const member LFSeparableModel::f must exist and the
//         output of \rho(x,y,z) must scale linearly with f. This variable
//         will be rescaled so that \rho(x,y,z) = f_user (where f_user is
//         the "f" value set by the user in the config file).
//
template<typename LFSeparableModel>
cuxTexture<float> load_lf(LFSeparableModel &m, const peyton::system::Config &cfg, std::string &lffile)
{
	cuxTexture<float> lf;

	// luminosity function
	if(cfg.count("lumfunc"))
	{
		lffile = cfg["lumfunc"];
		lf = load_resample_and_clip_texture_1D(lffile.c_str());
	}
	else
	{
		lffile = "";
		lf = load_constant_texture_1D(1.f, -100., 100.);
	}

	// location where the given luminosity function was sampled --
	// default is the Solar neighborhood
	float3 lfc;
	cfg.get(lfc.x, "lf_l", 0.f);	// Galactic longitude
	cfg.get(lfc.y, "lf_b", 0.f);	// Galactic latitude
	cfg.get(lfc.z, "lf_D", 0.f);	// Distance (in parsecs)
	lfc = direction(peyton::math::rad(lfc.x), peyton::math::rad(lfc.y)).xyz(lfc.z);

	// luminosity function renormalization
	m.f = 1.f;
	float norm = m.rho(lfc.x, lfc.y, lfc.z);	// Adjust f, so that future rho(x,y,z) give f at {lfc}
	assert(norm != 0.f && std::isnormal(norm));
	float f0;
	cfg.get(f0, "f", 1.f);				// normalization wrt the LF
	m.f = f0 / norm;

	// sanity check
	float norm2 = m.rho(lfc.x, lfc.y, lfc.z);
	ASSERT(fabsf(norm2/f0 - 1.f) < 1e-4);

	return lf;
}

inline void tex_get_bounds_1D(const cuxTexture<float> &tex, float &X0, float &X1)
{
	X0 = tex.coords[0].x;
	X1 = double(tex.width()-1)/tex.coords[0].y + X0;
}

inline bool absmag_adjust_from_lf(const cuxTexture<float> &tex, float &M0, float &M1)
{
	float X0, X1;
	tex_get_bounds_1D(tex, X0, X1);
	
	bool ret = false;
	if(X0 > M0) { ret = true; M0 = X0; }
	if(X1 < M1) { ret = true; M1 = X1; }
	return ret;
}

#endif
