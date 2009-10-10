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

#if 0
#include "config.h"

#include "model_J08.h"
#include "skyconfig_impl.h"

#include <astro/system/config.h>

void J08::prerun(host_state_t &hstate, bool draw)
{
	// bind the luminosity function texture to texture reference
	J08LFManager.bind(hstate.lf);
}

void J08::postrun(host_state_t &hstate, bool draw)
{
	// unbind LF texture reference
	J08LFManager.unbind();
}

void J08::load(host_state_t &hstate, const peyton::system::Config &cfg)
{
	// density distribution parameters
	rho0  = cfg.get("rho0");
	l     = cfg.get("l");
	h     = cfg.get("h");
	z0    = cfg.get("z0");
	f     = cfg.get("f");
	lt    = cfg.get("lt");
	ht    = cfg.get("ht");
	fh    = cfg.get("fh");
	q     = cfg.get("q");
	n     = cfg.get("n");

	// luminosity function
	if(cfg.count("lumfunc"))
	{
		hstate.lf = load_and_resample_texture_1D(cfg["lumfunc"].c_str());
	}
	else
	{
		hstate.lf = load_constant_texture_1D(1.f, -100., 100.);
	}

	// cutoff radius (default: 1Mpc)
	cfg.get(r_cut2,  "rcut",   1e6f);
	r_cut2 *= r_cut2;

	// component IDs
	cfg.get(comp_thin,  "comp_thin",  0);
	cfg.get(comp_thick, "comp_thick", 1);
	cfg.get(comp_halo,  "comp_halo",  2);
}

extern "C" skygenInterface *create_model_j08()
{
	return new skygenHost<J08>();
}

// Backwards compatibility
extern "C" skygenInterface *create_model_bahcallsoneira()
{
	return new skygenHost<J08>();
}
#endif
