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

#include "config.h"

#include "model_expDisk.h"
#include "skyconfig_impl.h"

#include "model_lib.h"

void expDisk::prerun(host_state_t &hstate, bool draw)
{
	// bind the luminosity function texture to texture reference
	expDiskLF.bind(hstate.lf);
}

void expDisk::postrun(host_state_t &hstate, bool draw)
{
	// unbind LF texture reference
	expDiskLF.unbind();
}

void expDisk::load(host_state_t &hstate, const peyton::system::Config &cfg)
{
	// density distribution parameters
	l     = cfg.get("l");
	h     = cfg.get("h");
	z0    = cfg.get("z0");
	Rg    = cfg.get("Rg");

	comp  = cfg.get("comp");

	// set this to 0. for now, to allow LF normalization
	// even beyond the user-specified model cutoff
	r_cut2 = 0.f;

	// Load luminosity function
	hstate.lf = load_lf(*this, cfg);

	// cutoff radius (default: no cutoff)
	cfg.get(r_cut2,  "rcut",   0.f);
	r_cut2 *= r_cut2;
}

extern "C" skygenInterface *create_model_exponentialdisk()
{
	return new skygenHost<expDisk>();
}
