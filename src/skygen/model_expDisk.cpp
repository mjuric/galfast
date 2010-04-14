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
#include "transform.h"

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
	l        = cfg.get("l");
	h        = cfg.get("h");

	assert(l > 0);
	assert(h > 0);
#if 0
	float z0 = cfg.get("z0");	// Solar offset from the Galactic plane
	Rg       = cfg.get("Rg");	// Distance to the Galactic center (assumed to be in l=0, b=0 direction)

	// compute the rotation of the Galactic plane, and cylindrical
	// coordinates of the Sun
	zsun = z0;
	rsun = sqrt(double(Rg*Rg) - double(z0*z0));
	asin = zsun / Rg;
	acos = rsun / Rg;
#else
	std::string ctr, orient;
	cfg.get(ctr,    "center",      "galplane");
	cfg.get(orient, "orientation", "galplane");
	load_transform(&T.x, M, ctr, orient, cfg);
#if 0
	std::cout << "translation = " << std::setprecision(10) << T.x << " " << T.y << " " << T.z << "\n";
	print_matrix(M);

	// Do some testing here...
	float3 v = { 0.f, 0.f, 0.f };
	v = transform(v, T, M);
	std::cout << std::setprecision(10) << v.x << " " << v.y << " " << v.z << "\n";
	abort();
#endif
#endif

	int userComp = cfg.get("comp");
	comp = componentMap.seqIdx(userComp);

	// set this to 0. for now, to allow LF normalization
	// even beyond the user-specified model cutoff
	r_cut2 = 0.f;

	// Load luminosity function
	std::string lffile;
	hstate.lf = load_lf(*this, cfg, lffile);

	// cutoff radius (default: no cutoff)
	cfg.get(r_cut2,  "rcut",   0.f);
	r_cut2 *= r_cut2;

	if(lffile.empty()) { lffile = "-"; }
	MLOG(verb1) << "Component " << userComp << " : " << "exponential disk {" << l << ", " << h << ", " << lffile << "}";
}

bool expDisk::hint_absmag(host_state_t &hstate, float &M0, float &M1) const
{
	return absmag_adjust_from_lf(hstate.lf, M0, M1);
}

extern "C" skygenInterface *create_model_exponentialdisk()
{
	return new skygenHost<expDisk>();
}
