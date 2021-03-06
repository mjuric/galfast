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

#include "galfast_config.h"

#include "model_powerLawEllipsoid.h"
#include "skyconfig_impl.h"
#include "model_lib.h"
#include "transform.h"

#include <astro/useall.h>

void powerLawEllipsoid::prerun(host_state_t &hstate, bool draw)
{
	// bind the luminosity function texture to texture reference
	powerLawEllipsoidLF.bind(hstate.lf);
}

void powerLawEllipsoid::postrun(host_state_t &hstate, bool draw)
{
	// unbind LF texture reference
	powerLawEllipsoidLF.unbind();
}

void powerLawEllipsoid::load_geometry(host_state_t &hstate, const peyton::system::Config &cfg)
{
#if 0
	// center of the ellipsoid
	cfg.get(c.x, "center_l", 0.f);	// Galactic longitude
	cfg.get(c.y, "center_b", 0.f);	// Galactic latitude
	cfg.get(c.z, "center_D", 0.f);	// Distance (in parsecs)
	c = direction(peyton::math::rad(c.x), peyton::math::rad(c.y)).xyz(c.z);

	std::string orientation; double M[3][3];
	cfg.get(orientation, "orientation", "");
	get_rotation_matrix(M, orientation, cfg);
	for(int i = 0; i != 3; i++)
		for(int j = 0; j != 3; j++)
			rot[i][j] = M[i][j];
#else
	std::string ctr, orient;
	cfg.get(ctr,    "center",      "galplane");
	cfg.get(orient, "orientation", "galplane");
	load_transform(&T.x, rot, ctr, orient, cfg);
#if 0
	std::cout << "translation = " << std::setprecision(10) << c.x << " " << c.y << " " << c.z << "\n";
	print_matrix(rot);

	// Do some testing here...
	float3 v = { 0.f, 100.f, -24.f };
	v = transform(v, c, rot);
	std::cout << std::setprecision(10) << v.x << " " << v.y << " " << v.z << "\n";
	abort();
#endif
#endif

	cfg.get(ca, "c_a", 1.f);	// ratio of c / a (z and x axes of the ellipsoid)
	cfg.get(ba, "b_a", 1.f);	// ratio of b / a (z and x axes of the ellipsoid)
}

void powerLawEllipsoid::load(host_state_t &hstate, const peyton::system::Config &cfg)
{
	// turn off density truncation
	rminSq = rmaxSq = 0.f;

	load_geometry(hstate, cfg);

	// radial functional dependence
	cfg.get(n, "n", -3.f);		// power law index

	// component ID
	int userComp = cfg.get("comp");
	comp = componentMap.seqIdx(userComp);

	// Load luminosity function
	std::string lffile;
	hstate.lf = load_lf(*this, cfg, lffile);

	// limits (\rho = 0 for d < rminSq || d > rmaxSq)
	// NOTE: This is intentionally after LF normalization computation, so
	// that normalizing to the values of the profile in the regions
	// excluded by this cut would still be possible.
	cfg.get(rminSq, "rmin", 0.f);	rminSq *= rminSq;
	cfg.get(rmaxSq, "rmax", 0.f);	rmaxSq *= rmaxSq;

	if(lffile.empty()) { lffile = "-"; }
	MLOG(verb1) << "Component " << userComp << " : " << "power law ellipsoid {" << n << ", " << ca << ", " << ba << ", " << lffile << "}";
}

bool powerLawEllipsoid::hint_absmag(host_state_t &hstate, float &M0, float &M1) const
{
	return absmag_adjust_from_lf(hstate.lf, M0, M1);
}

extern "C" skygenInterface *create_model_powerlawellipsoid()
{
	return new skygenHost<powerLawEllipsoid>();
}
