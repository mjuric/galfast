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

#include "model_powerLawEllipsoid.h"
#include "skyconfig_impl.h"

#include <astro/system/config.h>

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

void powerLawEllipsoid::load(host_state_t &hstate, const peyton::system::Config &cfg)
{
	rminSq = 0.f;
	rmaxSq = std::numeric_limits<float>::max();

	// center of the ellipsoid
	cfg.get(c[0], "x", 0.f);
	cfg.get(c[1], "y", 0.f);
	cfg.get(c[2], "z", 0.f);

	cfg.get(n, "n", 0.f);		// power law index
	cfg.get(ca, "c_a", 1.f);	// ratio of c / a (z and x axes of the ellipsoid)
	cfg.get(ba, "b_a", 1.f);	// ratio of b / a (z and x axes of the ellipsoid)

	// Euler rotation for a point in Galactic frame to the xyz frame of the ellipsoid
	//
	// If imagining this as the rotation of the _ellipsoid_, the following
	// rotation sequence occurs:
	//
	//	1) ccw around the z axis by an angle phi
	//	2) ccw around the x axis by theta
	//	3) ccw around the z acis by psi
	//
	// See http://mathworld.wolfram.com/EulerAngles.html, eqs 6-14
	//
	float phi, theta, psi;
	cfg.get(phi,   "phi",   0.f);	phi   = rad(phi);
	cfg.get(theta, "theta", 0.f);	theta = rad(theta);
	cfg.get(psi,   "psi",   0.f);	psi   = rad(psi);

	rot[0][0] =  cos(psi)   * cos(phi)   - cos(theta) * sin(phi) * sin(psi);
	rot[0][1] =  cos(psi)   * sin(phi)   + cos(theta) * cos(phi) * sin(psi);
	rot[0][2] =  sin(psi)   * sin(theta);
	rot[1][0] = -sin(psi)   * cos(phi)   - cos(theta) * sin(phi) * cos(psi);
	rot[1][1] = -sin(psi)   * sin(phi)   + cos(theta) * cos(phi) * cos(psi);
	rot[1][2] =  cos(psi)   * sin(theta);
	rot[2][0] =  sin(theta) * sin(phi);
	rot[2][1] = -sin(theta) * cos(phi);
	rot[2][2] =  cos(theta);

#if 0
	std::cerr << deg(phi) << " " << deg(theta) << " " << deg(psi) << "\n";
	FOR(0, 3)
	{
		FORj(j, 0, 3) { std::cerr << rot[i][j] << " "; }
		std::cerr << "\n";
	}

	// Do some testing here...
	float3 v = { 1.f, 0.f, 0.f };
	v = matmul3d(rot, v);
	std::cout << v.x << " " << v.y << " " << v.z << "\n";
	exit(0);
#endif

	// component ID
	cfg.get(comp, "comp", 0);

	// luminosity function
	// location where the given luminosity function was sampled
	// (this is usually the solar neighborhood, but since the
	// coordinates galactocentric, and the the position of the
	// Sun is user-settable, the default is set to the center of
	// the Galaxy)
	// NOTE: This default becomes invalid if the profile diverges
	// at the center.
	float lfc[3];
	cfg.get(lfc[0], "lfx", 0.f);
	cfg.get(lfc[1], "lfy", 0.f);
	cfg.get(lfc[2], "lfz", 0.f);

	cfg.get(f, "f", 1.f);				// normalization wrt the LF
	float norm = rho(lfc[0], lfc[1], lfc[2]);	// Adjust f, so that future rho(x,y,z) give f at {lfc}
	assert(norm != 0.f && std::isnormal(norm));
	f /= norm;

	// limits (\rho = 0 for d < rminSq || d > rmaxSq)
	// NOTE: This is intentionally after LF normalization computation, so
	// that normalizing to the values of the profile in the regions
	// excluded by this cut would still be possible.
	cfg.get(rminSq, "rmin", 0.f);	rminSq *= rminSq;
	cfg.get(rmaxSq, "rmax", 1e6f);	rmaxSq *= rmaxSq;

	if(cfg.count("lumfunc"))
	{
		hstate.lf = load_and_resample_1D_texture(cfg["lumfunc"].c_str());
	}
	else
	{
		hstate.lf = load_constant_texture(1.f);
	}
}

extern "C" skyConfigInterface *create_model_powerlawellipsoid()
{
	return new skyConfig<powerLawEllipsoid>();
}
