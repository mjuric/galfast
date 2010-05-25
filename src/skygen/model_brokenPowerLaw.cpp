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

#include "model_brokenPowerLaw.h"
#include "skyconfig_impl.h"
#include "model_lib.h"

#include <astro/math.h>

#include <astro/useall.h>
#include <fstream>

void brokenPowerLaw::load(host_state_t &hstate, const peyton::system::Config &cfg)
{
	// turn off density truncation
	rminSq = rmaxSq = 0.f;

	load_geometry(hstate, cfg);

	// radial functional dependence
	std::vector<float> n = cfg.get("n");		// power law index
	std::vector<float> rbreak;
	if(n.size() > 1)
	{
		rbreak = cfg.get("rbreak");	// power law index breaks (in parsecs)
	}

	ASSERT(n.size() < MAXPOWERLAWS);
	ASSERT(rbreak.size() == n.size()-1);

	FOR(0, n.size()) 	{ nArr[i]     = n[i];		}
	FOR(0, rbreak.size())	{ rbreakSq[i] = sqr(rbreak[i]); }
	nbreaks = rbreak.size();

	// component ID
	int userComp = cfg.get("comp");
	comp = componentMap.seqIdx(userComp);

	// compute piece-wise normalizations so that the density
	// profile remains continuous at power-law breaks
	fArr[0] = f = 1.f;
	FOR(0, nbreaks)
	{
		fArr[i+1] = fArr[i] * powf(rbreakSq[i], 0.5f*(nArr[i]-nArr[i+1]));

		// correctness tests
		float rho0 =   fArr[i] * powf(rbreakSq[i], 0.5f*nArr[i]);
		float rho1 = fArr[i+1] * powf(rbreakSq[i], 0.5f*nArr[i+1]);
		ASSERT(fabsf(rho0/rho1 - 1.f) < 1e-4);
	}

	// Load luminosity function (this will also compute the
	// overall normalization)
	std::string lffile;
	hstate.lf = load_lf(*this, cfg, lffile);

	// limits (\rho = 0 for d < rminSq || d > rmaxSq)
	// NOTE: This is intentionally after LF normalization computation, so
	// that normalizing to the values of the profile in the regions
	// excluded by this cut would still be possible.
	cfg.get(rminSq, "rmin", 0.f);	rminSq *= rminSq;
	cfg.get(rmaxSq, "rmax", 0.f);	rmaxSq *= rmaxSq;

#if 0
	std::ofstream out("res.txt");
	for(float x = -110000.; x < 110000.; x+= 500.)
	{
		float rho1 = rho(x, 0, 0);
		out << x << " " << rho1 << "\n";
	}
	out.close();
	exit(0);
#endif

	if(lffile.empty()) { lffile = "-"; }
	std::ostringstream ss;
	ss << nArr[0];
	FOR(0, nbreaks) { ss << "=(" << sqrt(rbreakSq[i])/1000 << "kpc)=" << nArr[i+1]; }
	MLOG(verb1) << "Component " << userComp << " : " << "broken power law {" << ss.str() << ", " << ca << ", " << ba << ", " << lffile << "}";
}


extern "C" skygenInterface *create_model_brokenpowerlaw()
{
	return new skygenHost<brokenPowerLaw>();
}
