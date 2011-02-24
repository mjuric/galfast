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

#include "model_LCBulge.h"
#include "skyconfig_impl.h"
#include "model_lib.h"
#include "transform.h"

#include <astro/math.h>

#include <astro/useall.h>
#include <fstream>

void LCBulge::load(host_state_t &hstate, const peyton::system::Config &cfg)
{
	// component ID
	int userComp = cfg.get("comp");
	comp = componentMap.seqIdx(userComp);

	// ellipsoid orientation
	std::string ctr, orient;
	cfg.get(ctr,    "center",      "galplane");
	cfg.get(orient, "orientation", "galplane");
	load_transform(&T.x, rot, ctr, orient, cfg);

	// exponential model parameters
	n = cfg.get("n");	// power law index of the exponent
	l = cfg.get("l");	// exponential scale length
	ba = cfg.get("b_a");	// b / a axis ratio
	ca = cfg.get("c_a");	// c / a axis ratio

	// Load and renormalize the luminosity function
	float f;
	cfg.get(f, "f", 1.f);
	std::string lffile;
	hstate.lf = load_lf_raw(*this, cfg, lffile);
	FOREACH(hstate.lf) { *i *= f; }

	//
	if(lffile.empty()) { lffile = "-"; }
	MLOG(verb1) << "Component " << userComp << " : " << " LCbulge {" << l << ", " << n << ", " << ca << ", " << ba << ", " << ", " << f << ", " << lffile << "}";
}

void LCBulge::prerun(host_state_t &hstate, bool draw)
{
	// bind the luminosity function texture to texture reference
	LCBulgeLF.bind(hstate.lf);
}

void LCBulge::postrun(host_state_t &hstate, bool draw)
{
	// unbind LF texture reference
	LCBulgeLF.unbind();
}

bool LCBulge::hint_absmag(host_state_t &hstate, float &M0, float &M1) const
{
	return absmag_adjust_from_lf(hstate.lf, M0, M1);
}

extern "C" skygenInterface *create_model_lcbulge()
{
	return new skygenHost<LCBulge>();
}
