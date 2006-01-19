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

#include "paralax.h"
#include "textstream.h"
#include "analysis.h"
#include "gslcc_min.h"
#include "gslcc_exceptions.h"
#include <fstream>

using namespace std;
//float plx_gri_locus::Mrc[4] = {4.6, 7.9, -3.0, 0.69};

void plx_gri_locus::ml_grri::setprior(const std::string &priorfile)
{
	ifstream ps(priorfile.c_str());
	itextstream pin(ps);
	
	vector<double> x, logp;
	load(pin, x, 0, logp, 1);
	double s = 0;
	FOREACH(logp) { s += *i; } // calculate sum for normalization
	FOREACH(logp) { *i = log(*i/s); }
	
	prior = gsl_spline_alloc(gsl_interp_cspline, x.size());
	gsl_spline_init(prior, &x[0], &logp[0], x.size());
	acc = gsl_interp_accel_alloc();
}

double paralax_with_prior(float ri, float gr, float sg, float sr, float si)
{
	static plx_gri_locus plx;
	if(plx.mlri.prior == NULL) { plx.mlri.setprior("ri_prior.txt"); }

	double RIp;
	try {
		RIp = plx.mlri(ri, gr, sg, sr, si);
	}
	catch(peyton::exceptions::EGSLMinimizer &e)
	{
		// the star is unfittable
		return -1;
	}

	return RIp;
}

double paralax_without_prior(float ri, float gr, float sg, float sr, float si)
{
	static plx_gri_locus plx;
	
	double RI;
	try {
		RI = plx.mlri(ri, gr, sg, sr, si);
	}
	catch(peyton::exceptions::EGSLMinimizer &e)
	{
		// the star is unfittable
		return -1;
	}

	return RI;
}
