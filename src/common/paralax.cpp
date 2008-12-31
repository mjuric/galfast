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

#include "paralax.h"
#include "textstream.h"
#include "analysis.h"
#include "gslcc_min.h"
#include "gslcc_exceptions.h"
#include <fstream>
#include <sstream>

#include <astro/system/fs.h>
#include <astro/useall.h>
using namespace std;
//float plx_gri_locus::Mrc[4] = {4.6, 7.9, -3.0, 0.69};

plx_gri_locus_ng paralax;

plx_gri_locus_ng::plx_gri_locus_ng()
: mlri(*this), plx_modified(false)
{
	// check if paralax relation was specified through an environment variable
	EnvVar plxvar("PARALAX");
	std::vector<double> tmp;
	std::string plx;
	if(plxvar)
	{
		DLOG(verb1) << "++++ Taking paralax from $PARALAX env. var.";
		plx = (std::string)plxvar;
	}
	else
	{
		DLOG(verb1) << "++++ Using default paralax coefficients.";
		plx = "3.2 13.30 -11.50 5.40 -0.65";
	}

	istringstream in(plx.c_str());
	double coef;
	while(in >> coef)
	{
		tmp.push_back(coef);
	}
	setParalaxCoefficients(tmp);
}

inline OSTREAM(const std::vector<double> &a)
{
	if(a.empty()) { return out; }
	out << a[0];
	FOR(1, a.size()) { out << " " << a[i]; }
	return out;
}

void plx_gri_locus_ng::setParalaxCoefficients(const std::vector<double> &Mcoef)
{
	Mrc = Mcoef;
	plx_modified = true;

	DLOG(verb1) << "color-absmag coef. = " << Mcoef;
}

void plx_gri_locus_ng::ml_grri::setprior(const std::string &priorfile)
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

double paralax_with_prior(float ri, float gr, float sg, float sr, float si, float *lnL)
{
	static plx_gri_locus_ng plx;
	if(plx.mlri.prior == NULL) { plx.mlri.setprior("ri_prior.txt"); }

	double RIp;
	try {
		RIp = plx.mlri(ri, gr, sg, sr, si, lnL);
	}
	catch(peyton::exceptions::EGSLMinimizer &e)
	{
		// the star is unfittable
		return -1;
	}

	return RIp;
}

double paralax_without_prior(float ri, float gr, float sg, float sr, float si, float *lnL)
{
//	static plx_gri_locus plx;
	
	double RI;
	try {
		RI = paralax.mlri(ri, gr, sg, sr, si, lnL);
	}
	catch(peyton::exceptions::EGSLMinimizer &e)
	{
		// the star is unfittable
		return -1;
	}

	return RI;
}

// Finds the distance of point (gr, ri) from the locus
struct distance_from_locus_t : public gsl::mmizer
{
	double y0, x0;

	distance_from_locus_t() : mmizer(-.5, 2.5, 0, .001) { }

	double operator()(float gr, float ri, float *ri_closest = NULL)
	{
		y0 = ri; x0 = gr;
		float dummy;
		if(ri_closest == NULL) { ri_closest = &dummy; }

		*ri_closest = evaluate(y0);
		return fn_to_minimize(*ri_closest);
	}

	double fn_to_minimize(double y)
	{
		double x = paralax.gr(y);
		return sqrt(sqr(x-x0)+sqr(y-y0));
	}
} distance_from_locus_obj;

double distance_from_locus(float gr, float ri, float *ri_closest)
{
	try {
		return distance_from_locus_obj(gr, ri, ri_closest);
	}
	catch(peyton::exceptions::EGSLMinimizer &e)
	{
		return -1;
	}
}
