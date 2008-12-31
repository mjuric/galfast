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

#include <astro/system/options.h>
#include <astro/math.h>
#include <astro/io/format.h>
#include <astro/useall.h>

#include "textstream.h"
#include "analysis.h"
#include "paralax.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <iostream>
#include <iomanip>

using namespace std;

float rri(float gr)
{
	double a = 4.9;
	double b = 2.45;
	double c = 1.68;
	double d = 0.049;
	double f = 1.39;

// 	define aa  ($b**2-3*$a*$c)
//         set a0 = 2*$b**3 - 9*$a*$b*$c + 27*$a**2*$d + 27*$a**2*ln(1-$1/$f)
//         set a1 = sqrt(-4*$aa**3 + a0**2)
//         set a2 = (-a0+a1)**(1/3.)
// 
//         set a3 = -2*$b + 2*2**(1/3.) * $aa / a2
//         set a4 = 2**(2/3.)*a2
// 
//         set $2 = 1/(6*$a)*(a3 + a4)
// 
//         if(!is_vector($1)) { define $2 ($2[0])\n }

	double aa = sqr(b) - 3*a*c;
	double a0 = 2*pow(b, 3) - 9*a*b*c + 27*sqr(a)*(d + std::log(1-gr/f));
	double a1 = sqrt(-4*pow(aa, 3) + sqr(a0));
	double a2 = pow(-a0+a1, 1./3.);

	double a3 = -2*b + 2*pow(2., (1./3.)) * aa / a2;
	double a4 = pow(2., 2./3.)*a2;

	return 1/(6*a)*(a3 + a4);
}

void locusCorrect(double gr, double ri, double sigma_g, double sigma_r, double sigma_i)
{
//	plx_gri_locus plx;
	double ml_ri = paralax.mlri(ri, gr, sigma_g, sigma_r, sigma_i);
	double ml_gr = paralax.gr(ml_ri);

	// calculate osculating error ellipse
	double x = ml_gr - gr;
	double y = ml_ri - ri;
	double z2 =
		peyton::sqr(x)/paralax.mlri.sx2
		- 2*paralax.mlri.sxy*(x)*(y)/(paralax.mlri.sx2*paralax.mlri.sy2)
		+ peyton::sqr(y)/paralax.mlri.sy2;

	cerr << "\tlocal set xp = " << gr << "\n";
	cerr << "\tlocal set yp = " << ri << "\n";
	cerr << "\tlocal define eg " << sigma_g << "\n";
	cerr << "\tlocal define er " << sigma_r << "\n";
	cerr << "\tlocal define ei " << sigma_i << "\n";
	cerr << "\tlocal define mlx " << ml_gr << "\n";
	cerr << "\tlocal define mly " << ml_ri << "\n";
	cerr << "\tlocal define d (" << z2 << "**.5)\n";

	cerr << "\n";
	cerr << "# sx2 = " << paralax.mlri.sx2 << "\n";
	cerr << "# sy2 = " << paralax.mlri.sy2 << "\n";
	cerr << "# sxy = " << paralax.mlri.sxy << "\n";
}

int main(int argc, char **argv)
{
// 	locusCorrect(0.935453, 0.407865, .05, .025, .02);
// 	locusCorrect(1.312,   1.15445, .07, .02, .04);
// 	locusCorrect(0.269789, 0.174989, .02, .03, .03);
// 	return 0;

try
{
	VERSION_DATETIME(version, "$Id: locus_demo.cpp,v 1.4 2006/07/13 23:27:30 mjuric Exp $");

	Options opts(
		argv[0],
		"This program has not been described",
		version,
		Authorship::majuric
	);

	//# add any arguments your program needs. eg:
	opts.argument("gr", "Location on the locus from which to draw source stars");
	opts.argument("sigma_g", "Add photometric errors to colors with given sigma");
	opts.argument("sigma_r", "Add photometric errors to colors with given sigma");
	opts.argument("sigma_i", "Add photometric errors to colors with given sigma");
	opts.argument("N", "Number of stars to draw");
	opts.argument("delri", "Halfwidth over which to draw gr colors (set to 0 for a point)");
	opts.argument("n", "Powerlaw index of the locus r-i projected density distribution");

	// add any options your program might need. eg:
	// opts.option("meshFactor", "meshFactor", 0, "--", Option::required, "4", "Resolution decrease between radial steps");

	parse_options(opts, argc, argv);

	/////// Start your application code here
	gsl_set_error_handler_off ();

	float gr = (float)opts["gr"];
	float sigma_g = (float)opts["sigma_g"];
	float sigma_r = (float)opts["sigma_r"];
	float sigma_i = (float)opts["sigma_i"];
	int N = (int)opts["N"];
	float delri = (float)opts["delri"];
	float n = (float)opts["n"];

//	plx_gri_locus plx;
	float ri = paralax.ri(gr);
	cerr << "back gr = " << paralax.gr(ri) << "\n";

	cout << "# gr = " << gr << "\n";
	cout << "# ri = " << ri << "\n";
	cout << "# sigma_g = " << sigma_g << "\n";
	cout << "# sigma_r = " << sigma_r << "\n";
	cout << "# sigma_i = " << sigma_i << "\n";
	cout << "# N = " << N << "\n";
	cout << "# delri = " << delri << "\n";
	cout << "# n = " << n << "\n";

	float ri0 = ri - delri;
	float ri1 = ri + delri;
	float ri0n1 = pow(ri0, n+1);
	float A = pow(ri1, n+1) - ri0n1;

	gsl_rng_env_setup();
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

	FOR(0, N)
	{
		double pri = gsl_rng_uniform(rng);
		pri = pow(A*pri + ri0n1, 1./(n+1.));
		double pgr = paralax.gr(pri);
		
		double dg = gsl_ran_gaussian(rng, sigma_g);
		double dr = gsl_ran_gaussian(rng, sigma_r);
		double di = gsl_ran_gaussian(rng, sigma_i);

		double dgr = dg - dr;
		double dri = dr - di;

		double ml_ri = paralax.mlri(pri + dri, pgr + dgr, sigma_g, sigma_r, sigma_i);
		double ml_gr = paralax.gr(ml_ri);

		cout << i << " " << pri + dri << " " << pgr + dgr << " " << ml_ri << " " << ml_gr << " ";
		cout << pri << " " << pgr << "\n";
	}
}
catch(EAny &e)
{
	e.print();
}
}
