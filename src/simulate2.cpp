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

/**
	Basic algorithm outline, definitions:
	  * (x, y) - Lambert map coordinates (transformable to l, b)
	  * m - r band magnitude

	Integrations:
	  * split (x, y) lambert space into 1deg^2 grid
	  * For each nonempty (x,y) pixel (== a beam in 3D space):
	     * split r-i into 0.01mag bins
	     * For each r-i bin, calculate:
	        * N(x, y, ri) = Integrate[rho(m; x,y,ri) dm, m_min, m_max]
		  where m is r band magnitude, m_min and m_max are flux
		  limits of the survey
	  * Calculate the sum - this is the total number of stars expected 
	    in survey volume. Normalize by this number to convert the
	    density to probability.
	  * Calculate marginal distributions:
	     * N(ri|x,y), N(y|x), N(x)
	  * Calculate spline approximations to integrals of marginalized
	    distributions for use when transforming uniform variates

	Sampling:
	  * Generating the sample of N stars, step 1
	     * 1: Sample p(x) to get x
	     * Sample p(y|x) to get y
	     * Reject and goto 1 if not in survey volume
	     * Sample p(ri|x,y) to get ri
	     * Store in array, goto 1 if number of stars != N
	  * Calculating magnitudes
	     * Sort the array of generated stars by x,y,ri
	     * For each star generated
	       * If not cached from previous star, calculate spline
	         approximation to p(m | x, y, ri)
	       * Draw a random variate from p(m | x, y, ri) - that is the
	         magnitude.
*/

#include "config.h"

#ifdef COMPILE_SIMULATE_X

#include "gsl/gsl_randist.h"

namespace lapack
{
	extern "C"
	{
		#include <g2c.h>
		#include "common/clapack.h"
	}
}

extern "C"
{
	#include "gpc/gpc.h"
}

//#include <nurbs++/nurbsS.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "simulate.h"
#include "projections.h"
#include "model.h"
#include "paralax.h"
#include "analysis.h"
#include "dm.h"

#include <vector>
#include <map>
#include <string>

#include <astro/io/binarystream.h>
#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

BLESS_POD(gpc_vertex);

// binary serialization routines for gpc_polygons and aux. data structures
BOSTREAM2(const gpc_vertex_list& p)
{
	out << p.num_vertices;
	FOR(0, p.num_vertices)
	{
		out << p.vertex[i];
	}
	return out;
}

BISTREAM2(gpc_vertex_list& p)
{
	in >> p.num_vertices;
	p.vertex = (gpc_vertex*)malloc(p.num_vertices*sizeof(gpc_vertex));
	FOR(0, p.num_vertices)
	{
		in >> p.vertex[i];
	}
	return in;
}

BOSTREAM2(const gpc_polygon& p)
{
	out << p.num_contours;
	FOR(0, p.num_contours)
	{
		out << p.hole[i];
		out << p.contour[i];
	}
	return out;
}

BISTREAM2(gpc_polygon& p)
{
	in >> p.num_contours;
	p.contour = (gpc_vertex_list*)malloc(p.num_contours*sizeof(gpc_vertex_list));
	p.hole = (int*)malloc(p.num_contours*sizeof(int));
	FOR(0, p.num_contours)
	{
		in >> p.hole[i];
		in >> p.contour[i];
	}
	return in;
}

void poly_bounding_box(double &x0, double &x1, double &y0, double &y1, const gpc_polygon &p);

// create a gpc_polygon rectangle
gpc_polygon poly_rect(double x0, double x1, double y0, double y1)
{
	static gpc_vertex v[4];
	static gpc_vertex_list vl = {4, v};
	static gpc_polygon rect = {1, NULL, &vl};

	v[0].x = x0; v[0].y = y0;
	v[1].x = x1; v[1].y = y0;
	v[2].x = x1; v[2].y = y1;
	v[3].x = x0; v[3].y = y1;

	return rect;
}

sim *sim::star::gsim = NULL;

// 	lambert proj(rad(90), rad(90));
// 	double dx = rad(1); // model grid angular resolution
// 	double dri = 0.01; // model CMD resolution
// 	double ri0 = 0.0, ri1 = 1.5;	// color limits
// 	double m0 = 14, m1 = 22;	// magnitude limits
// 	double dm = 0.1;	// model CMD magnitude resolution
// 	double x0, x1, y0, y1;	// lambert survey footprint bounding box
sim::sim(const std::string &pfx, galactic_model *model_)
: prefix(pfx), proj(rad(90), rad(90)), model(model_)
{
	sim::star::gsim = this;

	Config cfg(pfx + ".conf");

	cfg.get(dx,	"dx", 	rad(1));
	cfg.get(dri,	"dri",	0.01);
	cfg.get(ri0,	"ri0",	0.0);
	cfg.get(ri1,	"ri1",	1.5);
	cfg.get(m0,	"m0",	14.); 
	cfg.get(m1,	"m1",	22.);
	cfg.get(dm,	"dm",	0.1);
	
	LOG(app, verb1) << "dx  = " << deg(dx);
	LOG(app, verb1) << "dri = " << dri;
	LOG(app, verb1) << "ri0 = " << ri0;
	LOG(app, verb1) << "ri1 = " << ri1;
	LOG(app, verb1) << "m0  = " << m0;
	LOG(app, verb1) << "m1  = " << m1;
	LOG(app, verb1) << "dm  = " << dm;
}

bool operator <(const sim::star &a, const sim::star &b)
{
	int aX = a.X();
	int bX = b.X();
	int aY = a.Y();
	int bY = b.Y();
	return aX < bX || (aX == bX && (
		aY < bY || (aY == bY && 
		a.ri < b.ri)
		));
}

// precalculated sines and cosines of beam coordinates
struct pencil_beam
{
	Radians l, b;
	double cl, cb, sl, sb;

	pencil_beam(Radians l_, Radians b_, double d_ = 0.0)
	: l(l_), b(b_),
	  cl(cos(l_)), cb(cos(b_)), sl(sin(l_)), sb(sin(b_))
	{ }
};

// precalculated coordinates of a point in various coordinate systems
struct coord_pack
{
	const static double Rg = 8000.;

	double d, x, y, z, rcyl;

	coord_pack(const pencil_beam &p, const double d_)
	{
		d = d_;

		x = Rg - d*p.cl*p.cb;
		y = -d*p.sl*p.cb;
		z = d*p.sb;

		rcyl = sqrt(x*x + y*y);
	}
};

#define VARRAY(type, ptr, i) *((type*)(((void **)(ptr))[i]))
// this function returns the integrand rho(m) for integration
// of (rho(m) dm), where rho(m) is the distribution of stars
// in a given direction with a given magnitude m and color r-i
// in interval dm and solid angle dOmega
extern "C"
double model_fun_1(double m, void *param)
{
	// unpack the model and parameters
	pencil_beam &p = VARRAY(pencil_beam, param, 0);
	galactic_model &model = VARRAY(galactic_model, param, 1);
	double Mr = VARRAY(double, param, 2);
	double ri = VARRAY(double, param, 3);

	// calculate the 3D real space point corresponding to this
	// (l,b,m) triplet (given the Mr of the star)
	double d = stardist::D(m, Mr);
	coord_pack c(p, d);

	// rho is the number of stars at a point c, per pc^3, per dri
	double rho = model.rho(c.x, c.y, c.z, ri);
	static double ln10_over_5 = log(10.)/5.;
	// convert to the number of stars at point c, per dm, dOmega
	rho *= pow(d, 3.) * ln10_over_5;
	return rho;
}

int model_fun(double m, double const* y, double *dydx, void *param)
{
	double norm = VARRAY(double, param, 4);
	*dydx = model_fun_1(m, param) / norm;
	//std::cerr << "m = " << m << ", r = " << *dydx << "\n";
	return GSL_SUCCESS;
}

// these are defined in simulate.cpp
typedef int (*int_function)(double dd, double const* y, double *dydd, void *param);
bool sample_integral(const std::valarray<double> &xv, std::valarray<double> &yv, int_function f, void *param);

#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>

#if 0
double
integrate(double a, double b, int_function f, void *param, double abserr, double relerr)
{
	gsl_odeiv_step*    s = gsl_odeiv_step_alloc(gsl_odeiv_step_rk2, 1);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(abserr, relerr);
	gsl_odeiv_evolve*  e = gsl_odeiv_evolve_alloc(1);

	double h0 = (b - a) / 1000;	// initial step size
	double y = 0;			// initial integration value
	double x = a;			// initial integration point

	gsl_odeiv_system sys = {f, NULL, 1, param};

	// integration from 0 to dmax, but in ninterp steps
	// when integrating between dmin and dmax
	double h = h0;
	while (x < b)
	{
		gsl_odeiv_evolve_apply (e, c, s,
					&sys,
					&x, b, &h,
					&y);
		std::cerr << "# " << x << " " << y << " " << h << "\n";
	}
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);

	exit(-1);
	return y;
}
#else
/*bool
sample_integral2(const double a, const double b, double dx,
	std::vector<double> &xv, std::vector<double> &yv, int_function f, void *param,
	double abserr = 1e-6, double relerr = 0)*/
double
integrate(double a, double b, double dx, int_function f, void *param, double abserr, double relerr)
{
// 	abserr = 0;
// 	relerr = 1;

	gsl_odeiv_step*    s = gsl_odeiv_step_alloc(gsl_odeiv_step_rk2, 1);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(abserr, relerr);
	gsl_odeiv_evolve*  e = gsl_odeiv_evolve_alloc(1);

	double h = dx/2.;	// initial step size
	double y = 0;		// initial integration value
	double x = a;		// initial integration point

	ASSERT(VARRAY(double, param, 4) != 0);

	gsl_odeiv_system sys = {f, NULL, 1, param};

	double xto = x+dx;
	bool last = false;
	while(!last)
	{
		if(xto > b) { xto = b; last = true; }

		while (x < xto)
		{
			int status = gsl_odeiv_evolve_apply (e, c, s,
						&sys,
						&x, xto, &h,
						&y);
			if(status != GSL_SUCCESS)
			{
				std::cerr << "Error while integrating: " << gsl_strerror(status) << "\n";
				exit(-1);
			}
			double rho;
			GSL_ODEIV_FN_EVAL(&sys, x+h, &y, &rho);
			std::cout << "[" << VARRAY(double, param, 4) << "]: " << x << " " << y << " " << dx << " " << h << " " << rho << "\n";
/*			if(fabs(VARRAY(double, param, 3) - 1.27) < 0.01)
				std::cout << "[" << VARRAY(double, param, 3) << "]: " << x << " " << y << " " << dx << " " << h << "\n";*/
		}
		//std::cout << "[" << VARRAY(double, param, 4) << "]: " << x << " " << y << " " << dx << " " << h << "\n";

		xto += dx;
	}
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);

	return y;
}
#endif

//
// Populates the pdf vector with cumulative marginal PDFs of stars within a given
// field (x,y).
//
double sim::ri_mpdf(std::vector<double> &pdf, const double x, const double y)
{
	// deproject to (l,b)
	Radians l, b;
	proj.inverse(x, y, l, b);
	pencil_beam pb(l, b);

	const int Npts = 1000;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (Npts);
	pdf.clear();
	pdf.reserve((size_t)((ri1-ri0)/dri)+2);
	double sum = 0;

	bool act = false;
	// loop through spectral types (r-i slices), and calculate the expected
	// total number of stars with a given SpT (assuming the model returned by
	// model class), in a direction x, y
	for(double ri = ri0; ri <= ri1; ri += dri)
	{
		// get a model which will return the number of stars n in
		// direction (l,b) in the SpT interval [ri, ri+dri) (approximate
		// this by assuming a constant f=dN/dri in the short interval dri
		// and taking n = f*dri
		double ri_mid = ri + 0.5*dri;
		double Mr = model->absmag(ri_mid); // absolute magnitude for the stars with this SpT

		double norm = 1.;
		void *params[5] = { &pb, model, &Mr, &ri_mid, &norm };

		double n, error;
		gsl_function F = { &model_fun_1, params };
		#if 0
		ASSERT(gsl_integration_qag (&F, m0, m1, 0, 1e-7, Npts, GSL_INTEG_GAUSS61, w, &n, &error) == GSL_SUCCESS);
		#else
		n = integrate(m0, m1, dm, model_fun, params, 1e-10, 1e-5);
		#endif
		//ASSERT(n != 0) { std::cerr << "BU! " << x << " " << y << " " << ri_mid << " " << n << "\n"; }
		n *= sqr(dx)*dri; // multiply by solid angle. n is now the total number of stars in direction (l,b)
		                  // with colors in [ri, ri+dri)
		#if 0
		if(n != 0) { std::cerr << "BU!" << x << " " << y << " " << ri_mid << " " << n << "\n"; }
		if(std::abs(x + 0.0912734) < 0.0001 &&
		   std::abs(y + 0.0563668) < 0.0001 &&
		   std::abs(ri_mid - .215) < 0.0001)
		   {
 			std::cerr << "BU!" << x << " " << y << " " << ri_mid << "\n";
			std::cerr << "sum = " << sum << ", n = " << n << "\n";
			std::cerr << "size = " << pdf.size() << "\n";
			std::cerr << "err = " << error << "\n";
			act = true;
//			exit(-1);
 		   }
		#endif

		// store cumulative distribution (note: sum += n line _must_ come after push_back,
		// to have the cumulative pdf value for index I be the number of stars *less* (and
		// not .LE.) than I)
		pdf.push_back(sum);
		sum += n;
	}
	gsl_integration_workspace_free(w);

	// convert to probability
	FOREACH(pdf) { 
		//if(act && i - pdf.begin() == 22) { std::cerr << "XXX: " << *i << "\n"; }
		*i /= sum;
		//if(act && i - pdf.begin() == 22) { std::cerr << "XXX: " << *i << "\n"; }
	}

	//std::cerr << "Exit\n";
	return sum;
}

#define ITER(x) typeof((x).begin())

// serialize the object contained in boost::shared_ptr.
// useful for serializing of shared_ptr vectors or lists.
template <typename T>
	inline BOSTREAM2(const boost::shared_ptr<T> &p)
{
	if(p.get())
	{
		out << true << *p;
	} else {
		out << false;
	}
	return out;
}

template <typename T>
	inline BISTREAM2(boost::shared_ptr<T> &p)
{
	bool initialized;
	in >> initialized;
	if(initialized)
	{
		p.reset(new T);
		in >> *p;
	}
}

// cumulative marginal pdf in lambert y coodinate, for a given lambert x coordinate
// (Sy::pcum == probability that a star 
// within the observation footprint and given the model and x coordinate,
// will have the value of the lambert coordinate y *less* than Sx::Y)
struct Sy
{
	double pcum;
	int Y;
	std::vector<double> ripdf;
	Sy(double _pcum = 0, int _Y = 0) : pcum(_pcum), Y(_Y) {}
private:
	Sy(const Sy&);
	Sy& operator=(const Sy&);
};
bool operator<(const Sy& a, const Sy& b) { return a.pcum < b.pcum; }
inline BOSTREAM2(const Sy &p) { return out << p.pcum << p.Y << p.ripdf; }
inline BISTREAM2(Sy &p) { return in >> p.pcum >> p.Y >> p.ripdf; }

// cumulative marginal pdf in lambert x coodinate (Sx::pcum == probability that a star 
// within the observation footprint and given the model will have a 
// the value of the lambert coordinate x *less* than Sx::X)
struct Sx
{
	double pcum;
	int X;
	std::vector<boost::shared_ptr<Sy> > ypdf;
	Sx(double _pcum = 0, int _X = 0) : pcum(_pcum), X(_X) {}
private:
	Sx(const Sx&);
	Sx& operator=(const Sx&);
};
bool operator<(const Sx& a, const Sx& b) { return a.pcum < b.pcum; }
inline BOSTREAM2(const Sx &p) { return out << p.pcum << p.X << p.ypdf; }
inline BISTREAM2(Sx &p) { return in >> p.pcum >> p.X >> p.ypdf; }

// set of marginal pdfs in lambert x. It's the trunk of a tree
// of marginal PDFs, starting with mpdf in X, then Y, and
// finaly r-i. First branches of the tree are Sx structures,
// followed by Sy structures (stored in Sx::ypdf), followed by
// doubles of r-i mpdfs (in Sy::ripdf vector).
typedef std::vector<boost::shared_ptr<Sx> > XPDF;

OSTREAM(const XPDF &xpdf)
{
	// TODO: unfinished
	FOREACHj(xi, xpdf)
	{
		double x = (*xi)->X;
		FOREACHj(yi, (*xi)->ypdf)
		{
			double y = (*yi)->Y;
			FOREACHj(rii, (*yi)->ripdf)
			{
				//double ri = 
				out << x << " " << y << " " << *rii << "\n";
			}
		}
	}
}

gpc_polygon make_circle(double x0, double y0, double r, double dx);

// calculate marginal PDFs (``mpfds'') for stars within a given sky footprint,
// and with a given Galactic model. Store the mPDFs into an XPDF tree
// and serialize it to a file. These mpdfs are later used in generating
// Monte Carlo realisations of the sky
//
// This routine also calculates the 'skymap', for fast point-in-polygon
// lookups of the observed sky footprint.
void sim::precalculate_mpdf()
{
	using boost::shared_ptr;
	using std::cout;

	gpc_polygon sky;
#if 0
 	FILE *fp = fopen((prefix + ".gpc.txt").c_str(), "r");
 	gpc_read_polygon(fp, 1, &sky);
	fclose(fp);
#else
	sky = make_circle(0., 0., 0.1, dx);
	//sky = make_circle(-0.123308, 0.406791, 0.1, dx);
#endif

	std::map<std::pair<int, int>, gpc_polygon> skymap;	// a map of rectangular sections of the sky, for fast is-point-in-survey-area lookup
 	poly_bounding_box(x0, x1, y0, y1, sky);

// 	ri_mpdf(pdf[std::make_pair(X, Y)], factory, 0, 0);
// 	return;

	XPDF xpdf;
	int N = 0; // number of dx^2 cells processed
	double xsum = 0; // Grand total - number of stars in the whole field
	int X = 0;
	for(double x = x0; x < x1; x += dx, X++) // loop over all x values in the bounding rectangle
	{
		shared_ptr<Sx> sx(new Sx(xsum, X));

		double xa = x, xb = x+dx;
		int Y = 0;
		double xysum = 0.;
		for(double y = y0; y < y1; y += dx, Y++) // loop over all y values for a given x
		{
			double ya = y, yb = y+dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon poly;
			gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
			if(poly.num_contours == 0) continue; // if there are no observations in this direction

			skymap[std::make_pair(X, Y)] = poly; // store the polygon into a fast lookup map

			shared_ptr<Sy> sy(new Sy(xysum, Y));

			// calculate CMD integrals & marginal distributions
			double risum = ri_mpdf(sy->ripdf, (xa + xb)/2., (ya + yb)/2.);

			// store cumulative distribution
			sx->ypdf.push_back(sy);
			xysum += risum;

			++N;
//			cout << N << " " << X << " " << Y << ", " << risum << " stars\n";
		}
//		cout << "--- xysum = " << xysum << "\n";

		// convert starcounts in ypdf to marginal probabilities (prob. of finding a star w. y less then Y given x)
		if(xysum != 0.)
		{
			FOREACH(sx->ypdf) { (*i)->pcum /= xysum; }
			#if 1
			cout << "# N = " << N << "\n";
			FOREACH(sx->ypdf) { cout << X << " " << (*i)->Y << " " << (*i)->pcum << "\n"; }
			cout.flush();
			#endif
			xpdf.push_back(sx);
			xsum += xysum;
		}

//		static int kk = 0; if(++kk == 5) break;
	}
	cout << "# Total stars = " << xsum << " stars\n";

	// convert starcounts to cumulative probabilities
	FOREACH(xpdf) { (*i)->pcum /= xsum; }
	FOREACH(xpdf) { cout << (*i)->X << " -1 " << (*i)->pcum << "\n"; }

	// store all cumulative probability density maps to a binary file
	std::ofstream outf((prefix + ".pdf.dat").c_str());
	io::obstream out(outf);

	out << xsum;	// total number of stars (normalization) (see [1] below)
	out << xpdf;	// cumulative probability maps

	out << x0 << x1 << y0 << y1;	// sky footprint bounding rectangle
	out << skymap;			// skymap lookup table

	cout << io::binary::manifest << "\n";

	/* [1] The xsum number (and this whole mpdf tree) _is not_ the total number 
	of stars within the exact footprint. It's the total number of stars within 
	the "pixelized footprint" defined as the minimum set of pixels which completely 
	cover the whole polygonal sky footprint as given by skymap. This is intentional, 
	to make the monte-carlo drawings easier later, where a star can be drawn from
	a given pixel, it's coordinates within the pixels calculated, and then
	rejected if it's not within the footprint.
	If I calculated the mpdf for the exact footprint, we'd be in trouble when Monte
	Carlo drawing the exact coordinates of a star in a pixel containing an extremely
	small (and possibly very irregular!) survey polygon. It'd be unecesserily hard to randomly
	generate coordinates of a star within such a polygon.
	*/
}

void cumulative_dist::construct(const std::vector<double> &x, const std::vector<double> &y)
{
	ASSERT(y[0] == 0. && y[y.size()-1] == 1.)
	{
		std::cerr << y.size() << " " << y[0] << " " << y[y.size()-1] << "\n";
	}
	ASSERT(y.size() == x.size());
	hist.reserve(y.size());
	hist.clear();
	FOR(0, x.size())
	{
		hist.push_back(std::make_pair(x[i], y[i]));
	}
}

bool less_cdh(const double a, const cumulative_dist::value_type &b)
{
	return a < b.second;
}

double cumulative_dist::operator()(double prob)
{
	ASSERT(0 <= prob && prob <= 1);
	ITER(hist) i = upper_bound(hist.begin(), hist.end(), prob, less_cdh);
	
	if(i == hist.end()) { --i; }
	value_type b = *i; --i;
	value_type a = *i;
	
	double ret;
	if(b.second == a.second) { ret = (a.first + b.first) / 2.; }
	else { ret = a.first + (b.first - a.first)/(b.second - a.second) * (prob - a.second); }

	ASSERT(a.second <= prob && prob <= b.second) { std::cerr << a.second << " " << prob << " " << b.second << "\n"; }
	ASSERT(a.first <= ret && ret <= b.first) { std::cerr << prob << " " << a.first << " " << ret << " " << b.first << "\n"; }
	return ret;
}



typedef int (*int_function)(double dd, double const* y, double *dydd, void *param);

// integrate the function f from xv[0] to xv[xv.size()-1], storing the integral
// values at points xv in yv
bool
sample_integral2(const double a, const double b, double dx,
	std::vector<double> &xv, std::vector<double> &yv, int_function f, void *param,
	double abserr = 1e-6, double relerr = 0)
{
	using std::cout;

	gsl_odeiv_step*    s = gsl_odeiv_step_alloc(gsl_odeiv_step_rk2, 1);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(abserr, relerr);
	gsl_odeiv_evolve*  e = gsl_odeiv_evolve_alloc(1);

	double h = dx/2.;	// initial step size
	double y = 0;		// initial integration value
	double x = a;		// initial integration point

	xv.push_back(x);
	yv.push_back(y);

	ASSERT(VARRAY(double, param, 4) != 0);

	gsl_odeiv_system sys = {f, NULL, 1, param};

	// integration from 0 to dmax, but in ninterp steps
	// when integrating between dmin and dmax
	double xto = x+dx;
	bool last = false;
	while(!last)
	{
		if(xto > b) { xto = b; last = true; }

		while (x < xto)
		{
			//std::cout << "[" << VARRAY(double, param, 3) << "]: " << x << " " << y << " " << dx << " " << h << "\n";
			int status = gsl_odeiv_evolve_apply (e, c, s,
						&sys,
						&x, xto, &h,
						&y);
			if(status != GSL_SUCCESS)
			{
				std::cerr << "Error while integrating: " << gsl_strerror(status) << "\n";
				exit(-1);
			}

			xv.push_back(x);
			yv.push_back(y);
		}

		xto += dx;
	}
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);

	return true;
}


//
// Calculate a cumulative marginal pdf for a probability that a star has r-band magnitude
// less than r, for a given direction (x,y) and color r-i. Return the mpdf as a spline
// function in mspl object. Called from montecarlo().
//
#if 1
void sim::magnitude_mpdf(cumulative_dist &mspl, double x, double y, double ri)
{
	// deproject to (l,b)
	Radians l, b;
	proj.inverse(x, y, l, b);
	pencil_beam pb(l, b);

	// absolute magnitude for this spectral type
	double Mr = model->absmag(ri);

	double norm = 1.;
#if 0
	// calculate the total for normalization
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	const void *params[4] = { &pb, model, &Mr, &ri };
	gsl_function F = { &model_fun_1, params };
	double error;
	gsl_integration_qag(&F, m0, m1, 0, 1e-4, 1000, GSL_INTEG_GAUSS21, w, &norm, &error);
	gsl_integration_workspace_free(w);

	if(norm == 0) {
		std::cerr << "Norm = 0!: " << x << " " << y << " " << ri << "\n";
		exit(-1);
	}
#endif

	// estimate the number of magnitude bins
	int nm = (int)((m1 - m0) / dm) + 20;

	std::vector<double> m, f;
	m.reserve(nm); f.reserve(nm);
	m.push_back(m0); f.push_back(0);
	double dfmax = 0.01;

	const void *params2[5] = { &pb, model, &Mr, &ri, &norm };
	sample_integral2(m0, m1, dm, m, f, &model_fun, params2, 1e-10, 1e-5);

	// normalization - convert the integrated density to probability
	double fmax = f[f.size()-1];
	ASSERT(fmax > 0.);
#if 0
	ASSERT(fabs(fmax - 1) < 1e-1)
	{
		std::cerr << "norm = " << norm << " fmax=" << fmax << "\n";
		std::cerr << x << " " << y << " " << ri << "\n";
        	FOR(0, f.size()) { std::cout << m[i] << " " << f[i] << "\n"; }
		std::cout.flush();
	}
#endif
	FOREACH(f) { *i /= fmax; }
	//std::cerr << "total = " << total << " fmax=" << fmax << "\n";

// 	if(ri > 1.49)
// 	{
// 	FOR(0, f.size()) { std::cout << m[i] << " " << f[i] << "\n"; }
// 	exit(-1);
// 	}

	// construct spline approximation
//	std::cout << ri << " " << Mr << " " << m.size() << " " << f.size() << "\n";
	mspl.construct(m, f);

/*	FOR(0, f.size()) { std::cout << mspl(f[i]+0.001) << " " << f[i] << "\n"; }
	exit(-1);*/
/*	for(double f=0; f <= 1; f+=0.001) { std::cout << mspl(f) << " " << f << "\n"; std::cout.flush(); }
	exit(-1);*/
}
#else
void sim::magnitude_mpdf(cumulative_dist &mspl, double x, double y, double ri)
{
	// deproject to (l,b)
	Radians l, b;
	proj.inverse(x, y, l, b);
	pencil_beam pb(l, b);

	// absolute magnitude for this spectral type
	double Mr = model->absmag(ri);
	double dm = this->dm;

//	std::cout << x << " " << y << " " << l << " " << b << " " << Mr << "\n";

	// calculate the total for normalization
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	const void *params[4] = { &pb, model, &Mr, &ri };
	gsl_function F = { &model_fun_1, params };
	double total, error;
	gsl_integration_qag(&F, m0, m1, 0, 1e-4, 1000, GSL_INTEG_GAUSS21, w, &total, &error);
	//std::cout << "total = " << total << "\n";

	// estimate the number of magnitude bins
	int nm = (int)((m1 - m0) / dm) + 20;

	std::vector<double> m, f;
	m.reserve(nm); f.reserve(nm);
	m.push_back(m0); f.push_back(0);
	double dfmax = 0.01;
	double mcur = m0+dm;
	int k = 0;
	while(mcur < m1+dm)
	{
		if(mcur > m1) mcur = m1;

		const void *params[4] = { &pb, model, &Mr, &ri };
		gsl_function F = { &model_fun_1, params };

		// integrate the piece from mcur-dm to mcur, store to fcur
		double fcur;
		while(true)
		{
			gsl_integration_qag(&F, m.back(), mcur, 0, 1e-4, 1000, GSL_INTEG_GAUSS21, w, &fcur, &error);
			//std::cout << mcur-dm << " " << mcur << " " << fcur/total << "\n";
			if(fcur/total <= dfmax) break;
			dm /= 2;
			mcur -= dm;
			k = -3;
		}
		//std::cout << "OUT!\n";

		// store and cumulate
		m.push_back(mcur);
		f.push_back(f.back() + fcur);

		if(++k > 0) { dm *= 2; }
		mcur += dm;
	}
	gsl_integration_workspace_free(w);

	// normalization - convert the integrated density to probability
	double fmax = f[f.size()-1];
	FOREACH(f) { *i /= fmax; }
	ASSERT(abs(total - fmax)/total < 1e-5)
	{
		std::cerr << "total = " << total << " fmax=" << fmax << "\n";
	}
	//std::cerr << "total = " << total << " fmax=" << fmax << "\n";

//      	FOR(0, f.size()) { std::cout << m[i] << " " << f[i] << "\n"; }
//      	exit(-1);

	// construct spline approximation
	mspl.construct(m, f);
	
/*	FOR(0, f.size()) { std::cout << mspl(f[i]+0.001) << " " << f[i] << "\n"; }
	exit(-1);*/
/*	for(double f=0; f <= 1; f+=0.001) { std::cout << mspl(f) << " " << f << "\n"; std::cout.flush(); }
	exit(-1);*/
}
#endif

// ordering of Sx and Sy structures, by their cumulative marginal
// pdf values
template<typename S>
	struct less_S
	{
		bool operator()(const double& a, const boost::shared_ptr<S>& b)
			{ return a < b->pcum; }
	};

bool in_polygon(const gpc_vertex &t, const gpc_polygon &p);

// generate K stars in model-sky given in $prefix.pdf.dat
// file, precalculated with preprecalculate_mpdf
void sim::montecarlo(unsigned int K)
{
	std::map<std::pair<int, int>, gpc_polygon> skymap;	// a map of rectangular sections of the sky, for fast is-point-in-survey-area lookup
	XPDF xpdf; // marginal cumulative probability distribution function
	double N = 0; // Grand total - number of stars in the whole field

	// read the probability density maps and auxilliary data
	std::ifstream inf((prefix + ".pdf.dat").c_str());
	io::ibstream in(inf);
	in >> N;
	in >> xpdf;

	in >> x0 >> x1 >> y0 >> y1;
	in >> skymap;

	std::cerr << "# Generating " << K << " stars\n";

	// generate K stars in a two step process. First, draw (x,y) positions
	// from the distribution for all K stars. Then, sort the resulting array
	// by (x,y), and then for all stars in given (x,y) pixel draw ri and m.
	// This two-step process is necessary in order to avoid recalculating 
	// P(m|x,y,ri) distribution at each step (it's time-consuming to compute,
	// and there's not enough space to cache it for all pixels at once).
	int seed = 42;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, seed);
	std::vector<star> stars(K);
	FOR(0, K)
	{
		// pick x
		double ux = gsl_rng_uniform(rng);
		ITER(xpdf) ix = upper_bound(xpdf.begin(), xpdf.end(), ux, less_S<Sx>()); --ix;
		const Sx &X = **ix;

		// pick y, given an x
		double uy = gsl_rng_uniform(rng);
		ITER(X.ypdf) iy = upper_bound(X.ypdf.begin(), X.ypdf.end(), uy, less_S<Sy>()); --iy;
		const Sy &Y = **iy;

		// draw the exact location within the [x,x+dx], [y,y+dy] rectangle
		star &s = stars[i];
		s.x = x0 + dx*(X.X + gsl_rng_uniform(rng));
		s.y = y0 + dx*(Y.Y + gsl_rng_uniform(rng));

		// check that the star is inside survey footprint, reject if it's not
		gpc_polygon &poly = skymap[std::make_pair(X.X, Y.Y)];
		if(!in_polygon((gpc_vertex const&)std::make_pair(s.x, s.y), poly))
		{
			--i;
			continue;
		}

		// pick ri bin
		double uri = gsl_rng_uniform(rng);
		ITER(Y.ripdf) iri = upper_bound(Y.ripdf.begin(), Y.ripdf.end(), uri); --iri;
		int riidx = iri - Y.ripdf.begin();
		// pick ri within [ri, ri+dri)
		double ri = ri0 + dri*(riidx + gsl_rng_uniform(rng));
		ASSERT(ri0+riidx*dri <= ri && ri <= ri0+dri+riidx*dri);
		s.ri = ri;

//		ASSERT(Y.ripdf[riidx] != 0) { std::cerr << x0+dx*(s.X() + .5) << " " << y0+dx*(s.Y() + .5) << " " << ri0+dri*(s.RI()+.5) << "\n"; }
		if(std::abs(x0+dx*(s.X() + .5) + 0.0912734) < 0.0001 &&
		   std::abs(y0+dx*(s.Y() + .5) + 0.0563668) < 0.0001 &&
		   std::abs(ri0+dri*(s.RI() + .5) - .215) < 0.0001)
		   {
		    { std::cerr << x0+dx*(s.X() + .5) << " " << y0+dx*(s.Y() + .5) << " " << ri0+dri*(s.RI()+.5) << "\n"; }
		    std::cerr << riidx << " " << Y.ripdf[riidx] << " " << Y.ripdf[riidx+1] << "\n";
//		    exit(-1);
		   }
	}

	// sort stars by X, Y, ri
	std::sort(stars.begin(), stars.end());

	std::cerr << "# Assigning magnitudes \n";

	ticker tick(1000);
	// draw magnitudes for stars in each cell X,Y
	boost::tuple<int, int, int> cell = boost::make_tuple(-1000,-1000,-1000);
	cumulative_dist mspl;
	ITER(stars) it = stars.begin();
	while(it != stars.end())
	{
		//std::cerr << "+" << "\n";
		star &s = *it++;
		boost::tuple<int, int, int> cell2(s.X(), s.Y(), s.RI());
		// if we moved to a different (x, y, ri) bin, calculate
		// the cumulative marginal pdf in magnitude for that bin
		if(cell != cell2)
		{
			magnitude_mpdf(mspl, x0 + (s.X() + 0.5)*dx, y0 + (s.Y() + 0.5)*dx, ri0 + (s.RI() + 0.5)*dri);
			cell = cell2;
/*			double rix = ri0 + (s.RI() + 0.5)*dri;
			if(rix == 0.935) {
				static int cnt = 10;
				std::cerr << rix << "\t" << mspl(0) << "\t" << mspl(0.999) << "\n";
				if(--cnt == 0) {
					for(double f=0; f <= 1; f+=0.001) { std::cout << mspl(f) << " " << f << "\n"; std::cout.flush(); }
					exit(-1);
				}
			}*/
		}
		// draw the magnitude
		double um = gsl_rng_uniform(rng);
		s.m = mspl(um);
//		std::cout << um << " " << s.m << "\n";
		tick.tick();
	}

	// output
#if 1
	// simple ascii-text dump
	std::ofstream out((prefix + ".sim.txt").c_str());
	FOR(0, K)
	{
		star &s = stars[i];
		out << s.x << " " << s.y << " " << s.ri << " " << s.m << "\n";
	}
#else
	// NOTE: recheck this code!!
	unsigned int t = time(NULL), tstart = t;
	DMMArray<mobject> out;
	DMMArray<starmag> starmags;
	out.create("sim_stars.dmm");
	starmags.create("sim_starmags.dmm");
	std::vector<starmag> obsv;
	double Ar = .1;
	FORj(j, 0, K)
	{
		star &s = stars[j];
		int uniqId = out.size();
		ASSERT(uniqId == j);

		starmag sm;
		sm.fitsId = j;

		// calculate magnitudes and convert to flux
		float u, g, r, i, z;
		float uErr, gErr, rErr, iErr, zErr;
		u = z = 21.0;
		r = s.m;
		i = r - s.ri;
		double gr = factory.plx.gr(s.ri);
		g = r + gr;
		uErr = gErr = rErr = iErr = zErr = 0.02;
		// add extinction
/*		u += extinction(Ar, 0);
		g += extinction(Ar, 1);
		r += extinction(Ar, 2);
		i += extinction(Ar, 3);
		z += extinction(Ar, 4);*/
		// convert to flux
// 		phot::fluxfromluptitude(g, gErr, 1, sm.mag[0], sm.magErr[0], 1e9);
// 		phot::fluxfromluptitude(r, rErr, 2, sm.mag[1], sm.magErr[1], 1e9);
// 		phot::fluxfromluptitude(i, iErr, 3, sm.mag[2], sm.magErr[2], 1e9);
// 		phot::fluxfromluptitude(z, zErr, 4, sm.mag[3], sm.magErr[3], 1e9);
// 		phot::fluxfromluptitude(u, uErr, 0, sm.mag[4], sm.magErr[4], 1e9);
		sm.mag[0] = g;
		sm.mag[1] = r;
		sm.mag[2] = i;
		sm.mag[3] = z;
		sm.mag[4] = u;
		FOR(0, 5) { sm.mag[(i+4)%5] += extinction(Ar, i); }
		FOR(0, 5) { sm.magErr[i] = 0.2; }

		// store to database
		double l, b, ra, dec;
		proj.inverse(s.x, s.y, l, b);
		coordinates::galequ(l, b, ra, dec);
		ra = deg(ra); dec = deg(dec);

		using std::cout;
// 		cout << g << "\t" << sm.mag[0] << "\t" << sm.magErr[0] << "\n";
// 		cout << r << "\t" << sm.mag[1] << "\t" << sm.magErr[1] << "\n";
// 		cout << i << "\t" << sm.mag[2] << "\t" << sm.magErr[2] << "\n";
// 		cout << ra << "\t" << dec << "\n";
		
		obsv.clear();
		obsv.push_back(sm);
		mobject m = process_observations(starmags.size(), ra, dec, Ar, obsv);
		
// 		cout << g << "\t" << m.mag[0] << "\t" << m.magErr[0] << "\n";
// 		cout << r << "\t" << m.mag[1] << "\t" << m.magErr[1] << "\n";
// 		cout << i << "\t" << m.mag[2] << "\t" << m.magErr[2] << "\n";
// 		cout << m.ra << "\t" << m.dec << "\n";
//		print_mobject(cout, m);
		
		// store to database			
		FOREACH(obsv) { starmags.push_back(*i); }
		out.push_back(m);
		if((out.size() % 500000) == 0) { out.sync(); }

		// user interface stuff
		tick.tick();
		if((j % 10000) == 0)
		{
			unsigned int tnow = time(NULL);
			cout << tnow - t << "s [ " << (tnow - tstart)/60 << "min ] ";
			if(j > 20000) { cout << "[ " << io::format("% 5.2f") << ((double(K) / double(j) -1.) * (tnow - tstart)/3600.) << "hrs left]";}
			cout << "\n";
			t = tnow;
			cout.flush();
		}
	}
#endif
}

#if 0

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <valarray>

// template lapack_gb_matrix<typename T>
// 	class lapack_gb_matrix : public banded_matrix<T, boost::numeric::ublas::column_major>
// 	{
// 		
// 	}

namespace lapack {
// 	int dgbsv_(integer *n, integer *kl, integer *ku, integer *
// 	nrhs, doublereal *ab, integer *ldab, integer *ipiv, doublereal *b, 
// 	integer *ldb, integer *info);
	inline int solve(
		boost::numeric::ublas::vector<doublereal> &xy,
		boost::numeric::ublas::banded_matrix<doublereal> &m, 
		boost::numeric::ublas::vector<integer> &ipiv)
	{
		ASSERT(m.size1() == m.size2());
		ASSERT(xy.size() == m.size1());

		integer n = m.size1(), kl = m.lower(), ku = m.upper()-kl, nrhs = 1, ldab = 2*kl+ku+1;
/*		cout << n << " " << kl << " " << ku << " " << nrhs << " " << ldab << "\n";
		exit(-1);*/
		integer info;
		dgbsv_(&n, &kl, &ku, &nrhs, &m.data()[0], &ldab, &ipiv[0], &xy[0], &n, &info);
		return info;
	}

	inline int solve(
		boost::numeric::ublas::vector<doublereal> &xy,
		boost::numeric::ublas::banded_matrix<doublereal> &m)
	{
		boost::numeric::ublas::vector<integer> ipiv(m.size1());
		return solve(xy, m, ipiv);
	}
}

void prettyprint(double *v, int size)
{
	FOR(0, size)
	{
		if(i % 5 == 0) { cout << "\n"; }
		if(v[i] == 0) { cout << " * "; } else { cout << v[i] << " ";}
	}
	cout << "\n";
}

#include <iomanip>
void mprint(boost::numeric::ublas::banded_matrix<double> &m, boost::numeric::ublas::vector<double> *v = NULL)
{
	using namespace std;
	int n = m.size1(), kl = m.lower(), ku = m.upper();//-kl;
	FORj(r, 0, n)
	{
		FORj(c, 0, n)
		{
			r++; c++;
			if(std::max(1, c-ku) <= r && r <= std::min(n, c+kl))
			{
				cerr << setw(3) << m(r-1, c-1);
			} else {
				cerr << "  *";
			}
			r--; c--;
		}
		if(v != NULL)
		{
			cerr << setw(10) << (*v)[r];
		}
		cerr << "\n";
	}
	cerr << "\n";
}

void fit_int_spline();
void test_lapack()
{
 	fit_int_spline();
 	return;

	using namespace boost::numeric;
	using namespace lapack;
	using namespace std;

	integer n = 5, kl = 2, ku = 1, nrhs = 1;
	integer ldab = 2*kl+ku+1;
	ublas::banded_matrix<doublereal> m (n, n, kl, ku+kl);
	m(0, 0) = 11; m(0, 1) = 12;
	m(1, 0) = 21; m(1, 1) = 22; m(1, 2) = 23;
	m(2, 0) = 31; m(2, 1) = 32; m(2, 2) = 33; m(2, 3) = 34;
	              m(3, 1) = 42; m(3, 2) = 43; m(3, 3) = 44; m(3, 4) = 45;
	                            m(4, 2) = 53; m(4, 3) = 54; m(4, 4) = 55;
	double* v = &m.data()[0];
	prettyprint(v, m.data().size());
	mprint(m);
	ublas::vector<long int> ipiv(n);
	ublas::vector<double> rhs(n);
	FOR(0, 5) { rhs[i] = (double)i+1; }
	FOR(0, 5) { cerr << rhs[i] << "\t" << ipiv[i] << "\n"; }
	int info = lapack::solve(rhs, m, ipiv);

	cerr << "INFO = " << info << "\n";

	//cout << m << "\n";
	prettyprint(v, m.data().size()); cerr << "\n";
	mprint(m);
	FOR(0, 5) { cerr << rhs[i] << "\t" << ipiv[i] << "\n"; }
}

// double ib[] =
// {
// 	32, 27, 28, 15, 32, 8, 4, 1,
// 	1, 1, 1, 1, -1, 0, 0, 0,
// 	0, 1, 2, 3, 0, -1, 0, 0,
// 	0, 0, 2, 6, 0, 0, -2, 0
// };
// 
// double ib0[] =
// {
// 	32, 8, 4, 1,
// 	0, 1, 0, 0
// };
// 
// double ib1[] =
// {
// 	32, 27, 28, 15,
// 	0, 1, 2, 3
// };

void copy_block(boost::numeric::ublas::banded_matrix<double> &m, int r0, int c0, double *b, int w, int h)
{
	--b;
	FORj(i, 0, h)
		FORj(j, 0, w)
			if(*(++b) != 0.)
				m(r0+i, c0+j) = *b;
}

#if 0

double ib[] =
{
	1, 1, 1, 1, -1,  0,  0, 0,
	0, 1, 2, 3,  0, -1,  0, 0,
	0, 0, 1, 3,  0,  0, -1, 0,
	0, 0, 0, 0, 12,  6,  4, 3
};

double ib0[] =
{
	 0, 1, 0, 0,
	 0, 0, 1, 0,
	12, 6, 4, 3
};

double ib1[] =
{
	0, 1, 2, 3
};

#include <algorithm>

class ispline
{
public:
	double eval(double t, double *c)
	{
		return (((t*c[3]) + c[2])*t + c[1])*t + c[0];
	}
public:
	int nseg;
	std::valarray<double> x, h, c;
public:
	int find_x(double xv)
	{
		using namespace std;
		double *x = &this->x[0];
		int n = this->x.size();
		double *v = upper_bound(x, x+n, xv);
		--v;
		if(v-x >= n-1) { --v; } // if it's beyond the last segment, extrapolate
		return v-x;
	};

	double eval(double xv)
	{
		int i0 = find_x(xv);
		//cout << i0 << " ";
		double dx = x[i0+1] - x[i0];
		double t = (xv - x[i0])/dx;
		return eval(t, &c[4*i0]);
	};

	ispline(int np, double *vx, double *vy, double *coef)
		: nseg(np-1), x(vx, np), h(vy, np-1), c(coef, 4*(np-1))
	{}
};

void fit_int_spline()
{
	using namespace boost::numeric;
	using namespace lapack;
	using namespace std;

	valarray<double> x(6), h(x.size()-1);
	x[0] = 0; h[0] = 1;
	x[1] = 1; h[1] = 2;
	x[2] = 2; h[2] = 4;
	x[3] = 3; h[3] = 4;
	x[4] = 4; h[4] = 5;
	x[5] = 5;

	// construct and fill in the spline matrix and rhs vector
	int ns = h.size();	// number of segments	cout << "

	int n = 4*ns;
	ublas::vector<double> rhs(n); rhs *= 0.;
	ublas::banded_matrix<double> m(n, n, 3, 1 + 3);

	copy_block(m, 0, 0, ib0, 4, 3);
//	rhs[0] = 2.*(h[1]-h[0])/(x[2]-x[0]);	// D1
	rhs[2] = 12. * h[0] / (x[1] - x[0]);	// H1

	FOR(1, ns)
	{
		cerr << 4*i-1 << "\n";
		copy_block(m, 4*i-1, 4*(i-1), ib, 8, 4);
		rhs[4*i+2] = 12. * h[i] / (x[i] - x[i-1]);
	}

	copy_block(m, n-1, n-4, ib1, 4, 1);
//	rhs[n-1] = 2 * (h[n-1]-h[n-2]) / (x[n] - x[n-2]); // Dn
	mprint(m, &rhs);
	
	// solve spline matrix
 	int info = lapack::solve(rhs, m);
 	cerr << "LAPACK info = " << info << "\n\n";
 	mprint(m, &rhs);

	// repack to ispline data structure
	ispline isp(x.size(), &x[0], &h[0], &rhs[0]);
	for(double xx = 0; xx <= x[x.size()-1]+0.1; xx += 0.1)
	{
		cout << xx << "\t" << isp.eval(xx) << "\n";
	}
}

#else

double ib[] =
{
	1, 1, 1, -1, 0, 0,
	0, 1, 2, 0, -1, 0,
	0, 0, 0, 6,  3, 2
};

double ib0[] =
{
	 0, 1, 0,
	 6, 3, 1
};

double ib1[] =
{
	0, 1, 2
};

#include <algorithm>

class ispline
{
public:
	double eval(double t, double *c)
	{
		return ((c[2])*t + c[1])*t + c[0];
	}
public:
	int nseg;
	std::valarray<double> x, h, c;
public:
	int find_x(double xv)
	{
		using namespace std;
		double *x = &this->x[0];
		int n = this->x.size();
		double *v = upper_bound(x, x+n, xv);
		--v;
		if(v-x >= n-1) { --v; } // if it's beyond the last segment, extrapolate
		return v-x;
	};

	double eval(double xv)
	{
		int i0 = find_x(xv);
		//cout << i0 << " ";
		double dx = x[i0+1] - x[i0];
		double t = (xv - x[i0])/dx;
		return eval(t, &c[3*i0]);
	};

	ispline(int np, double *vx, double *vy, double *coef)
		: nseg(np-1), x(vx, np), h(vy, np-1), c(coef, 3*(np-1))
	{}
};

///////////////////////////////////

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <valarray>

class spline
{
private:
	gsl_interp *f;
	gsl_interp_accel *acc;
	std::valarray<double> xv, yv;
public:
	spline() : f(NULL), acc(NULL) {}
	spline(const double *x, const double *y, int n);
	void construct(const double *x, const double *y, int n);
	~spline();

 	double operator ()(double x)        { return gsl_interp_eval(f, &xv[0], &yv[0], x, acc); }
 	double deriv(double x)              { return gsl_interp_eval_deriv(f, &xv[0], &yv[0], x, acc); }
 	double deriv2(double x)             { return gsl_interp_eval_deriv2(f, &xv[0], &yv[0], x, acc); }
 	double integral(double a, double b) { return gsl_interp_eval_integ(f, &xv[0], &yv[0], a, b, acc); }

public:
	spline& operator= (const spline& a);
	spline(const spline& a) : f(NULL), acc(NULL) { *this = a; }
};

///////////////////////////////////

void fit_int_spline2(std::valarray<double> &x, std::valarray<double> &h)
{
	double sum = 0;
	FOR(0, h.size()-1)
	{
		h[i] = sum;
		sum += h[i+1];
		cerr << sum << "\n";
	}
	h[h.size()-1] = sum;
	cerr << x.size() << " " << h.size() << "\n";

	spline spl;
	spl.construct(&x[0], &h[0], h.size());
	for(double xx = 0; xx <= x[x.size()-1]+0.1; xx += 0.1)
	{
		cout << xx << "\t" << spl.deriv(xx) << "\n";
	}
	exit(0);
}

void fit_int_spline()
{
	using namespace boost::numeric;
	using namespace lapack;
	using namespace std;

	valarray<double> x(6), h(x.size()-1);
	x[0] = 0; h[0] = 1;
	x[1] = 1; h[1] = 2;
	x[2] = 2; h[2] = 4;
	x[3] = 3; h[3] = 4;
	x[4] = 4; h[4] = 5;
	x[5] = 5;

	fit_int_spline2(x, h);
	
	// construct and fill in the spline matrix and rhs vector
	int ns = h.size();	// number of segments	cout << "

	int n = 3*ns;
	ublas::vector<double> rhs(n); rhs *= 0.;
	ublas::banded_matrix<double> m(n, n, 2, 1 + 3);

	copy_block(m, 0, 0, ib0, 3, 2);
//	rhs[0] = 2.*(h[1]-h[0])/(x[2]-x[0]);	// D1
	rhs[1] = 6. * h[0] / (x[1] - x[0]);	// H1

	FOR(1, ns)
	{
		copy_block(m, 3*i-1, 3*(i-1), ib, 6, 3);
		rhs[3*i+1] = 6. * h[i] / (x[i] - x[i-1]);
	}

	copy_block(m, n-1, n-3, ib1, 3, 1);
//	rhs[n-1] = 2 * (h[n-1]-h[n-2]) / (x[n] - x[n-2]); // Dn
	mprint(m, &rhs);
	
	// solve spline matrix
 	int info = lapack::solve(rhs, m);
 	cerr << "LAPACK info = " << info << "\n\n";
 	mprint(m, &rhs);

	// repack to ispline data structure
	ispline isp(x.size(), &x[0], &h[0], &rhs[0]);
	for(double xx = 0; xx <= x[x.size()-1]+0.1; xx += 0.1)
	{
		cout << xx << "\t" << isp.eval(xx) << "\n";
	}
}

#endif
#endif

void test_lapack()
{
}

#endif // COMPILE_SIMULATE_X
