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

#if 0
namespace lapack
{
	extern "C"
	{
		#include <g2c.h>
		#include "common/clapack.h"
	}
}
#endif

//#include <nurbs++/nurbsS.h>
#include "gpc_cpp.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "simulate.h"
#include "projections.h"
#include "model.h"
#include "paralax.h"
#include "analysis.h"
#include "dm.h"

#include <vector>
#include <map>
#include <string>

#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

void poly_bounding_box(double &x0, double &x1, double &y0, double &y1, const gpc_polygon &p);

gpc_polygon poly_rect(double x0, double x1, double y0, double y1);

model_pdf::model_pdf()
{
}

// 	lambert proj(rad(90), rad(90));
// 	double dx = rad(1); // model grid angular resolution
// 	double dri = 0.01; // model CMD resolution
// 	double ri0 = 0.0, ri1 = 1.5;	// color limits
// 	double m0 = 14, m1 = 22;	// magnitude limits
// 	double dm = 0.1;	// model CMD magnitude resolution
// 	double x0, x1, y0, y1;	// lambert survey footprint bounding box
model_pdf::model_pdf(std::istream &in)
: proj(rad(90), rad(90))
{
	Config cfg; cfg.load(in);

	ASSERT(cfg.count("name"));
	pdfname = cfg["name"];

	ASSERT(cfg.count("footprint"));
	footprint = cfg["footprint"];

	cfg.get(skymap.dx,	"dx", 	1.);
	skymap.dx = rad(skymap.dx);
	ASSERT(cfg.count("dri")); cfg.get(dri,	"dri",	0.01);
	cfg.get(ri0,	"ri0",	0.0);
	cfg.get(ri1,	"ri1",	1.5);
	cfg.get(m0,	"r0",	14.);
	cfg.get(m1,	"r1",	22.);
	cfg.get(dm,	"dm",	0.1);

	std::string pole; double l0, b0;
	cfg.get(pole,	"projection",	std::string("90 90"));
	std::istringstream ss(pole);
	ss >> l0 >> b0;
	proj = lambert(rad(l0), rad(b0));

	// load the model
	ASSERT(cfg.count("model"));
	std::ifstream inconf(cfg["model"].c_str()); ASSERT(inconf);
	model.reset(galactic_model::load(inconf));
	ASSERT(model.get() != NULL);

	LOG(app, verb1) << "footprint  = " << footprint;
	LOG(app, verb1) << "dx  = " << deg(skymap.dx);
	LOG(app, verb1) << "dri = " << dri;
	LOG(app, verb1) << "ri0 = " << ri0;
	LOG(app, verb1) << "ri1 = " << ri1;
	LOG(app, verb1) << "r0  = " << m0;
	LOG(app, verb1) << "r1  = " << m1;
	LOG(app, verb1) << "dm  = " << dm;
	LOG(app, verb1) << "projection pole = " << l0 << " " << b0;
}

struct star_lessthan
{
	const model_pdf *gsim;
	star_lessthan(const model_pdf *gsim_) : gsim(gsim_) {}

	bool operator()(const model_pdf::star &a, const model_pdf::star &b) const
	{
		int aX = a.X(gsim->skymap);
		int bX = b.X(gsim->skymap);
		int aY = a.Y(gsim->skymap);
		int bY = b.Y(gsim->skymap);
		return aX < bX || (aX == bX && (
			aY < bY || (aY == bY && 
			a.ri < b.ri)
			));
	}
};

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
	
	//std::cerr << "FUN: " << Mr << " " << ri << " " << m << " " << d << " " << rho << "\n";
	
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
//			std::cerr << h << "\n";
//			std::cout << "[" << VARRAY(double, param, 4) << "]: " << x << " " << y << " " << dx << " " << h << " " << rho << "\n";
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

//
// Populates the pdf vector with cumulative marginal PDFs of stars within a given
// field (x,y).
//
double model_pdf::ri_mpdf(std::vector<double> &pdf, const double x, const double y)
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
	int nspt = (int)rint((ri1 - ri0)/dri);
//	for(double ri = ri0; ri <= ri1; ri += dri)
	FOR(0, nspt)
	{
		double ri = ri0 + i*dri;
		// get a model which will return the number of stars n in
		// direction (l,b) in the SpT interval [ri, ri+dri) (approximate
		// this by assuming a constant f=dN/dri in the short interval dri
		// and taking n = f*dri
		double ri_mid = ri + 0.5*dri;
		double Mr = model->absmag(ri_mid); // absolute magnitude for the stars with this SpT

		double norm = 1.;
		void *params[5] = { &pb, model.get(), &Mr, &ri_mid, &norm };

		double n, error;
		gsl_function F = { &model_fun_1, params };
		#if 0
		ASSERT(gsl_integration_qag (&F, m0, m1, 0, 1e-7, Npts, GSL_INTEG_GAUSS61, w, &n, &error) == GSL_SUCCESS);
		#else
		n = integrate(m0, m1, dm, model_fun, params, 1e-10, 1e-5);
		#endif
		//ASSERT(n != 0) { std::cerr << "BU! " << x << " " << y << " " << ri_mid << " " << n << "\n"; }
		n *= sqr(skymap.dx)*dri; // multiply by solid angle. n is now the total number of stars in direction (l,b)
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

BOSTREAM2(const Sy &p) { return out << p.pcum << p.Y << p.ripdf; }
BISTREAM2(Sy &p) { return in >> p.pcum >> p.Y >> p.ripdf; }
BOSTREAM2(const Sx &p) { return out << p.pcum << p.X << p.ypdf; }
BISTREAM2(Sx &p) { return in >> p.pcum >> p.X >> p.ypdf; }

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

void make_skymap(partitioned_skymap &m, Radians dx, const std::string &skypolyfn)
{
	gpc_polygon sky;
#if 1
	// north.gpc.txt
	std::cerr << "Reading " << skypolyfn << "... ";
 	FILE *fp = fopen(skypolyfn.c_str(), "r");
 	gpc_read_polygon(fp, 1, &sky);
	fclose(fp);
	std::cerr << "done.\n";
#else // while debugging - just a circular cutout at the pole
	sky = make_circle(0., 0., 0.1, dx);
	//sky = make_circle(-0.123308, 0.406791, 0.1, dx);
#endif
	m.skymap.clear();

 	poly_bounding_box(m.x0, m.x1, m.y0, m.y1, sky);
	m.dx = dx;

	std::cerr << "BBox: " << m.x0 << " " << m.y0 << "     " << m.x1 << " " << m.y1 << "\n";
	std::cerr << "dx: " << m.dx << "\n";
	std::cerr << "area: " << polygon_area(sky)*sqr(deg(1)) << "\n";

	double area = 0;
	int X = 0;
	ticker tick(1);
	for(double x = m.x0; x < m.x1; x += m.dx, X++) // loop over all x values in the bounding rectangle
	{
		double xa = x, xb = x+dx;
		int Y = 0;
		for(double y = m.y0; y < m.y1; y += m.dx, Y++) // loop over all y values for a given x
		{
			double ya = y, yb = y+dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon poly;
			gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
			if(poly.num_contours == 0) continue; // if there are no observations in this direction

			m.skymap[std::make_pair(X, Y)].poly = poly; // store the polygon into a fast lookup map
			m.skymap[std::make_pair(X, Y)].area = polygon_area(poly);
			area += m.skymap[std::make_pair(X, Y)].area;
			tick.tick();
		}
	}
	tick.close();

	std::cerr << "Total area: " << area*sqr(deg(1)) << " deg^2\n";
}

gpc_polygon make_circle(double x0, double y0, double r, double dx);

// calculate marginal PDFs (``mpfds'') for stars within a given sky footprint,
// and with a given Galactic model. Store the mPDFs into an XPDF tree
// and serialize it to a file. These mpdfs are later used in generating
// Monte Carlo realisations of the sky
//
// This routine also calculates the 'skymap', for fast point-in-polygon
// lookups of the observed sky footprint.
void model_pdf::precalculate_mpdf()
{
	using boost::shared_ptr;
	using std::cout;

	gpc_polygon sky;
#if 1
 	FILE *fp = fopen(footprint.c_str(), "r");
	ASSERT(fp != NULL) { std::cerr << "Could not open footprint file" << footprint << "\n"; }
	gpc_read_polygon(fp, 1, &sky);
	fclose(fp);
#else // while debugging - just a circular cutout at the pole
	sky = make_circle(0., 0., 0.1, dx);
	//sky = make_circle(-0.123308, 0.406791, 0.1, dx);
#endif

#if 0
	std::map<std::pair<int, int>, gpc_polygon> skymap;	// a map of rectangular sections of the sky, for fast is-point-in-survey-area lookup

	XPDF xpdf;
	int N = 0; // number of dx^2 cells processed
#else
	skymap.skymap.clear();
	xpdf.clear();
	int N = 0; // number of dx^2 cells processed
#endif
 	poly_bounding_box(skymap.x0, skymap.x1, skymap.y0, skymap.y1, sky);

	double xsum = 0; // Grand total - number of stars in the whole field
	int X = 0;
	for(double x = skymap.x0; x < skymap.x1; x += skymap.dx, X++) // loop over all x values in the bounding rectangle
	{
		std::cerr << "#";
		shared_ptr<Sx> sx(new Sx(xsum, X));

		double xa = x, xb = x+skymap.dx;
		int Y = 0;
		double xysum = 0.;
		for(double y = skymap.y0; y < skymap.y1; y += skymap.dx, Y++) // loop over all y values for a given x
		{
//			std::cerr << "-";
			double ya = y, yb = y+skymap.dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon poly;
			gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
			if(poly.num_contours == 0) continue; // if there are no observations in this direction

			// store the polygon into a fast lookup map
			partitioned_skymap::pixel_t &pix = skymap.skymap[std::make_pair(X, Y)];
			pix.poly = poly;
			pix.area = polygon_area(poly);

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
	cout << "# Total stars = " << xsum << " stars within " << skymap.skymap.size()*sqr(deg(skymap.dx)) << "deg^2 of pixelized area." << "\n";

	// convert starcounts to cumulative probabilities
	FOREACH(xpdf) { (*i)->pcum /= xsum; }
	FOREACH(xpdf) { cout << (*i)->X << " -1 " << (*i)->pcum << "\n"; }

	// store the normalization (it's unfortunate that the variable name in this method
	// is different from the one in the object :-( This needs fixing.).
	this->N = xsum;
}

io::obstream &model_pdf::serialize(io::obstream &out) const
{
	// store all cumulative probability density maps to a binary file
	//std::ofstream outf((prefix + ".pdf.dat").c_str());
	//io::obstream out(outf);

	out << N;	// total number of stars (normalization) (see [1] below)
	out << xpdf;	// cumulative probability maps

	out << proj;
	out << skymap;
	
	out << ri0 << ri1 << dri;
	out << m0 << m1 << dm;

	model->serialize(out);

	return out;
	
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
//			std::cout << "[" << VARRAY(double, param, 3) << "]: " << x << " " << y << " " << dx << " " << h << "\n";
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
void model_pdf::magnitude_mpdf(cumulative_dist &mspl, double x, double y, double ri)
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

	const void *params2[5] = { &pb, model.get(), &Mr, &ri, &norm };
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
void model_pdf::magnitude_mpdf(cumulative_dist &mspl, double x, double y, double ri)
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

peyton::io::ibstream &model_pdf::unserialize(peyton::io::ibstream &in)
{
	// read the probability density maps and auxilliary data
	in >> N;
	in >> xpdf;

	in >> proj;
	in >> skymap;

	in >> ri0 >> ri1 >> dri;
	in >> m0 >> m1 >> dm;

#if 1
	model.reset(galactic_model::unserialize(in));
#else
	// TODO: This is just a quick fix for debugging - the model should be serialized
	// together with the PDF
	std::ifstream inconf("galaxy.conf"); ASSERT(inconf);
	model.reset(galactic_model::load(inconf));
	ASSERT(model.get() != NULL);
#endif

	std::cerr << "LOAD: N = " << N << "\n";
	std::cerr << "LOAD: xpdf.size() = " << xpdf.size() << "\n";
	std::cerr << "LOAD: skymap.skymap.size() = " << skymap.skymap.size() << "\n";

	return in;
}

sky_generator::sky_generator(std::istream &in)
: rng(NULL), nstars(0), flags(0)
{
	Config cfg; cfg.load(in);

	// number of stars to read
	cfg.get(nstars,	"nstars", 	0);
	cfg.get(Ar,	"Ar", 		0.);	// extinction

	// random number generator
	int seed;
	cfg.get(seed,	"seed", 	42);
	rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, seed);

	// load pdfs
	ASSERT(cfg.count("pdfs"));
	std::istringstream ss(cfg["pdfs"].c_str());
	std::string pdfFn;
	while(ss >> pdfFn)
	{
		std::ifstream iin(pdfFn.c_str()); ASSERT(iin) { std::cerr << "Could not open " << pdfFn << "\n"; }
		io::ibstream in(iin);

		boost::shared_ptr<model_pdf> pdf(new model_pdf);
		if(!(in >> *pdf)) { ASSERT(0) { std::cerr << "Error loading " << pdfFn << "\n"; } break; }

		add_pdf(pdf);
	}

	// load magnitude error map
	if(cfg.count("photo_error_map"))
	{
		std::string fn = cfg["photo_error_map"];
		std::ifstream iin(fn.c_str()); ASSERT(iin) { std::cerr << "Could not open magnitude error map " << fn << "\n"; }
		io::ibstream in(iin);

		if(!(in >> magerrs_cvm)) { ASSERT(in >> magerrs); }
		magerrs.initialize(magerrs_cvm);
		
		cfg.get(magerrs.use_median_beam_for_all, "magerr_only_use_median", false);
		if(magerrs.use_median_beam_for_all) { std::cerr << "Using median error only.\n"; }
	}

	cfg.get(constant_photo_error, "constant_photo_error", 0.);
	cfg.get(paralax_dispersion, "paralax_dispersion", 0.);

	// load various flags
	int flag;
	cfg.get(flag, "apply_photo_errors", 0); 	if(flag) flags |= APPLY_PHOTO_ERRORS;
	
	// dump short status
	std::cerr << "constant_photo_error = " << constant_photo_error << "\n";
	std::cerr << "paralax_dispersion = " << paralax_dispersion << "\n";
	std::cerr << "magerrs.use_median_beam_for_all = " << magerrs.use_median_beam_for_all << "\n";
	std::cerr << "magerrs.empty() = " << magerrs.empty() << "\n";
	std::cerr << "Flags: "
		<< (flags & APPLY_PHOTO_ERRORS ? "APPLY_PHOTO_ERRORS" : "")
		<< "\n";
}

sky_generator::~sky_generator()
{
	if(rng != NULL) { gsl_rng_free(rng); }
}

struct null_deleter
{
	void operator()(void const *) const
	{
	}
};

// for adding static objects
void sky_generator::add_pdf(model_pdf &pdf)
{
	boost::shared_ptr<model_pdf> px(&pdf, null_deleter());
	return add_pdf(px);
}

void sky_generator::add_pdf(boost::shared_ptr<model_pdf> &pdf)
{
	pdfs.push_back(pdf);

	boost::shared_ptr<std::vector<model_pdf::star> > v(new std::vector<model_pdf::star>);
	stars.push_back(v);
}

// generate K stars in model-sky given in $prefix.pdf.dat
// file, precalculated with preprecalculate_mpdf
void sky_generator::montecarlo(star_output_function &out)
{
	ASSERT(pdfs.size() != 0);

	bool allowMisses = false;
	int Ktotal = nstars;

	if(Ktotal == 0)
	{
		double Kd = 0;
		std::cerr << "# " << (*pdfs.begin())->N << "\n";
		FOREACH(pdfs) { Kd += (*i)->N; }
		Ktotal = (int)rint(Kd);
		allowMisses = true;
	}

	std::cerr << "# Generating " << Ktotal << " stars";
	if(allowMisses) { std::cerr << " (approximately, based on model norm.)"; }
	std::cerr << ".\n";

	// calculate cumulative probability function of a star being a part
	// various pdfs. This will be used to decide for each star according
	// to which model_pdf are we going to draw it.
	std::vector<double> modelCPDF(pdfs.size());
	modelCPDF[0] = 0.;
	FOR(0, pdfs.size()-1) { modelCPDF[i+1] = modelCPDF[i] + pdfs[i]->N; }
	double norm = pdfs.back()->N + modelCPDF.back();
	FOR(1, pdfs.size()) { modelCPDF[i] /= norm; }
	
	// generate the data in smaller batches (to conserve memory)
	static const int Kbatch = 20000000;
//	static const int Kbatch = 10000;

	// prepare output
	sstruct::factory.useTag("lonlat[2]");
	sstruct::factory.useTag("color");
	sstruct::factory.useTag("mag");
	FOR(0, pdfs.size())
	{
		pdfs[i]->galmodel().setup_tags(sstruct::factory);
	}
	out.output_header(sstruct::factory);

	// do catalog generation
	int K = 0, Ngen = 0;
	while(K < Ktotal)
	{
		int Kstep = std::min(Kbatch, Ktotal - K);
		K += Kstep;
		Ngen += montecarlo_batch(out, Kstep, modelCPDF, allowMisses);
	}
	
	std::cerr << "# Generated " << Ngen << " stars.\n";
}

int sky_generator::montecarlo_batch(star_output_function &out, int K, const std::vector<double> &modelCPDF, bool allowMisses)
{
	// clear output vectors
	FOREACH(stars)
	{
		(*i)->clear();
	}

	// preallocate memory to avoid subsequent frequent reallocation
	FOR(0, modelCPDF.size())
	{
		double frac = i+1 != modelCPDF.size() ? modelCPDF[i+1] : 1;
		int n = (int)(K * 1.2 * (frac - modelCPDF[i]));
		stars[i]->reserve(n);
		std::cerr << "Preallocating space for " << n << " stars\n";
	}

	// generate K stars in a two step process. First, draw (x,y) positions
	// from the distribution for all K stars. Then, sort the resulting array
	// by (x,y), and then for all stars in given (x,y) pixel draw ri and m.
	// This two-step process is necessary in order to avoid recalculating 
	// P(m|x,y,ri) distribution at each step (it's time-consuming to compute,
	// and there's not enough space to cache it for all pixels at once).
	std::cerr << "# Generating...\n";
	ticker tick(10000);
	int Kgen = 0;
	FOR(0, K)
	{
		double u = gsl_rng_uniform(rng);
		ITER(modelCPDF) ix = upper_bound(modelCPDF.begin(), modelCPDF.end(), u); --ix;
		int idx = ix - modelCPDF.begin();
//		std::cerr << "Model index : " << idx << "\n";

		model_pdf::star s;
		bool succ = pdfs[idx]->draw_position(s, rng);
		if(!succ)
		{
			if(!allowMisses) { --i; }
			else { tick.tick(); }
			
			continue;
		} else {
			this->stars[idx]->push_back(s);
			Kgen++;
		}
		tick.tick();
	}
	tick.close();

	// 2nd stage - draw magnitudes and write output
	FOR(0, pdfs.size())
	{
		pdfs[i]->draw_magnitudes(*stars[i], rng);
//		observe(*stars[i], pdfs[i]->proj, out);
		observe2(*stars[i], pdfs[i]->galmodel(), pdfs[i]->proj, out);

		std::cerr << "model " << pdfs[i]->name() << " : " << stars[i]->size() << " stars (" << 
			100. * double(stars[i]->size()) / Kgen << "%)\n";
	}

	return Kgen;
}

bool model_pdf::draw_position(star &s, gsl_rng *rng)
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
	s.x = skymap.x0 + skymap.dx*(X.X + gsl_rng_uniform(rng));
	s.y = skymap.y0 + skymap.dx*(Y.Y + gsl_rng_uniform(rng));

	// check that the star is inside survey footprint, reject if it's not
	gpc_polygon &poly = skymap.skymap[std::make_pair(X.X, Y.Y)].poly;
	gpc_vertex vtmp = { s.x, s.y };
	if(!in_polygon(vtmp, poly))
//	if(!in_polygon((gpc_vertex const&)std::make_pair(s.x, s.y), poly))
	{
		return false;
	}

	// pick ri bin
	double uri = gsl_rng_uniform(rng);
	ITER(Y.ripdf) iri = upper_bound(Y.ripdf.begin(), Y.ripdf.end(), uri); --iri;
	int riidx = iri - Y.ripdf.begin();
	// pick ri within [ri, ri+dri)
	double ri = ri0 + dri*(riidx + gsl_rng_uniform(rng));
	ASSERT(ri0+riidx*dri <= ri && ri <= ri0+dri+riidx*dri);
	s.ri = ri;
	
	return true;
}

// assigns magnitudes to objects in stars vector
void model_pdf::draw_magnitudes(std::vector<model_pdf::star> &stars, gsl_rng *rng)
{
	std::cerr << "# Sorting... \n";

	// sort stars by X, Y, ri
	std::sort(stars.begin(), stars.end(), star_lessthan(this));

	std::cerr << "# Assigning magnitudes \n";

	ticker tick(10000);
	// draw magnitudes for stars in each cell X,Y
	boost::tuple<int, int, int> cell = boost::make_tuple(-1000,-1000,-1000);
	cumulative_dist mspl;
	ITER(stars) it = stars.begin();
	while(it != stars.end())
	{
//		std::cerr << "+" << "\n";
		star &s = *it++;
		boost::tuple<int, int, int> cell2(s.X(skymap), s.Y(skymap), s.RI(*this));
		// if we moved to a different (x, y, ri) bin, calculate
		// the cumulative marginal pdf in magnitude for that bin
		if(cell != cell2)
		{
			magnitude_mpdf(mspl, skymap.x0 + (s.X(skymap) + 0.5)*skymap.dx, skymap.y0 + (s.Y(skymap) + 0.5)*skymap.dx, ri0 + (s.RI(*this) + 0.5)*dri);
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
}

star_output_to_dmm::star_output_to_dmm(const std::string &objcat, const std::string &obscat, bool create)
{
	if(create)
	{
		out.create(objcat);
		starmags.create(obscat);
	} else {
		out.open(objcat, "rw");
		starmags.open(obscat, "rw");
	}
	out.setmaxwindows(2);
	starmags.setmaxwindows(2);
}

void star_output_to_textstream::output_header(const sstruct::factory_t &factory)
{
	// simple ascii-text dump
	//out << "# l    b    ri     r  " << factory << "\n";
	out << "# " << factory << "\n";
}

void star_output_to_textstream::output(sstruct &t)
{
	// simple ascii-text dump
//	out << deg(l) << " " << deg(b) << " " << ri << " " << r << " " << t << "\n";
	out << t << "\n";
}

void star_output_to_textstream::output(Radians ra, Radians dec, double Ar, std::vector<std::pair<observation, obsv_id> > &obsvs)
{
	// simple ascii-text dump
	double ri = obsvs[0].first.mag[1] - obsvs[0].first.mag[2];
	out << ra << " " << dec << " " << ri << " " << obsvs[0].first.mag[1] << "\n";
}

// Quick hack -- a position-independent luminosity function for binaries
void sky_generator::draw_companion(float &g, float &r, float &i, Radians l, Radians b, double dm /*distance modulus*/)
{
#if 1
	// HACK: just load a position-independent luminosity function from lumfun.txt
	static cumulative_dist lf; // luminosity function
	if(lf.empty())
	{
		// luminosity function
		text_input_or_die(lfin, "lumfun.txt")
		std::vector<double> ri, phi;
		load(lfin, ri, 0, phi, 1);
		// convert to cumulative distribution, from [0,1]
		FOR(1, phi.size()) { phi[i] += phi[i-1]; }
		FOR(0, phi.size()) { phi[i] = (phi[i] - phi.front()) / (phi.back() - phi.front()); }
		lf.construct(ri, phi);
	}
	float ri = lf(gsl_rng_uniform(rng));
#else
	float u = gsl_rng_uniform(rng);
	float ri = 0.1 + u*1.3; // HACK: draw r-i uniformly in 0.1-1.3 range
#endif

	r = paralax.Mr(ri) + dm;
	i = r - ri;
	g = r + paralax.gr(ri);
}

#define DBG_MAGCHECK 0

#if DBG_MAGCHECK
	float u, g, r, i, z;
#endif

void sky_generator::observe2(const std::vector<model_pdf::star> &stars, galactic_model &model, peyton::math::lambert &proj, star_output_function &sf)
{
	std::cerr << "Writing... ";
	ticker tick(10000);

	std::auto_ptr<sstruct> tagptr(sstruct::create());
	sstruct &t = *tagptr;

	const static double Rg = coord_pack::Rg;
	FORj(j, 0, stars.size())
	{
		const model_pdf::star &s = stars[j];

		// position -- convert to l,b
		Radians l, b;
		proj.inverse(s.x, s.y, l, b);
		double cl = cos(l), cb = cos(b), sl = sin(l), sb = sin(b);
		coordinates::galequ(l, b, l, b);

		// absolute magnitudes and distance
		const double Mr = paralax.Mr(s.ri);
		const double D = stardist::D(s.m, Mr);

		// galactocentric x,y,z
		double x = Rg - D*cl*cb;
		double y =    - D*sl*cb;
		double z =      D*sb;

		// draw model-dependent tags
		model.draw_tag(t, x, y, z, s.ri, rng);

		// write out this star and its tags
		t.lonlat()  = std::make_pair(deg(l), deg(b));
		t.color()   = s.ri;
		t.mag()     = s.m;
		sf.output(t);

		// user interface stuff
		tick.tick();
	}
	tick.close();
	std::cerr << "\n";
}

void test_tags()
{
	if(0) {
		sstruct::factory.useTag("star_name");
		sstruct::factory.useTag("xyz[3]");
		sstruct::factory.useTag("extinction.r");

		std::cerr << sstruct::factory.ivars[0] << "\n";
		std::cerr << sstruct::factory.ivars[1] << "\n";
		std::cerr << sstruct::factory.ivars[2] << "\n";
		std::cerr << sstruct::factory.ivars[3] << "\n";
		std::cerr << sstruct::factory.ivars[4] << "\n";

		if(1) {
			sstruct *tag = sstruct::create();
			sstruct::factory.useTag("extinction.r");

			float &Ar = tag->ext_r();
			float *xyz = tag->xyz();
			std::string &starname = tag->starname();
			starname = "Bla";

			sstruct *tag2 = sstruct::create();
			*tag2 = *tag;
			starname = "Newyyy";

			std::cerr << "Ar = " << Ar << "\n";
			std::cerr << "xyz = " << xyz[0] << " " << xyz[1] << " " << xyz[3] << "\n";
			std::cerr << "starname = " << starname << "\n";
			std::cerr << "starname2 = " << tag2->starname() << "\n";

			std::cerr << "Serializing tag list: "; sstruct::factory.serialize(std::cerr); std::cerr << "\n";
			std::cerr << "Serializing 1: ";  tag->serialize(std::cerr); std::cerr << "\n";
			std::cerr << "Serializing 2: "; tag2->serialize(std::cerr); std::cerr << "\n";

			delete tag;
			delete tag2;
		} else {
			sstruct *tags = sstruct::create(5);

			FOR(0, 5)
			{
				sstruct *tag = &tags[i];

				float &Ar = tag->ext_r();
				float *xyz = tag->xyz();
				tag->starname() = "Lipa moja";

				std::cerr << i << "  Ar = " << Ar << "\n";
				std::cerr << i << "  xyz = " << xyz[0] << " " << xyz[1] << " " << xyz[3] << "\n";
				std::cerr << "starname = " << tag->starname() << "\n";
			}
			delete [] tags;
		}
	}
	else
	{
		// Unserialize text file and print it out to the screen
		std::ifstream in("out.txt");
		sstruct::factory.unserialize(in);
		sstruct::factory.useTag("vel[3]");

		sstruct *tag = sstruct::create(100);
		size_t cnt = 0;
		std::cout << sstruct::factory << "\n";
		while(in >> tag[cnt])
		{
			std::cout << tag[cnt];
			std::cout << "\n";
			cnt++;
		}
		delete [] tag;
	}
}

void sky_generator::observe(const std::vector<model_pdf::star> &stars, peyton::math::lambert &proj, star_output_function &sf)
{
	std::cerr << "Observing... ";
	ticker tick(10000);
	std::vector<std::pair<observation, obsv_id> > obsvs;
	FORj(j, 0, stars.size())
	{
		model_pdf::star s = stars[j];

		obsv_id oid;
		observation obsv;

		#if !DBG_MAGCHECK
		float u, g, r, i, z;
		#endif

		// set extinction
		obsv.Ar = Ar;

		// position
		Radians l, b;
		proj.inverse(s.x, s.y, l, b);
		coordinates::galequ(l, b, obsv.ra, obsv.dec);
		obsv.ra = deg(obsv.ra); obsv.dec = deg(obsv.dec);

		// calculate magnitudes, absolute magnitudes and distance
		const double Mr = paralax.Mr(s.ri);
		const double D = stardist::D(s.m, Mr);
		u = s.m;
		g = s.m + paralax.gr(s.ri);
		r = s.m;
		i = s.m - s.ri;
		z = s.m;

		//
		// Physical things affecting the object's magnitudes and colors
		//

		float binaryFraction = 0.0;
		if(binaryFraction && (gsl_rng_uniform(rng) <= binaryFraction))
		{
			// this star has a binary companion -- add it

			// - draw the companion from the luminosity function
			float gb, rb, ib;
			draw_companion(gb, rb, ib, l, b, r - Mr);

#if DBG_PRINT_BINARIES
			std::cout << r - Mr << "   ";
			std::cout << g << " " << r << " " << i << "   ";
			std::cout << gb << " " << rb << " " << ib << "   ";
#endif
			// calculate joint magnitudes
			g = -2.5*log10(pow(10., -0.4*g) + pow(10., -0.4*gb));
			r = -2.5*log10(pow(10., -0.4*r) + pow(10., -0.4*rb));
			i = -2.5*log10(pow(10., -0.4*i) + pow(10., -0.4*ib));

#if DBG_PRINT_BINARIES
			std::cout << g << " " << r << " " << i << "\n";
#endif
		} else {
//			ASSERT(0);
		}
		// ignore stars fainter than flux limit (HACK)
//		if(r > 22.5) { continue; }
//		if(r > 23.0 && r-i < 0.3) { continue; }

		bool subdwarfs = false;
		if(subdwarfs)
		{
			ASSERT(0); // the logic of this code has to be rechecked, in light of algorithm changes here

			// calculate galactic coordinates
			V3 v; v.celestial(l, b, D);
			v.x = Rg - v.x;
			v.y *= -1;

			// calculate "metalicity factor" based on location in the Galaxy
			double f;
			if(v.z < 500) { f = 0; }
			else if(v.z < 2500) { f = (v.z - 500.) / 2000.; }
			else { f = 1; }

			// adjust the absolute magnitudes accordingly
			u += f;
			g += f;
			r += f;
			i += f;
			z += f;
		}

		if(paralax_dispersion)
		{
			// mix in the photometric paralax uncertancy - here is how this goes:
			// The star is at distance D. Give its' color, the normal absolute magnitude
			// it would have is Mr. However, due to age & metalicity, the star has
			// absolute magnitude Mr' = Mr + dMr, while _keeping_ the same color and
			// distance. The observer therefore sees the star to be dimmer/brighter than
			// a star with absmag Mr (==> change of magnitude).
			// Note: this is an oversimplification, as some colors do change (we should implement this)
			double dMr = gsl_ran_gaussian(rng, paralax_dispersion);
			u += dMr;
			g += dMr;
			r += dMr;
			i += dMr;
			z += dMr;
		}

		//
		// Things happening at observation ("in the telescope")
		//

		// add extinction and
		// store _extinction uncorrected_ magnitudes to obsv_mag
		obsv.mag[0] = g;
		obsv.mag[1] = r;
		obsv.mag[2] = i;
		obsv.mag[3] = z;
		obsv.mag[4] = u;
		FOR(0, 5) { obsv.mag[(i+4)%5] += extinction(Ar, i); }

		// calculate the magnitude of observational errors
		if(flags & APPLY_PHOTO_ERRORS)
		{
			if(!magerrs.empty())
			{
				FOR(0, 5) { obsv.magErr[i] = magerrs(l, b, obsv.mag[i]); }
			}

			if(constant_photo_error)
			{
				FOR(0, 5) { obsv.magErr[i] = constant_photo_error; }
			}

			// mix in magnitude errors
			FOR(0, 5) { obsv.mag[i] += gsl_ran_gaussian(rng, obsv.magErr[i]); }
		}
		else
		{
			FOR(0, 5) { obsv.magErr[i] = 0; }
		}

// 		std::cerr << "magerrs = ";
// 		FOR(0, 5) { std::cerr << obsv.magErr[i] << " "; }
// 		std::cerr << "\n";

		// "observe" the object and store it to object and observation database
		obsvs.clear();
		obsvs.push_back(std::make_pair(obsv, oid));
		sf.output(obsv.ra, obsv.dec, obsv.Ar, obsvs);

		// user interface stuff
		tick.tick();
	}
	tick.close();
	std::cerr << "\n";
}

void star_output_to_dmm::output(sstruct &t)
{
	std::cerr << "Output to DMM with galactic model tags not implemented. Aborting.";
	ASSERT(0);
}

void star_output_to_dmm::output(Radians ra, Radians dec, double Ar, std::vector<std::pair<observation, obsv_id> > &obsvs)
{
	mobject m = process_observations(starmags.size(), ra, dec, Ar, obsvs, po_info);

#if DBG_MAGCHECK
	// while debugging - if added no photometric error, the magnitudes in
	// mobject must be the same as the input magnitudes
	ASSERT(fabs(m.mag[0] - g) < 1e-4) { std::cerr << m.mag[0] << " " << g << "\n"; }
	ASSERT(fabs(m.mag[1] - r) < 1e-4) { std::cerr << m.mag[1] << " " << r << "\n"; }
	ASSERT(fabs(m.mag[2] - i) < 1e-4) { std::cerr << m.mag[2] << " " << i << "\n"; }
	ASSERT(fabs(m.mag[3] - z) < 1e-4) { std::cerr << m.mag[3] << " " << z << "\n"; }
	ASSERT(fabs(m.mag[4] - u) < 1e-4) { std::cerr << m.mag[4] << " " << u << "\n"; }

	ASSERT(fabs(m.ml_mag[0] - g) < 1e-4) { std::cerr << m.mag[0] << " " << g << "\n"; }
	ASSERT(fabs(m.ml_mag[1] - r) < 1e-4) { std::cerr << m.mag[1] << " " << r << "\n"; }
	ASSERT(fabs(m.ml_mag[2] - i) < 1e-4) { std::cerr << m.mag[2] << " " << i << "\n"; }
#endif

	FOREACH(obsvs) { starmags.push_back((*i).first); }
	out.push_back(m);

	if((out.size() % 500000) == 0) { out.sync(); }
}

void star_output_to_dmm::close()
{
	std::cerr << "# Observation processing statistics\n";
	std::cerr << po_info << "\n";

	out.close();
	starmags.close();
}

void pdfinfo(std::ostream &out, const std::string &pdffile)
{
	std::ifstream iin(pdffile.c_str()); ASSERT(iin) { std::cerr << "Could not open " << pdffile << "\n"; }
	io::ibstream in(iin);

	model_pdf pdf;
	if(!(in >> pdf)) { ASSERT(0) { std::cerr << "Error loading " << pdffile << "\n"; }; }

	out << "nstars\t" << pdf.N << "\n";
}

#endif // COMPILE_SIMULATE_X
