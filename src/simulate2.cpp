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

#include "gsl/gsl_randist.h"

#include "gpc_cpp.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/shared_ptr.hpp>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "simulate.h"
#include "projections.h"
#include "model.h"
#include "paralax.h"
#include "analysis.h"

#include <vector>
#include <map>
#include <string>

#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

void poly_bounding_box(double &x0, double &x1, double &y0, double &y1, const gpc_polygon &p);

gpc_polygon poly_rect(double x0, double x1, double y0, double y1);

model_pdf::model_pdf(const std::string &pdfname_)
	: pdfname(pdfname_)
{
}

model_pdf::model_pdf(std::istream &in, const std::string &pdfname_)
: proj(rad(90), rad(90)), pdfname(pdfname_)
{
	Config cfg; cfg.load(in);

	cfg.get(skymap.dx,	"dx", 	1.);
	skymap.dx = rad(skymap.dx);

	cfg.get(dri,	 "dri",    0.01);
	cfg.get(ri0,	 "ri0",	   0.0);
	cfg.get(ri1,	 "ri1",	   1.5);
	cfg.get(m0,	 "r0",	  14.);
	cfg.get(m1,	 "r1",	  22.);
	cfg.get(dm,	 "dm",	   0.1);
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
		}

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
	proj.deproject(l, b, x, y);
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
		n = integrate(m0, m1, dm, model_fun, params, 1e-10, 1e-5);
		n *= sqr(skymap.dx)*dri;	// multiply by solid angle. n is now the total number of stars in direction (l,b)
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
		*i /= sum;
	}

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


gpc_polygon make_circle(double x0, double y0, double r, double dx);

// calculate marginal PDFs (``mpfds'') for stars within a given sky footprint,
// and with a given Galactic model. Store the mPDFs into an XPDF tree
// and serialize it to a file. These mpdfs are later used in generating
// Monte Carlo realisations of the sky
//
// This routine also calculates the 'skymap', for fast point-in-polygon
// lookups of the observed sky footprint.
void model_pdf::construct_mpdf(const std::string &footfn, const std::string &modelfn)
{
	using boost::shared_ptr;
	using std::cout;

	// load the footprint
	gpc_polygon sky = xgpc_read(footfn, proj);

	// load the model
	std::ifstream inconf(modelfn.c_str());
	if(!inconf) { THROW(EFile, "Could not access " + modelfn); }
	model.reset(galactic_model::load(inconf));
	if(model.get() == NULL) { THROW(EAny, "Error loading galactic model from " + modelfn); }

	MLOG(verb1) << "footprint file  = " << footfn;
	MLOG(verb1) << "projection = " << proj;
	MLOG(verb1) << "dx   = " << deg(skymap.dx);
	MLOG(verb1) << "dri  = " << dri;
	MLOG(verb1) << "ri0  = " << ri0;
	MLOG(verb1) << "ri1  = " << ri1;
	MLOG(verb1) << "r0   = " << m0;
	MLOG(verb1) << "r1   = " << m1;
	MLOG(verb1) << "dm   = " << dm;
	MLOG(verb1) << "band = " << model->band();
	MLOG(verb1) << "color = " << model->color();

	skymap.skymap.clear();
	xpdf.clear();
	int N = 0; // number of dx^2 cells processed

 	poly_bounding_box(skymap.x0, skymap.x1, skymap.y0, skymap.y1, sky);

	double xsum = 0; // Grand total - number of stars in the whole field
	int X = 0;
	for(double x = skymap.x0; x < skymap.x1; x += skymap.dx, X++) // loop over all x values in the bounding rectangle
	{
		shared_ptr<Sx> sx(new Sx(xsum, X));

		double xa = x, xb = x+skymap.dx;
		int Y = 0;
		double xysum = 0.;
		for(double y = skymap.y0; y < skymap.y1; y += skymap.dx, Y++) // loop over all y values for a given x
		{
			double ya = y, yb = y+skymap.dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon poly;
			gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
			if(poly.num_contours == 0) continue; // if there are no observations in this direction

			// store the polygon into a fast lookup map
			partitioned_skymap::pixel_t &pix = skymap.skymap[std::make_pair(X, Y)];
			pix.poly = poly;
			pix.coveredArea = polygon_area(poly);
			pix.pixelArea = sqr(skymap.dx);

			shared_ptr<Sy> sy(new Sy(xysum, Y));

			// calculate CMD integrals & marginal distributions
			double risum = ri_mpdf(sy->ripdf, (xa + xb)/2., (ya + yb)/2.);

			// store cumulative distribution
			sx->ypdf.push_back(sy);
			xysum += risum;

			++N;
		}

		// convert starcounts in ypdf to marginal probabilities (prob. of finding a star w. y less then Y given x)
		if(xysum != 0.)
		{
			FOREACH(sx->ypdf) { (*i)->pcum /= xysum; }
			#if 1
			std::stringstream ss;
			FOREACH(sx->ypdf) { ss << "(" << X << " " << (*i)->Y << " " << (*i)->pcum << ")"; }
			MLOG(verb1) << "N = " << N << " " << ss.str();
			#endif
			xpdf.push_back(sx);
			xsum += xysum;
		}
	}
	MLOG(verb1) << "Total stars = " << xsum << " stars within " << skymap.skymap.size()*sqr(deg(skymap.dx)) << "deg^2 of pixelized area.";

	// convert starcounts to cumulative probabilities
	FOREACH(xpdf) { (*i)->pcum /= xsum; }

	std::stringstream ss;
	FOREACH(xpdf) { ss << "(" << (*i)->X << " -1 " << (*i)->pcum << ")"; }
	MLOG(verb1) << "X probabilities - " << ss.str();

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
void model_pdf::magnitude_mpdf(cumulative_dist &mspl, double x, double y, double ri)
{
	// deproject to (l,b)
	Radians l, b;
	proj.deproject(l, b, x, y);
	pencil_beam pb(l, b);

	// absolute magnitude for this spectral type
	double Mr = model->absmag(ri);

	double norm = 1.;

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
	FOREACH(f) { *i /= fmax; }

	// construct spline approximation
	mspl.construct(m, f);
}

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

	model.reset(galactic_model::unserialize(in));

	DLOG(verb1) << "N = " << N;
	DLOG(verb1) << "band = " << model->band();
	DLOG(verb1) << "color = " << model->color();
	DLOG(verb1) << "xpdf.size() = " << xpdf.size();
	DLOG(verb1) << "skymap.skymap.size() = " << skymap.skymap.size();

	return in;
}

sky_generator::sky_generator(std::istream &in, const std::string &pdfs)
: nstars(0)
{
	Config cfg; cfg.load(in);

	// number of stars to read
	cfg.get(nstars,	"nstars", 	0);

	// random number generator
	int seed;
	cfg.get(seed,	"seed", 	42);
	rng.reset(new rng_gsl_t(seed));

	// load pdfs
	std::istringstream ss(pdfs);
	std::string pdfFn;
	while(ss >> pdfFn)
	{
		std::ifstream iin(pdfFn.c_str());
		if(!iin) { THROW(EIOException, "Could not access " + pdfFn + "."); }
		io::ibstream in(iin);

		boost::shared_ptr<model_pdf> pdf(new model_pdf(pdfFn));
		if(!(in >> *pdf)) { THROW(EAny, "Error loading galactic model from " + pdfFn); }

		add_pdf(pdf);
	}

}

sky_generator::~sky_generator()
{
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
// file, precalculated with preconstruct_mpdf
void sky_generator::montecarlo(star_output_function &out)
{
	ASSERT(pdfs.size() != 0);

	bool allowMisses = false;
	int Ktotal = nstars;

	if(Ktotal == 0)
	{
		double Kd = 0;
		//std::cerr << "# " << (*pdfs.begin())->N << "\n";
		FOREACH(pdfs) { Kd += (*i)->N; }
		Ktotal = (int)rint(Kd);
		allowMisses = true;
	}

	MLOG(verb1) << "Generating " << Ktotal << " stars" << (allowMisses ? " (approximately, based on model norm.)" : "");

	// calculate cumulative probability function of a star being a part
	// various pdfs. This will be used to decide for each star according
	// to which model_pdf are we going to draw it.
	std::vector<double> modelCPDF(pdfs.size());
	modelCPDF[0] = 0.;
	FOR(0, pdfs.size()-1) { modelCPDF[i+1] = modelCPDF[i] + pdfs[i]->N; }
	double norm = pdfs.back()->N + modelCPDF.back();
	FOR(1, pdfs.size()) { modelCPDF[i] /= norm; }
	
	// generate the data in smaller batches (to conserve memory)
	static const int Kbatch = 10000000;

	// prepare output
	otable t(Kbatch);
	t.use_column("lb");

	std::string
		colorName = pdfs[0]->galmodel().color(),
		magName = pdfs[0]->galmodel().band(),
		absmagName = "abs" + pdfs[0]->galmodel().band();
	MLOG(verb1) << "color=" << colorName;
	MLOG(verb1) << "mag=" << magName;
	MLOG(verb1) << "absmag=" << absmagName;

	t.use_column(colorName);  t.alias_column(colorName, "color");
	t.use_column(absmagName); t.alias_column(absmagName, "absmag");
	t.use_column(magName);	  t.alias_column(magName, "mag");

	FOR(0, pdfs.size())
	{
		if(pdfs[i]->galmodel().band() != pdfs[0]->galmodel().band()) {
			THROW(EAny, "Bands in models of " + pdfs[0]->name() + " and " + pdfs[i]->name() + " differ (" + pdfs[0]->galmodel().band() + " vs. " + pdfs[i]->galmodel().band() + ")");
		}
		if(pdfs[i]->galmodel().color() != pdfs[0]->galmodel().color()) {
			THROW(EAny, "Colors in models of " + pdfs[0]->name() + " and " + pdfs[i]->name() + " differ (" + pdfs[0]->galmodel().color() + " vs. " + pdfs[i]->galmodel().color() + ")");
		}
	}

	// do catalog generation
	int K = 0, Ngen = 0;
	while(K < Ktotal)
	{
		t.clear();
		int Kstep = std::min(Kbatch, Ktotal - K);
		Ngen += montecarlo_batch(t, Kstep, modelCPDF, allowMisses);

		if(K == 0) { out.output_header(t); }
		out.output(t);

		K += Kstep;
	}

	MLOG(verb1) << "Generated " << Ngen << " stars.";
}

int sky_generator::montecarlo_batch(otable &t, int K, const std::vector<double> &modelCPDF, bool allowMisses)
{
	// clear output vectors
	FOREACH(stars) { (*i)->clear(); }

	// preallocate memory to avoid subsequent frequent reallocation
	FOR(0, modelCPDF.size())
	{
		double frac = i+1 != modelCPDF.size() ? modelCPDF[i+1] : 1;
		int n = (int)(K * 1.2 * (frac - modelCPDF[i]));
		stars[i]->reserve(n);
		DLOG(verb1) << "Preallocating space for " << n << " stars";
	}

	// generate K stars in a two step process. First, draw (x,y) positions
	// from the distribution for all K stars. Then, sort the resulting array
	// by (x,y), and then for all stars in given (x,y) pixel draw ri and m.
	// This two-step process is necessary in order to avoid recalculating 
	// P(m|x,y,ri) distribution at each step (it's time-consuming to compute,
	// and there's not enough space to cache it for all pixels at once).
	ticker tick("Generating", 10000);
	int Kgen = 0;
	FOR(0, K)
	{
		double u = rng->uniform();
		ITER(modelCPDF) ix = upper_bound(modelCPDF.begin(), modelCPDF.end(), u); --ix;
		int idx = ix - modelCPDF.begin();

		model_pdf::star s;
		bool succ = pdfs[idx]->draw_position(s, *rng);
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

	// 2nd stage - draw magnitudes and store output
	FOR(0, pdfs.size())
	{
		pdfs[i]->draw_magnitudes(*stars[i], *rng);
		draw_stars(*stars[i], pdfs[i]->galmodel(), pdfs[i]->proj, t);

		MLOG(verb1) << "model " << pdfs[i]->name() << " : " << stars[i]->size() << " stars (" <<
			100. * double(stars[i]->size()) / Kgen << "%)";
	}

	return Kgen;
}

bool model_pdf::draw_position(star &s, rng_t &rng)
{
	// pick x
	double ux = rng.uniform();
	ITER(xpdf) ix = upper_bound(xpdf.begin(), xpdf.end(), ux, less_S<Sx>()); --ix;
	const Sx &X = **ix;

	// pick y, given an x
	double uy = rng.uniform();
	ITER(X.ypdf) iy = upper_bound(X.ypdf.begin(), X.ypdf.end(), uy, less_S<Sy>()); --iy;
	const Sy &Y = **iy;

	// draw the exact location within the [x,x+dx], [y,y+dy] rectangle
	s.x = skymap.x0 + skymap.dx*(X.X + rng.uniform());
	s.y = skymap.y0 + skymap.dx*(Y.Y + rng.uniform());

	// check that the star is inside survey footprint, reject if it's not
	gpc_polygon &poly = skymap.skymap[std::make_pair(X.X, Y.Y)].poly;
	gpc_vertex vtmp = { s.x, s.y };
	if(!in_polygon(vtmp, poly))
	{
		return false;
	}

	// pick ri bin
	double uri = rng.uniform();
	ITER(Y.ripdf) iri = upper_bound(Y.ripdf.begin(), Y.ripdf.end(), uri); --iri;
	int riidx = iri - Y.ripdf.begin();
	// pick ri within [ri, ri+dri)
	double v;
	double ri = ri0 + dri*(riidx + (rng.uniform()));
	ASSERT(0.99999999*ri0+riidx*dri <= ri && ri <= ri0+dri+riidx*dri*1.00000001)
	{
		std::cerr << "ri0 = " << ri0 << " ri1=" << ri1 << " dri=" << dri << " ri=" << ri << "\n";
		std::cerr << "ri0+riidx*dri=" << ri0+riidx*dri << " ri0+dri+riidx*dri=" << ri0+dri+riidx*dri << "\n";
		std::cerr << "v=" << v << "\n";
	}
	s.ri = ri;
	
	return true;
}

// assigns magnitudes to objects in stars vector
void model_pdf::draw_magnitudes(std::vector<model_pdf::star> &stars, rng_t &rng)
{
	MLOG(verb1) << "Sorting...";

	// sort stars by X, Y, ri
	std::sort(stars.begin(), stars.end(), star_lessthan(this));

	ticker tick("Assigning magnitudes", 10000);
	// draw magnitudes for stars in each cell X,Y
	boost::tuple<int, int, int> cell = boost::make_tuple(-1000,-1000,-1000);
	cumulative_dist mspl;
	ITER(stars) it = stars.begin();
	while(it != stars.end())
	{
		star &s = *it++;
		boost::tuple<int, int, int> cell2(s.X(skymap), s.Y(skymap), s.RI(*this));
		// if we moved to a different (x, y, ri) bin, calculate
		// the cumulative marginal pdf in magnitude for that bin
		if(cell != cell2)
		{
			magnitude_mpdf(mspl, skymap.x0 + (s.X(skymap) + 0.5)*skymap.dx, skymap.y0 + (s.Y(skymap) + 0.5)*skymap.dx, ri0 + (s.RI(*this) + 0.5)*dri);
			cell = cell2;
		}
		// draw the magnitude
		double um = rng.uniform();
		s.m = mspl(um);
		tick.tick();
	}
}

void star_output_to_textstream::output_header(otable &t)
{
	// simple ascii-text dump
	out << "# ";
	t.serialize_header(out);
	out << "\n";
	out.flush();
}

void star_output_to_textstream::output(const otable &t)
{
	// simple ascii-text dump
	t.serialize_body(out);
}

void sky_generator::draw_stars(const std::vector<model_pdf::star> &stars, galactic_model &model, peyton::math::lambert &proj, otable &t)
{
	ticker tick("Writing", 10000);

	cdouble_t::host_t lb    = t.col<double>("lb");
	cfloat_t::host_t color  = t.col<float>("color");
	cfloat_t::host_t mag    = t.col<float>("mag");
	cfloat_t::host_t absmag = t.col<float>("absmag");
	cfloat_t::host_t XYZ    = t.col<float>("XYZ");	t.set_output("XYZ", true);
	cint_t::host_t comp     = t.col<int>("comp");	t.set_output("comp", true);

	const static double Rg = coord_pack::Rg;
	FORj(j, 0, stars.size())
	{
		const model_pdf::star &s = stars[j];

		// position -- convert to l,b
		Radians l, b;
		proj.deproject(l, b, s.x, s.y);
		double cl = cos(l), cb = cos(b), sl = sin(l), sb = sin(b);
		if(l < 0) { l = ctn::twopi + l; }

		// absolute magnitudes and distance
		const double Mr = model.absmag(s.ri);
		const double D = stardist::D(s.m, Mr);

		// galactocentric x,y,z
		double x = Rg - D*cl*cb;
		double y =    - D*sl*cb;
		double z =      D*sb;

		// store this star
		size_t row = t.add_row();
		lb(row, 0)  = deg(l);
		lb(row, 1)  = deg(b);
		color(row)  = s.ri;
		mag(row)    = s.m;
		absmag(row) = Mr;
		comp(row) = 0;
		XYZ(row, 0) = x;
		XYZ(row, 1) = y;
		XYZ(row, 2) = z;

		// user interface stuff
		tick.tick();
	}
	tick.close();

	// add any model-dependent details
	model.add_details(t, *rng);
}

void test_otable()
{
	return;

	otable result(10);
	std::string text =
		"# lb[2]{type=float;fmt=% .3f} SDSSu SDSSg SDSSr test[3]\n"
		"1.11 2.22	1 2 3	1.2 2.3 3.4\n"
		"2.11 4.22	2 5 8	2.2 3.3 4.4\n"
		"3.11 5.22	3 6 9	2.2 3.3 4.4\n"
		"4.11 6.22	4 7 0	2.2 3.3 4.4\n";

	std::istringstream in(text.c_str());
	result.unserialize_header(in);
	result.use_column("addition");
	result.unserialize_body(in);

//	result.set_output_all();
	std::cout << "# ";
	result.serialize_header(std::cout);
	std::cout << "\n";
	result.serialize_body(std::cout);

	exit(-1);
}

void pdfinfo(std::ostream &out, const std::string &pdffile)
{
	std::ifstream iin(pdffile.c_str());
	if(!iin) { THROW(EFile, "Could not open " + pdffile); }
	io::ibstream in(iin);

	model_pdf pdf;
	if(!(in >> pdf)) { THROW(EAny, "Error unserializing PDF from " + pdffile); }

	out << "nstars\t" << pdf.N << "\n";
}
