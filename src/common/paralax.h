#ifndef __mstars_h
#define __mstars_h

#include "xcat.h"
#include "gslcc_min.h"
#include <gsl/gsl_poly.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <astro/types.h>
#include <astro/math.h>
#include <astro/constants.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <valarray>

namespace stardist
{
	inline double logD(float m, double M) { return 1. + (m - M)/5.; }			///< log10 distance, based on absolute and apparent i magnitudes
	inline double D(float m, double M) { return std::pow(10.0, logD(m, M)); }	///< distance, based on absolute and apparent i magnitudes
	inline float m(double D, double M) { return M + 5*log10(D/10); }
};

double paralax_with_prior(float ri, float gr, float sg, float sr, float si);
double paralax_without_prior(float ri, float gr, float sg, float sr, float si);

class plx_gri_locus : public paralax_relation
{
public:
	static const double a = 4.9;	//	a = 5.0;
	static const double b = 2.45;//	b = 2.4;
	static const double c = 1.68;//	c = 1.6;
	static const double d = 0.035;//	d = 0.1;
	static const double f = 1.39;
public:
	float gr(float ri)
	{
		// old locus relation
		//float F = ((5.0*ri + 2.4)*ri + 1.6)*ri + 0.1;
		// new locus relation
		//float F = ((4.9*ri + 2.45)*ri + 1.68)*ri + 0.049;
		//float F = ((4.9*ri + 2.45)*ri + 1.68)*ri + 0.035;
		//float F = ((4.9*ri + 2.45)*ri + 1.68)*ri + 0.040;
		//float F = ((4.9*ri + 2.45)*ri + 1.68)*ri + 0.050;
		float F = ((a*ri + b)*ri + c)*ri + d;
		return f*(1-exp(-F));
	}

	#define sqr peyton::sqr
	float ri(float gr)
	{
		double aa = sqr(b) - 3*a*c;
		double a0 = 2*pow(b, 3) - 9*a*b*c + 27*sqr(a)*(d + std::log(1-gr/f));
		double a1 = sqrt(-4*pow(aa, 3) + sqr(a0));
		double a2 = pow(-a0+a1, 1./3.);

		double a3 = -2*b + 2*pow(2., (1./3.)) * aa / a2;
		double a4 = pow(2., 2./3.)*a2;

		return 1/(6*a)*(a3 + a4);
	}
	#undef sqr

public:
	// Maximum likelihood locus forcing
	struct ml_grri : public gsl::mmizer
	{
		plx_gri_locus &as;
		double y0, x0, sx2, sy2, sxy, denom, lnN;
		//double (*prior)(double);
		gsl_spline *prior;
		gsl_interp_accel *acc;

		ml_grri(plx_gri_locus &as_) : as(as_), mmizer(-.5, 2.5, 0, .001), prior(NULL), acc(NULL) { }
		void setprior(const std::string &priorfile);

		// finds r-i for which the likelihood is maximum
		double operator()(float ri, float gr, float sg, float sr, float si, float *lnL = NULL)
		{
			y0 = ri; x0 = gr;
			sx2 = peyton::sqr(sg) +  peyton::sqr(sr);	// var(g-r)
			sy2 = peyton::sqr(sr) +  peyton::sqr(si);	// var(r-i)
			sxy = -peyton::sqr(sr);						// cov(g-r,r-i)
			double rho2 = peyton::sqr(sxy) / (sx2*sy2);
			denom = 1.-rho2;
			double lnN = -log(2. * peyton::ctn::pi * sqrt(denom*sx2*sy2));
			double riml = evaluate(y0);
			if(lnL != NULL)
			{
				*lnL  = likelihood(as.gr(riml)-x0, riml-y0);
			}
#if 0
			double dx = 0.01; double sum = 0;
			for(double a = -.5; a < .5; a += dx)
			{
				for(double b = -.5; b < .5; b += dx)
				{
					sum += dx*dx * exp(likelihood(a, b) + lnN);
				}
			}
			ASSERT(abs(sum-1) < 1e-5) {
				std::cerr << "sum = " << sum << "\n";
				std::cerr << ri << " " << gr << " " << sg << " " << sr << " " << si << "\n";
			}
#endif
			return riml;
		}

		double likelihood(double x, double y)
		{
			double z =
				   peyton::sqr(x)/sx2
			     - 2*sxy*(x)*(y)/(sx2*sy2)
			     + peyton::sqr(y)/sy2;	
			
			double lnL = -z/(2.*denom);
			
			return lnL + lnN;
		}

		double fn_to_minimize(double y)
		{
			// Error ellipse equation
			double x = as.gr(y);
			double lnPrior = prior != NULL ? gsl_spline_eval(prior, y, acc) : 0;
			double lnL = likelihood(x-x0, y-y0);
			return -(lnL + lnPrior);
		}

		~ml_grri()
		{
			if(prior != NULL) { gsl_spline_free(prior); }
			if(acc != NULL) { gsl_interp_accel_free(acc); }
		}
	} mlri;

	void ml_magnitudes(float &gml, float &rml, float &iml,
		const float g, const float r, const float i,
		const float gErr, const float rErr, const float iErr,
		const float ml_ri, const float ml_gr)
	{
		float wg = 1./peyton::sqr(gErr);
		float wr = 1./peyton::sqr(rErr);
		float wi = 1./peyton::sqr(iErr);

		rml = (wr*r + wg*(g-ml_gr) + wi*(i+ml_ri)) / (wr + wi + wg);
		iml = rml - ml_ri;
		gml = rml + ml_gr;
	}

public:
/*	plx_gri_locus() : mlri(*this) { Mrc[0] = 4.6; Mrc[1] = 7.9; Mrc[2] = -3.0; Mrc[3] = 0.69; }*/
	std::valarray<double> Mrc;
	plx_gri_locus() : mlri(*this), Mrc(5)
		{ Mrc[0] = 3.2; Mrc[1] = 13.30; Mrc[2] = -11.50; Mrc[3] = 5.40; Mrc[4] = -0.65; } // this is ZI's "kinematic" relation
//		{ Mrc[0] = 4.0; Mrc[1] = 11.86; Mrc[2] = -10.74; Mrc[3] = 5.99; Mrc[4] = -1.2; } // this was the relation from astro-ph draft

	void distance_limits(double &Dmin, double &Dmax, float ri0, float ri1, float r_min, float r_max)
	{
		Dmin = stardist::D(r_min, Mr(ri0));
		Dmax = stardist::D(r_max, Mr(ri1));
	}

	//
	// return M_r absolute magnitude (from ZI's unified locus fits)
	//
	double Mr(const float ri)
	{
/*		if(ri < 0.0328631)	// assume these are horizontal branch stars if g-r < 0.2
		{
			// from Sirko et al., Mg ~ .69 for HB stars
			return .69 - gr(ri);
		}
		else
		{*/
			// return 4.5+(5.9-0.6*ri)*ri;
			// return 4.6+7.9*ri-3.0*pow(ri, 2)+0.69*pow(ri, 3);
			// return Mrc[0] + Mrc[1]*ri + Mrc[2]*pow(ri, 2) + Mrc[3]*pow(ri, 3);
			return gsl_poly_eval(&Mrc[0], Mrc.size(), ri);
/*		}*/
	}

	//
	// define abstracts from photometric_paralax class
	// calculate the distance and absolute magnitude of this star
	//
	virtual bool operator()(sdss_star &s)
	{
		try
		{
			// calculate maximum likelihood colors on the locus
			float ml_ri = mlri(s.ri(), s.gr(), s.gErr, s.rErr, s.iErr);
			float ml_gr = gr(ml_ri);

			// calculate ML magnitudes
			ml_magnitudes(s.ml_g, s.ml_r, s.ml_i, s.g, s.r, s.i, s.gErr, s.rErr, s.iErr, ml_ri, ml_gr);

			// absolute magnitude and distance for the max.like. r-i color
			s.Mr = Mr(ml_ri);
			s.earth.D = stardist::D(s.ml_r, s.Mr);

			ASSERT(fabs(s.ml_ri() - ml_ri) < 0.001);
			ASSERT(fabs(s.ml_gr() - ml_gr) < 0.001);

			return true;
		}
		catch(peyton::exceptions::EGSLMinimizer &e)
		{
			return false;
		}
	}

	double ml_r_band(float ri, float gr, float sg, float sr, float si, float *lnL = NULL)
	{
		double RI;
		try {
			RI = this->mlri(ri, gr, sg, sr, si, lnL);
		}
		catch(peyton::exceptions::EGSLMinimizer &e)
		{
			// the observed colors are inconsistent with stellar locus
			return -1;
		}
		return RI;
	}
};

#endif
