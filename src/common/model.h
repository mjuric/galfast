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

#ifndef model_h__
#define model_h__

#include <vector>
#include <map>
#include <iosfwd>
#include <cmath>
#include <string>
#include <valarray>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <astro/types.h>
#include <astro/system/config.h>
#include <astro/system/log.h>

#include "paralax.h"

struct rzpixel
{
	double r, rphi, z, N, V;
	double rho, sigma;
};

struct disk_model
{
	static const double Rg = 8000.;
	
	static const size_t nparams = 10;
	static const char *param_name[nparams];
	static const char *param_format[nparams];

	union {
		double p[nparams];
		struct { double rho0, l, h, z0, f, lt, ht, fh, q, n; };
	};

	disk_model() {}

	// Model functions
	double rho_thin(double r, double z)  const { return rho0 *     exp((Rg-r)/l  + (std::abs(z0) - std::abs(z + z0))/h); }
	double rho_thick(double r, double z) const { return rho0 * f * exp((Rg-r)/lt + (std::abs(z0) - std::abs(z + z0))/ht); }
	double rho_halo(double r, double z)  const { return rho0 * fh * pow(Rg/sqrt(halo_denom(r,z)),n); }
	double rho(double r, double z)       const { return rho_thin(r, z) + rho_thick(r, z) + rho_halo(r, z); }

	//double norm_at_Rg() const { return f*exp(Rg*(1./l - 1./lt)); }
	double norm_at_Rg() const { return rho_thick(Rg, 0)/rho_thin(Rg, 0); }

	// Derivatives of the model function
	double drho0(double r, double z, double rhom) const { return 1./rho0 * rhom; }
	double dl(double r, double z, double rhothin) const { return r/peyton::sqr(l) * rhothin; }
	double dh(double r, double z, double rhothin) const { return (-std::abs(z0)+std::abs(z+z0))/peyton::sqr(h) * rhothin; }
	double dz0(double r, double z, double rhothin, double rhothick) const { return (peyton::sgn(z0)-peyton::sgn(z+z0))*(rhothin/h + rhothick/ht); }
	double df(double r, double z, double rhothick) const { return 1./f * rhothick; }
	double dlt(double r, double z, double rhothick) const { return r/peyton::sqr(lt) * rhothick; }
	double dht(double r, double z, double rhothick) const { return (-std::abs(z0)+std::abs(z+z0))/peyton::sqr(ht) * rhothick; }
	// -- Halo derivatives assume z0 << z (which is why there's no halo component in dz0()
	double halo_denom(double r, double z) const { return peyton::sqr(r) + peyton::sqr(q)*peyton::sqr(z + z0); }
	double dfh(double r, double z, double rhoh) const { return 1./fh * rhoh; }
	double dq(double r, double z, double rhoh) const { return -n*q*peyton::sqr(z+z0)/halo_denom(r,z) * rhoh; }
	double dn(double r, double z, double rhoh) const { return log(Rg/sqrt(halo_denom(r,z))) * rhoh; }
};

struct model_fitter : public disk_model
{
	// model parameters
	std::vector<double> covar;
	std::vector<bool> fixed;
	double chi2_per_dof;
	double epsabs, epsrel;

	double variance(int i) { return covar.size() ? covar[i*nparams + i] : 0; }

	std::map<std::string, int> param_name_to_index;

	// Model data
	std::vector<rzpixel> *map, culled;
	std::pair<float, float> ri, r;
	std::pair<double, double> d;
	int ndata() { return map->size(); }

	model_fitter(const model_fitter& m)
		: disk_model(m),
		  covar(m.covar), fixed(m.fixed), chi2_per_dof(m.chi2_per_dof),
		  param_name_to_index(m.param_name_to_index),
		  map(m.map != NULL ? new std::vector<rzpixel>(*m.map) : NULL),
		  epsabs(1e-6), epsrel(1e-6)
		{
		}

	// Constructor
	model_fitter()
		: fixed(nparams, false)
	{
		for(int i = 0; i != nparams; i++) { param_name_to_index[param_name[i]] = i; }
	}

	// Generic model_fitter fitting functions
	int ndof() {
		int ndof = 0;
		FOR(0, fixed.size()) {
			if(!fixed[i]) ndof++;
		};
		return ndof;
		}

 	double &param(const std::string &name)
 	{
 		return p[param_name_to_index[name]];
	}
	std::vector<bool>::reference fix(const std::string &name)
	{
		return fixed[param_name_to_index[name]];
	}
 	void set_param(const std::string &name, double val, bool fixed)
	{
		param(name) = val;
		fix(name) = fixed;
	}

	void get_parameters(gsl_vector *x);
	void set_parameters(const gsl_vector *x);

	int fdf(gsl_vector *f, gsl_matrix *J);
	int fit(int cullIter=1, double nsigma = 3.);
	void cull(double nSigma);
	void residual_distribution(std::map<int, int> &hist, double binwidth);

	enum {PRETTY, HEADING, LINE};
	void print(std::ostream &out, int format = PRETTY);
};

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
	void construct(const std::valarray<double> &x, const std::valarray<double> &y)
		{ ASSERT(x.size() == y.size()); construct(&x[0], &y[0], x.size()); }
	void construct(const std::vector<double> &x, const std::vector<double> &y)
		{ ASSERT(x.size() == y.size()); construct(&x[0], &y[0], x.size()); }
	~spline();

 	double operator ()(double x)        { return gsl_interp_eval(f, &xv[0], &yv[0], x, acc); }
 	double deriv(double x)              { return gsl_interp_eval_deriv(f, &xv[0], &yv[0], x, acc); }
 	double deriv2(double x)             { return gsl_interp_eval_deriv2(f, &xv[0], &yv[0], x, acc); }
 	double integral(double a, double b) { return gsl_interp_eval_integ(f, &xv[0], &yv[0], a, b, acc); }

public:
	spline& operator= (const spline& a);
	spline(const spline& a) : f(NULL), acc(NULL) { *this = a; }
};

/*class model_factory
{
public:
	std::vector<std::pair<std::pair<float, float>, disk_model> > models;
	spline lf;		// luminosity function
	//plx_gri_locus plx;	// paralax relation
public:
	model_factory(const std::string &models = "");
	void load(const std::string &models);
	disk_model *get(float ri, double dri);
};
*/
class galactic_model
{
public:
	virtual double absmag(double ri) = 0;
	virtual double rho(double x, double y, double z, double ri) = 0;
	
	static galactic_model *load(std::istream &cfg);
};

class BahcallSoneira_model : public galactic_model
{
public:
	disk_model m;
public:
	BahcallSoneira_model(peyton::system::Config &cfg);
	BahcallSoneira_model(const std::string &prefix);
	virtual double absmag(double ri);
	virtual double rho(double x, double y, double z, double ri);
protected:
	void load(peyton::system::Config &cfg);
};

class toy_homogenious_model : public galactic_model
{
public:
	double rho0;
public:
	toy_homogenious_model(double rho0_) : rho0(rho0_) {}
	double absmag(double ri);
	double rho(double x, double y, double z, double ri);
};

// geocentric powerlaw model with a constant paralax relation
class toy_geocentric_powerlaw_model : public galactic_model
{
public:
	double rho0, alpha;
public:
	toy_geocentric_powerlaw_model(double rho0_, double alpha_) : rho0(rho0_), alpha(alpha_) {}
	double absmag(double ri);
	double rho(double x, double y, double z, double ri);
};

// geocentric powerlaw model with a polynomial Mr(ri) paralax relation
// reads its parameters from a config file, as keywords with prefix 'Mr_'
// for Mr_coef, and keywords rho0 and alpha for the powerlaw params.
class toy_geo_plaw_abspoly_model : public galactic_model
{
public:
	double rho0, alpha;
	std::vector<double> Mr_coef;
public:
	toy_geo_plaw_abspoly_model(const std::string &prefix);
	double rho(double x, double y, double z, double ri);
	double absmag(double ri);
};

#endif
