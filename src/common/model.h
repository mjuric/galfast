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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <astro/types.h>

struct rzpixel
{
	double r, z, N, V;
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

	double norm_at_Rg() const { return f*exp(Rg*(1./l - 1./lt)); }

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

	double variance(int i) { return covar[i*nparams + i]; }

	std::map<std::string, int> param_name_to_index;

	// Model data
	std::vector<rzpixel> *map;
	int ndata() { return map->size(); }

	model_fitter(const model_fitter& m)
		: disk_model(m),
		  covar(m.covar), fixed(m.fixed), chi2_per_dof(m.chi2_per_dof),
		  param_name_to_index(m.param_name_to_index),
		  map(m.map != NULL ? new std::vector<rzpixel>(*m.map) : NULL)
		{
		}

	// Constructor
	model_fitter()
		: ndof(0), fixed(false, nparams)
	{
		for(int i = 0; i != nparams; i++) { param_name_to_index[param_name[i]] = i; }
	}

	// Generic model_fitter fitting functions
	int ndof;

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

	enum {PRETTY, HEADING, LINE};
	void print(std::ostream &out, int format = PRETTY);
};

class model_factory
{
public:
	std::vector<std::pair<std::pair<float, float>, disk_model> > models;
public:
	model_factory(const std::string &models);
	disk_model *get(float ri);
};

#endif
