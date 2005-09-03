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

struct model
{
	static const double Rg = 8000.;

	// model parameters
	std::vector<double> p;
	std::vector<double> covar;
	std::vector<bool> fixed;
	double chi2_per_dof;

	double variance(int i) { return covar[i*p.size() + i]; }

	std::vector<std::string> param_name, param_format;
	std::map<std::string, int> param_name_to_index;

	// Model data
	std::vector<rzpixel> *map;
	int ndata() { return map->size(); }

	// Model function
	#define rho0	p[0]
	#define l	p[1]
	#define h	p[2]
	#define z0	p[3]
	#define f	p[4]
	#define lt	p[5]
	#define ht	p[6]
	#define fh	p[7]
	#define q	p[8]
	#define n	p[9]
	double rho_thin(double r, double z)  const { return rho0 *     exp((Rg-r)/l  + (std::abs(z0) - std::abs(z + z0))/h); }
	double rho_thick(double r, double z) const { return rho0 * f * exp((Rg-r)/lt + (std::abs(z0) - std::abs(z + z0))/ht); }
	double rho_halo(double r, double z)  const { return rho0 * fh * pow(Rg/sqrt(halo_denom(r,z)),n); }
	double rho(double r, double z)       const { return rho_thin(r, z) + rho_thick(r, z) + rho_halo(r, z); }

	double norm_at_Rg() const { return f*exp(8000.*(1./l - 1./lt)); }

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
	#undef rho0
	#undef l
	#undef h
	#undef z0
	#undef f
	#undef lt
	#undef ht
	#undef fh
	#undef q
	#undef n

	// Constructor
	model() : ndof(0) {}

	// Generic model fitting functions
	int ndof;

	void add_param(const std::string &name_ = "", double val = 0, bool fixed = false, const std::string &format = "%.9f");
	double &param(const std::string &name)
	{
		return p[param_name_to_index[name]];
	}

	void get_parameters(gsl_vector *x);
	void set_parameters(const gsl_vector *x);

	int fdf(gsl_vector *f, gsl_matrix *J);
	int fit(int cullIter=1, double nsigma = 3.);
	void cull(double nSigma);

	enum {PRETTY, HEADING, LINE};
	void print(std::ostream &out, int format = PRETTY);
};



#endif
