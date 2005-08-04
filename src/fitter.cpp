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

#include "analysis.h"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include <astro/system/options.h>
#include <astro/exceptions.h>
#include <astro/util.h>
#include <astro/math.h>
#include <astro/system/log.h>
#include <astro/io/format.h>
#include <astro/useall.h>

using namespace std;

const double Rg = 8000.;

struct rzpixel
{
	double r, z, N, V;
	double rho, sigma;
};

struct model
{
	// model parameters
	vector<double> p;
	vector<double> covar;
	vector<bool> fixed;
	double chi2_per_dof;

	double variance(int i) { return covar[i*p.size() + i]; }

	vector<string> param_name, param_format;
	std::map<string, int> param_name_to_index;

	// Model data
	vector<rzpixel> *map;
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
	double rho_thin(double r, double z)  const { return rho0 *     exp((Rg-r)/l  + (abs(z0) - abs(z + z0))/h); }
	double rho_thick(double r, double z) const { return rho0 * f * exp((Rg-r)/lt + (abs(z0) - abs(z + z0))/ht); }
	double rho_halo(double r, double z)  const { return rho0 * fh * pow(Rg/sqrt(halo_denom(r,z)),n); }
	double rho(double r, double z)       const { return rho_thin(r, z) + rho_thick(r, z) + rho_halo(r, z); }

	double norm_at_Rg() const { return f*exp(8000.*(1./l - 1./lt)); }

	// Derivatives of the model function
	double drho0(double r, double z, double rhom) const { return 1./rho0 * rhom; }
	double dl(double r, double z, double rhothin) const { return r/sqr(l) * rhothin; }
	double dh(double r, double z, double rhothin) const { return (-abs(z0)+abs(z+z0))/sqr(h) * rhothin; }
	double dz0(double r, double z, double rhothin, double rhothick) const { return (sgn(z0)-sgn(z+z0))*(rhothin/h + rhothick/ht); }
	double df(double r, double z, double rhothick) const { return 1./f * rhothick; }
	double dlt(double r, double z, double rhothick) const { return r/sqr(lt) * rhothick; }
	double dht(double r, double z, double rhothick) const { return (-abs(z0)+abs(z+z0))/sqr(ht) * rhothick; }
	// -- Halo derivatives assume z0 << z (which is why there's no halo component in dz0()
	double halo_denom(double r, double z) const { return sqr(r) + sqr(q)*sqr(z + z0); }
	double dfh(double r, double z, double rhoh) const { return 1./fh * rhoh; }
	double dq(double r, double z, double rhoh) const { return -n*q*sqr(z+z0)/halo_denom(r,z) * rhoh; }
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

	void add_param(const std::string &name_ = "", double val = 0, bool fixed = false, const std::string &format = "%.9f")
	{
		p.push_back(val);

		this->fixed.push_back(fixed);
		ndof += fixed == false;

		std::string name;
		name = name_.size() ? name_ : (io::format("param_%02d") << p.size()-1);
		param_name.push_back(name);

		param_name_to_index[name] = p.size() - 1;
		
		param_format.push_back(format);
	}
	double &param(const std::string &name)
	{
		return p[param_name_to_index[name]];
	}

	void get_parameters(gsl_vector *x)
	{
		int k = 0;
		FOR(0, p.size())
		{
			if(fixed[i]) continue;
			gsl_vector_set(x, k++, p[i]);
		}
	}
	void set_parameters(const gsl_vector *x)
	{
		int k = 0;
		FOR(0, p.size())
		{
			if(fixed[i]) continue;
			p[i] = gsl_vector_get(x, k++);
		}
	}

	int fdf(gsl_vector *f, gsl_matrix *J);
	int fit(int cullIter=1, double nsigma = 3.);
	void cull(double nSigma);

	enum {PRETTY, HEADING, LINE};
	void print(ostream &out, int format = PRETTY);
};

#define DFINIT int pcnt_ = 0, j_ = 0;
#define DFCALC(val) if(!fixed[pcnt_++]) gsl_matrix_set(J, i, j_++, (val)/x.sigma);
int model::fdf (gsl_vector * f, gsl_matrix * J)
{
	// calculate f_i values for all datapoints
	FOR(0, map->size())
	{
		const rzpixel &x = (*map)[i];
		double rhothin = rho_thin(x.r, x.z);
		double rhothick = rho_thick(x.r, x.z);
		double rhohalo = rho_halo(x.r, x.z);
		double rhom = rhothick + rhothin + rhohalo;

		if(f)
		{
			double df = rhom - x.rho;
			gsl_vector_set(f, i, df/x.sigma);
		}

		if(J)
		{
			DFINIT;
			DFCALC(drho0(x.r, x.z, rhom));
			DFCALC(dl(x.r, x.z, rhothin));
			DFCALC(dh(x.r, x.z, rhothin));
			DFCALC(dz0(x.r, x.z, rhothin, rhothick));
			DFCALC(df(x.r, x.z, rhothick));
			DFCALC(dlt(x.r, x.z, rhothick));
			DFCALC(dht(x.r, x.z, rhothick));
			DFCALC(dfh(x.r, x.z, rhohalo));
			DFCALC(dq(x.r, x.z, rhohalo));
			DFCALC(dn(x.r, x.z, rhohalo));
		}
	}

	return GSL_SUCCESS;
}
#undef DFCALC
#undef DFINIT

void model::cull(double nSigma)
{
	// calculate f_i values for all datapoints
	vector<rzpixel> newmap;
	FOR(0, map->size())
	{
		const rzpixel &x = (*map)[i];
		double rhom = rho(x.r, x.z);
		if(abs(x.rho - rhom) <= nSigma*x.sigma)
		{
			newmap.push_back(x);
		}
	}
	cerr << "Selected " << newmap.size() << " out of " << map->size() << " pixels\n";
	*map = newmap;
}

extern "C"
int model_fdf (const gsl_vector * v, void *params, gsl_vector * f, gsl_matrix * J)
{
	model *m = (model *)params;

	m->set_parameters(v);
	m->fdf(f, J);
}

extern "C"
int model_f (const gsl_vector * v, void *params, gsl_vector * f)
{
	return model_fdf(v, params, f, NULL);
}

extern "C"
int model_df (const gsl_vector * v, void *params, gsl_matrix * J)
{
	return model_fdf(v, params, NULL, J);
}

int print_state (size_t iter, gsl_multifit_fdfsolver * s, int dof)
{
	fprintf (stderr, "iter: %3u x = ", iter);
	FOR(0, dof)
	{
	  	fprintf(stderr, "%15.8f ", gsl_vector_get(s->x, i));
	}
	fprintf(stderr, "|f(x)| = %g\n", gsl_blas_dnrm2 (s->f));
}

int model::fit(int cullIter, double nsigma)
{
	//
	// Generic fit driver
	//
	int n = ndata();	// number of data points

	// allocate arrays, copy non-constant parameters into the array
	gsl_vector *v = gsl_vector_alloc(ndof);
	get_parameters(v);

	// initialize and set up the solver
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (gsl_multifit_fdfsolver_lmsder, n, ndof);

	gsl_multifit_function_fdf f;
	f.f = &model_f;
	f.df = &model_df;
	f.fdf = &model_fdf;
	f.n = n;
	f.p = ndof;
	f.params = (void *)this;
	gsl_multifit_fdfsolver_set (s, &f, v);

	// iterate
	int status = GSL_CONTINUE;
	for(int iter = 0; status == GSL_CONTINUE && iter < 500; iter++)
	{
		status = gsl_multifit_fdfsolver_iterate (s);

		print_state (iter, s, ndof);
		if (status && status != GSL_CONTINUE) { break; }; status = 0;
		//if (status) { break; }

		status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
	}
	fprintf (stderr, "status = %s\n", gsl_strerror (status));
	ASSERT(status == 0);

	// extract fitted parameters
	set_parameters(s->x);
	// expand the covariance matrix
	gsl_matrix *cov = gsl_matrix_alloc (ndof, ndof);
	gsl_multifit_covar (s->J, 0.0, cov);
	int np = p.size();
	covar.resize(np*np); fill(covar.begin(), covar.end(), 0.);
	int y = 0;
	FORj(r, 0, np) // rows
	{
		if(fixed[r]) continue;
		int x = 0;
		FORj(c, 0, np) // columns
		{
			if(fixed[c]) continue;
			covar[r*np + c] = gsl_matrix_get(cov, y, x);
			x++;
		}
		y++;
	}

	// chi^2/DOF
	double chi = gsl_blas_dnrm2(s->f);
	chi2_per_dof = sqr(chi)/ (n - ndof);

	// free the memory	
	gsl_matrix_free(cov);
	gsl_multifit_fdfsolver_free (s);

	if(cullIter--)
	{
		set_parameters(v);
		cull(nsigma);
		fit(cullIter, nsigma);
	}
	gsl_vector_free(v);
}

void model::print(ostream &out, int format)
{
	switch(format)
	{
	case PRETTY:
		out << io::format("%15s = %d") << "n(DOF)" << ndof << "\n";
		FOR(0, p.size())
		{
			out << io::format(std::string("%15s = ") + param_format[i]) << param_name[i] << p[i];
			out << " +- " << io::format(param_format[i]) << sqrt(variance(i));
			out << (fixed[i] ? " (const)" : " (var)");
			out << "\n";
		}
		out << io::format("%15s = %.5g") << "chi^2/dof" << chi2_per_dof << "\n";
		break;
	case HEADING:
		out << "# ";
		FOR(0, p.size())
		{
			out << param_name[i] << " ";
		}
		out << "\n# ";
		FOR(0, fixed.size())
		{
			out << (fixed[i] ? "const" : "var") << " ";
		}
		break;
	case LINE:
		FOR(0, p.size())
		{
			out << io::format(param_format[i]) << p[i];
			if(i != p.size()-1) { out << " "; }
		}
		out << "    ";
		FOR(0, p.size())
		{
			out << io::format(param_format[i]) << sqrt(variance(i));
			if(i != p.size()-1) { out << " "; }
		}
	}
}

void load_disk(vector<rzpixel> *data, const std::string &filename)
{
	data->clear();
	text_input_or_die(in, filename);

	rzpixel p;
	//bind(in, p.r, 0, p.z, 1, p.N, 2, p.V, 3); // if loading from cleanedrz.txt
	bind(in, p.r, 0, p.z, 1, p.N, 5, p.V, 6); // if loading from output of median3d.pl
	while(in.next())
	{
		p.rho = p.N / p.V;
		p.sigma = sqrt(p.N) / p.V;
		data->push_back(p);
	}
};

void clean_disk(vector<rzpixel> *data, const std::string &how)
{
	//
	// Remove the points near the plane of the galaxy
	// and close to the edges of the survey
	//
	vector<rzpixel> out;
	FOREACH(vector<rzpixel>::iterator, *data)
	{
		rzpixel &pix = *i;
		double phi = deg(atan2(pix.z,pix.r-8000));

 		if(-30 < phi && phi < 15) continue;
 		if(153 < phi && phi < 180) continue;
 		if(-180 < phi && phi < -150) continue;

		if(how == "thin")
		{
			if(!(abs(pix.z) >= 50 && abs(pix.z) <= 600)) continue;
		}
		else if(how == "thick")
		{
//			if(pix.z < 2500 && pix.z > 1000 && pix.r < 8000) continue;
//			if(pix.z < -2000 && pix.r > 8000) continue;

// /*			if(!(abs(pix.z) >= 1000 && abs(pix.z) <= 4000)) continue;
// 			if(!(abs(pix.z) >= 200 && abs(pix.z) <= 3000)) continue;*/
			
//			if(!(abs(pix.z) >= 200 && abs(pix.z) <= 1400)) continue;

			if(abs(pix.z) <= 500) continue;
			if(abs(pix.z) >= 2500) continue;
//			if(!(abs(pix.z) >= 200 && abs(pix.z) <= 3500)) continue;
			if(pix.z > 0 && deg(atan2(pix.z-500,pix.r-8000)) < 30. && pix.r < 10000) continue;
			if(pix.z > 0 && pix.r < 8000.) continue;
		}
		else if(how == "halo")
		{
			if(pix.z > 0 && deg(atan2(pix.z-500,pix.r-8000)) < 30. && pix.r < 10000) continue;
			if(pix.z > 0 && pix.z < 5000 && pix.r < 8000.) continue;

			if(pix.r < 5000) continue;
						
			if(!(abs(pix.z) >= 500)) continue;
		}

		out.push_back(pix);
	}

	cerr << data->size() - out.size() << " pixels rejected.\n";
	cerr << out.size() << " pixels accepted.\n";
	*data = out;
}

int main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version);

	Options opts(
		"Program for fitting r-z dataplanes of the Galaxy",
		version,
		Authorship::majuric
	);

	//# add any arguments your program needs. eg:
	opts.argument("binsfile", "File from which to read a list of rz-dataplanes to fit.");

	// add any options your program might need. eg:
	// opts.option("meshFactor", "meshFactor", 0, "--", Option::required, "4", "Resolution decrease between radial steps");

	try {
		opts.parse(argc, argv);
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}

	/////// Start your application code here
	std::string binsfile = opts["binsfile"];
	text_input_or_die(in, binsfile);

	std::string how = "halo";
//	std::string how = "thick";
//	std::string how = "thin";

	model m;
	double r0, r1; std::string rzfile, modelname;
	if(how == "thin")
	{
		m.add_param("rho0", 0, false);
		m.add_param("l", 0, false);
		m.add_param("h", 0, false);
		m.add_param("z0", 0, true);
		
		m.add_param("f", 0.0, true);
		m.add_param("lt", 3500, true);
		m.add_param("ht", 400, true);
		m.add_param("fh", 0.0, true);
		m.add_param("q", 1.5, true);
		m.add_param("n", 3, true);

		bind(in, modelname, 0, r0, 1, r1, 2, rzfile, 3,
			m.param("rho0"), 4,
			m.param("l"), 5,
			m.param("h"), 6,
			m.param("z0"), 7
			);
	}
	else if(how == "thick")
	{
		m.add_param("rho0", 0, true);
		m.add_param("l", 0, false);
		m.add_param("h", 0, false);
		m.add_param("z0", 0, true);

//		m.add_param("f", 0.012, true);
		m.add_param("f", 0.03, false);
		m.add_param("lt", 4800, false);
		m.add_param("ht", 1500, false);

		m.add_param("fh", 0.0, true);
		m.add_param("q", 1.5, true);
		m.add_param("n", 3, true);

		bind(in, modelname, 0, r0, 1, r1, 2, rzfile, 3,
			m.param("rho0"), 4,
			m.param("l"), 5,
			m.param("h"), 6,
			m.param("z0"), 7
			);
	}
	else if(how == "halo")
	{
		m.add_param("rho0", 0, false);
		m.add_param("l", 0, true);
		m.add_param("h", 0, true);
		m.add_param("z0", 0, true);
		
		m.add_param("f", 0.0647, true);
		m.add_param("lt", 3300, true);
		m.add_param("ht", 1100, true);

		m.add_param("fh", 0.0015, false);
		m.add_param("q", 2.3, false);
		m.add_param("n", 1.16, true);

		bind(in, modelname, 0, r0, 1, r1, 2, rzfile, 3,
			m.param("rho0"), 4,
			m.param("l"), 5,
			m.param("h"), 6,
			m.param("z0"), 7
			);
	}

	m.print(cout, model::HEADING);
	cout << "\n";
//	cout << "# name ri0 ri1 chi2/dof rho0 l h z0 err(rho0 l h z0)\n";

	vector<rzpixel> data;
	gsl_vector *v = gsl_vector_alloc(m.ndof);
	while(in.next())
	{
		m.get_parameters(v);

		load_disk(&data, rzfile);
		clean_disk(&data, how);
		m.map = &data;

//		m.param("l") = 1.;

		cerr << "Fitting " << rzfile << "\n";
		m.fit(0, 5);
		m.print(cerr);
		cerr << "norm_thick = " << m.norm_at_Rg() << "\n";
		cerr << "rho(R=8kpc,Z=5kpc) = " << m.rho(8000,5000) << "\n";
		cerr << "rho1(R=8kpc,Z=5kpc) = " << m.rho_thin(8000,5000) << "\n";
		cerr << "rho2(R=8kpc,Z=5kpc) = " << m.rho_thick(8000,5000) << "\n";
		cerr << "rhoh(R=8kpc,Z=5kpc) = " << m.rho_halo(8000,5000) << "\n";

		cout << setw(10) << modelname << setw(7) << r0 << setw(7) << r1 << setw(40) << rzfile;
		cout << " " << setw(10) << m.chi2_per_dof << " ";
		m.print(cout, model::LINE);
		cout << "\n";

		m.set_parameters(v);
	}
}
catch(EAny &e)
{
	e.print();
}
}
