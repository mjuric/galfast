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
 
#include "model.h"
#include "textstream.h"
#include "analysis.h"

#include <iostream>
#include <fstream>

#include <astro/system/options.h>
#include <astro/exceptions.h>
#include <astro/util.h>
#include <astro/math.h>
#include <astro/system/log.h>
#include <astro/io/format.h>
#include <astro/useall.h>

#include <gsl/gsl_poly.h>

using namespace std;

/////////////

const char *disk_model::param_name[disk_model::nparams] =
{
	"rho0", "l", "h", "z0", "f", "lt", "ht", "fh", "q", "n",
	"rho1", "rho2", "rho3", "rho4", "rho5", "rho6", "rho7", "rho8", "rho9", "rho10", 
};
const char *disk_model::param_format[disk_model::nparams] = 
{
	"%.5f", "%.0f", "%.0f", "%.2f", "%.3f", "%.0f", "%.0f", "%.5f", "%.2f", "%.2f",
	"%.5f", "%.5f", "%.5f", "%.5f", "%.5f", "%.5f", "%.5f", "%.5f", "%.5f", "%.5f"
};

/////////////

void model_fitter::get_parameters(gsl_vector *x)
{
	int k = 0;
	FOR(0, nparams)
	{
		if(fixed[i]) continue;
		gsl_vector_set(x, k++, p[i]);
	}
}

void model_fitter::set_parameters(const gsl_vector *x)
{
	int k = 0;
	FOR(0, nparams)
	{
		if(fixed[i]) continue;
		p[i] = gsl_vector_get(x, k++);
	}
}

#define DFINIT int pcnt_ = 0, j_ = 0;
#define DFCALC(val) if(!fixed[pcnt_++]) gsl_matrix_set(J, i, j_++, (val)/x.sigma);
int model_fitter::fdf (gsl_vector * f, gsl_matrix * J)
{
	// calculate f_i values for all datapoints
	FOR(0, map.size())
	{
		const rzpixel &x = map[i];

		int ri = x.ri_bin;

		double rhothin = rho_thin(x.r, x.z, ri);
		double rhothick = rho_thick(x.r, x.z, ri);
		double rhohalo = rho_halo(x.r, x.z, ri);
		double rhom = rhothick + rhothin + rhohalo;

		if(f)
		{
			double df = rhom - x.rho;
			gsl_vector_set(f, i, df/x.sigma);
		}

		if(J)
		{
			DFINIT;
			DFCALC(drho0(x.r, x.z, rhom, ri, 0));
			DFCALC(dl(x.r, x.z, rhothin));
			DFCALC(dh(x.r, x.z, rhothin));
			DFCALC(dz0(x.r, x.z, rhothin, rhothick));
			DFCALC(df(x.r, x.z, rhothick));
			DFCALC(dlt(x.r, x.z, rhothick));
			DFCALC(dht(x.r, x.z, rhothick));
			DFCALC(dfh(x.r, x.z, rhohalo));
			DFCALC(dq(x.r, x.z, rhohalo));
			DFCALC(dn(x.r, x.z, rhohalo));
			DFCALC(drho0(x.r, x.z, rhom, ri, 1));
			DFCALC(drho0(x.r, x.z, rhom, ri, 2));
			DFCALC(drho0(x.r, x.z, rhom, ri, 3));
			DFCALC(drho0(x.r, x.z, rhom, ri, 4));
			DFCALC(drho0(x.r, x.z, rhom, ri, 5));
			DFCALC(drho0(x.r, x.z, rhom, ri, 6));
			DFCALC(drho0(x.r, x.z, rhom, ri, 7));
			DFCALC(drho0(x.r, x.z, rhom, ri, 8));
			DFCALC(drho0(x.r, x.z, rhom, ri, 9));
			DFCALC(drho0(x.r, x.z, rhom, ri, 10));
		}

#if 0
		std::cerr << x.r << " " << x.z << " " << ri << " " << ri2idx(ri) << " " << rho0[ri2idx(ri)] << ":  ";
		std::cerr << rhothin << " ";
		std::cerr << rhothick << " ";
		std::cerr << rhohalo << " ";
		std::cerr << rhom << ": ";
		FORj(j, 0, ndof())
		{
			std::cerr << gsl_matrix_get(J, i, j) << " ";
		}
		std::cerr << "\n";
#endif
	}
//	exit(0);

	return GSL_SUCCESS;
}
#undef DFCALC
#undef DFINIT

void model_fitter::cull(double nSigma)
{
	// calculate f_i values for all datapoints
	map.clear();
	culled.clear();
	FOR(0, orig.size())
	{
		const rzpixel &x = orig[i];
		double rhom = rho(x.r, x.z, x.ri_bin);
		if(abs(x.rho - rhom) <= nSigma*x.sigma)
		{
			map.push_back(x);
		} else {
			culled.push_back(x);
		}
	}
	cerr << "Selected " << map.size() << " out of " << orig.size() << " pixels\n";
}

void model_fitter::residual_distribution(std::map<int, int> &hist, double binwidth)
{
	// calculate f_i values for all datapoints
	FOR(0, map.size())
	{
		const rzpixel &x = map[i];
		double rhom = rho(x.r, x.z, x.ri_bin);
		double r = (x.rho - rhom) / x.sigma;
		int ir = (int)floor((r+0.5*binwidth) / binwidth);
// 		std::cerr << r << " " << binwidth*ir << "\n";

		if(hist.find(ir) == hist.end()) hist[ir] = 1;
		else hist[ir]++;
	}
}

void model_fitter::print(ostream &out, int format, int ri_bin)
{
	int riidx = ri2idx(ri_bin);
	switch(format)
	{
	case PRETTY:
		out << io::format("%15s = (%.3f, %.3f)") << "ri" << ri[ri_bin].first << ri[ri_bin].second << "\n";
		out << io::format("%15s = %d") << "n(DOF)" << ndof() << "\n";
		out << io::format("%15s = %.5g") << "chi^2/dof" << chi2_per_dof << "\n";
		out << io::format("%15s = %.5g") << "eps{abs,rel}" << epsabs << " " << epsrel << "\n";
		FOR(0, nparams)
		{
			out << io::format(std::string("%15s = ") + param_format[i]) << param_name[i] << p[i];
			out << " +- " << io::format(param_format[i]) << sqrt(variance(i));
			out << (fixed[i] ? " (const)" : " (var)");
			out << "\n";
		}
		out << "\n";
		if(covar.size())
		{
			FORj(r, -1, nparams) // rows
			{
				if(r == -1) { out << io::format("%15s = ") << "corr. matrix"; }
				else { out << io::format("%15s = ") << param_name[r]; }
	
				FORj(c, 0, nparams) // columns
				{
					if(r == -1) { out << io::format(" %10s") << param_name[c]; continue; }
					double corr = fixed[c] || fixed[r] ? 0 : covar[r*nparams + c] / sqrt(variance(c)*variance(r));
					std::cerr << io::format(" %10.3g") << corr;
				}
				out << "\n";
			}
		}
		break;
	case HEADING:
		out << "# ";
		FOR(0, nparams)
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
		// parameters
		FORj(k, 0, nparams - nrho)
		{
			int i = k == 0 ? riidx : k;
			out << io::format(param_format[i]) << p[i];
			if(i != nparams-1) { out << " "; }
		}
		out << "       ";
		// errors
		FORj(k, 0, nparams - nrho)
		{
			int i = k == 0 ? riidx : k;
			out << io::format(param_format[i]) << sqrt(variance(i));
			if(i != nparams-1) { out << " "; }
		}
	}
}

/// model_factory class

#if 0
model_factory::model_factory(const std::string &modelsfile)
{
	if(modelsfile.size()) { load(modelsfile); }
}

void model_factory::load(const std::string &modelsfile)
{
	text_input_or_die(in, modelsfile);

	float ri0, ri1;	
	disk_model dm;
	bind(in, ri0, 0, ri1, 1);
	FOR(0, disk_model::nparams) { bind(in, dm.p[i], i+2); }

	while(in.next())
	{
		models.push_back(make_pair(make_pair(ri0, ri1), dm));
	}

	// luminosity function
	text_input_or_die(lfin, "lumfun.txt")
	vector<double> ri, phi;
	::load(lfin, ri, 0, phi, 1);
	lf.construct(ri, phi);
}

disk_model *model_factory::get(float ri, double dri)
{
#if 0
	FOR(0, models.size())
	{
		std::pair<float, float> &bin = models[i].first;
		if(bin.first <= ri && ri <= bin.second)
		{
			disk_model *dm = new disk_model(models[i].second);
/*			if(!(0.5 <= ri && ri <= 0.52)) { dm->rho0 /= 1000; }*/
			return dm;
		}
	}
#else
	std::auto_ptr<disk_model> dm(new disk_model(models[0].second));
	dm->rho0 = lf(ri);
	//dm->f = 0;
#endif
	return dm.release();
}

#endif

spline::spline(const double *x, const double *y, int n)
	: f(NULL), acc(NULL)
{
	construct(x, y, n);
}

spline& spline::operator= (const spline& a)
{
	if(a.f != NULL && a.acc != NULL)
	{
		construct(&a.xv[0], &a.yv[0], a.xv.size());
	} else {
		f = NULL; acc = NULL;
	}
	return *this;
}

void spline::construct(const double *x, const double *y, int n)
{
//	gsl_interp_accel_free(acc);
//	gsl_interp_free(f);

	// copy data
	xv.resize(n); yv.resize(n);
	copy(x, x+n, &xv[0]);
	copy(y, y+n, &yv[0]);
	
	construct_aux();
}

void spline::construct_aux()
{
	// construct spline
	f = gsl_interp_alloc(gsl_interp_linear, xv.size());
	//f = gsl_interp_alloc(gsl_interp_cspline, n);
	gsl_interp_init(f, &xv[0], &yv[0], xv.size());
	acc = gsl_interp_accel_alloc();
}

spline::~spline()
{
	if(acc != NULL) gsl_interp_accel_free(acc);
	if(f != NULL) gsl_interp_free(f);
}

BOSTREAM2(const spline &spl)
{
	return out << spl.xv << spl.yv;
}

BISTREAM2(spline &spl)
{
	if(!(in >> spl.xv >> spl.yv)) return in;
	ASSERT(spl.xv.size() == spl.yv.size());
	if(spl.xv.size()) { spl.construct_aux(); }
	return in;
}

////////////////////////////////////////////////////

double toy_homogenious_model::rho(double x, double y, double z, double ri)
{
	return rho0;
}

double toy_homogenious_model::absmag(double ri)
{
	//return plx.Mr(ri);
	return 4;
}

////////////////////////////////////////////////////

double toy_geocentric_powerlaw_model::rho(double x, double y, double z, double ri)
{
	x -= Rg;
	double d2 = sqr(x) + sqr(y) + sqr(z);
	return rho0 * pow(d2, alpha/2.);
}

double toy_geocentric_powerlaw_model::absmag(double ri)
{
	return 4+2*ri;
}

////////////////////////////////////////////////////

toy_geo_plaw_abspoly_model::toy_geo_plaw_abspoly_model(const std::string &prefix)
{
	Config conf(prefix + ".conf");

	// Get model parameters
	conf.get(alpha, "alpha", 0.);
	conf.get(rho0, "rho0", 1.);

	// Get absolute magnitude relation coefficients
	std::string key; ostringstream poly;
	double c;
	for(int i = 0; conf.get(c, key = std::string("Mr_") + str(i), 0.); i++)
	{
		Mr_coef.push_back(c);
		if(i == 0) { poly << c; }
		else { poly << " + " << c << "*x^" << i; }
	}

	for(double ri=0; ri < 1.5; ri += 0.1)
	{
		std::cerr << absmag(ri) << "\n";
	}

	LOG(app, verb1) << "rho(d) = " << rho0 << "*d^" << alpha;
	LOG(app, verb1) << "Mr(r-i) = " << poly.str();
}

double toy_geo_plaw_abspoly_model::rho(double x, double y, double z, double ri)
{
	// geocentric powerlaw distribution
#if 1
 	x -= Rg;
 	double d2 = sqr(x) + sqr(y) + sqr(z);
 	return rho0 * pow(d2, alpha/2.);
#else
	// geocentric shell distribution
	x -= Rg;
	double d2 = sqr(x) + sqr(y) + sqr(z);
	if(sqr(3000) < d2 && d2 < sqr(4000))
	{
		return rho0 * pow(d2, alpha/2.);
		//return pow(d2, -3./2.);
	}
	return 0;
	return 0.01*rho0 * pow(d2, alpha/2.);
#endif
}

double toy_geo_plaw_abspoly_model::absmag(double ri)
{
	// evaluate the polynomial
	ASSERT(0 <= ri && ri <= 1.5) { std::cerr << ri << "\n"; }
//	return 4;
//	return 4+2*ri-3;
	return gsl_poly_eval(&Mr_coef[0], Mr_coef.size(), ri);
}

void BahcallSoneira_model::load(peyton::system::Config &cfg)
{
	FOREACH(cfg) { std::cerr << "BS model: " << (*i).first << " = " << (*i).second << "\n"; }

	FOR(0, m.nparams - m.nrho)
	{
		std::string param = m.param_name[i];
		ASSERT(cfg.count(param)) { std::cerr << "Initial value for " << param << " not specified\n"; }

		m.p[i] = cfg[param];
	}
	
	// luminosity function
	if(cfg.count("lumfunc"))
	{
		ASSERT(cfg.count("rho0_ri"));
		rho0_ri = cfg["rho0_ri"];

		input_or_die(in, cfg["lumfunc"]);
		load_luminosity_function(in, rho0_ri);
	}
}

BahcallSoneira_model::BahcallSoneira_model()
{
}

BahcallSoneira_model::BahcallSoneira_model(peyton::system::Config &cfg)
{
	load(cfg);
}

BahcallSoneira_model::BahcallSoneira_model(const std::string &prefix)
{
	Config cfg(prefix);
	load(cfg);
}

double BahcallSoneira_model::absmag(double ri)
{
	return paralax.Mr(ri);
}

double BahcallSoneira_model::rho(double x, double y, double z, double ri)
{
	double r = sqrt(x*x + y*y);
	double norm = lf.empty() ? 1. : lf(ri);
//	norm = 1.;
	return norm * m.rho(r, z, 0);
}

void BahcallSoneira_model::load_luminosity_function(istream &in, std::pair<double, double> rho0_ri)
{
	// load the luminosity function and normalize to m.rho0.
	// rho0 is assumed to contain the number of stars per cubic parsec
	// per 1mag of r-i
	itextstream lfin(in);
	vector<double> ri, phi;
	::load(lfin, ri, 0, phi, 1);
	lf.construct(ri, phi);

	// make the LF dimensionless
	double dr = rho0_ri.second - rho0_ri.first;
	double stars_per_mag = lf.integral(rho0_ri.first, rho0_ri.second) / dr;
	FOREACH(phi) { *i /= stars_per_mag; };
	lf.construct(ri, phi);
	std::cerr << "Norm.: " << 1./stars_per_mag << "\n";
	std::cerr << "New int: " << lf.integral(rho0_ri.first, rho0_ri.second) / dr << "\n";
	std::cerr << "lf(1.0): " << lf(1.0) << "\n";
	std::cerr << "lf(1.1): " << lf(1.1) << "\n";
#if 0
	for(double ri=0; ri < 1.5; ri += 0.005)
	{
		std::cerr << ri << " " << lf(ri) << "\n";
	}
#endif
}

BLESS_POD(disk_model);
peyton::io::obstream& BahcallSoneira_model::serialize(peyton::io::obstream& out)
{
	std::string name("BahcallSoneira");
	out << name << m << lf;

	return out;
}
BISTREAM2(BahcallSoneira_model &m)
{
	return in >> m.m >> m.lf;
}

galactic_model *galactic_model::load(istream &cfgstrm)
{
	Config cfg;
	cfg.load(cfgstrm);

	FOREACH(cfg) { std::cerr << (*i).first << " = " << (*i).second << "\n"; }	
	if(cfg.count("model") == 0) { ASSERT(0); return NULL; }

	std::string model = cfg["model"];

	if(model == "BahcallSoneira") { return new BahcallSoneira_model(cfg); }

	ASSERT(0); return NULL;
}

peyton::io::obstream& galactic_model::serialize(peyton::io::obstream& out)
{
	ASSERT(0) { std::cerr << "This method must be overrided if you want to use serialization.\n"; }
	return out;
}

galactic_model *galactic_model::unserialize(peyton::io::ibstream &in)
{
	std::string model;
	in >> model;

	if(model == "BahcallSoneira")
	{
		std::auto_ptr<BahcallSoneira_model> m(new BahcallSoneira_model);
		in >> *m;
		return m.release();
	}

	ASSERT(0); return NULL;
}
