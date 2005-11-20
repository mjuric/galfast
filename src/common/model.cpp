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

using namespace std;

/////////////

const char *disk_model::param_name[disk_model::nparams] = { "rho0", "l", "h", "z0", "f", "lt", "ht", "fh", "q", "n" };
const char *disk_model::param_format[disk_model::nparams] = { "%.9f", "%.9f", "%.9f", "%.9f", "%.9f", "%.9f", "%.9f", "%.9f", "%.9f" };

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

void model_fitter::cull(double nSigma)
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

void model_fitter::print(ostream &out, int format)
{
	switch(format)
	{
	case PRETTY:
		out << io::format("%15s = %d") << "n(DOF)" << ndof << "\n";
		FOR(0, nparams)
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
		FOR(0, nparams)
		{
			out << io::format(param_format[i]) << p[i];
			if(i != nparams-1) { out << " "; }
		}
		out << "    ";
		FOR(0, nparams)
		{
			out << io::format(param_format[i]) << sqrt(variance(i));
			if(i != nparams-1) { out << " "; }
		}
	}
}

/// model_factory class

model_factory::model_factory(const std::string &modelsfile)
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
	load(lfin, ri, 0, phi, 1);
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
	dm->rho0 = lf(ri)*dri;
	//dm->f = 0;
#endif
	return dm.release();
}


spline::spline(const double *x, const double *y, int n)
	: f(NULL), acc(NULL)
{
	construct(x, y, n);
}

spline& spline::operator= (const spline& a)
{
	construct(&a.xv[0], &a.yv[0], a.xv.size());
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

	// construct spline
	f = gsl_interp_alloc(gsl_interp_cspline, n);
	gsl_interp_init(f, &xv[0], &yv[0], n);
	acc = gsl_interp_accel_alloc();
}

spline::~spline()
{
	gsl_interp_accel_free(acc);
	gsl_interp_free(f);
}

