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
#include <iostream>

#include <astro/system/options.h>
#include <astro/exceptions.h>
#include <astro/util.h>
#include <astro/math.h>
#include <astro/system/log.h>
#include <astro/io/format.h>
#include <astro/useall.h>

using namespace std;

void model::get_parameters(gsl_vector *x)
{
	int k = 0;
	FOR(0, p.size())
	{
		if(fixed[i]) continue;
		gsl_vector_set(x, k++, p[i]);
	}
}
void model::set_parameters(const gsl_vector *x)
{
	int k = 0;
	FOR(0, p.size())
	{
		if(fixed[i]) continue;
		p[i] = gsl_vector_get(x, k++);
	}
}

void model::add_param(const std::string &name_, double val, bool fixed, const std::string &format)
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

