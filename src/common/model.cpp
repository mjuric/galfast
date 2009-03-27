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
#define COMPILING_SIMULATE
 
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

#if HAVE_BOOST_REGEX
#include <boost/regex.hpp>
#endif

using namespace std;

// helpers for tag definitions
template<typename T>
sstruct::tagdef* createTag_aux(const std::string &tagName, sstruct::tagclass *tagClass, const size_t n = 1)
{
	if(n <= 1)
		return new sstruct::tagdefT<T>(tagName, tagClass);	// scalar
	else
		return new sstruct::tagdefTA<T>(tagName, n, tagClass);	// array
}

void parseAliases(std::vector<std::pair<std::string, size_t> > &aliases, const std::string &s)
{
	aliases.clear();

	// alias format: alias1,alias2,alias3=NNN,alias4...
	// where NNN is index into this tag's array (if the tag is an array)
	boost::regex e("(?:(\\w+)(?:=(\\d+))?)(?:,(\\w+)(?:=(\\d+))?)*");
	boost::smatch what;
	if(boost::regex_match(s, what, e, boost::match_extra))
	{
	#if 0
		unsigned i, j;
		std::cout << "** Match found **\n   Sub-Expressions:\n";
		std::cout << "Searched string: " << s << "\n";
		for(i = 0; i < what.size(); ++i)
			std::cout << "      $" << i << " = \"" << what[i] << "\"\n";
	#endif
		for(int i = 0; i < what.size(); i += 2)
		{
			if(i+1 >= what.size()) { break; }
			std::string alias = what[i].str();
			if(alias.size() == 0) { continue; }	// nothing was found

			std::string idxs = what[i].str();
			size_t idx = idxs.empty() ? 0 : atoi(idxs.c_str());
			aliases.push_back(make_pair(alias, idx));
		}
	}
}

sstruct::tagdef *sstruct::factory_t::createTag(const std::string &stringDef, size_t *offset_var)
{
	// parse stringDef and create a corresponding tag instance
	// the tag definition format is: name[N]{keyword1:xxx;keyword2:yyy;...}, where:
	//	[] part is optional -- if not specified, the quantity is a scalar
	//	N is the size of the array (integer)
	// The regex below matches: name, N, x as what[1],[2],[3], respectively
	//boost::regex e("(\\w+)(?:\\[(\\d+)\\])?(?:\\{([fgdis])(?:\\|(\\w+=[^,]+)(?:,(\\w+=[^,]+))*)\\})?");
	//boost::regex e("(\\w+)(?:\\[(\\d+)\\])?(?:\\{([fgdis])(?:\\|(?:\\s*(\\w+)\\s*=([^,]+))(?:,(?:\\s*(\\w+)\\s*=([^,]+)))*)?\\})?");
	boost::regex e("(\\w+)(?:\\[(\\d+)\\])?(?:\\{(?:(?:\\s*(\\w+)\\s*:([^;]+))(?:;(?:\\s*(\\w+)\\s*:([^;]+)))*)?\\})?");
	boost::smatch what;
	if(boost::regex_match(stringDef, what, e, boost::match_extra))
	{
	#if 0
		unsigned i, j;
		std::cout << "** Match found **\n   Sub-Expressions:\n";
		std::cout << "Searched string: " << stringDef << "\n";
		for(i = 0; i < what.size(); ++i)
			std::cout << "      $" << i << " = \"" << what[i] << "\"\n";
	#endif
		// slurp up the definition
		std::string tagName = what[1];
		size_t n = atoi(what[2].str().c_str());

		std::string className, fmt, type;
		std::vector<std::pair<std::string, size_t> > aliases;
		for(int i = 4; i < what.size(); i += 2)
		{
			if(i+1 >= what.size()) { break; }
			std::string key = what[i].str();
			if(key.size() == 0) { continue; }	// nothing was found

				if(key == "class") { className = what[i+1].str(); }
			else if(key == "fmt")   { fmt = what[i+1].str(); }
			else if(key == "type")  { type = what[i+1].str(); }
			else if(key == "alias") { parseAliases(aliases, what[i+1].str()); }
			else { THROW(EAny, "Unknown field definition key '" + key + "'."); }
		}

		MLOG(verb2) << "Defining tag from " << stringDef << ": type=" << type << ", class=" << className << ", fmt=" << fmt;

		// create a correctly subclassed tagdef
		if(type.empty()) { type = "float"; }
		tagclass *tagClass = getTagClass(className);
		tagdef *td = NULL;
		#define CREATETYPE(strtype, ctype) if(td == NULL && type == strtype) { td = createTag_aux<ctype>(tagName, tagClass, n); td->tagType = strtype; }
		CREATETYPE("int", int)
		CREATETYPE("double", double)
		CREATETYPE("string", std::string)
		CREATETYPE("float", float)
		#undef CREATETYPE
		if(td == NULL) { THROW(EAny, "Unknown tag data type '" + type + "'"); }

		// set overriden tag properties
		if(!fmt.empty()) { td->formatString = fmt; }

		// store the new tag definition
		definedTags[tagName] = td;
		if(offset_var) { td->offset_vars.push_back(offset_var); }

		// set aliases
		FOREACH(aliases)
		{
			aliasTag(tagName, i->first, i->second);
		}
	}
	else
	{
		MLOG(verb1) << "Failed to parse tag definition " << stringDef;
	}
}

sstruct::tagdef *sstruct::factory_t::aliasTag(const std::string &name, const std::string &alias, const size_t idx = -1)
{
	ASSERT(definedTags.count(name));
	tagdef *td = definedTags[name];
	ASSERT(td->offset != -1);		// the tag alias target has to be in use

	if(!definedTags.count(alias))
	{
		// autodefine tag by cloning the tag being aliased
		definedTags[alias] = idx < 0 ? td->clone(alias) : td->cloneScalar(alias);
	}

	tagdef *atd = definedTags[alias];
	ASSERT(atd->offset == -1);			// the alias must be unused
	ASSERT(atd->elementSize == td->elementSize);	// alias and the aliased tag must be of same element size
	ASSERT(idx < td->n);				// if aliasing a sub-element, must be within bounds

	atd->aliasingIndex = idx;
	atd->offset = td->offset + td->elementSize*idx;	// let the alias know where in the index it is
	FOREACH(atd->offset_vars) { **i = atd->offset; }	// notify aliases' listeners
	tagAliases[alias] = name;
	return atd;
}

// activation of tags that are in use
sstruct::tagdef *sstruct::factory_t::useTagRaw(const std::string &name, bool allowUndefined)
{
	if(!definedTags.count(name) && allowUndefined)
	{
		createTag(name);
	}

	if(!definedTags.count(name))
	{
		THROW(EAny, "Could not find/define tag " + name);
	}

	tagdef *td = definedTags[name];
	if(td->offset != -1) { return td; }	// if already in use
	return addTag(td);
}

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

double ToyHomogeneous_model::rho(double x, double y, double z, double ri)
{
	return rho0;
}

// double ToyHomogeneous_model::absmag(double ri)
// {
// 	return paralax.Mr(ri);
// }

ToyHomogeneous_model::ToyHomogeneous_model(peyton::system::Config &cfg)
	: galactic_model(cfg)
{
	cfg.get(rho0, "rho0", rho0);
	DLOG(verb1) << "rho0 = " << rho0;
}

peyton::io::obstream& ToyHomogeneous_model::serialize(peyton::io::obstream& out) const
{
	galactic_model::serialize(out);
	out << rho0;

	return out;
}

ToyHomogeneous_model::ToyHomogeneous_model(peyton::io::ibstream &in)
	: galactic_model(in)
{
	in >> rho0;
	ASSERT(in);
}

////////////////////////////////////////////////////

double ToyGeocentricPowerLaw_model::rho(double x, double y, double z, double ri)
{
	x -= Rg;
	double d2 = sqr(x) + sqr(y) + sqr(z);
	double norm = lf.empty() ? 1. : lf(ri);
	return norm * rho0 * pow(d2, n/2.);
}

// double ToyGeocentricPowerLaw_model::absmag(double ri)
// {
// 	return paralax.Mr(ri);
// }

ToyGeocentricPowerLaw_model::ToyGeocentricPowerLaw_model(peyton::system::Config &cfg)
	: galactic_model(cfg)
{
	cfg.get(rho0,  "rho0",  rho0);
	cfg.get(n, "n", n);

	if(cfg.count("lumfunc"))
	{
		text_input_or_die(in, cfg["lumfunc"]);
		vector<double> ri, phi;
		::load(in, ri, 0, phi, 1);
		lf.construct(ri, phi);
	}

	DLOG(verb1) << "rho0 = " << rho0 << ", n = " << n;
}

peyton::io::obstream& ToyGeocentricPowerLaw_model::serialize(peyton::io::obstream& out) const
{
	galactic_model::serialize(out);
	out << rho0 << n << lf;

	return out;
}
ToyGeocentricPowerLaw_model::ToyGeocentricPowerLaw_model(peyton::io::ibstream &in)
	: galactic_model(in)
{
	in >> rho0 >> n >> lf;
	ASSERT(in);
}

////////////////////////////////////////////////////
#if 0
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
#endif

void BahcallSoneira_model::load(peyton::system::Config &cfg)
{
	FOREACH(cfg) { MLOG(verb1) << (*i).first << " = " << (*i).second; }

	FOR(0, m.nparams - m.nrho)
	{
		std::string param = m.param_name[i];
		ASSERT(cfg.count(param)) { std::cerr << "Initial value for " << param << " not specified\n"; }

		m.p[i] = cfg[param];
	}
	
	// luminosity function
	if(cfg.count("lumfunc"))
	{
		if(cfg.count("rho0_ri") == 0)
		{
			rho0_ri = make_pair(0., 0.);
		} else {
			rho0_ri = cfg["rho0_ri"];
		}

		input_or_die(in, cfg["lumfunc"]);
		load_luminosity_function(in, rho0_ri);
	}

	// cutoff radius (default: 100kpc)
	cfg.get(r_cut2,  "rcut",   1e5);
	r_cut2 *= r_cut2;
}

// BahcallSoneira_model::BahcallSoneira_model()
// {
// }

BahcallSoneira_model::BahcallSoneira_model(peyton::system::Config &cfg)
	: galactic_model(cfg)
{
	load(cfg);
}
/*
BahcallSoneira_model::BahcallSoneira_model(const std::string &prefix)
{
	Config cfg(prefix);
	load(cfg);
}*/

// double BahcallSoneira_model::absmag(double ri)
// {
// 	return paralax.Mr(ri);
// }

sstruct::factory_t sstruct::factory;
std::map<sstruct *, char *> sstruct::owner;

bool BahcallSoneira_model::draw_tag(sstruct &t, double x, double y, double z, double ri, gsl_rng *rng)
{
	galactic_model::draw_tag(t, x, y, z, ri, rng);

	double r = sqrt(x*x + y*y);

	double thin = m.rho_thin(r, z, 0);
	double thick = m.rho_thick(r, z, 0);
	double halo = m.rho_halo(r, z, 0);
	double rho = thin+thick+halo;

	double pthin  = thin / rho;
	double pthick = (thin + thick) / rho;

	double u = gsl_rng_uniform(rng);
	if(u < pthin) { t.component() = THIN; }
	else if(u < pthick) { t.component() = THICK; }
	else { t.component() = HALO; }

//	float *f = t.XYZ(); f[0] = x; f[1] = y; f[2] = z;
//	std::cerr << r << " " << z << " : " << pthin << " " << pthick << " -> " << u << " " << t.comp << "\n";

	return true;
}

double BahcallSoneira_model::rho(double x, double y, double z, double ri)
{
	double r = sqrt(x*x + y*y);
	double norm = lf.empty() ? 1. : lf(ri);
//	norm = 1.;
	double rho = norm * m.rho(r, z, 0);

#if 1
	// Galactocentric cutoff: model it as a smooth transition, so that the integrator driver doesn't barf
	// The exponential is an analytic approximation of a step function
	double rc = (x*x + y*y + z*z) / r_cut2 - 1;

	double f;
	     if(rc < -0.01) { f = 1.; }
	else if(rc > 0.01)  { f = 0.; }
	else                { f = 1. / (1. + exp(1000 * rc)); };

	rho = rho * f;
#else
	// simple cutoff
	if(x*x + y*y + z*z > r_cut2) { return 0; }
#endif
	return rho;
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
	if(dr > 0) {
		double stars_per_mag = lf.integral(rho0_ri.first, rho0_ri.second) / dr;
		FOREACH(phi) { *i /= stars_per_mag; };
	}
	lf.construct(ri, phi);
#if 0
	std::cerr << "Norm.: " << 1./stars_per_mag << "\n";
	std::cerr << "New int: " << lf.integral(rho0_ri.first, rho0_ri.second) / dr << "\n";
	std::cerr << "lf(1.0): " << lf(1.0) << "\n";
	std::cerr << "lf(1.1): " << lf(1.1) << "\n";
#endif
#if 0
	for(double ri=0; ri < 1.5; ri += 0.005)
	{
		std::cerr << ri << " " << lf(ri) << "\n";
	}
#endif
}

BLESS_POD(disk_model);
peyton::io::obstream& BahcallSoneira_model::serialize(peyton::io::obstream& out) const
{
	galactic_model::serialize(out);
	out << m << lf << r_cut2;

	return out;
}
BahcallSoneira_model::BahcallSoneira_model(peyton::io::ibstream &in)
	: galactic_model(in)
{
	in >> m >> lf >> r_cut2;
	ASSERT(in);
}


galactic_model *galactic_model::load(istream &cfgstrm)
{
	Config cfg;
	cfg.load(cfgstrm);

	FOREACH(cfg) { DLOG(verb1) << (*i).first << " = " << (*i).second; }
	if(cfg.count("model") == 0) { ASSERT(0); return NULL; }

	std::string model = cfg["model"];

	if(model == "BahcallSoneira") { return new BahcallSoneira_model(cfg); }
	if(model == "ToyHomogeneous") { return new ToyHomogeneous_model(cfg); }
	if(model == "ToyGeocentricPowerLaw") { return new ToyGeocentricPowerLaw_model(cfg); }

	ASSERT(0); return NULL;
}

galactic_model *galactic_model::unserialize(peyton::io::ibstream &in)
{
	std::string model;
	in >> model;

	if(model == "BahcallSoneira")
	{
		std::auto_ptr<BahcallSoneira_model> m(new BahcallSoneira_model(in));
		return m.release();
	}
	if(model == "ToyHomogeneous")
	{
		std::auto_ptr<ToyHomogeneous_model> m(new ToyHomogeneous_model(in));
		return m.release();
	}
	if(model == "ToyGeocentricPowerLaw")
	{
		std::auto_ptr<ToyGeocentricPowerLaw_model> m(new ToyGeocentricPowerLaw_model(in));
		return m.release();
	}

	ASSERT(0); return NULL;
}

peyton::io::obstream& galactic_model::serialize(peyton::io::obstream& out) const
{
	out << name();
	out << m_band << m_color;
	out << paralax_loaded << paralax;
	return out;
}

galactic_model::galactic_model(peyton::io::ibstream &in)
{
	in >> m_band >> m_color;
	in >> paralax_loaded >> paralax;
	ASSERT(in);
}

galactic_model::galactic_model(peyton::system::Config &cfg)
{
	cfg.get(m_band,  "band",   "mag");
	cfg.get(m_color, "color", "color");

	// load paralax coefficients, if there are any
	if(cfg.count("col2absmag.poly"))
	{
		std::vector<double> coeff = cfg["col2absmag.poly"];
		paralax.setParalaxCoefficients(coeff);
		paralax_loaded = true;
	}
}

bool galactic_model::setup_tags(sstruct::factory_t &factory)
{
	factory.useTag("comp");
	factory.useTag("XYZ[3]");
	return true;
}

bool galactic_model::draw_tag(sstruct &t, double x, double y, double z, double ri, gsl_rng *rng)
{
	t.component() = 0;
	float *f = t.XYZ(); f[0] = x; f[1] = y; f[2] = z;

	return true;
}
