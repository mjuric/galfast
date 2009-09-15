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
	This file contains the bits and pieces of the model fitter
	used for Juric et al. (2008) paper. Refer to these when
	rewriting the fitter to use the new codebase.
*/

#if 0
struct rzpixel
{
	double r, rphi, z, N, V;
	double rho, sigma;
	int ri_bin;
};

struct disk_model
{
	static const double Rg = 8000.;
	
	static const size_t nrho = 10;
	static const size_t nparams = 10 + nrho;
	static const char *param_name[nparams];
	static const char *param_format[nparams];

	//int disk;
/*	void set_ri_bin(int k)
	{
		disk = ri2idx(k);
	}*/
	int ri2idx(int k) const
	{
		return k == 0 ? k : (nparams - nrho) + (k-1);
	}

	union {
		double p[nparams];
		double rho0[nparams];
		struct
		{
			double rho0x, l, h, z0, f, lt, ht, fh, q, n;
			double rho1, rho2, rho3, rho4, rho5, rho6, rho7, rho8, rho9, rho10;
		};
	};

	disk_model() {}

// 	// Model functions
	double rho_thin(double r, double z, int ri)  const { return rho0[ri2idx(ri)] *     exp((Rg-r)/l  + (std::abs(z0) - std::abs(z + z0))/h); }
	double rho_thick(double r, double z, int ri) const { return rho0[ri2idx(ri)] * f * exp((Rg-r)/lt + (std::abs(z0) - std::abs(z + z0))/ht); }
	double rho_halo(double r, double z, int ri)  const { return rho0[ri2idx(ri)] * fh * pow(Rg/sqrt(halo_denom(r,z)),n); }
	double rho(double r, double z, int ri)       const { return rho_thin(r, z, ri) + rho_thick(r, z, ri) + rho_halo(r, z, ri); }

	//double norm_at_Rg() const { return f*exp(Rg*(1./l - 1./lt)); }
	double norm_at_Rg(int ri) const { return rho_thick(Rg, 0, ri)/rho_thin(Rg, 0, ri); }

	// Derivatives of the model function
	double drho0(double r, double z, double rhom, int ri, int rij) const {
		double tmp = ri == rij ? 1./rho0[ri2idx(ri)] * rhom : 0.;
		//std::cerr << ri_bin << " " << tmp << "\n";
		return tmp;
	}
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
	std::vector<rzpixel> orig;		/// all input data
	std::vector<rzpixel> map;		/// data used in last fit
	std::vector<rzpixel> culled;		/// culled data (culled = orig - map)
public:
	void setdata(const std::vector<rzpixel> &data) { orig = map = data; }
	std::vector<std::pair<float, float> > ri, r;
	std::vector<std::pair<double, double> > d;
	int ndata() { return map.size(); }

	model_fitter(const model_fitter& m)
		: disk_model(m),
		  covar(m.covar), fixed(m.fixed), chi2_per_dof(m.chi2_per_dof),
		  param_name_to_index(m.param_name_to_index),
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
	int fit(int cullIter, const std::vector<double> &nsigma);
	void cull(double nSigma);
	void residual_distribution(std::map<int, int> &hist, double binwidth);

	enum {PRETTY, HEADING, LINE};
	void print(std::ostream &out, int format = PRETTY, int ri_bin = 0);
};

#endif

#if 0
class galactic_model
{
protected:
	std::string m_band;	// apparent/absolute magnitude band (e.g., "sdss_r")
	std::string m_color;	// the name of the color in ri -- note: the "color" here can be the absolute magnitude

	plx_gri_locus_ng paralax;	// polynomial converting between the "color" and the absolute magnitude
	bool paralax_loaded;		// whether polynomial coefficients were loaded
public:
	const std::string &band() const { return m_band; }
	const std::string &color() const { return m_color; }

public:
	virtual bool add_details(otable &out, rng_t &rng) {};	// does nothing by default

public:
	virtual double absmag(double ri) { ASSERT(paralax_loaded); return paralax.Mr(ri); }
	virtual double rho(double x, double y, double z, double ri) = 0;

	virtual peyton::io::obstream& serialize(peyton::io::obstream& out) const;	// needed for serialization
	virtual const std::string &name() const = 0;

public:
	static galactic_model *load(std::istream &cfg);
	static galactic_model *unserialize(peyton::io::ibstream &in);

protected:
	galactic_model() : paralax_loaded(false) {};
	galactic_model(peyton::system::Config &cfg);				// config file load constructor
	galactic_model(peyton::io::ibstream &in);				// unserialization constructor
};

class BahcallSoneira_model : public galactic_model
{
public:
	disk_model m;
	std::pair<double, double> rho0_ri;	/// interval accross which m.rho0 was calculated

	spline lf;		/// dimensionless local luminosity function

	double r_cut2;		/// Galactocentric radius squared of density cutoff (rho beyond r_cut is 0)
public:
	static const int THIN = 0, THICK = 1, HALO = 2;

	virtual const std::string &name() const { static std::string s = "BahcallSoneira"; return s; }
	virtual bool add_details(otable &out, rng_t &rng);	// by default, adds comp=0 and XYZ tags
public:
	BahcallSoneira_model();
	BahcallSoneira_model(peyton::system::Config &cfg);
	BahcallSoneira_model(peyton::io::ibstream &in);

	virtual double rho(double x, double y, double z, double ri);
	virtual peyton::io::obstream& serialize(peyton::io::obstream& out) const;
protected:
	void load(peyton::system::Config &cfg);
	void load_luminosity_function(std::istream &in, std::pair<double, double> rho0_ri);
};

class ToyHomogeneous_model : public galactic_model
{
public:
	double rho0;
public:
	ToyHomogeneous_model(double rho0_ = 1.) : rho0(rho0_) {}
	ToyHomogeneous_model(peyton::system::Config &cfg);
	ToyHomogeneous_model(peyton::io::ibstream &in);

	virtual const std::string &name() const { static std::string s = "ToyHomogeneous"; return s; }
public:
	virtual double rho(double x, double y, double z, double ri);
	virtual peyton::io::obstream& serialize(peyton::io::obstream& out) const;
};

// geocentric powerlaw model with a constant paralax relation
class ToyGeocentricPowerLaw_model : public galactic_model
{
public:
	double rho0, n;
	spline lf;		/// local luminosity function (if given)
public:
	ToyGeocentricPowerLaw_model(double rho0_ = 1., double n_ = -3.) : rho0(rho0_), n(n_) {}
	ToyGeocentricPowerLaw_model(peyton::system::Config &cfg);
	ToyGeocentricPowerLaw_model(peyton::io::ibstream &in);

	virtual const std::string &name() const { static std::string s = "ToyGeocentricPowerLaw"; return s; }
public:
	double rho(double x, double y, double z, double ri);
	virtual peyton::io::obstream& serialize(peyton::io::obstream& out) const;
};
#endif

#if 0
////////////////////////////////////////////////////

double ToyHomogeneous_model::rho(double x, double y, double z, double ri)
{
	return rho0;
}

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
	x -= Rg();
	double d2 = sqr(x) + sqr(y) + sqr(z);
	double norm = lf.empty() ? 1. : lf(ri);
	return norm * rho0 * pow(d2, n/2.);
}

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

BahcallSoneira_model::BahcallSoneira_model(peyton::system::Config &cfg)
	: galactic_model(cfg)
{
	load(cfg);
}

bool BahcallSoneira_model::add_details(otable &t, rng_t &rng)
{
	cfloat_t::host_t XYZ = t.col<float>("XYZ");
	cint_t::host_t  comp = t.col<int>("comp");

	for(size_t row = 0; row != t.size(); row++)
	{
		float x = XYZ(row, 0);
		float y = XYZ(row, 1);
		float z = XYZ(row, 2);
		float r = sqrt(x*x + y*y);

		float thin = m.rho_thin(r, z, 0);
		float thick = m.rho_thick(r, z, 0);
		float halo = m.rho_halo(r, z, 0);
		float rho = thin+thick+halo;

		float pthin  = thin / rho;
		float pthick = (thin + thick) / rho;

		float u = rng.uniform();
		if(u < pthin) { comp(row) = THIN; }
		else if(u < pthick) { comp(row) = THICK; }
		else { comp(row) = HALO; }
	}

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
#endif

/////////////

#if 0
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

	}

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

#endif
