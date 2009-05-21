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
#include <sstream>
#include <valarray>
#include <list>

#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <astro/types.h>
#include <astro/system/config.h>
#include <astro/system/log.h>
#include <astro/io/binarystream.h>
#include <astro/io/format.h>
#include <astro/exceptions.h>

#include "paralax.h"
#include "binarystream.h"
#include "column.h"

const std::string &datadir(); // return the path to built-in datafiles (TODO: move it to someplace where it belongs)

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

class spline
{
public:
	gsl_interp *f;
	gsl_interp_accel *acc;
	std::valarray<double> xv, yv;

	friend BOSTREAM2(const spline &spl);
	friend BISTREAM2(spline &spl);
protected:
	void construct_aux();
public:
	spline() : f(NULL), acc(NULL) {}
	spline(const double *x, const double *y, int n);
	void construct(const double *x, const double *y, int n);
	void construct(const std::valarray<double> &x, const std::valarray<double> &y)
		{ ASSERT(x.size() == y.size()); construct(&x[0], &y[0], x.size()); }
	void construct(const std::vector<double> &x, const std::vector<double> &y)
		{ ASSERT(x.size() == y.size()); construct(&x[0], &y[0], x.size()); }
	~spline();

 	double operator ()(double x)        const { return gsl_interp_eval(f, &xv[0], &yv[0], x, acc); }
	double deriv(double x)              const { return gsl_interp_eval_deriv(f, &xv[0], &yv[0], x, acc); }
	double deriv2(double x)             const { return gsl_interp_eval_deriv2(f, &xv[0], &yv[0], x, acc); }
	double integral(double a, double b) const { return gsl_interp_eval_integ(f, &xv[0], &yv[0], a, b, acc); }

	bool empty() const { return xv.size() == 0; }
public:
	spline& operator= (const spline& a);
	spline(const spline& a) : f(NULL), acc(NULL) { *this = a; }
};
BOSTREAM2(const spline &spl);
BISTREAM2(spline &spl);

inline OSTREAM(const float x[3]) { return out << x[0] << " " << x[1] << " " << x[2]; }
inline ISTREAM(float x[3]) { return in >> x[0] >> x[1] >> x[2]; }

template<typename T, std::size_t N>
		inline std::ostream &operator <<(std::ostream &out, const boost::array<T, N> &x)
		{
			out << x[0];
			FOR(1, N) { out << " " << x[i]; }
			return out;
		}

template<typename T, std::size_t N>
		inline std::istream &operator >>(std::istream &in, boost::array<T, N> &x)
		{
			FOR(0, N) { in >> x[i]; }
			return in;
		}

template<typename T, std::size_t N>
		inline BOSTREAM2(const boost::array<T, N> &x)
		{
			FOR(0, N) { out << x[i]; }
			return out;
		}

template<typename T, std::size_t N>
		inline BISTREAM2(boost::array<T, N> &x)
		{
			FOR(1, N) { in >> x[i]; }
			return in;
		}

struct rng_gsl_t : public rng_t
{
	bool own;
	gsl_rng *rng;

	rng_gsl_t(gsl_rng *r, bool own_ = false) : rng(r), own(own_) {}
	rng_gsl_t(unsigned long int seed)
		: own(true), rng(NULL)
	{
		rng = gsl_rng_alloc(gsl_rng_default);
		gsl_rng_set(rng, seed);
	}

	virtual ~rng_gsl_t() { if(own && rng) gsl_rng_free(rng); }

	virtual float uniform()
	{
		return (float)gsl_rng_uniform(rng);
	}
	virtual float gaussian(const float sigma) { return (float)gsl_ran_gaussian(rng, sigma); }
};

class fmtout
{
protected:
	static const size_t BUFMAX = 20000;
	char buf[BUFMAX+1];	// line buffer
	size_t pos;
public:
	fmtout() : pos(0) {}
	const char *c_str() const { return buf; }

	size_t prep_buf()
	{
		if(pos == BUFMAX)
		{
			// This should really never happen ...
			buf[BUFMAX] = 0;
			THROW(peyton::exceptions::EAny, "Line buffer exhausted");
		}

		// Spaces between fields
		if(pos != 0) { buf[pos] = ' '; pos++; }
		return BUFMAX - pos;
	}

	template<typename T>
	int printf_aux(char *dest, size_t len, const char *fmt, const T &v)
	{
		// Default action: explode, because this has to be overloaded for
		// every printf-legal type
		THROW(peyton::exceptions::EAny, "Internal error");
	}

	int printf_aux(char *dest, size_t maxlen, const char *fmt, const double &v) 	{ pos += snprintf(dest, maxlen, fmt, v); }
	int printf_aux(char *dest, size_t maxlen, const char *fmt, const float &v) 	{ pos += snprintf(dest, maxlen, fmt, v); }
	int printf_aux(char *dest, size_t maxlen, const char *fmt, const int &v) 	{ pos += snprintf(dest, maxlen, fmt, v); }

	template<typename T>
	int printf(const std::string &fmt, const T &v)
	{
		size_t len = prep_buf();

		if(!fmt.size())
		{
			// No format specification -- revert to iostreams
			std::ostringstream ss;
			ss << v;
			std::string out = ss.str();
			strncpy(buf+pos, out.c_str(), std::min(len, out.size()));
			pos += out.size();
		}
		else
		{
			// sprintf format specified, use it
			printf_aux(buf+pos, len, fmt.c_str(), v);
		}
	}
};

struct column_type_traits
{
	const std::string typeName;
	const size_t elementSize;

	virtual void  serialize(fmtout &out, const std::string &format, const void *val) const = 0;
	virtual void  unserialize(void *val, std::istream &in) const = 0;
	virtual void* constructor(void *p) const = 0;
	virtual void  destructor(void *val) const = 0;

	static const column_type_traits *get(const std::string &datatype);
	template<typename T> static const column_type_traits *get() { ASSERT(0); }
protected:
	static std::map<std::string, boost::shared_ptr<column_type_traits> > defined_types;
	column_type_traits(const std::string &name, const size_t size) : typeName(name), elementSize(size) {}
};
// These are C type->traits mappings. Specialize them for each datatype supported by column_type_traits::get
template<> const column_type_traits *column_type_traits::get<float>();
template<> const column_type_traits *column_type_traits::get<int>();
template<> const column_type_traits *column_type_traits::get<double>();
template<> const column_type_traits *column_type_traits::get<char>();

class otable
{
protected:
	class kv
	{
	public:
		std::string what;
	protected:
		friend class otable;

		virtual void set_property(const std::string &key, const std::string &value) = 0;
		virtual void serialize_def(std::ostream &out) const = 0;

		kv(const std::string &what_) : what(what_) {}
	};

	class columnclass : public kv
	{
	protected:
		friend class otable;

		otable &parent;
		std::string className;			// primary name of this class (e.g., "photometry", "color", "astrometry", ...)

		std::string formatString;		// default format of fields of this class
		const column_type_traits *typeProxy;	// default datatype of field of this class

		columnclass(otable &parent_);

		virtual void set_property(const std::string &key, const std::string &value);
		virtual void serialize_def(std::ostream &out) const
		{
			out << className << "{";

			// aliases
			FOREACH(parent.cclasses)
			{
				if(i->second.get() != this) { continue; }
				if(i->first == className) { continue; }
				out << "alias=" << i->first << ";";
			}

			// keywords
			if(!formatString.empty()) { out << "fmt=" << formatString << ";"; }

			out << "}";
		}
	};

protected:
	struct columndef : public kv
	{
	protected:
		friend class otable;
		friend struct save_column_default;

		std::string columnName;				// primary name of the column

		const columnclass *columnClass;			// class of this column (note: it's never NULL)
		const column_type_traits *typeProxy;		// a proxy for type's serialization/construction/properties (note: must be accessed through type())
		std::string formatString;			// io::formatter format string of the column (note: must be accessed through getFormatString())

		otable &parent;					// parent table of this column

	protected:
		boost::shared_ptr<columndef> clone(const std::string &newColumnName) const
		{
			boost::shared_ptr<columndef> c(new columndef(parent));
			c->columnName = newColumnName;

// 			c->ptr.base = ptr.base;
// 			c->ptr.base = 0;
			c->ptr.reshape(ptr);

			c->columnClass = columnClass;
			c->formatString = formatString;
			c->typeProxy = typeProxy;

			return c;
		}

		const std::string &getFormatString() const
		{
			if(!formatString.empty()) { return formatString; }
			return columnClass->formatString;
		}

		const column_type_traits *type() const
		{
			return typeProxy ? typeProxy : columnClass->typeProxy;
		}

		/*
			Note: The data is stored in 'structure of arrays format'. E.g., if
			this column is a vector of double[3], and the length of the table is 4,
			the memory layout is:
				11 21 31 41 xx
				12 22 32 42 xx
				13 23 33 43 xx
			where xx is some possible padding to ensure proper word-alignment
			of adjacent rows.

			The memory location of element j in row i, use:
				Aij = data + pitch*i + elementSize*j
		*/
/*		void *data;			// the actual data
		size_t width;			// number of data elements (1 scalar, >1 if array)
		size_t length;			// length of the column (the number of rows)
		size_t pitch;			// the actuall row-length in bytes (may include some padding for proper memory alignment)*/
 		column<char>	ptr;

		friend struct cmp_in;
		friend struct cmp_out;

	protected:
		columndef(otable &parent_);

		void alloc(const size_t len);
		void dealloc();

		virtual void set_property(const std::string &key, const std::string &value);
		template<typename T> column<T> &dataptr()
		{
			ASSERT(column_type_traits::get<T>() == type())
			{
				std::cerr << "Attempting to access a " << type()->typeName << " column as " << column_type_traits::get<T>()->typeName << "\n";
			}
			//xptr<T> ptr(sizeof(T), length, width, pitch, (T*)data);
			//return column<T>(ptr);
			return (column<T> &)ptr;
			//return column<T>(data, pitch);
		}
	public:
		~columndef();
		columndef &add_alias(const std::string &name)		// add an additional name (alias) to this column
		{
			set_property("alias", name);
		}
		void serialize(fmtout &line, const size_t row) const;	// write out the element at row row
		void unserialize(std::istream &in, const size_t row);	// read in an element into row row
		virtual void serialize_def(std::ostream &out) const;	// write out the definition of this column
		size_t capacity() const { return ptr.nrows(); }
	};
	friend struct cmp_in;
	friend struct cmp_out;
	friend struct record_loaded_columns;
	friend struct save_column_default;

	std::map<std::string, boost::shared_ptr<columnclass> > cclasses;
	std::map<std::string, boost::shared_ptr<columndef> > columns;
	size_t length;	// maximum number of rows in the table
	size_t nrows;	// rows actually in the table
	std::vector<std::string> colInput, colOutput;

public:
	size_t size() const { return nrows; }
	size_t capacity() const { return length; }
	void clear() { nrows = 0; }
	size_t add_row()
	{
		size_t tmp = nrows++;
		if(nrows > capacity()) { THROW(peyton::exceptions::EAny, "Maximum number of rows in otable reached"); }
		return tmp;
	}

public:
	struct parse_callback { virtual bool operator()(kv *kvobj) = 0; };
	
protected:
	void init();
	kv *parse(const std::string &defs, parse_callback *cback = NULL);

	columndef &getColumn(const std::string &name);
	void getColumnsForOutput(std::vector<const columndef*> &cols) const;	// aux helper
	void getColumnsForInput(std::vector<columndef*> &cols);			// aux helper

public:
	// column lookup by name
	template<typename T>
	column<T> &col(const std::string &name)
	{
		return getColumn(name).dataptr<T>();
	}

	columndef &use_column(const std::string &coldef, bool setOutput = true);
	size_t get_used_columns(std::set<std::string> &cols) const;		// returns the list of columns in use
	void alias_column(const std::string &column, const std::string &alias)
	{
		getColumn(column).add_alias(alias);
	}
//	bool del_column(const std::string &name);

	otable(const size_t len)
	{
		length = len; nrows = 0;
		init();
	}

	// serialization/unserialization routines
	std::ostream& serialize_header(std::ostream &out) const;
	std::istream& unserialize_header(std::istream &in, std::set<std::string> *columns = NULL);
	std::ostream& serialize_body(std::ostream& out, size_t from = 0, size_t to = -1) const;
	std::istream& unserialize_body(std::istream& in);
	size_t set_output(const std::string &colname, bool output);
	size_t set_output_all(bool output = true);
};

class osink;
class opipeline_stage
{
	protected:
		std::set<std::string> prov, req;
		std::string uniqueId;
		stopwatch swatch;			// times how long it takes to process() this stage

		osink *nextlink;
	public:
		void chain(osink *nl) { nextlink = nl; }
		virtual size_t run(otable &t, rng_t &rng) = 0;
		float getProcessingTime() { return swatch.getTime(); }
		void setUniqueId(const std::string &uid) { uniqueId = uid; }
		const std::string &getUniqueId() const { return uniqueId; }

	public:
		static boost::shared_ptr<opipeline_stage> create(const std::string &name);
		virtual const std::string &name() const = 0;
		virtual const std::string &type() const { static std::string s("stage"); return s; }

	public:
		virtual bool init(const peyton::system::Config &cfg, otable &t) = 0;
		virtual bool prerun(const std::list<opipeline_stage *> &pipeline, otable &t);
		const std::set<std::string> &requires() const { return req; }

		bool satisfied_with(const std::set<std::string> &haves);

		static const int PRIORITY_INPUT      = -10000;
		static const int PRIORITY_STAR       =      0;
		static const int PRIORITY_INSTRUMENT =   1000;
		static const int PRIORITY_OUTPUT     =  10000;
		virtual int priority() { return PRIORITY_STAR; }

	public:
		bool inits(const std::string &cfgstring, otable &t) { return inits(cfgstring.c_str(), t); }
		bool inits(const char *cfgstring, otable &t)
		{
			std::istringstream ss(cfgstring);
			peyton::system::Config cfg;
			cfg.load(ss);
			return init(cfg, t);
		}

		opipeline_stage() : nextlink(NULL)
		{
		}
};

class osink : public opipeline_stage
{
	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng) = 0;

	public:
		virtual size_t run(otable &t, rng_t &rng) { THROW(peyton::exceptions::ENotImplemented, "We should have never gotten here"); } // we should never reach this place

	public:
		osink() : opipeline_stage()
		{
//			req.insert("_source");
		}
};


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
