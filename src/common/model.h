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

#include <boost/array.hpp>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>

#include <astro/types.h>
#include <astro/system/config.h>
#include <astro/system/log.h>
#include <astro/io/binarystream.h>

#include "paralax.h"

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

 	double operator ()(double x)        { return gsl_interp_eval(f, &xv[0], &yv[0], x, acc); }
 	double deriv(double x)              { return gsl_interp_eval_deriv(f, &xv[0], &yv[0], x, acc); }
 	double deriv2(double x)             { return gsl_interp_eval_deriv2(f, &xv[0], &yv[0], x, acc); }
 	double integral(double a, double b) { return gsl_interp_eval_integ(f, &xv[0], &yv[0], a, b, acc); }

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


class sstruct	// "Smart struct" -- a structure with variable number (in runtime) of predefined members
{
	protected:
		char *tags;

	public:
		struct tagdef
		{
			const std::string tagName;		// unique name of the tag
			const size_t size;			// size of the tag data (in bytes)
			size_t offset;				// the offset of this tag, if active, -1 otherwise
			std::vector<size_t*> offset_vars;	// variable to update with tag offset, if/when this tag gets activated

			virtual void  serialize1(const void *, peyton::io::obstream &) const = 0;
			virtual void  serialize2(const void *, std::ostream &) const = 0;
			virtual void  unserialize1(void *, peyton::io::ibstream &) = 0;
			virtual void  unserialize2(void *, std::istream &) = 0;
			virtual void* constructor(void *) = 0;
			virtual void  destructor(void *val) = 0;
			virtual void  copy(void *dest, void *src) = 0;

			tagdef(const std::string &tid, const size_t s) : tagName(tid), size(s), offset(-1) {}
		protected:
			tagdef(const tagdef &);
			tagdef &operator =(const tagdef &);
		};
		template<typename T> struct tagdefT : public tagdef
		{
			virtual void  serialize1(const void *val, peyton::io::obstream &out) const { const T *v = reinterpret_cast<const T*>(val); out << *v; }
			virtual void  serialize2(const void *val, std::ostream &out) const { const T *v = reinterpret_cast<const T*>(val); out << *v; }
			virtual void  unserialize1(void *val, peyton::io::ibstream &in)  { T *v = reinterpret_cast<T*>(val); in >> *v; }
			virtual void  unserialize2(void *val, std::istream &in) { T *v = reinterpret_cast<T*>(val); in >> *v; }
			virtual void* constructor(void *p)  { return new (p) T(); }
			virtual void  destructor(void *val) { reinterpret_cast<T*>(val)->~T(); }
			virtual void  copy(void *dest, void *src) { *reinterpret_cast<T*>(dest) = *reinterpret_cast<T*>(src); }

			tagdefT(const std::string &tid) : tagdef(tid, sizeof(T)) {}
		};
		// simple array of types T
		template<typename T> struct tagdefTA : public tagdef
		{
			size_t n;

			virtual void  serialize1(const void *val, peyton::io::obstream &out) const { const T *v = reinterpret_cast<const T*>(val); FOR(0,n) { out << v[i]; } }
			virtual void  serialize2(const void *val, std::ostream &out) const { const T *v = reinterpret_cast<const T*>(val); FOR(0,n) { out << (i ? " " : "") << v[i]; } }
			virtual void  unserialize1(void *val, peyton::io::ibstream &in)  { T *v = reinterpret_cast<T*>(val); FOR(0,n) { in >> v[i]; } }
			virtual void  unserialize2(void *val, std::istream &in) { T *v = reinterpret_cast<T*>(val); FOR(0,n) { in >> v[i]; } }
			virtual void* constructor(void *p)  { FOR(0,n) { new (reinterpret_cast<T*>(p)+i) T(); } }
			virtual void  destructor(void *val) { FOR(0,n) { reinterpret_cast<T*>(val)[i].~T(); } }
			virtual void  copy(void *dest, void *src) { FOR(0,n) { reinterpret_cast<T*>(dest)[i] = reinterpret_cast<T*>(src)[i];} }

			tagdefTA(const std::string &tid, size_t n_) : tagdef(tid, sizeof(T)*n_), n(n_) {}
		};

		struct factory_t	// singleton used to initialize arrays
		{
			std::map<size_t,      tagdef*> usedTags;	// tags currently in use (offset -> tagdef map)
			std::map<std::string, tagdef*> definedTags;	// all defined tags that can be activated with useTag()
			std::map<std::string, std::string> tagAliases;	// tag alias -> tag name

			size_t nextOffset;			// next available offset for a tag activated with useTag()
			size_t tagSize;				// size of all active tags, once the structure has been frozen

			std::vector<tagdef *> streamTags;	// list of tags to be unserialized

			// helpers for tag definitions
			template<typename T> void defineScalarTag(const std::string &tagName, size_t *offset_var = NULL)
			{
				tagdef *td = new tagdefT<T>(tagName);
				if(offset_var) { td->offset_vars.push_back(offset_var); }

				definedTags[tagName] = td;
			}
			template<typename T> void defineArrayTag(const std::string &tagName, const size_t n, size_t *offset_var = NULL)
			{
				if(n == 0) { return defineScalarTag<T>(tagName, offset_var); }

				tagdef *td = new tagdefTA<T>(tagName, n);
				if(offset_var) { td->offset_vars.push_back(offset_var); }

				definedTags[tagName] = td;
			}

			// activation of tags that are in use
			tagdef *useTagRaw(const std::string &name, bool allowUndefined = false);
			size_t useTag(const std::string &name, bool allowUndefined = false)
			{
				return useTagRaw(name, allowUndefined)->offset;
			}
			tagdef *addTag(tagdef *td)
			{
				die_if_frozen();

				td->offset = nextOffset;
				nextOffset += td->size;

				usedTags[td->offset] = td;
				FOREACH(td->offset_vars) { **i = td->offset; }

				return td;
			}
			tagdef *aliasTag(const std::string &name, const std::string &alias)
			{
				ASSERT(definedTags.count(name));
				ASSERT(definedTags.count(alias));
				
				tagdef *td = definedTags[name];
				ASSERT(td->offset != -1);		// the tag alias target has to be in use

				tagdef *atd = definedTags[alias];
				ASSERT(atd->offset == -1);		// the alias must be unused
				ASSERT(atd->size == td->size);		// alias and the aliased tag must be of same size
				
				atd->offset = td->offset;		// let the alias know where in the index it is
				FOREACH(atd->offset_vars) { **i = atd->offset; }	// notify aliases' listeners
				tagAliases[alias] = name;
				return atd;
			}
			size_t getOffset(const std::string &name)
			{
				ASSERT(definedTags.count(name));
				tagdef *td = definedTags[name];
				return td->offset;
			}

			// in-use tags serialization/unserialization
			peyton::io::ibstream& unserialize(peyton::io::ibstream& in)
			{
				std::string tagName;
				size_t size;
				FOREACH(usedTags)
				{
					in >> tagName >> size;
					if(!definedTags.count(tagName)) { ASSERT(0) { std::cerr << "tagName = " << tagName << " not registered."; } }
					tagdef *td = definedTags[tagName];
					ASSERT(size == td->size) { std::cerr << "tagName = " << tagName << "\n"; }
					addTag(td);
				}
				return in;
			}
			peyton::io::obstream& serialize(peyton::io::obstream& out)
			{
				FOREACH(usedTags) { out << *i; }
				return out;
			}
			size_t gettags(std::set<std::string> &tags) const
			{
				FOREACH(usedTags) { tags.insert(i->second->tagName); }
				FOREACH(tagAliases) { tags.insert(i->first); }
			}
			std::ostream& serialize(std::ostream& out) const
			{
				bool first = true;
				FOREACH(usedTags)
				{
					if(!first) out << " ";
					tagdef *td = i->second;
					out << td->tagName;
					first = false;
				}
				if(tagAliases.size())
				{
					out << "  |  ";
					first = true;
					FOREACH(tagAliases)
					{
						if(!first) out << " ";
						out << i->first << "=" << i->second; // alias=name pairs
						first = false;
					}
				}
				return out;
			}
			std::istream& unserialize(std::istream& in)
			{
				std::string tagName, line;
				size_t size;

				std::getline(in, line);
				ASSERT(in);
				std::istringstream ss(line.c_str());

				// gobble up any "#" characters
				do { ss >> tagName; } while(tagName == "#" && ss);
				ASSERT(ss);

				// split tagName to tagName and size
				streamTags.clear();
				do {
					//std::cerr << "Tag: " << tagName << "\n";
					streamTags.push_back(useTagRaw(tagName, true));
				} while(ss >> tagName && tagName != "|");

				// check if there are alias definitions following the | sign
				if(!ss) { return in; }

				std::string aliasName;
				ss >> aliasName;
				ASSERT(ss);
				do {
					// parse alias=tagName pairs
					size_t idx = aliasName.find('=');
					if(idx == std::string::npos) { ASSERT(0); continue; }
					tagName = aliasName.substr(idx+1);
					aliasName = aliasName.substr(0, idx);
					aliasTag(tagName, aliasName);
				} while(ss >> aliasName);
				return in;
			}
	
			void freeze_tags()
			{
				if(tagSize != -1) { return; }
				tagSize = nextOffset;
			}
			void die_if_frozen()
			{
				if(tagSize != -1)
				{
					ASSERT(0) { std::cerr << "Tags have been frozen!\n"; }
					abort();
				}
			}

			static const size_t max_ovars = 1000;
			size_t ovars[max_ovars]; // offset variables

			// some special (shared) ovars
			static const size_t IVAR_COLOR =  max_ovars-1;	// offset variable with "color" information
			static const size_t IVAR_MAG =    max_ovars-2;	// offset variable with magnitude information
			static const size_t IVAR_ABSMAG = max_ovars-3;	// offset variable with absolute magnitude information

			static const size_t SDSS_BASE	= 20;		// offset variable from which SDSS-related stuff begins
			static const size_t DEBUG_BASE	= 100;		// offset variable from which debugging/test-related stuff begins

			factory_t()
				: nextOffset(0), tagSize(-1)
			{
				//std::cerr << "Factory: defining tags\n";
				defineScalarTag<int>("comp", &ovars[0]);
				defineScalarTag<float>("extinction.r", &ovars[1]);
				defineArrayTag<double>("lb[2]", 2, &ovars[3]);
				defineArrayTag<float>("XYZ[3]", 3, &ovars[4]);
				defineScalarTag<float>("FeH", &ovars[5]);
				defineArrayTag<float>("vcyl[3]", 3, &ovars[6]);
				defineScalarTag<int>("photoFlags", &ovars[7]);
				defineScalarTag<std::string>("star_name", &ovars[DEBUG_BASE+0]);		// test thingee

				// SDSS
				defineScalarTag<float>("absSDSSr", &ovars[SDSS_BASE+0]);		// absolute magnitude
				defineScalarTag<float>("SDSSr", &ovars[SDSS_BASE+1]);			// SDSS r band
				defineScalarTag<float>("SDSSri", &ovars[SDSS_BASE+2]);		// LF color
				defineArrayTag<float>("SDSSugriz[5]", 5, &ovars[SDSS_BASE+3]);	// SDSS ugriz colors

				// LSST
//				defineScalarTag<float>("absLSSTr", &ovars[IVAR_COLOR]);	// absolute magnitude
//				defineScalarTag<float>("LSSTr", &ovars[IVAR_MAG]);	// LSST r band

				// built-in generics
				defineScalarTag<float>("color", &ovars[IVAR_COLOR]);	// generic
				defineScalarTag<float>("mag",   &ovars[IVAR_MAG]);	// generic
				defineScalarTag<float>("absmag", &ovars[IVAR_ABSMAG]); // absolute magnitude in IVAR_MAG's band
			}
			~factory_t()
			{
				FOREACH(definedTags) { delete i->second; }
			}
		};

		// tag accessors -- WARNING: The indices here MUST match the indices in factory_t::factory_t()
		float &color()		{ return get<float>(factory.ovars[factory_t::IVAR_COLOR]); }
		float &mag()		{ return get<float>(factory.ovars[factory_t::IVAR_MAG]); }
		float &absmag()		{ return get<float>(factory.ovars[factory_t::IVAR_ABSMAG]); }

		int &component()	{ return get<int>(factory.ovars[0]); }
		float &ext_r()		{ return get<float>(factory.ovars[1]); }
//		float *vel()		{ return get<float[3]>(factory.ovars[2]); }
		std::pair<double, double> &lb()
					{ return get<std::pair<double,double> >(factory.ovars[3]); }
		float *XYZ()		{ return get<float[3]>(factory.ovars[4]); }
		float &FeH()		{ return get<float>(factory.ovars[5]); }
		float *vcyl()		{ return get<float[3]>(factory.ovars[6]); }

		int &photoFlags()	{ return get<int>(factory.ovars[7]); }
	
		std::string &starname()	{ return get<std::string>(factory.ovars[factory_t::DEBUG_BASE+0]); }

		// SDSS
		float &absSDSSr()	{ return get<float>(factory.ovars[factory_t::SDSS_BASE+0]); }
		float &SDSSr()		{ return get<float>(factory.ovars[factory_t::SDSS_BASE+1]); }
		float *SDSSugriz()	{ return get<float[5]>(factory.ovars[factory_t::SDSS_BASE+3]); }

		static factory_t factory;			// factory singleton
		static std::map<sstruct *, char* > owner;	// list of objects that own their tags pointer

		// tag lookup by offset
		template<typename T> T& get(const size_t offset)
		{
			ASSERT(factory.usedTags.count(offset));
			ASSERT(factory.usedTags[offset]->size == sizeof(T));
			return *reinterpret_cast<T*>(tags + offset);
		}
		template<typename T> T* getptr(const size_t offset)
		{
			ASSERT(factory.usedTags.count(offset));
			return reinterpret_cast<T*>(tags + offset);
		}
		// tag lookup by name (slow, should almost never be used)
		template<typename T> T& get(const std::string &name)
		{
			ASSERT( factory.definedTags.count(name) );		// ensure this tag exists
			const tagdef *td = factory.definedTags[name];
			ASSERT( td->offset != -1 && td->size == sizeof(T) );	// ensure this tag is active
			return *reinterpret_cast<T*>(tags + td->offset);
		}

	public:
		std::ostream& serialize(std::ostream& out) const
		{
			bool first = true;
			FOREACH(factory.usedTags)
			{
				if(!first) { out << " "; }
				i->second->serialize2(tags + i->second->offset, out);
				first = false;
			}
			return out;
		};
		std::istream& unserialize(std::istream& in)
		{
			FOREACH(factory.streamTags)
			{
				(*i)->unserialize2(tags + (*i)->offset, in);
			}
			return in;
		};

		peyton::io::obstream& serialize(peyton::io::obstream& out)
		{
			out << (char)factory.usedTags.size();
			FOREACH(factory.usedTags)
			{
				i->second->serialize1(tags + i->second->offset, out);
			}
			return out;
		};
		peyton::io::ibstream& unserialize(peyton::io::ibstream& in)
		{
			char count;
			in >> count;
			FOREACH(factory.usedTags)
			{
				i->second->unserialize1(tags + i->second->offset, in);
				count--;
				if(!count) { break; }
			}
			return in;
		};

	public:
		sstruct &operator=(const sstruct &s)
		{
			FOREACHj(j, factory.usedTags)
			{
				j->second->copy(tags + j->second->offset, s.tags + j->second->offset);
			}
			return *this;
		}

	public:
		static sstruct* create()
		{
			factory.freeze_tags();

			sstruct *t = new sstruct;
			// optimization when no tags are defined
			if(factory.tagSize == 0) { t->tags = NULL; return t; }

			t->tags = new char[factory.tagSize];
			//std::cerr << "Allocated " << (void*)t->tags << " as array (size=" << factory.tagSize << ").\n";
			owner[t] = t->tags;
			// tag construction
			FOREACH(factory.usedTags)
			{
				i->second->constructor(t->tags + i->second->offset);
			}
			return t;
		}
		static sstruct* create(size_t n)
		{
			factory.freeze_tags();

			sstruct *t = new sstruct[n];
			if(n == 0) { return t; }

			// optimization when no tags are defined
			if(factory.tagSize == 0) { FOR(0, n) { t[i].tags = NULL; }; return t; }

			char *tags = new char[factory.tagSize*n];
			//std::cerr << "Allocated " << (void*)tags << " as array (size=" << factory.tagSize*n << ").\n";
			for(int i=0; i != n; i++)
			{
				t[i].tags = tags + factory.tagSize*i;
				FOREACHj(j, factory.usedTags)
				{
					j->second->constructor(t[i].tags + j->second->offset);
				}
			}
			owner[t] = tags;
			//std::cerr << "Array " << (void*)tags << " bound to " << t << "\n";
			return t;
		}
		~sstruct()
		{
			// see if we're the owner of the tags memory, delete if we are
			if(tags && owner.count(this))
			{
				//std::cerr << "Calling destructors on " << (void*)tags << "\n";
				FOREACHj(j, factory.usedTags)
				{
					j->second->destructor(tags + j->second->offset);
				}
				//std::cerr << "Deleting " << (void*)tags << " as array.\n";
				delete [] tags;
				owner.erase(this);
			}
		}
	protected:
		sstruct()
		{
			//std::cerr << "In constructor for " << this << "\n";
		};
		friend class galactic_model;
};
inline OSTREAM(const sstruct &ss) { return ss.serialize(out); }
inline ISTREAM(sstruct &ss) { return ss.unserialize(in); }
inline OSTREAM(const sstruct::factory_t &ss) { return ss.serialize(out); }
inline ISTREAM(sstruct::factory_t &ss) { return ss.unserialize(in); }

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
	virtual bool draw_tag(sstruct &t, double x, double y, double z, double ri, gsl_rng *rng);	// by default, adds comp=0 and XYZ tags
	virtual bool setup_tags(sstruct::factory_t &factory);	// by default, sets up comp and XYZ[3] tags

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

	virtual bool draw_tag(sstruct &t, double x, double y, double z, double ri, gsl_rng *rng);
	virtual const std::string &name() const { static std::string s = "BahcallSoneira"; return s; }
public:
	BahcallSoneira_model();
	BahcallSoneira_model(peyton::system::Config &cfg);
	BahcallSoneira_model(peyton::io::ibstream &in);

//	virtual double absmag(double ri);
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
//	virtual double absmag(double ri);
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
//	double absmag(double ri);
	double rho(double x, double y, double z, double ri);
	virtual peyton::io::obstream& serialize(peyton::io::obstream& out) const;
};

#if 0
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
#endif
