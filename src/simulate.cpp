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

#if 0
 
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <locale>
#include <astro/util.h>

#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

inline bool isspace(const char c) { return c == ' ' || c == '\t'; }

#include <cstdlib>
template<typename F>
	bool parsenum(F& var, const char *s0, int& at)
	{
	#if 0
		char *endptr;
		F ret = std::strtod(s + at, &endptr);
		at = endptr - s;
		return ret;
	#else
		const char *s = s0 + at;

		while(isspace(*s)) { ++s; };
		if(*s == 0) { return false; }

		F num = 0.;
		bool neg = *s == '-';
		if((neg || *s == '+') && s[1] == 0) { return false; }

		// before decimal point
		char c;
		while(isdigit(c = *(++s)))
		{
			num *= 10; num += c - '0';
		}

		if(c == '.')
		{
			// after decimal point
			F div = 10.;
			while(isdigit(c = *(++s)))
			{
				num += (c - '0') / div;
				div *= 10.;
			}
		}
		
		if(c == 'e' || c == 'E')
		{
			// sign
			bool eneg;
			TODO!!

			// exponent
			int e = 0;
			while(isdigit(c = *(++s)))
			{
				e *= 10; e += c - '0';
			}
			num *= pow(10., e);
		}

		if(neg) { num *= -1.; }

//		cout << "\t=>" << num << "\n";
		var = num;
		at = s - s0;
//		std::cout << "at=" << at << "\n";
		return true;
	#endif
	}

template<typename F>
	std::istream& parsenum(F &var, std::istream& in)
	{
		using namespace std;

		in >> ws;
		if(!in) { return in; }

		F num = 0.; char c;
		in.get(c);
		bool neg = c == '-';
		if(!neg) { in.putback(c); }

		// before decimal point
		while(in.get(c) && isdigit(c))
		{
			num *= 10; num += c - '0';
		}

		if(c == '.' && in)
		{
			// after decimal point
			F div = 10.;
			while(in.get(c) && isdigit(c))
			{
				num += double(c - '0') / div;
				div *= 10.;
			}
		}
		
		if((c == 'e' || c == 'E') && in)
		{
			// exponent
			int e = 0;
			while(in.get(c) && isdigit(c))
			{
				e *= 10; e += c - '0';
			}
			num *= pow(10., e);
		}

		if(neg) { num *= -1.; }

//		cout << "\t=>" << num << "\n";		
		var = num;
		return in;
	}

#include <memory>

class itstream
{
protected:
	static const int maxline = 1000;
	char buffer[maxline];

	int ptr;
	std::istringstream helper_stream;
	bool helper_stream_needs_update, helper_stream_used;

	std::istream &stream()
	{
		//std::cout << "using stream()\n";
		if(helper_stream_needs_update)
		{
			helper_stream.str(buffer);
			helper_stream_needs_update = false;
		}
		else //if((int)helper_stream.tellg() != ptr)
		{
			//std::cout << ptr << "\n";
			helper_stream.seekg(ptr);
		}

		helper_stream_used = true;
		return helper_stream;
	}
protected:
	std::istream &in;

	typedef bool (itstream::*parser_function)(void *varptr);
	struct field
	{
		parser_function parser;
		void *varptr;
		int column;
	};

 	typedef std::map<int, field> field_map;
//	typedef std::vector<field> field_map;
 	field_map fields;

public:
	int nread, at;
	bool is_comment, return_comments;
	std::string comment;
	
public:
	itstream(std::istream &f_) : in(f_), nread(0), return_comments(false), helper_stream_needs_update(true) {}

	template<typename T>
		bool parse(void *varptr)
		{
			T& var = *static_cast<T*>(varptr);
			return stream() >> var;
		}

	template<typename T>
		void bind(T &var, int position);

	bool returncomments(bool rc) { return_comments = rc; }

	itstream &skip_lines(int n = 1)
	{
		std::string line;
		while(n > 0)
		{
			getline(in, line);
			n--;

			if(in.eof()) return *this;
			if(return_comments && line[0] == '#') { n++; }
		}
	}

	bool ignore_fields(int n);
	itstream &next();

	operator bool() { return nread == fields.size() || is_comment; }
	bool iscomment() { return is_comment; }
};

template<>
	bool itstream::parse<std::string>(void *varptr)
	{
		std::istream &in = stream();

		// read in a string, optionally enclosed in single or double quotes
		char q; bool quoted = false;
		in >> q;
		if(q != '\'' && q != '"') { in.putback(q); }
		else { quoted = true; }
		
		std::string &var = *static_cast<std::string*>(varptr);
		if(!quoted) { return in >> var; }

		// reading of quoted strings
		char buf[itstream::maxline];
		std::string bstr;
		while(in)
		{
			in.getline(buf, itstream::maxline, q);
			bstr.append(buf, in.gcount());
			
			if(bstr.find('\n', bstr.size() - in.gcount()) != std::string::npos)
			{
				// runaway quoted string
				in.setstate(std::ios::failbit);
				return false;
			}
			
			if(bstr.size() && *bstr.rbegin() == '\\')
			{
				// this quotation mark has been escaped
				*bstr.rbegin() = q;
				continue;
			}

			break;
		}

		if(in) { var = bstr; }
		return true;
	}

template<>
	bool itstream::parse<double>(void *varptr)
	{
		return parsenum(*static_cast<double *>(varptr), buffer, ptr);
	}

template<>
	bool itstream::parse<float>(void *varptr)
	{
		return parsenum(*static_cast<float *>(varptr), buffer, ptr);
	}

bool itstream::ignore_fields(int n)
{
#if 0
	std::string dummy;
	FOR(0, n) { parse<std::string>(static_cast<void*>(&dummy), in); }
	at += n;
	return in;
#else
	const char *s = buffer + ptr;
	char q;
	while(n)
	{
		s--;

		while(isspace(*(++s)));
		if((q = *s) == 0) break;

		if(q != '\'' && q != '"')
		{
			while(!isspace(*(++s)));
		}
		else
		{
			while(*(++s) != q && *s);
			if(*s != q) break;
		}
		
		--n;
	}
	at += n;
	ptr = s - buffer;
	return n == 0;
#endif
}

itstream &itstream::next()
{
	// check for comments
	while(in && in.peek() == '#')
	{
		if(return_comments)
		{
			getline(in, comment);
			is_comment = true;
		}
		else
		{
			in.ignore(1000000, '\n');
		}
	}

	// advance to the next line and read in all of the fields on the current line
	is_comment = false;
	nread = 0; at = 1; ptr = 0;
	helper_stream_needs_update = true;
	
	in.getline(buffer, maxline);
	if(in)
	{
		FOREACH(field_map::iterator, fields)
		{
			int id = (*i).first;
			field &f = (*i).second;
			helper_stream_used = false;

			if(id != at && !ignore_fields(id - at)) break;
			if(!CALL_MEMBER_FN(*this,f.parser)(f.varptr)) break;
			if(helper_stream_used) { ptr = helper_stream.tellg(); }

			++at;
			++nread;
		}
	}

	return *this;
}

template<typename T>
	void itstream::bind(T &var, int column)
	{
		field f =
			{
				&itstream::parse<T>,
				static_cast<void*>(&var),
				column
			};
		fields[column] = f;
/*		fields.push_back(f);*/
	}

#include <fstream>
#include "textstream.h"
#include <sys/time.h>

double seconds()
{
	timeval tv;
	int ret = gettimeofday(&tv, NULL);
	assert(ret == 0);
	return double(tv.tv_sec) + double(tv.tv_usec)/1e6;
}

using namespace std;

double parse_itext(const std::string &buf)
{
	istringstream f(buf);
	double begin = seconds();

	float a; double b;
	double b3, b4, b5, b6, b7, b8, b9;

	itstream tf(f);

	tf.bind(a, 1);
	tf.bind(b, 2);
	tf.bind(b3, 3);
	tf.bind(b4, 4);
	tf.bind(b5, 5);
	tf.bind(b6, 6);
	tf.bind(b7, 7);
	tf.bind(b8, 8);
	tf.bind(b9, 9);

	while(tf.next())
	{
//		cout << a << " " << b << " " << b9 << "\n";
		static int i = 0; if(i == 5) { cout << a << " " << b << " " << b9 << "\n"; }; i++;
	}

	return seconds() - begin;
}

double parse_stdio(const std::string &buf)
{
	istringstream f(buf);
	double begin = seconds();

	float a; double b;
	double b3, b4, b5, b6, b7, b8, b9;

	char bufx[1000];
	while(f.getline(bufx, 1000))
	{
		//sscanf(bufx, "%f %*f %*f %*f %*f %*f %*f %f", &a, &b);
		//sscanf(bufx, "%f %f", &a, &b);
		sscanf(bufx, "%f %lf %lf %lf %lf %lf %lf %lf %lf", 
			&a, &b, &b3, &b4, &b5, &b6, &b7, &b8, &b9);
		//static int i = 0; if(i == 5) { cout << a << " " << b << " " << b9 << "\n"; }; i++;
	}
	return seconds() - begin;
}

double parse_iostr(const std::string &buf)
{
	istringstream f(buf);
	double begin = seconds();

	float a; double b;
	double b3, b4, b5, b6, b7, b8, b9;

	while(!f.eof())
	{
		f >> a >> b >> b3 >> b4 >> b5 >> b6 >> b7 >> b8 >> b9;
		f.ignore(100000, '\n');
		//static int i = 0; if(i == 5) { cout << a << " " << b << " " << b9 << "\n"; }; i++;
	}
	return seconds() - begin;
}

#if 1
double parse_hand(const std::string &buf)
{
	istringstream f(buf);
	double begin = seconds();

	float a; double b;
	double b3, b4, b5, b6, b7, b8, b9;

	char bufx[1000];
	int n = 0;
	while(f.getline(bufx, 1000))
	{
		++n;
		int at = 0;
		parsenum<typeof(a)>(a, bufx, at);
		parsenum<typeof(b)>(b, bufx, at);
		parsenum<typeof(b3)>(b3, bufx, at);
		parsenum<typeof(b4)>(b4, bufx, at);
		parsenum<typeof(b5)>(b5, bufx, at);
		parsenum<typeof(b6)>(b6, bufx, at);
		parsenum<typeof(b7)>(b7, bufx, at);
		parsenum<typeof(b8)>(b8, bufx, at);
		parsenum<typeof(b9)>(b9, bufx, at);

		//static int i = 0; if(i == 5) { cout << a << " " << b << " " << b9 << "\n"; }; i++;
	}
	//cerr << n << " lines read\n";
	return seconds() - begin;
}
#else
double parse_hand(const std::string &buf)
{
	istringstream f(buf);
	double begin = seconds();

	float a; double b;
	double b3, b4, b5, b6, b7, b8, b9;

	char bufx[1000];
	int n = 0;
	while(f)
	{
		++n;
		int at = 0;
		parsenum<typeof(a)>(a, f);
		parsenum<typeof(b)>(b, f);
		parsenum<typeof(b3)>(b3, f);
		parsenum<typeof(b4)>(b4, f);
		parsenum<typeof(b5)>(b5, f);
		parsenum<typeof(b6)>(b6, f);
		parsenum<typeof(b7)>(b7, f);
		parsenum<typeof(b8)>(b8, f);
		parsenum<typeof(b9)>(b9, f);
		f.ignore(100000, '\n');

		//static int i = 0; if(i == 5) { cout << a << " " << b << " " << b9 << "\n"; }; i++;
	}
	cerr << n << " lines read\n";
	return seconds() - begin;
}
#endif

int main()
{
	double begin = seconds();

	ifstream af("test23.txt");
	stringbuf ss(ios::out);
	af >> &ss;
	string buf = ss.str();
	cerr << "loading done in " << seconds() - begin << " seconds\n";

	double sitext = 0., sstdio = 0., siostr = 0., shand = 0.;
	int N = 5;
	FOR(0, N)
	{
		double titext = parse_itext(buf); sitext += titext;
		double tstdio = parse_stdio(buf); sstdio += tstdio;
		double tiostr = parse_iostr(buf); siostr += tiostr;
		double thand  = parse_hand (buf); shand  += thand;
		cout << titext << " " << tstdio << " " << tiostr << " " << thand << "\n";
		cout.flush();
	}
	sitext /= N; sstdio /= N; siostr /= N; shand /= N;
	cout << "\n" << sitext << " " << sstdio << " " << siostr << " " << shand << "\n";
	cout << sitext/sstdio << " " << sstdio/sstdio << " " << siostr/sstdio << " " << shand /sstdio << "\n";

	return 0;
}

#else

#include "raytrace.h"
#include "interval_arithmetic.h"
#include "dm.h"
#include "projections.h"
#include "paralax.h"
#include "textstream.h"

#include <astro/math/vector.h>
#include <astro/coordinates.h>
#include <astro/math.h>
#include <astro/util.h>
#include <astro/sdss/rungeometry.h>
#include <astro/io/fits.h>
#include <astro/system/options.h>

#include <astro/useall.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>

#include "ximage.h"
#include "floodfill.h"
#include "analysis.h"
#include "binarystream.h"

#include <nag.h>
#include <nage02.h>

using namespace std;

#if 0
#if 0
struct mobject
{
	int obs_offset;			// offset in starmags array to where the obsv. of this object begin
	int n;				// total number of observations available

	float Ar;

	float mag[5];			// corrected averaged magnitudes
	float magErr[5];		// note: magnitude ordering is grizu
	float N[5];			// number of observations of this object, which were used to calculate mag[]
	short flags;			// how exactly we calculated mag[] and magErr[]

	float ml_mag[3];		// magnitudes deduced by max.likelihood fitting to locus
	float D;			// geocentric distance

	double ra, dec;
#endif

struct galaxy_model
{
	double den0, den1, l0, l1, h0, h1, z0;	// thin and thick disk parameters
	double denhalo, r0, alpha;	// halo parameters
	float ri0, ri1;

	double dx;			// scale for this model

	double density(const V3 &p);
};

double galaxy_model::density(const V3 &p)
{
	double rho = p.rho(),
	       r = abs(p),
	       z = p.z;

	double den = den0*exp(rho/l0 + (z-z0)/h0) +
	             den1*exp(rho/l1 + (z-z0)/h1) +
	             denhalo*pow(r/r0, alpha);
	return den;
}

int galaxy_model::number(double den, double V)
{
	double n = den*V;
	int N = (int)n;

	// add Poisson noise
	// TODO: draw a number from N(mu=n, sigma=sqrt(n)) and add it to n

	// the fraction is the probability of having a particle there
	n -= N;
	float prob = float(rand())/RAND_MAX;
	if(prob < n)
	{
		N++;
	}

	return N;
}

int galaxy_model::star()
{
	// generate stellar colors
	float ri = ri0 + (ri1-ri0)*rnd();
	// find the absolute magnitude:
	float Mr = paralax.mr(ri);
	// assuming we're on the main sequence:
	float gr = paralax.gr(ri);
	// R magnitude from distance
	float r = stardist::m(D, s.Mr);
}

struct real_star
{
	float ri, gr, r;
};

class generate_model : public bin3d_pixelize_volume
{
public:
	map<float, galaxy_model> ri_models;	// disk and halo models for stars of a given ri color (spectral type)
public:
	cylindrical_pixelize_volume(const std::string &volfn, float dx, int n, pair<float, float> r, pair<float, float> ri, Radians phi0_ = 0)
		: bin3d_pixelize_volume(volfn, dx, n, r, ri), phi0(phi0_)
		{}

	double 	rho,	// galactocentric radius (set in action())
		D;	// earthcentric distance (set in action())
	mobject m;
	sdss_star ss;

	bool observe(real_star s)
	{
		// find r magnitude
		m.mag[0] = s.r;
		m.mag[1] = s.gr + s.r;
		m.mag[2] = s.r - s.ri;

		// add photometric errors, based on some braindead model of photometric errors
		// first approximation - all are drawn from some gaussian distribution, irrespective of magnitude
		m.magErr[0] = m.magErr[1] = m.magErr[2] = .02;
		FOR(0, 3)
		{
			m.mag[i] += gsl_ran_gaussian(rng, m.magErr[i]);
			ss.mag[i+1] = m.mag[i]; ss.magErr[i+1] = m.magErr[i]; // for ML color calculation
		}

		if(paralax(s))	// calculate the absolute magnitude and distances
		{
			m.D = s.earth.D;
			m.ml_mag[0] = s.ml_g;
			m.ml_mag[1] = s.ml_r;
			m.ml_mag[2] = s.ml_i;
		} else {
			m.D = 0;
		}

		// calculate position
		ecenequ(V3(x, y, z), D, m.ra, m.dec);
		DEG(m.ra); DEG(m.dec);

		return true;
	}

	virtual double action(float ddv)
	{
		// generate stars for each model, in this pixel
		rho = sqrt(sqr(x+Rg)+sqr(y));	// galactocentric radius
		D = sqrt(sqr(x)+sqr(y)+sqr(z));	// earthcentric distance
		m.ra = -1;
		FOREACH(ri_models)
		{
			const galaxy_model &model = (*i).second;
			double den = model.density(rho, z);
			int N = model.number(den, ddv);
			while(N)
			{
				real_star &s = model.star();
				if(observe(m, s))
				{
					// << write the star to DMM
				}
			}
		}
	}
};

cat_generator::generate(std::string &fn)
{
	binned_run br;
	binary_input_or_die(in, fn);
	in >> br;
	brs.dx = br.dx;
	FOREACH(br.pixels)
	{
		V3 v = (*i).first;
		v *= b
		binned_run::pixel &p = (*i).second;
		p.volume;
		
		double den = density(p);
		int N = number(den, V);
		if(N == 0) { return 0; }

		brs.pixels[k].uniqueVolume = p.volume;
	}
}

class cat_generator
{
public:
	double dx;	// volume pixel dx, used when generating star coordinates

protected:
	int N;
	double V;
	V3 p;

	int number(double den);
public:
	int init(const V3 &p, double v);
//	bool next(mobject &m);
};

int den_generator::init(const V3 &p_, double v)
{
	srand(time(NULL));

	V = v;
	p = p_;

	double den = density(p);
	N = number(den);
	
	return N;
}

#if 0
bool den_generator::next(mobject &m)
{
	if(N == 0) return false;
	--N;
	
	// assign a random coordinate within the pixel
	V3 pos;
	rand(pos, -dx/2, +dx/2);
	pos += p;

	// convert to equatorial coordinates
	ecenequ(pos, m.D, m.ra, m.dec);
	
	return true;
}
#endif

/////////////////////////////
#endif

#include "model.h"

struct direction
{
	Radians l, b;
	double cl, cb, sl, sb;
	
	double d, d2,
		x, y, z, r2, r, rcyl;

	direction(Radians l_, Radians b_, double d_ = 0.0)
	: l(l_), b(b_),
	  cl(cos(l_)), cb(cos(b_)), sl(sin(l_)), sb(sin(b_))
	{
		set_dist(d_);
	}

	void set_dist(const double dnew)
	{
		d = dnew;
		d2 = d*d;
		x = Rg - d*cl*cb;
		y = -d*sl*cb;
		z = d*sb;
		r2 = x*x + y*y + z*z;
		r = sqrt(r2);
		rcyl = sqrt(r2 - z*z);
	}
};

//
// Calculate differential of dr wrt. to dd
//
inline double dr_dd(const direction &d)
{
	return (d.d - Rg*d.cl*d.cb) / d.r;
}

//
// Number of stars per cubic parsec at point d (the Galactic model)
//
model_fitter m;
double model_rho(const direction &d)
{
	double rho = m.rho(d.rcyl, d.z);
	return rho;
}

//
// Calculate the number of stars per unit solid angle,
// per unit d.d, at a point d
//
inline double drho_dd(const direction &d)
{
/*	double drdd = dr_dd(d);
	double rho = model_rho(d) * d.r2 * drdd;*/
//	double rho = model_rho(d) * d.d2;
	double rho = model_rho(d) * d.d2 * d.d;
//	double rho = model_rho(d);
	return rho;
}

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <valarray>

///////////////////////////////////

#if 1
typedef int (*int_function)(double dd, double const* y, double *dydd, void *param);

bool
sample_integral(const std::valarray<double> &xv, std::valarray<double> &yv, int_function f, void *param)
{
	gsl_odeiv_step*    s = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 1);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(1e-6, 0.0);
	gsl_odeiv_evolve*  e = gsl_odeiv_evolve_alloc(1);

	const double a = xv[0];
	const double b = xv[xv.size()-1];

	double h = (b - a) / xv.size();	// initial step size
	double y = 0;			// initial integration value
	double x = a;			// initial integration point

	gsl_odeiv_system sys = {f, NULL, 1, param};

	// integration from 0 to dmax, but in ninterp steps
	// when integrating between dmin and dmax
	FOR(0, xv.size())
	{
		while (x < xv[i])
		{
			gsl_odeiv_evolve_apply (e, c, s,
						&sys,
						&x, xv[i], &h,
						&y);
		}
		yv[i] = y;
	}
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);

	return true;
}

double
integrate(double a, double b, int_function f, void *param)
{
	gsl_odeiv_step*    s = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 1);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(1e-6, 0.0);
	gsl_odeiv_evolve*  e = gsl_odeiv_evolve_alloc(1);

	double h = (b - a) / 1000;	// initial step size
	double y = 0;			// initial integration value
	double x = a;			// initial integration point

	gsl_odeiv_system sys = {f, NULL, 1, param};

	// integration from 0 to dmax, but in ninterp steps
	// when integrating between dmin and dmax
	while (x < b)
	{
		gsl_odeiv_evolve_apply (e, c, s,
					&sys,
					&x, b, &h,
					&y);
	}
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);

	return y;
}

#endif

///////////////////////////////////

//
// Calculate the integral of probability from which we can
// draw stars
//
struct distanceSampler
{
/*	gsl_spline *f;
	gsl_interp_accel *acc;*/
	spline f;
	
	double N; // Total number of stars in [dmin,dmax]

	distanceSampler(const direction &d, const double dmin, const double dmax);
	~distanceSampler();

/*	inline double distance(double u) { return gsl_spline_eval(f, u, acc); }*/
	inline double distance(double u) { return f(u); }
};

distanceSampler::~distanceSampler()
{
/*	gsl_spline_free(f);
	gsl_interp_accel_free(acc);*/
}

int
func(double dd, double const* y, double *dydd, void *param)
{
	direction &d = *(direction *)param;
	d.set_dist(dd);

	*dydd = drho_dd(d);
	return GSL_SUCCESS;
}

#if 1
distanceSampler::distanceSampler(const direction &d0, const double dmin, const double dmax)
{
	std::valarray<double> xv(1001), yv(1001);
	FOR(0, xv.size()) { xv[i] = double(i) / (xv.size()-1.); }
	xv = dmin + xv*(dmax-dmin);

	direction d = d0;
	sample_integral(xv, yv, func, (void *)&d);

	// normalization - convert the integrated density to probability
	N = yv[yv.size()-1];
	yv /= N;

	// construct spline approximation
	f.construct(&yv[0], &xv[0], yv.size());
}
#else
distanceSampler::distanceSampler(const direction &d0, const double dmin, const double dmax)
{
	gsl_odeiv_step* s = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 1);
	gsl_odeiv_control* c = gsl_odeiv_control_y_new(1e-6, 0.0);
	gsl_odeiv_evolve* e = gsl_odeiv_evolve_alloc(1);

	double h = 10;
	double y = 0;
	double dist = dmin;
	direction d = d0;

	gsl_odeiv_system sys = {func, NULL, 1, (void *)&d};

	// integration from 0 to dmax, but in ninterp steps
	// when integrating between dmin and dmax
	const int ninterp = 1000;
	std::valarray<double> yv(ninterp+1), dv(ninterp+1);
	for(int i = 0; i < yv.size(); i++)
	{
		double di = dmin + i * (dmax - dmin) / ninterp;

		while (dist < di)
		{
			gsl_odeiv_evolve_apply (e, c, s,
									&sys,
									&dist, di, &h,
									&y);
		}
		//cout << "(dist, y) = " << dist << " " << y << "\n";
		yv[i] = y;
		dv[i] = di;
	}
	gsl_odeiv_evolve_free(e);
	gsl_odeiv_control_free(c);
	gsl_odeiv_step_free(s);

	// normalization - convert the integrated density to probability
	N = yv[yv.size()-1];
	yv /= N;

/*	FOR(0, yv.size())
	{
		cout << "(d, y) = " << dv[i] << " " << yv[i] << "\n";
	}*/
	
	// construction of a spline spanning the dmin, dmax interval
	// This will, for range [0, 1] return distances from [dmin,dmax]
// 	f = gsl_spline_alloc(gsl_interp_cspline, yv.size());
// 	gsl_spline_init(f, &yv[0], &dv[0], yv.size());
// 	acc = gsl_interp_accel_alloc();
	f.construct(&yv[0], &dv[0], yv.size());

/*	FOR(0, 1001)
	{
		double u = double(i)/1000.;
		cout << u << " " << distance(u) << "\n";
	}
	exit(0);*/
}
#endif

class spline_func
{
protected:
	gsl_spline *f;
	gsl_interp_accel *acc;
public:
	spline_func(const std::string &fn, int x_idx = 0, int y_idx = 1);
	~spline_func();
	double operator()(const double x) const { return gsl_spline_eval(f, x, acc); }
};

spline_func::~spline_func()
{
	gsl_spline_free(f);
	gsl_interp_accel_free(acc);
}

spline_func::spline_func(const std::string &fn, int x_idx, int y_idx)
	: f(NULL), acc(NULL)
{
	text_input_or_die(in, fn);
	vector<double> x, y;
	load(in, x, x_idx, y, y_idx);
	ASSERT(x.size() == y.size());

	f = gsl_spline_alloc(gsl_interp_cspline, x.size());
	gsl_spline_init(f, &x[0], &y[0], x.size());
	acc = gsl_interp_accel_alloc();
}

void hdr(ostream &out, float ri, float Mr, float rmin, float rmax, double dmin, double dmax, double N)
{
	out << "# (ri, Mr)     = " << ri << " " << Mr << "\n";
	out << "# (rmin, rmax) = " << rmin << " " << rmax << "\n";
	out << "# (dmin, dmax) = " << dmin << " " << dmax << "\n";
	out << "# Ntot = " << N << "\n";
}

Nag_Spline fit_nag_spline(std::valarray<double> &x, std::valarray<double> &y, std::valarray<double> &w)
{
	int m = y.size();
	Nag_Spline spline;
	Nag_Comm warmstart;
	double s = 2;
	double fp;
	cerr << "SF" << x.size() << " " << y.size() << " " << w.size() << " " << m << "\n";
#if 0
	double a[] = {1., 2., 3., 4., 5., 6.};
	double b[] = {1.1, 2.3, 3.2, 4.0, 5.1, 6.4};
	double wt[] = {1., 1., 1., 1., 1., 1.};
	nag_1d_spline_fit(Nag_Cold, 6, &a[0], &b[0], &wt[0], s, 10,
		&fp, &warmstart, &spline, NAGERR_DEFAULT);
#endif
	FOR(0, m) { cerr << x[i] << " " << y[i] << " " << w[i] << "\n"; }
	nag_1d_spline_fit(Nag_Cold, m, &x[0], &y[0], &w[0], s, m/2.,
		&fp, &warmstart, &spline, NAGERR_DEFAULT);
	cerr << "EF\n";

	return spline;
}

void histogram(std::valarray<double> &h, 
	const std::valarray<double> &b, 
	const std::valarray<double> &x,
	double *nless = NULL, double *nmore = NULL)
{
	const double *b0 = &b[0], *b1 = &b[b.size()];
	h.resize(b.size()-1);
	h *= 0.;
	double dummy;
	if(nless == NULL) { nless = &dummy; }
	if(nmore == NULL) { nmore = &dummy; }
	*nmore = *nless = 0.;
	FOR(0, x.size())
	{
		double v = x[i];
		const double *bin = upper_bound(b0, b1, v);
		if(bin == b0) {(*nless)++; continue; }	// lower than lowest
		if(bin == b1) {(*nmore)++; continue; }	// higher than highest bin
		long idx = bin - b0 - 1;
		h[idx]++;
	}
	cerr << *nmore << " " << *nless << "\n";
}

void midpoints(std::valarray<double> &bh, const std::valarray<double> &b)
{
	bh.resize(b.size()-1);
	FOR(0, bh.size())
	{
		bh[i] = (b[i] + b[i+1]) / 2.;
	}
}

void make_bins(std::valarray<double> &b, double b0, double b1, double dx)
{
	int nm = (int)((b1 - b0) / dx);
	if(gsl_fcmp(dx*nm, b1 - b0, 2*GSL_DBL_EPSILON) < 0)
	{
		++nm;
	}
	++nm;

	b.resize(nm);
	FOR(0, b.size())
	{
		b[i] = b0 + dx*double(i);
	}
	b[b.size()-1] = b1;
}

void print_hist(std::ostream &out, 
	const std::valarray<double> &bins, const std::valarray<double> &h,
	double *hl = NULL, double *hh = NULL)
{
	if(hl != NULL) { out << bins[0] << "\t" << *hl << "\n"; }
	FOR(0, bins.size()-1)
	{
		double b = (bins[i] + bins[i+1]) / 2.;
		out << b << "\t" << h[i] << "\n";
	}
	if(hh != NULL) { out << bins[bins.size()-1] << "\t" << *hh << "\n"; }
}

void simulate(int n, Radians l, Radians b, float ri, float rmin, float rmax)
{
	plx_gri_locus plx;
	double dmin, dmax;
	plx.distance_limits(dmin, dmax, ri, ri, rmin, rmax);
	double Mr = plx.Mr(ri);

	direction d0(l, b);
	distanceSampler ds(d0, dmin, dmax);
	hdr(cout, ri, Mr, rmin, rmax, dmin, dmax, ds.N);

	spline_func photoerr("sigma_mean.txt");

	int seed = 42;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, seed);

	ofstream denf("mq_den.txt");
	hdr(denf, ri, Mr, rmin, rmax, dmin, dmax, ds.N);
	model_fitter mold = m;
	FOR(0, 1001)
	{
		double dd = dmin + (dmax-dmin)*double(i)/1000.;
		d0.set_dist(dd);
		denf << dd;
		m = mold;
		denf << " " << drho_dd(d0);
/*		m.param("ht") = 510;
		m.param("h") = 260;*/
/*		m.param("rho0") = 0.95;
		m.param("f") = 0.062;
		m.param("ht") = 800;
		m.param("h") = 260;*/
		m.rho0 = 1.;
		m.f = 0.04;
		m.ht = 1000;
		m.h = 260;
/*		for(int i = 240; i != 300; i+=10)
		{
			m.param("h") = i;
			denf << " " << drho_dd(d0);
		}*/
		denf << " " << drho_dd(d0);
		denf << "\n";
	}
	m = mold;

	ofstream simf("mq_sim.txt");
	hdr(simf, ri, Mr, rmin, rmax, dmin, dmax, ds.N);
	valarray<double> x(n), x0(n);
	FOR(0, n)
	{
		double u = gsl_rng_uniform(rng);
		double d = ds.distance(u);

		float m_r = stardist::m(d, Mr);
		double s_mean = photoerr(m_r);
		double err = gsl_ran_gaussian(rng, s_mean);
		//err = sqrt(sqr(err) + sqr(0.1));
		err = 0.3;
		m_r += err;
		double d_with_err = stardist::D(m_r, Mr);

		x[i] = d_with_err;
		x0[i] = d;
		simf << u << " " << d << " " << d_with_err << " " << err << "\n";
	}
	gsl_rng_free(rng);

	// construct histogram
	valarray<double> h, bh, bins, w;
	double hl, hm;
	make_bins(bins, 0, 1800, 10);
	midpoints(bh, bins);
	histogram(h, bins, x, &hl, &hm);
	w.resize(h.size());
	FOR(0, w.size()) { w[i] = h[i] > 0. ? 1./h[i] : 1.; }

	// fit spline
	Nag_Spline spl = fit_nag_spline(bh, h, w);
	//cerr << "Spline fitted.\n";
	FOR(0, bh.size())
	{
		double y;
		nag_1d_spline_evaluate(bh[i], &y, &spl, NAGERR_DEFAULT);
		//cout << i << " " << bh[i] << " " << h[i] << " " << y << "\n";
	}

	// dump histogram
	ofstream hist("mq_hist.txt");
	print_hist(hist, bins, h, &hl, &hm);

	// apply the correction
	valarray<double> xc(x.size());
	FOR(0, x.size())
	{
		xc[0] = -1;
		double s = x[i];
		double y = log(s);
		double delta2 = sqr(0.46*.5);
		double dd = exp(y + delta2);
		if(dd > bh[bh.size()-1]) { continue; }
		if(s > bh[bh.size()-1]) { continue; }
		double ydd, yd;
		nag_1d_spline_evaluate(dd, &ydd, &spl, NAGERR_DEFAULT);
		nag_1d_spline_evaluate(s, &yd, &spl, NAGERR_DEFAULT);
		double ratio = ydd/yd;
		double r = s*exp(.5*delta2)*ratio;
//		cout << s << " " << r << " " << ratio << "\n";
		xc[i] = r;
	}

	valarray<double> hc, h0;
	histogram(hc, bins, xc, &hl, &hm);
	histogram(h0, bins, x0, &hl, &hm);
	FOR(0, bh.size())
	{
		double y;
		nag_1d_spline_evaluate(bh[i], &y, &spl, NAGERR_DEFAULT);
		cout << bh[i] << " " << y << " " << hc[i] << " " << h0[i] << "\n";
	}
}

namespace stlops
{
	template<typename T>
		inline OSTREAM(const std::vector<T> &a)
		{
			out << "[";
			FOREACH(a)
			{
				if(i != a.begin()) { out << ", "; }
				out << *i;
			}
			out << "]";
			return out;
		}
}

void simulate()
{
// 	using namespace stlops;
// 	vector<int> ii; ii.push_back(10); ii.push_back(20);
// 	cout << ii << "\n";
// 	return;

// 	m.add_param("rho0", 1);
// 	m.add_param("l", 3000);
// 	m.add_param("h", 270);
// 	m.add_param("z0", 25);
// 
// 	m.add_param("f", 0.04);
// 	m.add_param("lt", 3500);
// //	m.add_param("ht", 400);
// 	m.add_param("ht", 950);
// 
// 	m.add_param("fh", 0.0001);
// 	m.add_param("q", 1.5);
// 	m.add_param("n", 3);

	m.rho0 = 1;
	m.l = 3000;
	m.h = 270;
	m.z0 = 25;

	m.f = 0.04;
	m.lt = 3500;
	m.ht = 1400;

	m.fh = 0.0001;
	m.q = 1.5;
	m.n = 3;

	//simulate(100000, 0., rad(90.), 0., 15000.);
	//simulate(100000, rad(33.), rad(15.), 0., 15000.);
//	simulate(100000, rad(0), rad(90.), 1.1, 15.0, 21.5);
//	simulate(500000, rad(180), rad(30.), 1.1, 15.0, 21.5);
	simulate(500000, rad(0), rad(90.), 1.1, 15.0, 21.5);
//	simulate(500000, rad(0), rad(10.), 0.1, 15.0, 21.5);
}

extern "C"
{
	#include "gpc/gpc.h"
}

gpc_polygon make_polygon(const RunGeometry &geom, const lambert &proj, double dx)
{
	Mask mask(geom);
	Radians mu0 = geom.muStart;
	Radians mu1 = geom.muStart + mask.length();

	vector<double> ex, ey;

	gpc_polygon p = {0, 0, NULL};

	FORj(col, 0, 6)
	{
		Radians nu0 = mask.lo(col);
		Radians nu1 = mask.hi(col);
		Radians l, b, mu, nu;

		// bottom edge
		for(mu = mu0; mu < mu1; mu += dx)
		{
			coordinates::gcsgal(geom.node, geom.inc, mu, nu0, l, b);
			ex.push_back(l);
			ey.push_back(b);
		}
		// right edge
		for(nu = nu0; nu < nu1; nu += dx)
		{
			coordinates::gcsgal(geom.node, geom.inc, mu1, nu, l, b);
			ex.push_back(l);
			ey.push_back(b);
		}

		// top edge
		for(mu = mu1; mu > mu0; mu -= dx)
		{
			coordinates::gcsgal(geom.node, geom.inc, mu, nu1, l, b);
			ex.push_back(l);
			ey.push_back(b);
		}

		// left edge
		for(nu = nu1; nu > nu0; nu -= dx)
		{
			coordinates::gcsgal(geom.node, geom.inc, mu0, nu, l, b);
			ex.push_back(l);
			ey.push_back(b);
		}

		// convert to contour in lambert coordinates
		gpc_vertex v[ex.size()];
		FOR(0, ex.size())
		{
			proj.convert(ex[i], ey[i], v[i].x, v[i].y);
		}

		// add contour to polygon
		gpc_vertex_list c = { ex.size(), v };
		gpc_add_contour(&p, &c, false);

		ex.clear();
		ey.clear();
	}
	return p;
}

//
// Simple crossing number algorithm, valid for non-self intersecting polygons.
// See http://softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm for details.
//
bool in_contour(const gpc_vertex &t, const gpc_vertex_list &vl)
{
	int cn = 0;
	const int n = vl.num_vertices;
	for(int i = 0; i != n; i++)
	{
		gpc_vertex &a = vl.vertex[i];
		gpc_vertex &b = (i + 1 == n) ? vl.vertex[0] : vl.vertex[i+1];

		if   (((a.y <= t.y) && (b.y > t.y))    		// an upward crossing
		   || ((a.y > t.y)  && (b.y <= t.y)))  		// a downward crossing
		{
			// compute the actual edge-ray intersect x-coordinate
			double vt = (float)(t.y - a.y) / (b.y - a.y);
			if (t.x < a.x + vt * (b.x - a.x)) 	// t.x < intersect
                		++cn;				// a valid crossing of y=t.y right of t.x
		}
	}

	return (cn&1);    // 0 if even (out), and 1 if odd (in)	
}

bool in_polygon(const gpc_vertex &t, const gpc_polygon &p)
{
	bool in = false;
	FOR(0, p.num_contours)
	{
		bool inc = in_contour(t, p.contour[i]);
		in = (in != inc);
//		if(inc) { cerr << "+"; }
	}
//	cerr << "\n";
	return in;
}

class striped_polygon
{
public:
	int n;
	double dx;
	vector<gpc_polygon> stripes;
	double x0, x1;
public:
	striped_polygon(gpc_polygon &p, int n);
	~striped_polygon();
};

void sm_write(const std::string &fn, const gpc_polygon &p);

void poly_bounding_box(double &x0, double &x1, double &y0, double &y1,
	const gpc_polygon &p)
{
	// find minimum and maximum vertex
	x0 = x1 = p.contour[0].vertex[0].x;
	y0 = y1 = p.contour[0].vertex[0].y;
	FOR(0, p.num_contours)
	{
		const int n = p.contour[i].num_vertices;
		FORj(j, 0, n)
		{
			const double x = p.contour[i].vertex[j].x;
			const double y = p.contour[i].vertex[j].y;
			x0 = std::min(x0, x);
			x1 = std::max(x1, x);
			y0 = std::min(y0, y);
			y1 = std::max(y1, y);
		}
	}
}

striped_polygon::striped_polygon(gpc_polygon &p, int n_)
	: n(n_)
{
	// split the polygon p into n stripes, for faster lookup
	// when doing poin-in-polygon searches.
	double y0, y1;
	poly_bounding_box(x0, x1, y0, y1, p);
	y0 += .01*y0; y1 += .01*y1;

	// split into n polygons
	dx = (x1-x0) / n;
	stripes.resize(n);
	FOR(0, n)
	{
		double xa = x0 +     i * dx;
		double xb = x0 + (i+1) * dx;
		gpc_vertex v[] = {{xa, y0}, {xb, y0}, {xb, y1}, {xa, y1}};
		gpc_vertex_list vl = {4, v};
		gpc_polygon rect = {1, NULL, &vl};

		gpc_polygon_clip(GPC_INT, &p, &rect, &stripes[i]);
//		cerr << i << " " << xa << " " << xb << " " << i * (x1-x0) / n << " " << (i+1) * (x1-x0) / n << "\n";

/*		sm_write("stripe.gpc.txt", stripes[i]);
		exit(0);*/
	}
}

striped_polygon::~striped_polygon()
{
	FOR(0, stripes.size()) { gpc_free_polygon(&stripes[i]); }
}

bool in_polygon(const gpc_vertex& t, striped_polygon &p)
{
	// find the right poly stripe
	int s = (int)floor((t.x - p.x0) / p.dx );
	ASSERT(s < p.n)
	{
		cerr << t.x << " " << t.y << " " << s << " " << p.n << "\n";
	}
	return in_polygon(t, p.stripes[s]);
}

#include <sys/time.h>
double seconds()
{
	timeval tv;
	int ret = gettimeofday(&tv, NULL);
	assert(ret == 0);
	return double(tv.tv_sec) + double(tv.tv_usec)/1e6;
}

void testIntersection(const std::string &prefix)
{
	lambert proj(rad(90), rad(90));
	double l = rad(300);
	double b = rad(60);

	gpc_polygon allsky;
	FILE *fp = fopen((prefix + ".gpc.txt").c_str(), "r");
	gpc_read_polygon(fp, 1, &allsky);
	fclose(fp);
	striped_polygon sky(allsky, 100);
//	gpc_polygon &sky = allsky;

	// specific point test
	double x, y;
	proj.convert(l, b, x, y);
	x = 0.272865; y = -0.387265;
	//x = 0.280607; y = -0.387462;
	x = 1.05179; y = 0.807653;

	cerr << "Inside: " << in_polygon((gpc_vertex const&)make_pair(x, y), sky) << "\n";
//exit(0);
	// monte-carlo tests
	double begin = seconds();

	srand(23);
	RunGeometryDB db;
	for(int k = 0; k != 500; k++) {
	FOREACH(db.db)
	{
		const RunGeometry &geom = (*i).second;
		Mask mask(geom);

		Radians mu0 = geom.muStart;
		Radians mu1 = geom.muStart + mask.length();

		FORj(col, 0, 6)
		{
			Radians nu0 = mask.lo(col);
			Radians nu1 = mask.hi(col);
			Radians l, b, mu, nu;
			mu = math::rnd(mu0, mu1);
			nu = math::rnd(nu0, nu1);
			coordinates::gcsgal(geom.node, geom.inc, mu, nu, l, b);
			if(b < 0) { continue; }

			proj.convert(l, b, x, y);
			if(!in_polygon((gpc_vertex const&)make_pair(x, y), sky))
			{
				cerr << "-";
// 				cerr << "Test failed for:\n";
// 				cerr << " l, b = " << deg(l) << " " << deg(b) << "\n";
// 				cerr << " mu, nu = " << deg(mu) << " " << deg(nu) << "\n";
// 				cerr << " x, y = " << x << " " << y << "\n";
// 				cerr << " run, col = " << geom.run << " " << col << "\n";
// 				cerr << " mu0, mu1 = " << deg(mu0) << " " << deg(mu1) << "\n";
// 				cerr << " nu0, nu1 = " << deg(nu0) << " " << deg(nu1) << "\n";
			}
		}
	}
		cerr << ".";
	}
	cerr << "\n";
	cerr << "time: " << seconds() - begin << "\n";
	gpc_free_polygon(&allsky);
}

gpc_polygon make_circle(double x0, double y0, double r, double dx)
{
	gpc_polygon p = {0, 0, NULL};
	int n = (int)(2*ctn::pi*r / dx);
	gpc_vertex v[n];
	gpc_vertex_list c = { n, v };
	FOR(0, n)
	{
		v[i].x = x0+r*cos(i*dx/r);
		v[i].y = y0+r*sin(i*dx/r);
	}
	gpc_add_contour(&p, &c, false);
	return p;
}

double polygon_area(const gpc_polygon &p)
{
	double A = 0;
	FOR(0, p.num_contours)
	{
		int n = p.contour[i].num_vertices;
		gpc_vertex *v = p.contour[i].vertex;
		double cA = 0;
		FORj(j, 0, n)
		{
			gpc_vertex &a = v[j];
			gpc_vertex &b = (j + 1 == n) ? v[0] : v[j+1];

			cA += a.x*b.y - b.x*a.y;
		}
		cA = abs(cA) * (p.hole[i] ? -1 : 1);
		A += cA;
	}
	A *= 0.5;
	return abs(A);
}

void sm_write(const std::string &fn, const gpc_polygon &p)
{
	// dump polygons in SM compatible format
	ofstream out(fn.c_str());
	FOR(0, p.num_contours)
	{
		gpc_vertex *v = p.contour[i].vertex;
		FORj(j, 0, p.contour[i].num_vertices)
		{
			out << v[j].x << " " << v[j].y << " " << i << "\n";
		}
		out << "#\n";
	}
}

//lambert lnorth(rad(90), rad(90)), lsouth(rad(90), rad(-90));
void makeSkyMap(const std::string &prefix, const lambert &proj)
{
	Radians dx = rad(.25); /* polygon sampling resolution in radians */
//	Radians dx = rad(1); /* polygon sampling resolution - good for footprint plots */
	cerr << "hemisphere = " << prefix << ", dx = " << dx << " radians\n";

	RunGeometryDB db;
	double x, y;
	gpc_polygon sky = {0, 0, NULL};
/*	proj.convert(rad(0.), rad(-30.), x, y);
	cerr << sqrt(x*x+y*y) << "\n";*/
	gpc_polygon circle = make_circle(0., 0., sqrt(2.), dx);

	text_input_or_die(in, "/home/scratch/projects/galaxy/workspace/catalogs/runs.txt");
	std::set<int> runs;
	load(in, runs, 0);
	cerr << "Processing " << runs.size() << " runs.\n";

	int k = 0;
	FOREACH(runs)
	{
		const RunGeometry &geom = db.getGeometry(*i);
// 		cerr << geom.run << " ";
 		cerr << ".";

		//if(geom.run != 752) continue;
		//if(geom.run != 752 && geom.run != 756) continue;

		gpc_polygon rpoly = make_polygon(geom, proj, dx);
		gpc_polygon_clip(GPC_INT, &rpoly, &circle, &rpoly);
		gpc_polygon_clip(GPC_UNION, &sky, &rpoly, &sky);

// 		double A = polygon_area(rpoly);
// 		gpc_free_polygon(&rpoly);
// 
// 		int nvert = 0;
// 		FOR(0, sky.num_contours) { nvert += sky.contour[i].num_vertices; }
// 		cerr << " [" << A*sqr(deg(1)) << "] [" << sky.num_contours << " contours, " << nvert << " vertices]\n";
	}
	cerr << "\n";

	int nvert = 0;
	FOR(0, sky.num_contours) { nvert += sky.contour[i].num_vertices; }
	cerr << "total [" << polygon_area(sky)*sqr(deg(1)) << "deg2 area, "
	     << sky.num_contours << " contours, " << nvert << " vertices]\n";

	// store the footprint polygon
	sm_write(prefix + ".foot.txt", sky);
	FILE *ofp = fopen((prefix + ".gpc.txt").c_str(), "w");
	gpc_write_polygon(ofp, 1, &sky);
	fclose(ofp);

	// free memory
	gpc_free_polygon(&sky);
	gpc_free_polygon(&circle);
}

void test_nurbs();
void test_lapack();
namespace sim
{
	void calculate_pdf(const std::string &prefix);
	void montecarlo(const std::string &prefix, unsigned int K);
}

int main(int argc, char **argv)
{
//	simulate();
/*	makeSkyMap("north", lambert(rad(90), rad(90)));
 	makeSkyMap("south", lambert(rad(-90), rad(-90)));*/
//	testIntersection("north");
//	test_nurbs();
//	test_lapack();

 	sim::calculate_pdf("north");
	sim::montecarlo("north", 100);

	return 0;
}
#if 0

void binBitImage(Radians ndx, XImage &img, const BitImage &b)
{
	Radians dx = b.dx;
	int I0 = img.x() / 2;
	int J0 = img.y() / 2;
	double A = sqr(dx);
	FORj(i, 0, b.w)
	{
		FORj(j, 0, b.h)
		{
			if(!b.get(i-b.x0, j-b.y0)) continue;
			Radians x = (i - b.x0)*dx;
			Radians y = (j - b.y0)*dx;

			int I = (int)floor((x - dx/2.) / ndx + .5) + I0;
			int J = (int)floor((y - dx/2.) / ndx + .5) + J0;

			img(I, J) += A;
		}
	}
}

struct gmake
{
	XImage img;
	Radians dx;
};

struct bixel
{
	Radians x[4], y[4], l[4], b[4];	// lambert and galactic coordinates of bixel side centers
	Radians dx, dy, dm;	// bixel size
	double m[2];	// magnitude space

	void init(double d0, double d1, Radians x0, Radians y0, Radians dx, Radians dy)
	{
		m[0] = m0; m[1] = m1; dm = m[1] - m[0];

		this->dx = dx; this->dy = dy;
		x[0] = x0; x[1] = x0-dx/2; x[2] = x0; x[3] = x0+dx/2;
		y[3] = y0; y[0] = y0-dy/2; y[1] = y0; y[2] = y0+dy/2;

		inverse();
	}

	void inverse()
	{
		FOR(0, 4) { proj.inverse(x[i], y[i], l[i], b[i]); }
	}

	void split(int what, bixel &a, bixel &b)
	{
		switch(what)
		{
		case 0:
			a.init(x[0]-dx/4., y[3], dx/2, dy);	// left
			b.init(x[0]+dx/4., y[3], dx/2, dy);	// right
			break;
		case 1:
			a.init(x[0], y[3]-dx/4., dx, dy/2);	// top
			b.init(x[0], y[3]+dx/4., dx, dy/2);	// bottom
			break;
		case 2:
			a = *this; a.m[1] = m[0] + dm/2;	// closer
			b = *this; a.m[0] = m[0] + dm/2;	// farther
			break;
		}
	}
};

class galaxy_model
{
public:
	virtual bool delta(bool split[], const bixel &b) = 0;
	virtual bool generate_stars(const bixel &b) = 0;
	virtual bool next_star() = 0;
};

class disk_model : public galaxy_model
{
public:
	double den0, h, l, z0;		// disk parameters
	float ri0, ri1;
	
	double eps;
public:
	double density(const V3 &p);

	virtual bool delta(bool split[], const bixel &b);
};

double disk_model::density(const double r, const double z)
{
	return den0 * exp(r/l + (z-z0)/h);
}

plx_gri_locus paralax;
double disk_model::delta(bool split[], const bixel &b)
{
	// convert bixel side points to galactic coordinates
	double d0 = paralax.
}

int galaxy_model::number(double den, double V)
{
	double n = den*V;
	int N = (int)n;

	// add Poisson noise
	// TODO: draw a number from N(mu=n, sigma=sqrt(n)) and add it to n

	// the fraction is the probability of having a particle there
	n -= N;
	float prob = float(rand())/RAND_MAX;
	if(prob < n)
	{
		N++;
	}

	return N;
}

void ()
{
	DMMArray<mobject> out;
	out.create("simulation.dmm");

	int I0 = img.x() / 2;
	int J0 = img.y() / 2;

	lambert proj(rad(90), rad(90));

	m0 = 14; m1 = 23; dm = 1;
	epsilon = 0.01;

	// for each pixel in the image
	FOREACH(img)
	{
		Radians x = dx*(i.x - I0);
		Radians y = dx*(i.y - J0);

		// process the beam
		stack<bixel> bixels;
		bixel b;
		for(double m = m0; m <= m1; m += md)
		{
			b.init(m, m+md, x, y, dx, dx);

			// for each model
			FOREACHj(mod, models)
			{
				galaxy_model &model = *mod;
				bixels.push_back(b);
	
				// process the bixel
				bool delta[3];
				vector<bixels> tmpbix;
				while(bixels.size())
				{
					bixel b = bixels.top();
					bixels.pop();

					// check for tolerances
					if(model.delta(delta, b))
					{
						// split the bixel
						tmpbix.push_back(b);
						int at0 = 0, at = 0;
						bixel a, b;
						FORj(j, 0, 3)
						{
							if(!delta[j]) continue;

							at = tmpbix.size();
							FOR(at0, at)
							{
								tmpbix[i].split(j, a, b);
								tmpbix.push_back(a);
								tmpbix.push_back(b);
							}
							at0 = at;
						}

						// push split pixels to stack
						FOR(at0, at) { bixels.push_back(tmpbix[i]); }
						tmpbix.clear();
		
						continue;
					}

					// generate stars
					int N = model.generate_stars(b);
					mobject &obj;
					FORj(n, 0, N)
					{
						model.next_star(obj);
						out.push_back(obj);
					}
				}
			}
		}
	}
}

#endif
#endif
