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
model m;
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
	double rho = model_rho(d) * d.d2;
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

//
// Calculate the integral of probability from which we can
// draw stars
//
struct distanceSampler
{
	gsl_spline *f;
	gsl_interp_accel *acc;
	
	double N; // Total number of stars in [dmin,dmax]

	distanceSampler(const direction &d, const double dmin, const double dmax);
	~distanceSampler();

	inline double distance(double u) { return gsl_spline_eval(f, u, acc); }
};

int
func(double dd, double const* y, double *dydd, void *param)
{
	direction &d = *(direction *)param;
	d.set_dist(dd);

	*dydd = drho_dd(d);
	return GSL_SUCCESS;
}


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
	f = gsl_spline_alloc(gsl_interp_cspline, yv.size());
	gsl_spline_init(f, &yv[0], &dv[0], yv.size());
	acc = gsl_interp_accel_alloc();

/*	FOR(0, 1001)
	{
		double u = double(i)/1000.;
		cout << u << " " << distance(u) << "\n";
	}
	exit(0);*/
}

distanceSampler::~distanceSampler()
{
	gsl_spline_free(f);
	gsl_interp_accel_free(acc);
}

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
	model mold = m;
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
		m.param("rho0") = 1.;
		m.param("f") = 0.04;
		m.param("ht") = 1000;
		m.param("h") = 260;
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
	FOR(0, n)
	{
		double u = gsl_rng_uniform(rng);
		double d = ds.distance(u);

		float m_r = stardist::m(d, Mr);
		double s_mean = photoerr(m_r);
		double err = gsl_ran_gaussian(rng, s_mean);
		err = sqrt(sqr(err) + sqr(0.1));
		m_r += err;
		double d_with_err = stardist::D(m_r, Mr);

		simf << u << " " << d << " " << d_with_err << " " << err << "\n";
	}

	gsl_rng_free(rng);
}

namespace stlops
{
	template<typename T>
		inline OSTREAM(const std::vector<T> &a)
		{
			out << "[";
			FOREACH(typename std::vector<T>::const_iterator, a)
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

	m.add_param("rho0", 1);
	m.add_param("l", 3000);
	m.add_param("h", 270);
	m.add_param("z0", 25);

	m.add_param("f", 0.04);
	m.add_param("lt", 3500);
	m.add_param("ht", 1400);

	m.add_param("fh", 0.0001);
	m.add_param("q", 1.5);
	m.add_param("n", 3);
	
	//simulate(100000, 0., rad(90.), 0., 15000.);
	//simulate(100000, rad(33.), rad(15.), 0., 15000.);
//	simulate(100000, rad(0), rad(90.), 1.1, 15.0, 21.5);
//	simulate(500000, rad(180), rad(30.), 1.1, 15.0, 21.5);
	simulate(500000, rad(0), rad(90.), 1.1, 15.0, 21.5);
//	simulate(500000, rad(0), rad(10.), 0.1, 15.0, 21.5);
}



struct BitImage
{
	const int nbits;

	int w, h;
	int x0, y0;
	double dx;

	valarray<int> data;
	
	BitImage(double dx_, int w_, int h_)
		: nbits(sizeof(int)*8), dx(dx_), w(w_), h(h_), x0(w/2), y0(h/2), data(((long long)w)*h/nbits + 1)
	{
		cerr << "w*h: " << ((long long)w)*h << "\n";
		cerr << "w*h/nbits + 1: " << ((long long)w)*h/nbits + 1 << "\n";
		cerr << "nbits: " << nbits << "\n";
		cerr << "Size: " << data.size() << "\n";
	}

	void index(int x, int y, int &byte, int &bit) const
	{
		long long ind = ((long long)(y+y0))*h + (x+x0);
		byte = ind / nbits;
		bit = ind % nbits;

		if(byte < 0 || byte >= data.size())
		{
			cerr << "Error: x = " << x << ", y = " << y << ", byte = " << byte << ", data.size() = " << data.size() << "\n";
			ASSERT(!(ind < 0 || ind >= data.size()));
		}
	}

	bool get(int x, int y) const
	{
		int byte, bit;
		index(x, y, byte, bit);

		int v = data[byte];
		return v & (1 << bit);
	}
	
	void set(int x, int y)
	{
		int byte, bit;
		index(x, y, byte, bit);

		int &v = data[byte];
		v |= (1 << bit);
	}
};

void makeSkyMap()
{
	Radians dx = rad(.125 * 1/60.); /* sampling resolution in radians */
	cerr << "dx = " << dx << "\n";
	int l = int(2*sqrt(2.)/dx + 1) + 2;
	cerr << "l = " << l << "\n";
	BitImage north(dx, l, l),
		south(dx, 1, 1);
	lambert lnorth(rad(90), rad(90)), lsouth(rad(90), rad(-90));

	RunGeometryDB db;
	double x, y;
	FOREACH(RunGeometryDB::db_t::iterator, db.db)
	{
		const RunGeometry &geom = (*i).second;
		Mask mask(geom);
		Radians muEnd = geom.muStart + mask.length() + dx;
		cerr << geom.run << " ";
		FORj(col, 0, 6)
		{
			for(Radians nu = mask.lo(col) + dx/2; nu < mask.hi(col) + dx; nu += dx)
			{
				for(Radians mu = geom.muStart + dx/2; mu < muEnd; mu += dx)
				{
					Radians ra, dec;
					coordinates::gcsequ(geom.node, geom.inc, mu, nu, ra, dec);

					dec *= -1;
					if(dec <= 0) continue;

					BitImage &img = dec > 0 ? north : south;
					lambert &proj = dec > 0 ? lnorth : lsouth;

					proj.convert(ra, dec, x, y);
//					proj.convert(0, 0, x, y);
					int I = (int)floor(x/dx + 0.5);
					int J = (int)floor(y/dx + 0.5);
					if(abs(I) > img.x0 || abs(J) > img.y0)
					{
						cerr << "Problem: " << I << ", " << J;
						cerr << " " << deg(ra) << ", " << deg(dec);
						cerr << " " << x << ", " << y << "\n";
						exit(-1);
					}
					img.set(I, J);
				}
			}
			cerr << col+1;
		}
		cerr << "\n";
	}

	int sum = 0;
	FOR(0, north.data.size()) {
		int k = north.data[i];
		int m = 1;
		FORj(j, 0, 32) {
			sum += (k & m) != 0;
			m <<= 1;
		}
	}
	cerr << "North area: " << sum/*sqr(deg(dx))*/ << " square degrees";
}

int main(int argc, char **argv)
{
	simulate();
	//makeSkyMap();
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
