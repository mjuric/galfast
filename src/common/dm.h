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

#ifndef __dm_h
#define __dm_h

#include "fitsloader.h"
#include "analysis.h"
#include "xcat.h"
#include "paralax.h"

#include <astro/constants.h>
#include <astro/coordinates.h>
#include <astro/sdss/rungeometry.h>
#include <astro/sdss/photometry.h>
#include <astro/math.h>
#include <astro/math/vector.h>
#include <astro/io/format.h>
#include <astro/system/memorymap.h>
#include <astro/system/fs.h>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <cmath>
#include <map>
#include <fstream>
#include <valarray>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <iomanip>

class packed_id
{
protected:
	long long pack;
	static const long long RUN;
	static const long long COL;
	static const long long FIELD;
	static const long long OBJID;
public:
	packed_id(int run_, int col_, int field_, int objid_)
	{
		pack =	(objid_ << 0) |
			(field_ << 14) | 
			(col_ << 26) |
			(((long long)run_) << 29);

		ASSERT(run() == run_);
		ASSERT(col() == col_);	
		ASSERT(field() == field_);
		ASSERT(objid() == objid_);
	}
	packed_id() : pack(0) {}
	
	int run() const { return (pack & RUN) >> 29; }
	int col() const { return (pack & COL) >> 26; }
	int field() const { return (pack & FIELD) >> 14; }
	int objid() const { return (pack & OBJID) >> 0; }
	
	long long id() const { return pack; }
	operator long long&() { return pack; }
};

inline OSTREAM(const packed_id &pid)
{
	return out << "[" << pid.run() << " " << pid.col() << " " << pid.field() << " " << pid.objid() << "]";
}

///

// magnitudes
#define M_ML_G			(10 + 0)
#define M_ML_R			(10 + 1)
#define M_ML_I			(10 + 2)

struct sdss_color : public std::pair<int, int>
{
public:
	#define COLOR_P1S	-1
	#define COLOR_P2S	-2
	#define COLOR_MR	-3
	#define COLOR_MV	-4

	std::string name;

	#define TYPE_COLORMAG	0
	#define TYPE_SIGMA	-1
	#define TYPE_FIELD	-2	// general accessor for any field of sdss_star
	int type;

public:
	static int bandIdx(const char c)
	{
		switch(c)
		{
			case 'g' : return 0;
			case 'r' : return 1;
			case 'i' : return 2;
			case 'z' : return 3;
			case 'u' : return 4;
			case 'G' : return M_ML_G;
			case 'R' : return M_ML_R;
			case 'I' : return M_ML_I;
			default:
				std::cerr << "Unknown band (" << c << ").\n";
				ASSERT(0);
		}
	}

	static char bandName(int i)
	{
		switch(i)
		{
			case 0 : return 'g';
			case 1 : return 'r';
			case 2 : return 'i';
			case 3 : return 'z';
			case 4 : return 'u';
			case M_ML_G : return 'G';
			case M_ML_R : return 'R';
			case M_ML_I : return 'I';
			default:
				std::cerr << "Unknown band index (" << i << ").\n";
				ASSERT(0);
		}
	}

	bool getColor(const std::string &color)
	{
		if(color.size() == 1)
		{
			first = second = bandIdx(color[0]);
		}
		else if(color == "p1s") { first = COLOR_P1S; second = 0; }
		else if(color == "p2s") { first = COLOR_P2S; second = 0; }
		else if(color == "M_r") { first = COLOR_MR; second = 0; }
		else if(color == "M_V") { first = COLOR_MV; second = 0; }
		else if(color.size() == 2) {
			// assume this is a color
			first = bandIdx(color[0]);
			second = bandIdx(color[1]);
		} else {
			return false;
		}

		return true;
	}

	inline void set(const std::string &color); /* defined below */

	sdss_color(const std::string &color = "")
	{
		if(color.size()) { set(color); }
	}
};

/////////////////

struct starid
{
	int		fitsId;
	packed_id	sloanId;
	int		uniqId;
};

struct starmag
{
	int fitsId;
	float mag[5];			// measured magnitudes (uncorrected for extinction)
	float magErr[5];		// note: magnitude ordering is grizu
};
OSTREAM(const starmag &sm);

struct star : public starmag
{
	double ra, dec;
	float Ar;			// extinction
};

#define MAG_MEAN		1
#define MAG_MEDIAN		2
#define MAG_CLIPPED_MEAN	3

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

public:
	// accessors
	inline float gr() const { return mag[0] - mag[1]; }
	inline float ri() const { return mag[1] - mag[2]; }

	inline float u() const { return mag[4]; }
	inline float g() const { return mag[0]; }
	inline float r() const { return mag[1]; }
	inline float i() const { return mag[2]; }
	inline float z() const { return mag[3]; }

	inline float ml_gr() const { return ml_mag[0] - ml_mag[1]; }
	inline float ml_ri() const { return ml_mag[1] - ml_mag[2]; }
	
	inline float ml_g() const { return ml_mag[0]; }
	inline float ml_r() const { return ml_mag[1]; }
	inline float ml_i() const { return ml_mag[2]; }

	inline float sigma(const sdss_color &c) const
	{
		if(c.first == c.second) { return magErr[c.second]; }

		if(c.first >= 0) { return std::sqrt(peyton::sqr(magErr[c.first]) + peyton::sqr(magErr[c.second])); }

		ASSERT(0);
	}

	enum {FIELD_UNKNOWN = 0, FIELD_RA, FIELD_DEC};
	static int field_id(const std::string &name)
	{
		if(name == "ra") { return FIELD_RA; }
		if(name == "dec") { return FIELD_DEC; }
		return FIELD_UNKNOWN;
	}

	inline double field(const sdss_color &c) const
	{
		switch(c.type) {
		case TYPE_SIGMA:    return sigma(c);
		case TYPE_COLORMAG: return color(c);
		case TYPE_FIELD:
			switch(c.first) {
			case FIELD_RA: return ra;
			case FIELD_DEC: return dec;
			default: ASSERT(0);
			}
		default: ASSERT(0);
		}
	}

	inline float color(const sdss_color &c) const
	{
		// this is a magnitude request
		if(c.first == c.second) { return magnitude(c.second); }

		// this is a true color request
		if(c.first >= 0) { return magnitude(c.first) - magnitude(c.second); }

		// this is a special color/field request
		switch(c.first)
		{
			case COLOR_P1S: return 0.910*u() - 0.495*g() - 0.415*r() - 1.280;
			case COLOR_P2S: return -0.249*u() +0.794*g() - 0.555*r() + 0.234;
			case COLOR_MR: return ml_r() - 5*log10(D/10);
			case COLOR_MV: return (ml_r() - 5*log10(D/10)) + 0.44*ml_gr() - 0.02;
		}
		ASSERT(0);
	}

	inline float magnitude(const int m) const
	{
		if(m < 10)
		{
			return mag[m];
		}
		else
		{
			// ML magnitude
			return ml_mag[m-10];
		}
	}
};
OSTREAM(const mobject &m);

void sdss_color::set(const std::string &color)
{
	name = color;

	if(first = mobject::field_id(color))
	{
		type = TYPE_FIELD;
	}
	else if(color.compare(0, 6, "sigma_") == 0)
	{
		type = TYPE_SIGMA;
		ASSERT(getColor(color.substr(6))) {
			std::cerr << "color = " << color << "\n";
		}
	}
	else
	{
		type = TYPE_COLORMAG;
		ASSERT(getColor(color)) {
			std::cerr << "color = " << color << "\n";
		}
	}
}

//////////////////////

/*
* 	Data structures for binner
*/

typedef std::pair<float, float> ribin;

struct less_S3 {
	bool operator()(const peyton::math::S3 &v1, const peyton::math::S3 &v2) const
	{
#if 0
		if(v1.x < v2.x) return true;
		if(v1.x > v2.x) return false;
	
		if(v1.y < v2.y) return true;
		if(v1.y > v2.y) return false;
		
		if(v1.z < v2.z) return true;
		return false;
#endif
		if(v1.z < v2.z) return true;
		if(v1.z > v2.z) return false;
	
		if(v1.x < v2.x) return true;
		if(v1.x > v2.x) return false;
		
		if(v1.y < v2.y) return true;
		return false;
	}
};

struct binned_run
{
public:
	struct pixel
	{
		std::set<int> stars; /// all stars in given pixel, by their uniqueIds
		double volume;

		pixel() : volume(0) {}
		~pixel() { }
	};

protected:
	bool sniffOnly;
public:
	double dx;		/// pixel lengthscale
	int run;		/// which run is this
	ribin colorbin;		/// color range included in the binned data
	double Dmin, Dmax;	/// volume limits
	
	typedef std::map<peyton::math::S3, pixel, less_S3> pixelmap;
	pixelmap pixels;  /// pixel data

	~binned_run() { }
	
	std::string filename;
	bool loaded;

	binned_run() : loaded(false), sniffOnly(false) {}
	binned_run(const std::string &fn) : filename(fn), loaded(false), sniffOnly(false) { load(true); }

	bool load(bool sniffOnly = false);

	bool unload()
	{
		if(!loaded) return false;
		pixels.clear();
		return true;
	}

	friend inline ibinarystream &operator >>(ibinarystream &in, binned_run &br);
};

inline OSTREAM(const binned_run &brs)
{
	int nstars = 0, nobs = 0;
	FOREACH(binned_run::pixelmap::const_iterator, brs.pixels)
	{
		const binned_run::pixel &p = (*i).second;
		nstars += p.stars.size();
		nobs += p.stars.size();
	}

	out << "# dx = " << brs.dx << "\n";
	out << "# runs = {" << 1 << " : " << brs.run << " }\n";
	out << "# colorbins = {"; out << " [" << brs.colorbin.first << ", " << brs.colorbin.second << ")"; out << " }\n";
	out << "# pixels = " << brs.pixels.size() << "\n";
	out << "# stars = " << nstars << "\n";
	out << "# observations = " << nobs << "\n";
	out << "#\n";
	out << "#    x     y     z  obs  volume  unq Nrun runs[Nrun]\n";

	FOREACH(binned_run::pixelmap::const_iterator, brs.pixels)
	{
		const peyton::math::S3 &idx = (*i).first;
		const binned_run::pixel &p = (*i).second;

		out << std::setw(8) << brs.dx*idx.x << std::setw(8) << brs.dx*idx.y << std::setw(8) << brs.dx*idx.z;
		out << std::setw(7) << p.stars.size() << " " << std::setw(11) << p.volume << std::setw(7) << p.stars.size();
		out << std::setw(4) << 1;
		out << std::setw(6) << brs.run;

		out << "\n";
	}

	return out;
}

inline obinarystream &operator <<(obinarystream &out, const binned_run::pixel &p)
{
	out << p.volume;
	out << p.stars.size();
	FOREACH(std::set<int>::const_iterator, p.stars) { out << *i; }
	return out;
}

inline ibinarystream &operator >>(ibinarystream &in, binned_run::pixel &p)
{
	size_t size;

	p.stars.clear();

	in >> p.volume;
	in >> size;
	FOR(0, size) { int v; in >> v; p.stars.insert(v); }
	return in;
}

inline obinarystream &operator <<(obinarystream &out, const binned_run &br)
{
	ASSERT(br.loaded);

	out << br.dx << br.run << br.colorbin << br.Dmin << br.Dmax;
	out << br.pixels.size();
	FOREACH(binned_run::pixelmap::const_iterator, br.pixels) { out << (*i).first << (*i).second; }
	return out;
}

inline ibinarystream &operator >>(ibinarystream &in, binned_run &br)
{
	in >> br.dx >> br.run >> br.colorbin >> br.Dmin >> br.Dmax;

	if(br.sniffOnly) { return in; }
	
	peyton::math::S3 k;
	size_t size;
	
	br.pixels.clear();
	
	in >> size;
	FOR(0, size)
	{
		in >> k;
		in >> br.pixels[k];
	}

	br.loaded = true;
	return in;
}

inline bool binned_run::load(bool sniff)
{
	if(loaded) return false;
	binary_input_or_die(in, filename);
	
	sniffOnly = sniff;
	in >> *this;
	sniffOnly = false;

	return true;
}

struct less_S2 {
	bool operator()(const peyton::math::S2 &v1, const peyton::math::S2 &v2) const
	{
		if(v1.x < v2.x) return true;
		if(v1.x > v2.x) return false;
	
		if(v1.y < v2.y) return true;
		return false;
	}
};

struct binned_runset
{
public:
	struct pixel
	{
		int N;	// number of star observations in pixel (non-unique)
		int uniqueN; 	// number of unique stars in pixel (for poisson error estimation)
		double volume;	// total volume observed in this pixel
		double uniqueVolume;	// unique volume observed in this pixel
		std::set<short> runs; // runs contributing to this pixel
		
		pixel() : N(0), uniqueN(0), volume(0), uniqueVolume(0) {}
	};
public:
	double dx;	/// pixel lengthscale
	typedef std::map<peyton::math::S3, pixel, less_S3> pixelmap;
	pixelmap pixels; /// pixel data
	std::set<int> runs;	/// runs included in this skymap
	std::set<ribin> colorbins;	/// color bins included in this skymap
};

OSTREAM(const binned_runset &brs);

#endif
