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
#include <astro/io/binarystream.h>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <cmath>
#include <map>
#include <fstream>
#include <valarray>
#include <numeric>
#include <set>
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
		ASSERT(col() == col_) {
			std::cerr << " .. " << run_ << " " << col_ << " " << field_ << " " << objid_ << "\n";
			std::cerr << " .. " << "[" << run() << " " << col() << " " << field() << " " << objid() << "]" << "\n";
		}
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
#define M_ML_G			(100 + 0)
#define M_ML_R			(100 + 1)
#define M_ML_I			(100 + 2)

struct sdss_color : public std::pair<int, int>
{
public:
	#define COLOR_P1S	-1
	#define COLOR_P2S	-2
	#define COLOR_MR	-3
	#define COLOR_MV	-4
	#define COLOR_DRIP	-5
	#define COLOR_RDM	-6	// distance modulus, r band

	std::string name;

	#define TYPE_COLORMAG	0
	#define TYPE_SIGMA	-1
	#define TYPE_FIELD	-2	// general accessor for any field of sdss_star
	#define TYPE_UNDEFINED	-1000
	int type;

public:
	sdss_color() : type(TYPE_UNDEFINED) {}
	operator bool() const { return type != TYPE_UNDEFINED; }

	static int bandIdx(const char c)
	{
		switch(c)
		{
			case 'g' : return 0;
			case 'r' : return 1;
			case 'i' : return 2;
			case 'z' : return 3;
			case 'u' : return 4;
			case 'A' : return 5;
			case '1' : return 10;
			case '2' : return 11;
			case '3' : return 12;
			case '4' : return 13;
			case '5' : return 14;
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
			case 5 : return 'A' ;
			case 10 : return '1'; // u uncorrected for ext.
			case 11 : return '2'; // g uncorrected for ext.
			case 12 : return '3'; // r uncorrected for ext.
			case 13 : return '4'; // i uncorrected for ext.
			case 14 : return '5'; // z uncorrected for ext.
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
		else if(color == "rdm") { first = COLOR_RDM; second = 0; }
		else if(color == "p2s") { first = COLOR_P2S; second = 0; }
		else if(color == "M_r") { first = COLOR_MR; second = 0; }
		else if(color == "M_V") { first = COLOR_MV; second = 0; }
		else if(color == "dRIp") { first = COLOR_DRIP; second = 0; }
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

	sdss_color(const std::string &color)
	: type(TYPE_UNDEFINED) 
	{
		if(color.size()) { set(color); }
	}

	sdss_color(const char *color)
	: type(TYPE_UNDEFINED) 
	{
		if(color != NULL) { set(color); }
	}
};
BOSTREAM2(const sdss_color &c);
BISTREAM2(sdss_color &c);

inline ISTREAM(sdss_color &c)
{
	std::string color;
	in >> color;
	c.set(color);
	return in;
}

inline OSTREAM(const sdss_color &c)
{
	return out << c.name;
}

/////////////////

// various identifiers an observation can have.
struct obsv_id
{
	int		fitsId;		// source catalog index
	packed_id	sloanId;	// SDSS identification
	int		uniqId;		// unique object to which this observation belongs (mobject DMM file index)
};

// observation photometry
struct obsv_mag
{
	int fitsId;			// source catalog index
	float mag[5];			// measured magnitudes (uncorrected for extinction)
	float magErr[5];		// note: magnitude ordering is grizu
};
OSTREAM(const obsv_mag &sm);

#define PROPER_MOTION_SUPPORT 0

// observation with photometry, astrometry and extinction information
struct observation : public obsv_mag
{
	double ra, dec;			// ra & dec (degrees)
	float Ar;			// extinction
	
#if PROPER_MOTION_SUPPORT
	float colc;			// frame column
#endif
};

#define MAG_MEAN		1
#define MAG_MEDIAN		2
#define MAG_CLIPPED_MEAN	3

class mobject;
void print_mobject(std::ostream &out, const mobject &m);

// supporting fields that are not physically stored on disk
// but are sometimes calculated during processing and passed
// along with a mobject to various routines
struct mobject_details
{
	peyton::Radians lon, lat;	// transformed longitude, latitude (see T_COORD in selector.cpp)
};

struct mobject
{
	int obs_offset;			// offset in obsv_mags array to where the obsv. of this object begin
	int n;				// total number of observations available

	float Ar;

	float mag[5];			// corrected averaged magnitudes
	float magErr[5];		// note: magnitude ordering is grizu
	float N[5];			// number of observations of this object, which were used to calculate mag[]
	short flags;			// how exactly we calculated mag[] and magErr[]

	float ml_mag[3];		// magnitudes deduced by max.likelihood fitting to locus (g, r, i)
	float D;			// geocentric distance, calculated through phot. paralax

	double ra, dec;			// object coordinates (degrees)

#if PROPER_MOTION_SUPPORT
	struct pm_t
	{
		// variability information
		float variance[5], chi2[5], magMin[5], magMax[5]; // grizu ordering

		// proper motion information
		struct { double a, b; double cov00, cov01, cov11, r; } vra, vdec;
		peyton::MJD t0, t1; // first/last epoch
		
		// minimum distance of any observation from the border of a frame
		float colcdist;
	} pm;
#endif

public:
	// accessors
	inline float gr() const { return mag[0] - mag[1]; }
	inline float ri() const { return mag[1] - mag[2]; }

	inline float u() const { return mag[4]; }
	inline float g() const { return mag[0]; }
	inline float r() const { return mag[1]; }
	inline float i() const { return mag[2]; }
	inline float z() const { return mag[3]; }

	inline float uErr() const { return magErr[4]; }
	inline float gErr() const { return magErr[0]; }
	inline float rErr() const { return magErr[1]; }
	inline float iErr() const { return magErr[2]; }
	inline float zErr() const { return magErr[3]; }
	
	inline float ml_gr() const { return ml_mag[0] - ml_mag[1]; }
	inline float ml_ri() const { return ml_mag[1] - ml_mag[2]; }
	
	inline float ml_g() const { return ml_mag[0]; }
	inline float ml_r() const { return ml_mag[1]; }
	inline float ml_i() const { return ml_mag[2]; }

	// non-const accessors
	inline float& u() { return mag[4]; }
	inline float& g() { return mag[0]; }
	inline float& r() { return mag[1]; }
	inline float& i() { return mag[2]; }
	inline float& z() { return mag[3]; }

	inline float& uErr() { return magErr[4]; }
	inline float& gErr() { return magErr[0]; }
	inline float& rErr() { return magErr[1]; }
	inline float& iErr() { return magErr[2]; }
	inline float& zErr() { return magErr[3]; }
	
	inline float& ml_g() { return ml_mag[0]; }
	inline float& ml_r() { return ml_mag[1]; }
	inline float& ml_i() { return ml_mag[2]; }

		
	inline float sigma(const sdss_color &c) const
	{
		if(c.first == c.second) { return magErr[c.second]; }

		if(c.first >= 0) { return std::sqrt(peyton::sqr(magErr[c.first]) + peyton::sqr(magErr[c.second])); }

		ASSERT(0);
	}

	enum {FIELD_UNKNOWN = 0, FIELD_RA, FIELD_DEC, FIELD_LOCUSDIST, FIELD_COLORLOGL, FIELD_DISTANCE, FIELD_LON, FIELD_LAT};
	static int field_id(const std::string &name)
	{
		if(name == "ra") { return FIELD_RA; }
		if(name == "dec") { return FIELD_DEC; }
		if(name == "lon") { return FIELD_LON; }
		if(name == "lat") { return FIELD_LAT; }
		if(name == "locdist") { return FIELD_LOCUSDIST; }
		if(name == "colorlogL") { return FIELD_COLORLOGL; }
		if(name == "distance") { return FIELD_DISTANCE; }
		return FIELD_UNKNOWN;
	}

	inline double field(const sdss_color &c, const mobject_details *mi = NULL) const
	{
		switch(c.type) {
		case TYPE_SIGMA:    return sigma(c);
		case TYPE_COLORMAG: return color(c);
		case TYPE_FIELD:
			switch(c.first) {
			case FIELD_RA: return ra;
			case FIELD_DEC: return dec;
			case FIELD_LON: ASSERT(mi != NULL); return peyton::math::deg(mi->lon);
			case FIELD_LAT: ASSERT(mi != NULL); return peyton::math::deg(mi->lat);
#if 1
			case FIELD_LOCUSDIST: return distance_from_locus(gr(), ri());
#else
			case FIELD_LOCUSDIST: {
				double d = distance_from_locus(gr(), ri());
				if(d == -1 && ml_mag[0] != 0.) {
//					print_mobject(std::cout, *this);
					std::cerr << "gr, ri, d(ml,real) = " << ml_gr() << " " << ml_ri()
						<< " " << sqrt(peyton::sqr(ml_gr() - gr()) + peyton::sqr(ml_ri() - ri()))
						<< "\n";
				}
				return d;
			}
#endif
			case FIELD_COLORLOGL: 
			{
				float lnL = 1;
				paralax_without_prior(ri(), gr(), magErr[0], magErr[1], magErr[2], &lnL);
				return lnL;
			}
			case FIELD_DISTANCE: return D;
			default: ASSERT(0) { std::cerr << "Unknown field = " << c.first << "\n"; }
			}
		default: ASSERT(0) { std::cerr << "Unknown color/field type = " << c.type << " [name = " << c.name << "]\n"; }
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
			case COLOR_RDM: return 5*log10(D/10);
			case COLOR_MR: return ml_r() - 5*log10(D/10);
			case COLOR_MV: return (ml_r() - 5*log10(D/10)) + 0.44*ml_gr() - 0.02;
			case COLOR_DRIP: 
				double RIp =    paralax_with_prior(ri(), gr(), magErr[0], magErr[1], magErr[2]);
				double RI = paralax_without_prior(ri(), gr(), magErr[0], magErr[1], magErr[2]);
				//std::cerr << RIp - RI << " " << RIp << " " << RI << "\n";
				//std::cout << ml_gr() << " " << ml_ri() << "\n";
				return RIp-RI;
			break;
		}
		ASSERT(0);
	}

	inline float magnitude(const int m) const
	{
		if(m < 5)
		{
			return mag[m];
		}
		else
		{
/*			if(m >= 10 && m <= 14)
				std::cerr << g() << " " << mag[((m-10)+4)%5] << " " << extinction(Ar, m-10) << " " << Ar << "\n";*/
			switch(m)
			{
				case 5: return Ar;
				case 10:
				case 11:
				case 12:
				case 13:
				case 14: return mag[((m-10)+4)%5] + extinction(Ar, m-10);
				default: return ml_mag[m-100]; // ML magnitude
			}
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

	//ASSERT(0) { std::cerr << "Unknown color/field " << color << "\n"; };
}

struct proc_obs_info
{
	std::map<float, zero_init<double> > lnLhist;
	int lf_failed, magerr_too_small;
	proc_obs_info() { lf_failed = magerr_too_small = 0; }
};
OSTREAM(const proc_obs_info &inf);
mobject process_observations(int obs_offset, double ra, double dec, float Ar, std::vector<std::pair<observation, obsv_id> > &obsv, proc_obs_info &inf);
void loadRuns(std::set<int> &runs, const std::string &runfile = "");
class catalog_streamer;
void makelookup(
	catalog_streamer &cat,		// input catalog of observations (FITS files, text files, etc.)
	const std::string &select, 	// name of the object catalog (indexed by uniq ID) (DMM file)
	const std::string &selindex,// name of the fitsID -> uniq ID map (DMM file)
	const std::string &runindex // name of the sloanID -> uniqID (DMM file)
	);
void make_run_index_offset_map(std::ostream &out, const std::string &runidxfn);

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
	FOREACH(brs.pixels)
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

	FOREACH(brs.pixels)
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
	FOREACH(p.stars) { out << *i; }
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
	FOREACH(br.pixels) { out << (*i).first << (*i).second; }
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
