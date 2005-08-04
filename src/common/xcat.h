#ifndef __xcat_h
#define __xcat_h

#include "binarystream.h"

#include <astro/types.h>
#include <astro/util.h>
#include <astro/math.h>
#include <astro/system/log.h>
#include <astro/system/memorymap.h>
#include <astro/coordinates.h>

#include <map>
#include <vector>
#include <algorithm>

// abstract class for encapsulating different paralax relations
class sdss_star;
class paralax_relation
{
public:
	virtual bool operator()(sdss_star &s) = 0;
};

const double Rg = 8000;

inline float extinction(float Ar, int band)
{
	switch(band)
	{
		case 0: return Ar*5.155/2.751;
		case 1: return Ar*3.793/2.751;
		case 2: return Ar*1;
		case 3: return Ar*2.086/2.751;
		case 4: return Ar*1.479/2.751;
		default: ASSERT(0);
	};
}

struct sdss_star
{
	////////---------------- data

	int id;

	// things read in from the catalog
	short run, col, field, objid;
	peyton::Radians ra, dec;
	peyton::Radians l, b;				// note: l and b are calculated, not read in from files
	float colc;

	float Ar;					// r-band extinction [mag]

	// magnitudes (already extinction corrected)
	union
	{
		float mag[5];
		struct { float u, g, r, i, z; };
	};

	union
	{
		float magErr[5];
		struct { float uErr, gErr, rErr, iErr, zErr; };
	};

	// distance related stuff
	float Mr;			// absolute magnitude (r band), calculated from photometric paralax (see paralax_relation class)
	float ml_g, ml_r, ml_i;		// ML magnitudes determined by gr-ri locus fitting (see plx_gri_locus class)

	// this is calculated using Mr value above
	struct {
		float D;		// distance from Earth
		float x, y, z;		// rectangular earth centered coordinates
	} earth;
	struct {
		float rho, phi, z;	// galactocentric cylindrical coordinates
	} gc;

	// these are calculated by unique.x
	short Nobservations;		// number of observations of this star
	bool primaryobservation;	// is this star counted as the primary
	int obsptr;			// pointer to an avgmag object in a catalog generated using unique.x
					// if this is not the primary, the id of the primary

	/////////---------------  methods



	// colors
	inline float ug() const { return u - g; }
	inline float gr() const { return g - r; }
	inline float ri() const { return r - i; }
	inline float iz() const { return i - z; }
	inline float ur() const { return u - r; }
	inline float uz() const { return u - z; }
	inline float rz() const { return r - z; }
	inline float gi() const { return g - i; }

	inline float ml_gr() const { return ml_g - ml_r; }
	inline float ml_ri() const { return ml_r - ml_i; }

	inline float sqr(float x) const { return x*x; }
	inline float grErr() const { return sqrt(sqr(gErr)+sqr(rErr)); }
	inline float giErr() const { return sqrt(sqr(gErr)+sqr(iErr)); }
	inline float riErr() const { return sqrt(sqr(rErr)+sqr(iErr)); }
	inline float ugErr() const { return sqrt(sqr(uErr)+sqr(gErr)); }
	inline float izErr() const { return sqrt(sqr(iErr)+sqr(zErr)); }
	inline float urErr() const { return sqrt(sqr(uErr)+sqr(rErr)); }
	inline float uzErr() const { return sqrt(sqr(uErr)+sqr(zErr)); }
	inline float rzErr() const { return sqrt(sqr(rErr)+sqr(zErr)); }

	float extinction(int band) const
	{
		switch(band)
		{
			case 0: return Ar*5.155/2.751;
			case 1: return Ar*3.793/2.751;
			case 2: return Ar*1;
			case 3: return Ar*2.086/2.751;
			case 4: return Ar*1.479/2.751;
			default: ASSERT(0);
		};
	}

	static std::string binary(int x)
	{
		char s[33] = {0};
		FOR(0, 32) { s[i] = x & (0x1 << i) ? '1' : '0'; }
		return s;
	}
	
	static std::string binary(long long x)
	{
		char s[65] = {0};
		FOR(0, 64) { s[i] = x & (((long long)0x1) << i) ? '1' : '0'; }
		return s;
	}
	
	#define deg peyton::math::deg
	void print(std::ostream &out) const
	{
		out << "id  = " << id << "\n";
		out << "pos = [ " << run << " " << col << " " << field << " " << objid << "]\n";
		out << "equ = [ " << deg(ra) << " " << deg(dec) << "]\n";
		out << "gal = [ " << deg(l) << " " << deg(b) << "]\n";
		out << "extinction = ["; FORj(k, 0, 5) { out << " " << extinction(k); } out << "]\n";
		out << "mag = [ " << u << " " << g << " " << r << " " << i << " " << z << "]\n";
		out << "err = [ " << uErr << " " << gErr << " " << rErr << " " << iErr << " " << zErr << "]\n";
//		out << "flags1 = ["; FORj(k, 0, 5) { out << " " << binary(flags[k]); } out << "]\n";
//		out << "flags2 = ["; FORj(k, 0, 5) { out << " " << binary(flags2[k]); } out << "]\n";
		out << "colors [ug, gr, ri, iz, ur, uz, ri, gi] :\n";
		out << "        " << ug() << " " << gr() << " " << ri() << " " << iz() << " " << ur() << " " << uz() << " " << ri() << " " << gi() << "\n";
		out << "colErr  " << ugErr() << " " << grErr() << " " << riErr() << " " << izErr() << " " << urErr() << " " << uzErr() << " " << riErr() << " " << giErr() << "\n";
		out << "ml_mag =    [ " << ml_g << " "  << ml_r << " " << ml_i << " ]\n";
		out << "ml_colors = [ " << ml_gr() << " " << ml_ri() << " ]\n";
		out << "Mr  = " << Mr << "\n";
		out << "distance from Earth [x, y, z, D] :\n";
		out << "        " << earth.x << " " << earth.y << " " << earth.z << " " << earth.D << "\n";
		out << "distance from gal. center [rho, phi, z] :\n";
		out << "        " << gc.rho << " " << deg(gc.phi) << " " << gc.z << "\n";
	}
	#undef deg

	bool calculate_distances(paralax_relation &paralax, bool calc_paralax = true)
	{
		if(calc_paralax)
		{
			if(!paralax(*this)) return false;	// calculate the absolute magnitude and distances
		}

		float &D = earth.D;
		float x, y, z;

		// switching to rectangular geocentric coordinates
		float Re = D*cos(b);			// geocentric cylindrical coord. sys. distance
		earth.x = x = Re*cos(l);		// x distance from Earth
		earth.y = y = Re*sin(l);		// y distance from Earth
		earth.z = z = D*sin(b);			// z distance from Earth
		earth.x *= -1; earth.y *= -1;		// rotate the coordinate system to point _away_ from the galactic center (180deg rotation)

		// switching to galactocentric cylindrical
		x = Rg - x;				// switching to galactocentric coordinates
		y = -y;					// keeping the coordinate system righthanded
		gc.rho = sqrt(sqr(x)+sqr(y));		// galactocentric cylindric sys. distance
		gc.phi = atan2(y, x);			// polar angle
		gc.z = z;				// height

		return true;
	}

};

//
// Used for simple transformations of data, when importing stars from ZI's
// or Schlegel's files: apply extinction, calculate galactic coordinates,
// etc...
//
inline void preprocess_sdss_star(sdss_star &s)
{
	s.col--;			// make the columns to start from 0

	RAD(s.ra); RAD(s.dec);
	peyton::coordinates::equgal(s.ra, s.dec, s.l, s.b);

	FOR(0, 5) { s.mag[i] -= s.extinction(i); }
}

inline OSTREAM(const sdss_star &s)
{
	s.print(out);
	return out;
}

// StrictWeakOrdering implementation for comparing stars by run
struct sdss_star_less
{
	bool operator()(const sdss_star &a, const sdss_star &b)
	{
		if(a.run == b.run) {
			if(a.col == b.col)
			{
				if(a.field == b.field)
				{
					return a.objid < b.objid;
				}
				return a.field < b.field;
			}
			return a.col < b.col;
		}
		return a.run < b.run;
	}
};


#ifndef NO_SDSS_STAR_CAT

//
//
//	Binary star catalog classes and methods
//
//

#include <valarray>

template <typename K, typename V, typename L>
inline obinarystream &operator <<(obinarystream &out, const std::map<K, V, L> &m)
{
	out << m.size();
	typedef std::map<K, V, L> M;
	FOREACH(typename M::const_iterator, m) { out << (*i).first << (*i).second; }
	return out;
}

template <typename K, typename V, typename L>
inline ibinarystream &operator >>(ibinarystream &in, std::map<K, V, L> &m)
{
	size_t size; K k; V v;
	m.clear();

	in >> size;
	FOR(0, size) { in >> k >> v; m[k] = v; }
	return in;
}


template <typename V>
inline obinarystream &operator <<(obinarystream &out, const std::vector<V> &a)
{
	out << a.size();
	FOREACH(typename std::vector<V>::const_iterator, a) { out << *i; }
	return out;
}

template <typename V>
inline ibinarystream &operator >>(ibinarystream &in, std::vector<V> &a)
{
	size_t size; V v;

	in >> size;
	a.clear();
	a.reserve(size);

	FOR(0, size) { in >> v; a.push_back(v); }
	return in;
}

template <typename V>
inline obinarystream &operator <<(obinarystream &out, const std::valarray<V> &a)
{
	out << a.size();
	FOR(0, a.size()) { out << a[i]; }
	return out;
}

template <typename V>
inline ibinarystream &operator >>(ibinarystream &in, std::valarray<V> &a)
{
	size_t size; V v;

	in >> size;
	a.resize(size);

	FOR(0, size) { in >> a[i]; }
	return in;
}

typedef std::vector<sdss_star> starvector;

#include <astro/system/memorymap.h>
#include <set>

struct sdss_star_cat : peyton::system::DMMArray<sdss_star>
{
public:
	struct subset
	{
		typedef sdss_star_cat::iterator iterator;
		iterator beg, en;

	public:
		iterator &begin() { return beg; }
		iterator begin() const { return beg; }
		iterator &end() { return en; }
		iterator end() const { return en; }

		subset(const iterator &b = NULL, const iterator &e = NULL) : beg(b), en(e) {}

		sdss_star &operator[](int i) { return beg[i]; }
		sdss_star operator[](int i) const { return beg[i]; }

		int size() { return en - beg; }
	};
public:
	std::string dir;
	std::set<int> runs;

	std::map<int, subset> runsubsets;
public:
	sdss_star_cat(std::string dir_);

	subset getrun(int run);
};

#endif

#endif
