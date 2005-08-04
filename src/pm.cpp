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

#define NO_SDSS_STAR_CAT
 
#include "fitsloader.h"
#include "analysis.h"
#include "xcat.h"

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

#include <cmath>
#include <map>
#include <fstream>
#include <valarray>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <time.h>

#include <astro/useall.h>
using namespace std;

/////////////

//////////////////////////
//////////////////////////

/////

struct matcher
{
public:
	struct object
	{
		Radians ra, dec;
		union {
			int parent;		// the object this object is linked to, if positive (that is, parent > 0)
			int minusN;		// number of observations, if this object is the primary
		};
		int next;			// next object equivalent to this object, -1 if this is the last such object
						// - while matching: the next object in this bucket
	};

	struct index
	{
		int x, y, z;
		index(int x_, int y_, int z_) : x(x_), y(y_), z(z_) {}

		bool operator <(const index &a) const
		{
			if(a.x < x) return true;
			if(a.x > x) return false;

			if(a.y < y) return true;
			if(a.y > y) return false;

			if(a.z < z) return true;
			return false;
		}
		bool operator ==(const index &a) const
		{
			return x == a.x && y == a.y && z == a.z;
		}
	};
	struct hash_index { size_t operator()(const index &v) const { return (v.x << 22) + (v.y << 11) + v.z; } };

	struct bucket
	{
		int first, last;
		bucket() : first(-1), last(-1) {}
	};

public:
	typedef map<index, bucket> map_type;
	map_type buckets;
	Radians dmax, dphi;
	double r;
	DMMArray<object> objects;

	matcher(Radians dmax_, Radians dphi_) : dmax(dmax_), dphi(dphi_), r(1./dphi_) {}

public:
	void addtobucket(bucket &b, int oid)
	{
		if(b.last != -1) { objects[b.last].next = oid; b.last = oid; }
		else { b.first = b.last = oid; }
	}

	bool nearest(bucket &b, int oid)
	{
		object &o = objects[oid];

		int curid = b.first;
		while(curid >= 0)
		{
			object &cur = objects[curid];
			Radians d = distance(o.ra, o.dec, cur.ra, cur.dec);

			if(d < dmax)
			{
				o.parent = curid;
				return true;
			}
			curid = cur.next;
		}

		return false;
	}

	int match()
	{
		cout << "Opened memory map\n";
		cout << "Length : " << objects.size() << "\n";

		V3 p;
		map_type::iterator it;
		ticker tick(10000);
		int n = 0, nlinked = 0;
		FORj(oid, 0, objects.size())
		{
			object &o = objects[oid];
			p.celestial(r, o.ra, o.dec);
			o.parent = -1;
			o.next = -1;

			index idx0((int)p.x, (int)p.y, (int)p.z);
			index idx(idx0.x - 1, idx0.y - 1, idx0.z - 1);

			FORj(x, -1, 2)
			{
				FORj(y, -1, 2)
				{
					FORj(z, -1, 2)
					{
						it = buckets.find(idx);
						if(it != buckets.end())
						{
//							cout << "Found bucket !:" << &(*it).second << "\n";
							if(nearest((*it).second, oid))
							{
//								cout << "+\n";
								nlinked++;
								goto linked;
							}
						}
						idx.z++;
					}
					idx.z -= 3;
					idx.y++;
				}
				idx.y -= 3;
				idx.x++;
			}

		notlinked:
			addtobucket(buckets[idx0], oid);

		linked:
			// if it's linked, don't add it to bucket
			tick.tick();
			n++;
		}
		cout << "n = " << n << "\n";
		return nlinked;
	}

	void init(const std::string &lookupfn)
	{
		objects.open(lookupfn, "rw");
	}

	void makelookup(std::set<int> &runs, const std::string &lookupfn)
	{
#if 1
		sdss_star s;
		valarray<int> f1(5), f2(5);
		fits_set_streamer cat(runs, s, f1, f2);

		objects.create(lookupfn);

		cout << "sizeof(object): " << sizeof(object) << "\n";

		ticker tick(10000);
		int id = 0;
		while(cat.next())
		{
			RAD(s.ra); RAD(s.dec);
			object o = { s.ra, s.dec, -1, -1 };
			objects[id] = o;

			id++;
			tick.tick();
			if(id % 1000000 == 0) { objects.sync(); }
		}
		cout << "\n";

		objects.close();
#endif
	}
	
	void makerings()
	{
		FOR(0, objects.size())
		{
			object &o = objects[i];
			if(o.parent < 0)
			{
				o.next = -1;
				continue;
			}
			else
			{
				// insert this child at the start of it's parent's next chain
				object &p = objects[o.parent];
				o.next = p.next;
				p.next = i;
				p.parent--;	// we're storing the number of observations of this objects into p.parent, as a negative number
			}
		}
	}

	void makelinear()
	{
		binary_output_or_die(out, "match_groups.lut");

		// a map from linear rerun ID's to match_groups file
		DMMArray<int> idx;
		idx.create("match_index.dmm");

		int uniq = 0;
		ticker tick(10000);
		out << uniq;
		FOR(0, objects.size())
		{
			object &o = objects[i];
			if(o.parent >= 0) continue;

			int n = 0, N = -o.minusN;
			int grouploc = out.f.tellp();
			out << N;// << o.ra << o.dec;
			for(int k = i; k != -1; k = objects[k].next, n++)
			{
				out << k;
				idx[k] = grouploc;
			}
			uniq++;
			tick.tick();
			ASSERT(n == N);
		}
		out.f.seekp(0);
		out << uniq;
		cout << "Unique objects: " << uniq << "\n";
	}

	void stats(map<int, int> &hist)
	{
		ticker tick(10000);
		FOREACH(DMMArray<object>::iterator, objects)
		{
			tick.tick();
			if((*i).parent >= 0) continue;

			if(!hist.count(-(*i).parent))
			{
				hist[-(*i).parent] = 1;
			}
			else
			{
				hist[-(*i).parent]++;
			}
		}
	}
};

#include "phObjc.h"

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

#if 0
		cout << sdss_star::binary(pack) << "\n";
		cout << sdss_star::binary(RUN) << "\n";
		cout << sdss_star::binary(COL) << "\n";
		cout << sdss_star::binary(FIELD) << "\n";
		cout << sdss_star::binary(OBJID) << "\n";

		cout << run_ << " " << col_ << " " << field_ << " " << objid_ << "\n";
		cout << run() << " " << col() << " " << field() << " " << objid() << "\n";

		cout << sdss_star::binary(run()) << "\n";
		cout << sdss_star::binary(col()) << "\n";
		cout << sdss_star::binary(field()) << "\n";
		cout << sdss_star::binary(objid()) << "\n";
#endif

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
const long long packed_id::RUN = ((1 << 14) - 1) << 29;
const long long packed_id::COL = ((1 << 3) - 1) << 26;
const long long packed_id::FIELD = ((1 << 12) - 1) << 14;
const long long packed_id::OBJID = ((1 << 14) - 1) << 0;

OSTREAM(const packed_id &pid)
{
	return out << "[" << pid.run() << " " << pid.col() << " " << pid.field() << " " << pid.objid() << "]";
}

class selector
{
public:
	struct star
	{
		int id;				// linear rerun 137 catalog ID - index to FITS files
		packed_id sloanid;		// packed SDSS identification
		float colc;
		double ra, dec;
		float mag[5], magErr[5];
		float Ar;
	};
public:
	bool filter(sdss_star &s, const valarray<int> &flags, const valarray<int> &flags2)
	{
		//
		// Flags cuts on magnitudes
		//
		const int f1req = OBJECT1_BRIGHT | OBJECT1_SATUR;
		const int f2req = OBJECT2_DEBLENDED_AS_MOVING;

		FOR(0, 5)
		{
			if(flags[i]  & f1req) return false;
			if(flags2[i] & f2req) return false;

			// if there's a detection with sigma < 0.5 in at least one band, accept the star
			if(finite(s.mag[i]) && finite(s.magErr[i]) && s.magErr[i] < 0.5) return true;
		}

		return true;
	}

	void makelookup(std::set<int> &runs, const std::string &select, const std::string &selindex)
	{
		sdss_star s;
		valarray<int> f1(5), f2(5);
		fits_set_streamer cat(runs, s, f1, f2);

		DMMArray<star> sel;
		DMMArray<int> selidx;

		sel.create(select);
		selidx.create(selindex);

		ticker tick(10000);
		int id = -1;		// FITS catalog index
		star ss;
		while(cat.next())
		{
			id++; tick.tick();
			if(id % 1000000 == 0) { sel.sync(); selidx.sync(); }

			bool accept = filter(s, f1, f2);
			if(!accept) { selidx[id] = -1; continue; }

			ss.id = id;
			ss.ra = s.ra;
			ss.dec = s.dec;
			ss.sloanid = packed_id(s.run, s.col, s.field, s.objid);
			ss.colc = s.colc;
			FOR(0, 5) { ss.mag[i] = s.mag[i]; ss.magErr[i] = s.magErr[i]; }
			ss.Ar = s.Ar;
			
			selidx[id] = sel.size();
			sel.push_back(ss);
		}
		cout << "Total stars in FITS files: " << id+1 << "\n";
		cout << "Total stars accepted     : " << sel.size() << "\n";
		cout << "\n";
	}

	/*
		match_groups file has complete catalog indices in it. Create a new match_groups
		file which will point to objects in the _selected_ catalog.
	*/
	void reindex_match_groups(const string &matchgroups, const string &selectindex, const string &mgsel)
	{
		binary_input_or_die(in, matchgroups);
		binary_output_or_die(out, mgsel);

		// open selected catalog index
		DMMArray<int> selidx;
		selidx.open(selectindex, "r");

		//
		int N, id, size;

		in >> size;

		ticker tick(10000);
		unsigned int t = time(NULL), tstart = t;
		vector<int> obsv;
		int g = 0;
		out << g;
		FORj(k, 0, size)
		{
			in >> N;
			tick.tick();
			
			obsv.clear();
			FOR(0, N)
			{
				in >> id;
				int sid = selidx[id];
				if(sid == -1) continue;
				obsv.push_back(sid);
			}
			if(obsv.size() < 2) continue;		// ignore groups with less than two selected observations

			out << obsv.size();
			FOREACH(vector<int>::iterator, obsv) { out << *i; }
			g++;
		}
		out.f.seekp(0);
		out << g;
		cout << "Number of groups selected: " << g << "\n";
	}

	void groupinfo(int fitsID)
	{
		DMMArray<int> groupindex("match_index.dmm");	// mapping from fitsID to group file offset
		DMMArray<int> pmindex("pm_catalog_index.dmm");	// mapping from fitsID to star entry in PM source catalog
		DMMArray<star> catalog("pm_catalog.dmm");

		int unique;
		binary_input_or_die(groups, "match_groups.lut");
		groups >> unique;

		cout << groupindex.size() << " stars in catalog\n";
		if(fitsID < 0 || fitsID >= groupindex.size()) { cerr << "id out of catalog range\n"; }

		int N;
		int ptr = groupindex[fitsID];
		groups.f.seekg(ptr);
		groups >> N;
		cout << "Object " << fitsID << " has " << N << " observations:\n\n";

		float wt[N], f[N];
		float mean = 0, sigma = 0;

		int k;
		cout << "ra	dec	r      rsigma	identification\n";
		int n = 0;
		FOR(0, N)
		{
			groups >> k;
			k = pmindex[k];
			f[i] = -1;
			if(k == -1) { cout << "-- ignored\n"; continue; }

			star s = catalog[k];
			if(s.magErr[2] < 0.01)
			{
				s.magErr[2] = sqrt(1./2.*(sqr(0.01) + sqr(s.magErr[2])));
			}
			cout << s.ra << "\t" << s.dec << "\t" << s.mag[2] << "\t" << s.magErr[2] << "\t" << s.sloanid << "\n";

				
			phot::fluxfromluptitude(s.mag[2], s.magErr[2], 2, f[i], wt[i], 1e9);
			wt[i] = 1./sqr(wt[i]);
			
			mean += f[i]*wt[i];
			sigma += wt[i];
			n++;
		}
		mean /= sigma;
		float chi2 = 0;
		FOR(0, N)
		{
			if(f[i] < 0) continue;
			chi2 += sqr(f[i] - mean)*wt[i];
		}
		chi2 /= n;
		sigma = 1./sqrt(sigma);
		phot::luptitude(mean, sigma, 2, mean, sigma, 1e9);

		cout << "mean = " << mean << ", sigma = " << sigma << ", chi2/N = " << chi2 << "\n";
		
		cout << "\n";
	}
};

#if 1

class avgmag
{
public:
	struct mobject
	{
		int id;
		double ra, dec;
		int npos;
		float mag[5], wt[5], chi2[5];		
		float min[5], max[5];
		unsigned char N[5];
		struct { double a, b; double cov00, cov01, cov11, r; } vra, vdec;
		MJD t0, t1;
		float colcdist;
		float Ar;
	};
public:
	avgmag() {}

	void getfluxes(const sdss_star &s, float *f, float *w, float *mag)
	{
		float fErr;
		FOR(0, 5)
		{
			mag[i] = s.mag[i]+s.extinction(i);
			phot::fluxfromluptitude(mag[i], s.magErr[i], i, f[i], fErr, 1e9);
			w[i] = 1./sqr(fErr);
		}
	}
	

	#define MAXRUNS 10000
	#define MJD_OF_J2000_NOON    51544.5
	double runepochs[MAXRUNS];	// time of runs since J2000, in years
	RunGeometryDB geomDB;
	void initrunepochs(set<int> runs)
	{
		FOR(0, MAXRUNS) { runepochs[i] = 0; }
		FOREACH(set<int>::iterator, runs)
		{
			runepochs[*i] = (geomDB.getGeometry(*i).tstart - MJD_OF_J2000_NOON)
				/ (365.);
		}
	}
	MJD epochtomjd(double epochTime) { return epochTime*365. + MJD_OF_J2000_NOON; }
	
	float dist_from_edge(const float colc) const { return min(colc, 2048 - colc); }
	
	mobject process_observations(vector<selector::star> &obsv)
	{
		mobject m;
		int N = obsv.size();

		// initialize
		memset(&m, 0, sizeof(m));
		selector::star &primary = obsv.front();
		m.id = primary.id;
		m.npos = N;
		m.Ar = primary.Ar;
		m.ra = primary.ra;
		m.dec = primary.dec;
		m.colcdist = dist_from_edge(primary.colc);
		FOR(0, 5) { m.max[i] = m.min[i] = primary.mag[i]; }

		// magnitudes
		FOR(0, 5)
		{
			// sum things up
			float f[N], w[N], mag[N];
			int n = 0;
			FORj(k, 0, N)
			{
				selector::star &s = obsv[k];

				// throw out bad measurements
				if(!finite(s.mag[i]) || !finite(s.magErr[i]))
				{
					int kk = 2;
					continue;
				}

				// fix rerun 137 magnitude error bugs
				if(s.magErr[i] < 0.01)
				{
					s.magErr[i] = sqrt(1./2.*(sqr(0.01) + sqr(s.magErr[i])));
				}
				n++;

				// calculate the magnitude, flux and inverse variance of flux
				mag[k] = s.mag[i];
				phot::fluxfromluptitude(mag[k], s.magErr[i], i, f[k], w[k], 1e9);
				w[k] = 1./sqr(w[k]);

				m.mag[i] += f[k]*w[k];
				m.wt[i] += w[k];
	
				m.min[i] = std::min(m.min[i], mag[k]);
				m.max[i] = std::max(m.max[i], mag[k]);
			}
			m.N[i] = n;

			if(n > 0)	// averages make sense if we have at least one valid observation in this band
			{
				// average
				m.mag[i] /= m.wt[i];
	
				// chi squared
				n = 0;
				FORj(k, 0, N)
				{
					selector::star &s = obsv[k];

					// throw out bad measurements
					if(!finite(s.mag[i]) || !finite(s.magErr[i])) continue;
					n++;

					m.chi2[i] += sqr(f[k] - m.mag[i]) * w[k];
				}
				m.chi2[i] /= m.N[i];
				ASSERT(m.N[i] == n);

				// convert to luptitudes and stdevs
				m.wt[i] = 1./sqrt(m.wt[i]);
				phot::luptitude(m.mag[i], m.wt[i], i, m.mag[i], m.wt[i], 1e9);
			} else {
				// leave everything at 0 (to what we initialzed it)
			}
		}

		// velocities
		double ra[N], dec[N], t[N];	// ra and dec are in milli arcsecs, time in years since J2000
		m.t0 = 1e+10;			// earliest time of observation
		m.t1 = -1e-10;			// latest time of observation
		FORj(k, 0, N)
		{
			selector::star &s = obsv[k];

			t[k] = runepochs[s.sloanid.run()];
 			ra[k] = (s.ra - primary.ra)*3600.*1000.;
			dec[k] = (s.dec - primary.dec)*3600.*1000.;

			m.t0 = min(m.t0, epochtomjd(t[k]));
			m.t1 = max(m.t1, epochtomjd(t[k]));

			// distance from the edge
			m.colcdist = std::min(m.colcdist, dist_from_edge(s.colc));

			// doublecheck, just in case...
#if 0
			Radians dist = distance(rad(s.ra), rad(s.dec), rad(primary.ra), rad(primary.dec));
			if(dist >= rad(1./3600.))
			{
				cout << s.ra << s.dec << " !~ " << primary.ra << " " << primary.dec << " | " << deg(dist)*3600. << "\n";
			}
#endif
			ASSERT(distance(rad(s.ra), rad(s.dec), rad(primary.ra), rad(primary.dec)) < rad(1./3600.)*2.);
		}

		if(N > 1)	// velocities make sense only if N > 1
		{
			double sumsq;
			gsl_fit_linear(t, 1, ra, 1, N, &m.vra.b, &m.vra.a, &m.vra.cov00, &m.vra.cov01, &m.vra.cov11, &sumsq);
			gsl_fit_linear(t, 1, dec, 1, N, &m.vdec.b, &m.vdec.a, &m.vdec.cov00, &m.vdec.cov01, &m.vdec.cov11, &sumsq);
			
			if(N == 2)
			{
				// if we have only two points, errors are "perfect"
				m.vra.cov00 = m.vra.cov01 = m.vra.cov11 = m.vra.r = 0.;
				m.vdec.cov00 = m.vdec.cov01 = m.vdec.cov11 = m.vdec.r = 0.;
			}
			else
			{
				double sxx = gsl_stats_variance(t, 1, N);
	
				double syy = gsl_stats_variance(ra, 1, N);
				double sxy = gsl_stats_covariance(t, 1, ra, 1, N);
				m.vra.r = sxy / sqrt(sxx*syy);
				
				syy = gsl_stats_variance(dec, 1, N);
				sxy = gsl_stats_covariance(t, 1, dec, 1, N);
				m.vdec.r = sxy / sqrt(sxx*syy);
	
				//cout << m.vra.cov00 << " " << m.vra.cov11 << " " << m.vra.cov01 << " " << m.vra.r << "\n";
				//cout << m.vdec.cov00 << " " << m.vdec.cov11 << " " << m.vdec.cov01 << " " << m.vdec.r << "\n";
				//exit(0);
			}
		}

		return m;
	}

	void multiepoch(set<int> &runs, const string &matchgroupssel, const string &select)
	{
		binary_input_or_die(in, matchgroupssel);
		initrunepochs(runs);

		// open catalogs
		DMMArray<selector::star> sel;
		sel.setmaxwindows(40);
		sel.setwindowsize(1024*1024 / 5);
		sel.open(select, "r");

		//
		int N, id, size;
		//double ra, dec;
		vector<selector::star> obsv;

		in >> size;
		DMMArray<mobject> out;
		out.create("multiepoch.dmm");

		ticker tick(10000);
		unsigned int t = time(NULL), tstart = t;
		FORj(k, 0, size)
		{
			in >> N;// >> ra >> dec;
			obsv.clear();
			for(int i = 0; i != N; i++)
			{
				in >> id;
				selector::star &s = sel[id];
				obsv.push_back(s);
			}

			ASSERT(obsv.size() >= 2);

			mobject m = process_observations(obsv);
			out.push_back(m);

			if((out.size() % 500000) == 0) { out.sync(); }

			// user interface stuff
			tick.tick();
			if((k % 10000) == 0)
			{
				unsigned int tnow = time(NULL);
				cout << sel.winopenstat << " " << tnow - t << "s [ " << (tnow - tstart)/60 << "min ] ";
				if(k > 20000) { cout << "[ " << io::format("% 5.2f") << ((double(size) / double(k) -1.) * (tnow - tstart)/3600.) << "hrs left]";}
				cout << "\n";
				t = tnow;
				cout.flush();
			}
		}
	}

	void dumpascii(const string &matchgroupssel)
	{
		// output all of this as text, if the object is unique, to files sorted by declination
		valarray<ofstream *> files((ofstream *)NULL, 180);

		DMMArray<mobject> mags;
		mags.open("multiepoch.dmm", "r");

		ticker tick(10000);
		int n = 0, npm = 0;
		FOR(0, mags.size())
		{
			mobject m = mags[i];
			tick.tick();
			
			// this is a quick fix - these object should be kicked out upstream
			if(m.t0 == m.t1) continue;

			int decbin = (int)floor(m.dec);
			ASSERT(decbin > -90 && decbin < 90);

if(decbin == -1) {
			ofstream *&out = files[decbin + 90];
			if(out == NULL)
			{
				// open the file if it has not yet been opened
				string fn = io::format("matched-%+03d-%+03d.txt") << decbin << decbin+1;
				out = new ofstream(fn.c_str());
				(*out) << "# id ra dec coldist t0 t1 (6)\n";
				(*out) << "#\n";
				(*out) << "# Ar (7)\n";
				(*out) << "# Nobsv(u) <u> sigma(<u>) chi2(u) min(u) max(u) (13)\n";
				(*out) << "# Nobsv(g) <g> sigma(<g>) chi2(g) min(g) max(g) (19)\n";
				(*out) << "# Nobsv(r) <r> sigma(<r>) chi2(r) min(r) max(r) (25)\n";
				(*out) << "# Nobsv(i) <i> sigma(<i>) chi2(i) min(i) max(i) (31)\n";
				(*out) << "# Nobsv(z) <z> sigma(<z>) chi2(z) min(z) max(z) (37)\n";
				(*out) << "#\n";
				(*out) << "# N(pos) (38)\n";
				(*out) << "# vra  sigma(vra)  correlation(vra)  del_ra(J2000)  sigma(del_ra(J2000)) (43)\n";
				(*out) << "# vdec sigma(vdec) correlation(vdec) del_dec(J2000) sigma(del_dec(J2000)) (48)\n";
				(*out) << "#\n";
			}

			char buf[1000];
			// identification & photom. aux. info.
			//*out << i << "\t" << N << "\t" << deg(ra) << "\t" << deg(dec) << "\t" << m.colcdist;
			sprintf(buf, "%9i %12.8f % 12.8f %7.2f %12.6f %12.6f", m.id, m.ra, m.dec, m.colcdist, m.t0, m.t1);
//			sprintf(buf, "%8i %12.8f % 12.8f %7.2f", i, 0., 0., m.colcdist);
			*out << buf;

			//*out << "     " << m.Ar;
			sprintf(buf, "%5.3f", m.Ar);
			*out << "     " << buf;
			// photometry
			FORj(k, 0, 5)
			{
				if(finite(m.mag[k]))
				{
					//*out << "\t" << m.mag[k] << "\t" << m.wt[k] << "\t" << m.chi2[k] << "\t" << m.min[k] << "\t" << m.max[k];
					sprintf(buf, "%3i %5.2f %6.3f %8.3f %5.2f %5.2f", (int)m.N[k], m.mag[k], m.wt[k], m.chi2[k], m.min[k], m.max[k]);
				} else {
					//*out << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << 0;
					sprintf(buf, "%3i %5.2f %6.3f %8.3f %5.2f %5.2f", 0, 0., 0., 0., 0., 0.);
				}
				*out << "     " << buf;
			}

			// TODO: temporary fix which should also be implemented upstream
			// for objects observed in short period of time, the proper motion is nonsense
			ASSERT(m.t1 - m.t0 > 0);
			if(m.t1 - m.t0 < 30)
			{
				npm++;
//				cerr << "+" << (m.t1 - m.t0);
				memset(&m.vra, 0, sizeof(m.vra));
				memset(&m.vdec, 0, sizeof(m.vdec));
			}
			
			// proper motions
			sprintf(buf, "%3i", (int)m.npos);
			*out << "     " << buf;
			//*out << "\t\t" << m.vra.a << "\t" << m.vra.b << "\t" << r_ra << "\t" << m.vra.sumsq << "\t" << sqrt(m.vra.cov11) << "\t" << m.vra.cov01 << "\t" << sqrt(m.vra.cov00);
			sprintf(buf, "% 10.2f %10.2f % 5.3f  % 10.2f %10.2f", m.vra.a, sqrt(m.vra.cov11), m.vra.r, m.vra.b, sqrt(m.vra.cov00));
			*out << "     " << buf;
			//*out << "\t\t" << m.vdec.a << "\t" << m.vdec.b << "\t" << r_dec << "\t" << m.vdec.sumsq << "\t" << sqrt(m.vdec.cov11) << "\t" << m.vdec.cov01 << "\t" << sqrt(m.vdec.cov00);
			sprintf(buf, "% 10.2f %10.2f % 5.3f  % 10.2f %10.2f", m.vdec.a, sqrt(m.vdec.cov11), m.vdec.r, m.vdec.b, sqrt(m.vdec.cov00));
			*out << "     " << buf;
			
			*out << "\n";
			(*out).flush();
}
			n++;
		}
		tick.close();
		cout << "Number of unique stars dumped: " << n << "\n";
		cout << "Number of stars with too short pm interval: " << npm << "\n";

		FOR(0, files.size()) { if(files[i]) { delete files[i]; } }
	}
	
#if 0
	void starinfo(int id, const string &catalog, const string &lookupfn)
	{
		// open catalog
		sdss_star_cat cat(catalog);

		// open matches LUT
		DMMArray<matcher::object> lut;
 		lut.addfile(lookupfn);
		lut.open("r");

		// open magnitude averages
		DMMArray<mobject> mags;
		mags.setwindowsize(1*1024*1024);
		mags.setmaxwindows(200);
 		mags.addfiles("mags", cat.size());
		mags.open("r");

		// let's hit it rolling
		matcher::object *o = &lut[id];
		if(o->parent > 0)
		{
			// this is not a primary - redirect to primary
			cout << "This object is a child of " << o->parent << "\n";
			id = o->parent;
		}
		
		cout << "Primary: id = " << id << "\n";
		cout << "\t" << deg(o->ra) << "\t" << deg(o->dec) << "\n";

		cout << "Magnitude average record:\n";
		mobject &m = mags[id];
		cout << "\t" << m.N << "\n";
		char gri[] = {'g', 'r', 'i'};
		FOR(0, 3)
		{
			cout << "\t" << gri[i] << "   " << m.mag[i] << " " << m.wt[i] << " " << m.chi2[i] << " " << m.min[i] << " " << m.max[i] << "\n";
		}

		cout << "Number of children: " << -o->parent << "\n";
		cout << "Children:\n";
		int i = id;
		while(i >= 0)
		{
			sdss_star &s = cat[i];
			cout << "\t" << i << "\t" << deg(s.ra) << "\t" << deg(s.dec) << "\t" << s.r << "\t" << s.rErr << "\n";
			i = lut[i].next;
		}
	}
#endif
};
#endif

#include <sys/types.h>
#include <sys/stat.h>

void loadRuns(set<int> &runs, const std::string &runfile = "")
{
	string fn(runfile.size() ? runfile : "catalogs/runs.txt");
	cout << "Loading runs from: " << runfile << "\n";
	text_input_or_die (in, runfile);
	load(in, runs, 0);
	in_stream.close();
}

///////////////////

struct catline
{
	int run;
	int fitsIDstarts[7];
	
	bool operator <(const catline &c) const { return run < c.run; }
};

struct fitscat
{
	set<catline> cat;
	std::string fitsdir;

	fitscat();
	void open(const std::string &catfile)
	{
		text_input_or_die(cf, catfile);
		cf_stream >> fitsdir;

		catline cl;
		bind(cf, cl.run, 0, cl.fitsIDstarts[0], 1, cl.fitsIDstarts[1], 2, cl.fitsIDstarts[2], 3, cl.fitsIDstarts[3], 4, cl.fitsIDstarts[4], 5, cl.fitsIDstarts[5], 6, cl.fitsIDstarts[6], 7);
		
		while(cf.next())
		{
			cat.insert(cl);
		}
	}
	
	std::string fitsfile(int run, int col)
	{
		return io::format(fitsdir + "/calibObj-%06d-%d-star.fits.gz") << run << (col+1);
	}

	int firstID(int run, int col)
	{
		catline cl; cl.run = run;
		ASSERT(cat.count(cl) != 0);
		(*cat.find(cl)).fitsIDstarts[col];
	}
};

/*
	Creates a fitsID lookup table for FITS catalog.
*/
void importCatalog(const std::string &catfile, const std::string fitsdir)
{
	// figure out all of the runs
	set<int> runs;
#if 0
	string pattern = fitsdir + "/calibObj-*-1-star.fits.gz";
	string sspattern = fitsdir + "/calibObj-%d-1-star.fits.gz";
	dir inputfiles(pattern);
	FOREACH(inputfiles)
	{
		int run;
		sscanf((*i).c_str(), sspattern.c_str(), &run);
		runs.insert(run);
	}
#else
	loadRuns(runs, fitsdir + "/runs.txt");
#endif
	output_or_die(cf, catfile);
	cf << fitsdir << "\n";

	// make linear fitsID maps
	int id = 0;
	set<catline> cat;
	sdss_star s;
	valarray<int> flags1(5), flags2(5);
	FOREACHj(set<int>::iterator, run, runs)
	{
		catline c;
		c.run = *run;
		c.fitsIDstarts[0] = id;
		FOR(1, 7)
		{
			string file = io::format(fitsdir + "/calibObj-%06d-%d-star.fits.gz") << *run << i;
			cerr << "\t" << file << " ";
			
			fits_loader in(file, s, flags1, flags2, true);
			id += in.rows;

			c.fitsIDstarts[i] = id;
			cerr << id << "\n";
		}
		cat.insert(c);
		cout << "\n";
	}
	
	// store linear fitsID maps together with run map
	FOREACH(set<catline>::iterator, cat)
	{
		cf << io::format("%6d") << (*i).run;
		FORj(j, 0, 7)
		{
			cf << io::format("%10d") << (*i).fitsIDstarts[j];
		}
		cf << "\n";
	}
}

///////////////////

int main(int argc, char **argv)
{
try
{
	set<int> runs;
	loadRuns(runs, "/home/mjuric/projects/galaxy/workspace/catalogs/runs.txt");
	
	selector sel;
	matcher m(rad(1/3600.), rad(1./60.));
	avgmag am;

//	packed_id pid(745, 3, 11, 344);
//	return 0;

#if 1
	importCatalog("recalib.cat.txt", "data/schlegel");
	return 0;
#endif

#if 0
#if 1	// Index creation
	cerr << "Making lookup catalog\n";
	m.makelookup(runs, "catlookup.dmm");
#endif

	m.init("catlookup.dmm");

#if 1	// matching
	cerr << "Matching\n";
	int nlinked = m.match();
	cout << "nlinked = " << nlinked << "\n";
	
	cout << "Making rings...\n";
	m.makerings();
#endif

#if 1	// final output
	cerr << "Making groups and group index files.\n";
	m.makelinear();
#endif

#if 1	// statistics
	cout << "Statistics:\n";
	map<int, int> hist;
	m.stats(hist);
	int total = 0, stars = 0;
	FOREACH(hist)
	{
		cout << (*i).first << " " << (*i).second << "\n";
		total += (*i).first * (*i).second;
		stars += (*i).second;
	}
	cout << "total observations: " << total << "\n";
	cout << "total stars:        " << stars << "\n";
#endif

#endif

#if 0
	//sel.makelookup(runs, "pm_catalog.dmm", "pm_catalog_index.dmm");
	//sel.reindex_match_groups("match_groups.lut", "pm_catalog_index.dmm", "pm_groups.lut");
	sel.groupinfo(atoi(argv[1]));
#endif

#if 0
	//am.go("catalogs", "catlookup");
	//am.starinfo(302726, "catalogs", "catlookup");
	//m.groupinfo(830039);
	
	//am.multiepoch(runs, "pm_groups.lut", "pm_catalog.dmm");
	//am.dumpascii("pm_groups.lut");
#endif


} catch (peyton::exceptions::EAny &e)
{
	e.print();
} catch (...)
{
	cout << "Uncaught exception!\n";
}

}
