#include <astro/constants.h>
#include <astro/coordinates.h>
#include <astro/system/log.h>
#include <astro/util.h>
#include <astro/sdss/rungeometry.h>
#include <astro/sdss/photometry.h>
#include <astro/math.h>
#include <astro/system/options.h>
#include <astro/io/gzstream/fstream.h>
#include <astro/io/compress.h>
#include <astro/io/magick.h>
#include <astro/system/config.h>
#include <astro/util/varexpand.h>
#include <astro/image/indexers.h>
#include <astro/io/fits.h>
#include <astro/io/format.h>
#include <astro/system/memorymap.h>
#include <astro/system/fs.h>

#include "textstream.h"
#include "analysis.h"

#include "vrml.h"
#include "gslcc.h"
#include "ximage.h"
#include "integerbin.h"

#include <astro/useall.h>

#include <cmath>
#include <map>
#include <fstream>
#include <string>
#include <set>
#include <list>
#include <valarray>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_fit.h>
#include <time.h>

using namespace std;

/////////////


#include "xcat.h"

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
//	typedef __gnu_cxx::hash_map<index, bucket, hash_index> map_type;
	typedef map<index, bucket> map_type;
	map_type buckets;
	Radians dmax, dphi;
	double r;
//	MemoryMapVector<object> objects;
	DMMArray<object> objects;

	matcher(Radians dmax_, Radians dphi_) : dmax(dmax_), dphi(dphi_), r(1./dphi_) {}

public:
	void addtobucket(bucket &b, int oid)
	{
		if(b.last != -1) { objects[b.last].next = oid; b.last = oid; }
		else { b.first = b.last = oid; }
//		cout << "Bucket: " << &b << " " << b.first << " " << b.last << " " << k << "\n";
	}

	bool nearest(bucket &b, int oid)
	{
		object &o = objects[oid];
//		cout << "nearest Bucket: " << &b << " " << b.first << " " << b.last << "\n";

		int curid = b.first;
		while(curid >= 0)
		{
			object &cur = objects[curid];
//			cout << "nearest cur: " << cur << " " << cur.parent << " " << cur.next << "\n";
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
//			cout << n << " : " << deg(o.ra) << " " << deg(o.dec) << " :: ";

			index idx0((int)p.x, (int)p.y, (int)p.z);
			index idx(idx0.x - 1, idx0.y - 1, idx0.z - 1);

//			cout << dphi << " " << r << " " << p << " " << idx0.x << " " << idx0.y << " " << idx0.z << "\n";

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
#if 0
			if(n % 1000 == 0) {
				cout << n << " :: Buckets: " << buckets.size() << "\n";
				//cout << n << " : " << deg(o.ra) << " " << deg(o.dec) << " :: " << idx0.x << " " << idx0.y << " " << idx0.z << "\n";
				FOREACHj(j, buckets)
				{
					cout << "\t" << (*j).first.x << " " << (*j).first.y << " " << (*j).first.z << "\n";
				}
			}
#endif
		}
		cout << "n = " << n << "\n";

#if 0
		FOREACHj(j, buckets)
		{
			int k = 1;
			bucket &b = (*j).second;
			object *cur = b.first;
			while(cur)
			{
				k++;
				cur = cur->next;
			}
			cout << "\t" << (*j).first.x << " " << (*j).first.y << " " << (*j).first.z << "| k = " << k << "\n";
		}
		FOREACH(objects)
		{
			object &o = *i;
			if(o.parent)
			{
				object &p = *o.parent;
				cout << deg(o.ra) << " " << deg(o.dec) << " == " << deg(p.ra) << " " << deg(p.dec) << " [ " << arcsec(o.ra - p.ra) << " " << arcsec(o.dec - p.dec) << "]\n";
			}
		}
#endif
		return nlinked;
	}

	void init(const std::string &lookupfn)
	{
		objects.open(lookupfn, "rw");
	}

	void makelookup(const std::string &catalog, const std::string &lookupfn)
	{
		sdss_star_cat cat(catalog);
		
		objects.create(lookupfn);

		cout << "Pagesize:       " << MemoryMap::pagesize << "\n";
		cout << "sizeof(object): " << sizeof(object) << "\n";
		cout << "Catalog size:   " << cat.size() << " objects\n";
		cout << "LUT file size:  " << cat.size() * sizeof(object) / (1024*1024) << "MB\n";

		ticker tick(10000);
		FOREACH(sdss_star_cat::iterator, cat)
		{
			sdss_star &s = *i;
			object o = { s.ra, s.dec, -1, -1 };
			objects[s.id] = o;
			tick.tick();
		}
		cout << "\n";

		objects.close();
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

		// create new DMM set for storage
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
			out << N << o.ra << o.dec;
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

	void groupinfo(int id)
	{
		int unique;
		binary_input_or_die(groups, "match_groups.lut");
		groups >> unique;

		// open the match index table
		DMMArray<int> idx;
		idx.open("match_index.dmm", "r");

		cout << idx.size() << " stars in catalog\n";
		if(id < 0 || id >= idx.size()) { cerr << "id out of catalog range\n"; }

		int N; double ra, dec;
		int ptr = idx[id];
		groups.f.seekg(ptr);
		groups >> N >> ra >> dec;
		cout << "Object " << id << " has " << N << " observations:\n\t";
		
		int k;
		FOR(0, N)
		{
			groups >> k;
			cout << k << " ";
		}
		cout << "\n";
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

#if 1

class avgmag
{
public:
	struct mobject
	{
		int id;
		float mag[5], wt[5], chi2[5];
		float min[5], max[5];
		unsigned char N[5];
		struct { double a, b; double cov00, cov01, cov11, sumsq; } vra, vdec;
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
	
#if 0
	void go(const string &catalog, const string &lookupfn)
	{
		int size = Filename(lookupfn).size() / sizeof(matcher::object);
		cout << "NStars: " << size << "\n";
		cout << "Size: " << size*sizeof(mobject) << "\n";

		DMMArray<mobject> mags;
		mags.setwindowsize(1*1024*1024);
		mags.setmaxwindows(200);
 		mags.addfiles("mags", size);
		mags.open("rw");

		binary_input_or_die(lookup, lookupfn);
		sdss_star_cat cat(catalog);

		matcher::object o;

		// make magnitude sums
		int n = 0;

		ticker tick(10000);
#if 1
		// First pass: sum the magnitudes, calculate total weights, find min/max
		// Streaming in: cat, lookup
		// Random access: mags
		float f[5], w[5], mag[5];
		FOREACHj(ss, cat)
		{
			sdss_star &s = *ss;
			lookup >> o;
			ASSERT(o.ra != s.ra);

			// calculate observed flux and error in flux
			getfluxes(s, f, w, mag);
			float coldist = min(s.colc, 2048 - s.colc);

			// do processing
			int id;
			if(o.parent < 0)
			{
				// primary observation
				mobject &m = mags[n];
				memset(&m, 0, sizeof(m));

				// initialize min/max to first value
				FOR(0, 5) { m.min[i] = m.min[i] = mag[i]; }
				m.colcdist = coldist;

				// leave this object as the current one
				id = n;
			}
			else
			{
				// a child - set a pointer back to parent
				mobject &m = mags[n];
				m.minusParent = -o.parent;

				// store the fluxes into mag array so that
				// we don't have to iterate through sdss_star_cat again
				// when calculating chi2
				FOR(0, 5) { m.mag[i] = f[i]; m.wt[i] = w[i]; }

				// and change the id to point to parent
				id = o.parent;
			}

			// calculating averages
			ASSERT(id <= n);
			mobject &m = mags[id];

			FOR(0, 5)
			{
				m.mag[i] += f[i]*w[i];
				m.wt[i] += w[i];

				m.min[i] = std::min(m.min[i], mag[i]);
				m.max[i] = std::max(m.max[i], mag[i]);
			}
			m.colcdist = std::min(m.colcdist, coldist);
			m.Ar = s.Ar;
			m.N++;

			n++;
			tick.tick();
		}
		tick.close();
		mags.sync();

		// Second pass:
		// calculate mean magnitude
		// calculate chi squared sums
		// **Random access: mags, cat
		FORj(n, 0, size)
		{
			float f[5], w[5];
			mobject &m0 = mags[n];
			if(m0.N > 0)
			{
				// primary - use the opportunity to calculate mean mag
				FOR(0, 5) { m0.mag[i] /= m0.wt[i]; }
				
				// get fluxes
				float mag[5];
				getfluxes(cat[n], f, w, mag);
			} else {
				// child
				FOR(0, 5) { f[i] = m0.mag[i]; w[i] = m0.wt[i]; }
			}

			// chi squared of flux
			mobject &m = m0.N > 0 ? m0 : mags[-m0.minusParent];
			FOR(0, 5) { m.chi2[i] += sqr(f[i] - m.mag[i]) * w[i]; }

			tick.tick();
		}
		mags.sync();
		
		// Third pass:
		// convert fluxes to luptitudes (aka asinh magnitudes)
		// convert chi2 to chi2/N
		// **Streaming: mags
		FOREACH(mags)
		{
			mobject &m = *i;
			if(m.N < 1) continue;
			FOR(0, 5)
			{
				m.wt[i] = sqrt(1/m.wt[i]);
				phot::luptitude(m.mag[i], m.wt[i], i, m.mag[i], m.wt[i], 1e9);
				m.chi2[i] /= m.N;
			}
		}
#endif

#if 1
		// Third pass:
		// output all of this as text, if the object is unique, to files sorted by declination
		valarray<ofstream *> files((ofstream *)NULL, 180);
		text_output_or_die(out, "mags.txt");
		DMMArray<matcher::object> lut;
 		lut.addfile(lookupfn);
		lut.open("r");
		tick.close();
		n = 0;
		FOR(0, size)
		{
			if(mags[i].N <= 1) continue;	// we want just multi-epoch observations (with N > 1)

			matcher::object &o = lut[i];
			int decbin = (int)floor(deg(o.dec));
			ASSERT(decbin > -90 && decbin < 90);

			ofstream *&out = files[decbin + 90];
			if(out == NULL)
			{
				// open the file if it has not yet been opened
				string fn = io::format("matched-%+03d-%+03d.txt") << decbin << decbin+1;
				out = new ofstream(fn.c_str());
				(*out) << "# id Nobsv ra dec Ar coldist   <g> sigma(<g>) chi2(g) min(g) max(g)   <r> sigma(<r>) chi2(r) min(r) max(r)   <i> sigma(<i>) chi2(i) min(i) max(i)" << "\n";
			}

			mobject &m = mags[i];
			*out << i << "\t" << m.N << "\t" << deg(o.ra) << "\t" << deg(o.dec) << "\t" << m.Ar << "\t" << m.colcdist;
			FORj(k, 0, 5)
			{
				*out << "\t" << m.mag[k] << "\t" << m.wt[k] << "\t" << m.chi2[k] << "\t" << m.min[k] << "\t" << m.max[k];
			}
			*out << "\n";

			n++;
			tick.tick();
		}
		tick.close();
		cout << "Number of unique tars dumped: " << n << "\n";
		
		FOR(0, files.size()) { if(files[i]) { delete files[i]; } }
#endif
	}
#endif

	#define MAXRUNS 10000
	#define MJD_OF_J2000_NOON    51544.5
	MJD runepochs[MAXRUNS];	// time of runs since J2000, in years
	RunGeometryDB geomDB;
	void initrunepochs(sdss_star_cat &cat)
	{
		FOR(0, MAXRUNS) { runepochs[i] = 0; }
		FOREACH(std::set<int>::iterator, cat.runs)
		{
			runepochs[*i] = (geomDB.getGeometry(*i).tstart - MJD_OF_J2000_NOON)
				/ (365.);
		}
	}
	MJD epochtomjd(MJD epochTime) { return epochTime + MJD_OF_J2000_NOON; }
	
	float dist_from_edge(const float colc) const { return min(colc, 2048 - colc); }
	
	mobject process_observations(vector<sdss_star> &obsv)
	{
		mobject m;
		int N = obsv.size();

		// initialize
		memset(&m, 0, sizeof(m));
		sdss_star &primary = obsv.front();
		m.id = primary.id;
		m.Ar = primary.Ar;
		m.colcdist = dist_from_edge(primary.colc);
		FOR(0, 5) { m.max[i] = m.min[i] = primary.mag[i]+primary.extinction(i); }

		// magnitudes
		FOR(0, 5)
		{
			// sum things up
			float f[N], w[N], mag[N];
			int n = 0;
			FORj(k, 0, N)
			{
				sdss_star &s = obsv[k];
				
				// throw out bad measurements
				if(!finite(s.mag[i]) || !finite(s.magErr[i])) continue;
				n++;

				// calculate the uncorrected magnitude, flux and inverse variance of flux
				mag[k] = s.mag[i]+s.extinction(i);
				phot::fluxfromluptitude(mag[k], s.magErr[i], i, f[k], w[k], 1e9);
				w[k] = 1./sqr(w[k]);

				m.mag[i] += f[k]*w[k];
				m.wt[i] += w[k];
	
				m.min[i] = std::min(m.min[i], mag[k]);
				m.max[i] = std::max(m.max[i], mag[k]);
			}
			m.N[i] = n;

			if(n > 0)	// averages make sense if we have at least one observation in this band
			{
				// average
				m.mag[i] /= m.wt[i];
	
				// chi squared
				n = 0;
				FORj(k, 0, N)
				{
					sdss_star &s = obsv[k];

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
			}
		}

		// velocities
		double ra[N], dec[N], t[N];	// ra and dec are in milli arcsecs, time in years since J2000
		m.t0 = 1e+10;			// earliest time of observation
		m.t1 = -1e-10;			// latest time of observation
		FORj(k, 0, N)
		{
			sdss_star &s = obsv[k];

			t[k] = runepochs[s.run];
 			ra[k] = deg(s.ra)*3600.*1000.;
			dec[k] = deg(s.dec)*3600.*1000.;

			m.t0 = min(m.t0, epochtomjd(t[k]));
			m.t1 = max(m.t1, epochtomjd(t[k]));

			// distance from the edge
			m.colcdist = std::min(m.colcdist, dist_from_edge(s.colc));

			// doublecheck, just in case...
			ASSERT(distance(s.ra, s.dec, primary.ra, primary.dec) < rad(1./3600.));
		}

		if(N > 1)	// velocities make sense only if N > 1
		{
			gsl_fit_linear(t, 1, ra, 1, N, &m.vra.b, &m.vra.a, &m.vra.cov00, &m.vra.cov01, &m.vra.cov11, &m.vra.sumsq);
			gsl_fit_linear(t, 1, dec, 1, N, &m.vdec.b, &m.vdec.a, &m.vdec.cov00, &m.vdec.cov01, &m.vdec.cov11, &m.vdec.sumsq);
			if(N == 2)
			{
			}
		}

		return m;
	}
	
	void multiepoch(const string &catalog, const string &matchgroups)
	{
		binary_input_or_die(in, matchgroups);
		sdss_star_cat cat(catalog);
		initrunepochs(cat);

		int N, id, size;
		double ra, dec;
		vector<sdss_star> obsv;

		in >> size;
		DMMArray<mobject> out;
		out.create("multiepoch.dmm");

		ticker tick(10000);
		unsigned int t = time(NULL), tstart = t;
		FORj(k, 0, size)
		{
			in >> N >> ra >> dec;
			obsv.clear();
			for(int i = 0; i != N; i++)
			{
				in >> id;
				obsv.push_back(cat[id]);
			}
			mobject m = process_observations(obsv);
			out[k] = m;
			tick.tick();
			if(!(k % 10000))
			{
				unsigned int tnow = time(NULL);
				cout << cat.winopenstat << " " << tnow - t << "s [ " << (tnow - tstart)/60 << "min ] ";
				if(k > 20000) { cout << "[ " << io::format("% 5.2f") << ((double(size) / double(k) -1.) * (tnow - tstart)/3600.) << "hrs left]";}
				cout << "\n";
				t = tnow;
			}
			if(k == 500000) { out.sync(); }
		}
	}

	void dumpascii(const string &matchgroups)
	{
		// output all of this as text, if the object is unique, to files sorted by declination
		valarray<ofstream *> files((ofstream *)NULL, 180);
		
		binary_input_or_die(in, matchgroups);
		int size;
		in >> size;

		DMMArray<mobject> mags;
		mags.open("multiepoch.dmm", "r");

		ticker tick(10000);
		int n = 0, i = -1, N;
		double ra, dec;
		FORj(k, 0, size)
		{
			tick.tick();

			in >> N >> ra >> dec; i++;
			in.f.ignore(N*sizeof(int));
			if(N == 1) continue;

			int decbin = (int)floor(deg(dec));
			ASSERT(decbin > -90 && decbin < 90);

			ofstream *&out = files[decbin + 90];
			if(out == NULL)
			{
				// open the file if it has not yet been opened
				string fn = io::format("matched-%+03d-%+03d.txt") << decbin << decbin+1;
				out = new ofstream(fn.c_str());
				(*out) << "# id Nobsv ra dec coldist\n";
				(*out) << "#\n";
				(*out) << "# Ar\n";
				(*out) << "# Nobsv(u) <u> sigma(<u>) chi2(u) min(u) max(u)\n";
				(*out) << "# Nobsv(g) <g> sigma(<g>) chi2(g) min(g) max(g)\n";
				(*out) << "# Nobsv(r) <r> sigma(<r>) chi2(r) min(r) max(r)\n";
				(*out) << "# Nobsv(i) <i> sigma(<i>) chi2(i) min(i) max(i)\n";
				(*out) << "# Nobsv(z) <z> sigma(<z>) chi2(z) min(z) max(z)\n";
				(*out) << "#\n";
				(*out) << "# vra  sigma(vra)  correlation(vra)  ra(J2000)  sigma(ra(J2000))\n";
				(*out) << "# vdec sigma(vdec) correlation(vdec) dec(J2000) sigma(dec(J2000))\n";
				(*out) << "\n";
			}

			mobject m = mags[i];

			char buf[1000];
			// identification & photom. aux. info.
			//*out << i << "\t" << N << "\t" << deg(ra) << "\t" << deg(dec) << "\t" << m.colcdist;
			sprintf(buf, "%8i %12.8f % 12.8f %7.2f", i, deg(ra), deg(dec), m.colcdist);
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

			// proper motions
			double r_ra = sqr(m.vra.cov01) / (m.vra.cov00*m.vra.cov11);
			double r_dec = sqr(m.vdec.cov01) / (m.vdec.cov00*m.vdec.cov11);
			if(N == 2) {
				r_ra = r_dec = 
				m.vra.cov00 = m.vra.cov11 = m.vra.cov01 = m.vra.sumsq =
				m.vdec.cov00 = m.vdec.cov11 = m.vdec.cov01 = m.vdec.sumsq = 0;
			}
			//*out << "\t\t" << m.vra.a << "\t" << m.vra.b << "\t" << r_ra << "\t" << m.vra.sumsq << "\t" << sqrt(m.vra.cov11) << "\t" << m.vra.cov01 << "\t" << sqrt(m.vra.cov00);
			sprintf(buf, "% 10.2f %10.2f %5.3f  %.8e %.8e", m.vra.a, sqrt(m.vra.cov11), r_ra, m.vra.b / (3600.*1000.), sqrt(m.vra.cov00) / (3600.*1000.));
			*out << "     " << buf;
			//*out << "\t\t" << m.vdec.a << "\t" << m.vdec.b << "\t" << r_dec << "\t" << m.vdec.sumsq << "\t" << sqrt(m.vdec.cov11) << "\t" << m.vdec.cov01 << "\t" << sqrt(m.vdec.cov00);
			sprintf(buf, "% 10.2f %10.2f %5.3f  %.8e %.8e", m.vdec.a, sqrt(m.vdec.cov11), r_ra, m.vdec.b / (3600.*1000.), sqrt(m.vdec.cov00) / (3600.*1000.));
			*out << "     " << buf;
			
			*out << "\n";
			(*out).flush();
			
			n++;
		}
		tick.close();
		cout << "Number of unique stars dumped: " << n << "\n";

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

int main(int argc, char **argv)
{
try
{
	matcher m(rad(1/3600.), rad(1./60.));
#if 1
//	m.makelookup("catalogs", "catlookup");

	m.init("catlookup");

//	int nlinked = m.match();
//	cout << "nlinked = " << nlinked << "\n";

//	cout << "Making rings...\n";
//	m.makerings();
	m.makelinear();

#if 0
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
#else
	avgmag am;
	//am.go("catalogs", "catlookup");
	//am.starinfo(302726, "catalogs", "catlookup");
	//am.multiepoch("catalogs", "match_groups.lut");
	//am.dumpascii("match_groups");
	m.groupinfo(830039);
#endif

} catch (peyton::exceptions::EAny &e)
{
	e.print();
} catch (...)
{
	cout << "Uncaught exception!\n";
}

}
