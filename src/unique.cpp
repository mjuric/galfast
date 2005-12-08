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
#include "dm.h"

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

#define TRACE_TIME 1
#define RUN_CHECK 1
#define FILTER 1

#include "xcat.h"

//////////////////////////
//////////////////////////

/////

static int npass = 0;
struct matcher
{
public:
	struct object
	{
		Radians ra, dec;
		int nextObs;	// next observation of this object, -1 if this is the last observation
		int nextObj;	// while matching: the next object in this matcher::bucket, or -1 if this is the last object
						// if equal to -2, this matcher::object is a linked observation, and not a primary
		int run;
		//char col; float colc;
		#define LINKED_OBSERVATION -2
		
		bool operator<(const object &b) const { return run < b.run; }
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
		int nextObj;
		bucket() : nextObj(-1) {}
	};
	static bucket *sbucket;

	std::string printbucket(bucket *b)
	{
		std::string s;
		int curid = b->nextObj;
		while(curid != -1)
		{
			s += " ";
			s += str(curid);
			curid = objects[curid].nextObj;
		}
		return s;
	}
	
public:
//	typedef __gnu_cxx::hash_map<index, bucket, hash_index> map_type;
	typedef map<index, bucket> map_type;
	map_type buckets;
	Radians dmax, dphi;
	double r;
	DMMArray<object> objects;

	// dmax -> maximum distance of an observation to adjacent observation to be
	//			considered as an observation of the same object
	// dphi -> resolution of subdivision of the celestial sphere to bitmap lookup 
	//			table. Has to be greater than dmax.
	matcher(Radians dmax_, Radians dphi_) : dmax(dmax_), dphi(dphi_), r(1./dphi_) {}

public:
	void init(const std::string &lookupfn, const char *access = "rw")
	{
		objects.open(lookupfn, access);
		cerr << "Opened observation lookup table [" << lookupfn << "], " << access << " mode, " << objects.size() << " observations.\n";
	}

	void addtobucket(bucket &b, int oid)
	{
		// insert to the front of the bucket
		objects[oid].nextObj = b.nextObj;
		b.nextObj = oid;
		if(oid == 192581 || oid == 62) {
			cerr << "**bucketadd** oid = " << oid << "\n";
			cerr << "62(j,s) = " << objects[62].nextObj << " " << objects[62].nextObs << "\n";
		}
		if(oid == 62) {
			sbucket = &b;
			cerr << "Xbucket: " << printbucket(&b) << "\n";
		}
//		cout << "Bucket: " << &b << " " << b.first << " " << b.last << " " << k << "\n";
	}

	bool link_to_nearest(bucket &b, int oid)
	{
		if(&b == sbucket) { cerr << "oid = " << oid << " bucket: " << printbucket(&b) << "\n"; }
		object o = objects[oid];
//		cout << "nearest Bucket: " << &b << " " << b.first << " " << b.last << "\n";

		if(oid == 192581 || oid == 62) {
			cerr << "**** oid = " << oid << "\n";
			cerr << "62(j,s) = " << objects[62].nextObj << " " << objects[62].nextObs << "\n";
		}

		int closest_id_ptr = -1;
		Radians closest_dist = dmax;
		for(int curid = b.nextObj, id_ptr = -2; curid >= 0; id_ptr = curid, curid = objects[curid].nextObj)
		{
			object &cur = objects[curid];
			// special handling of rejects
			if(o.ra < 0)
			{
				if(cur.ra >= 0) { continue; }

				closest_dist = 0;
				closest_id_ptr = id_ptr;
				break;
			}

			#if RUN_CHECK
				if(cur.run == o.run) { continue; }
				#if !TRACE_TIME
				if(cur.nextObs != -1 && objects[cur.nextObs].run == o.run) { continue; }
				#endif
			#endif

			if(oid == 192581 || oid == 62) {
				cerr << "curid = " << curid << "\n";
				cerr << "cur = " << cur.run << " " << deg(cur.ra) << " " << deg(cur.dec) << "\n";
				cerr << "o   = " <<   o.run << " " << deg(  o.ra) << " " << deg(  o.dec) << "\n";
				cerr << "62(j,s) = " << objects[62].nextObj << " " << objects[62].nextObs << "\n";
			}

//			cout << "nearest cur: " << cur << " " << cur.parent << " " << cur.next << "\n";
			npass++;
			Radians d = distance(o.ra, o.dec, cur.ra, cur.dec);
			if(oid == 192581 || oid == 62) {
				cerr << " d = " << arcsec(d) << "  >?  dmax = " << arcsec(dmax) << "\n";
			}

			if(d < closest_dist)
			{
				closest_dist = d;
				closest_id_ptr = id_ptr;
				if(oid == 192581 || oid == 62) {
					cerr << " LINKED!\n";
				}
			}
		}

		if(closest_id_ptr != -1)
		{
			object &o = objects[oid];
			// insert this observation as the second from
			// the start of the observation chain
			int *root = closest_id_ptr == -2 ? &b.nextObj : &objects[closest_id_ptr].nextObj;

			if(oid == 192567) {
				cerr << "curid = " << *root << "\n";
				//cerr << "cur = " << cur.run << " " << deg(cur.ra) << " " << deg(cur.dec) << "\n";
				//cerr << "o   = " <<   o.run << " " << deg(  o.ra) << " " << deg(  o.dec) << "\n";
				cerr << "oid = " << oid << " bucket: " << printbucket(&b) << "\n";
			}
	
			#if TRACE_TIME
				int *curid = root;
				if(*curid != -1)
				{
					o.nextObj = objects[*curid].nextObj;
					objects[*curid].nextObj = LINKED_OBSERVATION;
				}
			#else
				int *curid = &objects[*root].nextObs;
				o.nextObj = LINKED_OBSERVATION;
			#endif

			o.nextObs = *curid;
			*curid = oid;

			if(oid == 192567) {
				cerr << "oid = " << oid << " bucket: " << printbucket(&b) << "\n";
			}

			if(o.ra >= 0) {
			#if 0
			cout << "Linking: " << deg(d)*3600 << " " 
				<< o.run << "," << (int)o.col << "," << o.colc
				<< " -> "
				<< cur.run << "," << (int)cur.col << "," << cur.colc
				<< "\n";
			#endif
			}
		}

		return closest_id_ptr != -1;
	}

	int match(const std::string &lookupfn)
	{
		init(lookupfn);

		V3 p;
		map_type::iterator it;
		ticker tick(10000);
		int n = 0, nlinked = 0;
//		ofstream fff("st.txt");
		FORj(oid, 0, objects.size())
		{
			object &o = objects[oid];
			p.celestial(r, o.ra, o.dec);
			o.nextObs = -1;
			o.nextObj = -1;
#if 0
			char buf[1000];
			sprintf(buf, "%.10f %.10f", deg(o.ra), deg(o.dec));
			fff << o.run << " " << (int)o.col << " " << buf << "\n";
			tick.tick();
			continue;
#endif
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
							if(link_to_nearest((*it).second, oid))
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
//		cout << "n = " << n << "\n";

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
		objects.close();
		return nlinked;
	}

	void makelookup(catalog_streamer &cat, const std::string &lookupfn)
	{
		sdss_star &s = cat.record();
		objects.create(lookupfn);

		cerr << "Pagesize:       " << MemoryMap::pagesize << "\n";
		cerr << "sizeof(object): " << sizeof(object) << "\n";

		ticker tick(10000);
		int ncat = 0;
		while(cat.next())
		{
			#if FILTER
			if(s.r >= 22) { s.ra = -1; s.dec = 0; }
			#endif

			object o = { rad(s.ra), rad(s.dec), -1, -1, s.run}; //, s.col, s.colc };
//			if(s.r >= 22) { cout << deg(o.ra) <<" " << deg(o.dec) << "\n"; }
			objects[s.id] = o;
			tick.tick();
			++ncat;
		}
		tick.close();

		cerr << "Catalog size:   " << ncat << " observations\n";
		cerr << "LUT file size:  " << ncat * sizeof(object) / (1024*1024) << "MB\n";
		cerr << "\n";

		objects.close();
	}

	void makelinear(const std::string &lookupfn, const std::string &matchgroupsfn, const std::string &matchindexfn)
	{
		init(lookupfn);
		binary_output_or_die(out, matchgroupsfn);
		//ofstream dump("dump.txt");

		// create new DMM set for storage
		DMMArray<int> idx;
		idx.create(matchindexfn);

		int uniq = 0;
		ticker tick(10000);
		out << uniq;
		vector<int> obsv;
		FOR(0, objects.size())
		{
			object obj = objects[i];
			if(obj.nextObj == LINKED_OBSERVATION) continue;

			obsv.clear();
			int grouploc = out.f.tellp();
			if(obj.ra >= 0) // do not store the indices for rejected objects
			{
				for(int k = i; k != -1; k = objects[k].nextObs)
				{
					obsv.push_back(k);
				}
				sort(obsv.begin(), obsv.end());
				int nobsv = obsv.size();

				out << nobsv << obj.ra << obj.dec;
				//dump << " " << nobsv;
				FOREACH(obsv)
				{
					out << *i;
					//dump << " " << *i;
					idx[*i] = grouploc;
				}
				//dump << "\n";

				// sanity check
				object o = objects[obsv.front()];
				std::map<int, object> oruns;
				FOREACH(obsv)
				{
					object cur = objects[*i];
					Radians d = distance(o.ra, o.dec, cur.ra, cur.dec);
					ASSERT(d < dmax) {
						cerr << "obj = " << cur.run << " " << deg(obj.ra) << " " << deg(obj.dec) << "\n";
						cerr << "cur = " << cur.run << " " << deg(cur.ra) << " " << deg(cur.dec) << "\n";
						cerr << "o   = " <<   o.run << " " << deg(  o.ra) << " " << deg(  o.dec) << "\n";
						cerr << " d = " << arcsec(d) << "  >  dmax = " << arcsec(dmax) << "\n";
					}
					ASSERT(!RUN_CHECK || oruns.count(cur.run) == 0) {
						//FOREACH(obsv) { cerr << objects[*i].run << " "; }
						cerr << ", " << obsv.size() << " observations \n";
						cerr << "obj = " << cur.run << " " << deg(obj.ra) << " " << deg(obj.dec) << "\n";
						o = oruns[cur.run];
						cerr << "cur = " << cur.run << " " << deg(cur.ra) << " " << deg(cur.dec) << "\n";
						cerr << "o   = " <<   o.run << " " << deg(  o.ra) << " " << deg(  o.dec) << "\n";
					}
					oruns[cur.run] = cur;

					#if TRACE_TIME
						o = cur;
					#endif
				}
			} else {
				int nobsv = 0;
				for(int k = i; k != -1; k = objects[k].nextObs)
				{
					nobsv++;
				}
				out << (int)0 << obj.ra << (double)nobsv;
			}
			uniq++;
			tick.tick();
		}
		out.f.seekp(0);
		out << uniq;

		tick.close();
		cerr << "Unique objects: " << uniq << "\n";
#if 0
		FOR(0, objects.size())
		{
			ASSERT(idx[i] != 0) {
				cerr << "Bug: " << i << " missing\n";
				object &o = objects[i];
				cerr << "\t" << deg(o.ra) << " " << deg(o.dec) << " " << o.nextObj << " " << o.nextObs << "\n";
			}
		}
#endif
		objects.close();
	}
#if 0
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
	#endif

	//
	// Statistics of match catalog
	//
	void stats(map<int, int> &hist, const std::string &matchgroupsfn)
	{
		binary_input_or_die(in, matchgroupsfn);
		hist.clear();

		ticker tick(10000);
		int uniq = 0;
		in >> uniq;
		vector<int> obsv;

		int nobsv; Radians ra, dec;
		FOR(0, uniq)
		{
			in >> nobsv >> ra >> dec;
			ASSERT(in.f.good());
			in.f.seekg(nobsv*sizeof(int), ios::cur);

			if(deg(ra) < 0)
			{
				ASSERT(!hist.count(0));
				hist[0] = (int)dec;
				continue;
			}

			if(!hist.count(nobsv))
			{
				hist[nobsv] = 1;
			}
			else
			{
				hist[nobsv]++;
			}

			tick.tick();
		}
	}
};
matcher::bucket *matcher::sbucket = NULL;

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
		FOREACH(cat.runs)
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

template<typename C>
class container_aux
{
	C &c;
	typedef typename C::value_type value_type;

	class aux
	{
		adapt_c<C> c;
	public:
		aux(C &c_) : c(c_) { c.clear(); }
		aux insert(const value_type &v) { c.push_back(v); return *this; }
		aux operator ,(const value_type &v) { return insert(v); }
	};

public:
	container_aux(C &c_) : c(c_) {}
	aux operator=(const value_type &v) { return aux(c).insert(v); }
};
template<typename C>
container_aux<C> container(C &c) { return container_aux<C>(c); }

template<typename C>
std::string join(const std::string &separator, const C& c)
{
	ostringstream ss;
	FOREACH(c)
	{
		if(i != c.begin()) { ss << separator; }
		ss << *i;
	}
	return ss.str();
}

int main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version);
	Options opts(
		"Identify observations of unique objects, and generate a lookup table. The code "
		"assumes the observations in the catalog are sorted in time-ascending order.",
		version, Authorship::majuric
	);

	double match_radius = 1;
	double map_resolution = 60;
	std::string
		inputCatalog = "/home/mjuric/projects/galaxy/workspace/catalogs/test.txt",
		stage = "stats",
		output_dir("."),
		lutfn("match_groups.lut"),
		indexfn("match_index.dmm"),
		cattype("fits"),
		tmp_dir("."),
		tmp_fn("catlookup.dmm");
	bool filter = false;

	std::vector<std::string> stages;
	container(stages) = "mklookup" , "match" , "mkgroups" , "stats";
	std::string allstages = join(", ", stages);
	stages.push_back("all");

	{ // setup arguments and options
		using namespace peyton::system::opt;

		opts.argument("inputCatalog", "Input catalog file", Option::optional);

		opts.option("stage", binding(stage), desc("Execute a stage of matching algorithm. Can be one of " + allstages + ". Special value 'all' causes all of the stages to be executed in sequence."));
		opts.option("radius", binding(match_radius), shortname('r'), desc("Match radius (arcsec)"));
		opts.option("map-resolution", binding(map_resolution), shortname('m'), desc("Size of bins into which the sky is divided while matching (arcsec). Has to be greater than 2*radius."));
		opts.option("type", binding(cattype), desc("Type of the input catalog (valid choices are 'fits' and 'text'"));
		opts.option('f', binding(filter), name("filter"), desc("Turn on filtering of input catalog to reduce the number of input stars"));
		opts.option("output-dir", binding(output_dir), desc("Directory where the group table and index will be stored"));
		opts.option("group-table-file", binding(lutfn), desc("Name of the group table file"));
		opts.option("group-index-file", binding(indexfn), desc("Name of the group index file"));
		opts.option("tmp-dir", binding(tmp_dir), desc("Directory where the temporary lookup file, used while matchin, will be stored."));
		opts.option("tmp-lookup-file", binding(tmp_fn), desc("Name of the temporary lookup file. This file can be deleted after the matching is complete."));
	}

	try {
		opts.parse(argc, argv);
		if(opts.found("inputCatalog")) { inputCatalog = opts["inputCatalog"]; }

		if(find(stages.begin(), stages.end(), stage) == stages.end())
			{ THROW(EOptions, "Argument to option --stage must be one of " + join(", ", stages) + "."); }
		if(!Filename(inputCatalog).exists())
			{ THROW(EOptions, "Input catalog file '" + inputCatalog + "' does not exist or is inaccessible"); }
		if(2*match_radius > map_resolution)
			{ THROW(EOptions, "2*match-radius must be less than map-resolution"); }
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}

#if 0
	cout << "stage = " << stage << "\n";
	cout << "match_radius = " << match_radius << "\n";
	cout << "map_resolution = " << map_resolution << "\n";
	cout << "inputCatalog = " << inputCatalog << "\n";
	cout << "output_dir = " << output_dir << "\n";
	cout << "lutfn = " << lutfn  << "\n";
	cout << "indexfn = " << indexfn << "\n";
	cout << "cattype = " << cattype << "\n";
	cout << "filter = " << filter << "\n";
	cout << "tmp_dir = " << tmp_dir << "\n";
	cout << "tmp_fn = " << tmp_fn << "\n";
#endif

	/////// Start your application code here
	matcher m(rad(match_radius/3600.), rad(map_resolution/3600.));
#if 1
	std::string tmp_fookup = tmp_dir + '/' + tmp_fn;
	std::string match_lut_path = output_dir + '/' + lutfn;
	std::string match_idx_path = output_dir + '/' + indexfn;

	if(stage == "mklookup" || stage == "all")
	{
		if(cattype == "fits")
		{
			#if 1
			cerr << "\n[1/3] Creating matching lookup table\n";
			std::set<int> runs;
			loadRuns(runs, inputCatalog);
			fits_set_streamer fitscat(runs);
			m.makelookup(fitscat, tmp_fookup);
			#endif
		} else {
			ASSERT(cattype == "fits");
		}
	}

	if(stage == "match" || stage == "all")
	{
		#if 1
		cerr << "\n[2/3] Matching observations\n";
		int nlinked = m.match(tmp_fookup);
		cerr << "nlinked = " << nlinked << "\n";
		cerr << "nbuckets = " << m.buckets.size() << "\n";
		#endif
	}

	if(stage == "mkgroups" || stage == "all")
	{
		#if 1
		cerr << "\n[3/3] Creating list of groups and match index\n";
		m.makelinear(tmp_fookup, match_lut_path, match_idx_path);
		#endif
	}

	if(stage == "stats" || stage == "all")
	{
		#if 1
		cerr << "\nCalculating statistics:\n";
		map<int, int> hist;
		m.stats(hist, match_lut_path);
		int total = 0, stars = 0;
		FOREACH(hist)
		{
			cerr << (*i).first << " " << (*i).second << "\n";
			if((*i).first != 0)
			{
				total += (*i).first * (*i).second;
				stars += (*i).second;
			} else {
				total += (*i).second;
			}
		}
		cerr << "total observations (counting rejected obsv.): " << total << "\n";
		cerr << "total stars (not counting rejected obsv.):    " << stars << "\n";
		#endif
	}

	cerr << "TRACE_TIME = " << TRACE_TIME << "\n";
	cerr << "npass = " << npass << "\n";
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
	cerr << "Uncaught exception!\n";
}

}
