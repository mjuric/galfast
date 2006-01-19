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

#include "config.h"
//#ifdef HAVE_LIBCCFITS

#define NO_SDSS_STAR_CAT

#include "dm.h" 
#ifdef HAVE_LIBCCFITS
#include "fitsloader.h"
#endif
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

#include <astro/useall.h>
using namespace std;

/////////////

//////////////////////////
//////////////////////////

#include "../phObjc.h"

////////////////////////////////////////////////////////////////

plx_gri_locus paralax;

struct filter_info
{
	int no_r_band, r_too_dim, no_two_bands, run_rejected;
	filter_info() : no_r_band(0), r_too_dim(0), no_two_bands(0), run_rejected(0) {}
};
static filter_info definf;
OSTREAM(const filter_info &inf)
{
	out << "# No r-band observation                : " << inf.no_r_band << " observations.\n";
	out << "# No observation in two bands          : " << inf.no_two_bands << " observations.\n";
	out << "# r-band observation dimmer than limit : " << inf.r_too_dim << " observations.\n";
	out << "# Run explicitly rejected              : " << inf.run_rejected << " observations.\n";
	return out;
}

//
// The first filter, applied when importing observations from source catalog
// to DMM arrays.
//
bool filter(sdss_star &s, const valarray<int> &flags, const valarray<int> &flags2, filter_info &inf = definf)
{
	//
	// Flags cuts on magnitudes
	//
	const int f1req = OBJECT1_BRIGHT | OBJECT1_SATUR;
	const int f2req = OBJECT2_DEBLENDED_AS_MOVING;

	int nValidBands = 0;
	FOR(1, 4)	// check flags for g, r, i magnitudes
	{
		if(flags[i]  & f1req) return false;
		if(flags2[i] & f2req) return false;

		if(finite(s.mag[i]) && finite(s.magErr[i])) { nValidBands++; }
	}

	//
	// Apply magnitude cuts
	//
	if(!finite(s.mag[2])) { inf.no_r_band++; return false; }	// must have r magnitude
	if(nValidBands < 2) { inf.no_two_bands++; return false; }	// must have at least one color
	if(s.r >= 22) { inf.r_too_dim++; return false; }			// r magnitude cutoff

	// Other case-by-case rejections
	//if(s.run == 5224) { inf.run_rejected++; return false; }
	
	return true;
}

//
// Converts the source observation catalog to easily searchable DMM arrays.
// Also filters the source catalog, according to the given filter.
//
// Three arrays are created: select, selindex, runindex
//    select is the DMM array of selected obsv. from the source catalog
//    selindex is the DMM array which maps the indices of objects in orig. catalog. to indices in selected catalog
//    runindex is the DMM array which maps the indices of selected obsv. catalog, to all other indices (source cat, unique object, etc.)
//
void makelookup(
	catalog_streamer &cat,		// input catalog of observations (FITS files, text files, etc.)

	const std::string &select, 	// name of the object catalog (indexed by uniq ID) (DMM file)
	const std::string &selindex,// name of the fitsID -> uniq ID map (DMM file)
	const std::string &runindex // name of the sloanID -> uniqID (DMM file)
	)
{
//#ifdef HAVE_LIBCCFITS
	// setup SDSS FITS file catalog streamer
	sdss_star &s = cat.record();
	valarray<int> 
		&f1 = cat.flags1(),
		&f2 = cat.flags2();	// photometry flags (FLAGS and FLAGS2 FITS fields)

	// open arrays for catalog and index
	DMMArray<observation> sel;	// catalog of selected observations (name of the index: selIndex)
	DMMArray<int> selidx;	// fitsID -> selIndex map. May be sparse, because not all FITS records will be accepted
	DMMArray<obsv_id> runidx;// selIndex -> uniqID map - just creating a placeholder array here

	sel.create(select);
	selidx.create(selindex);
	runidx.create(runindex);

	ticker tick(10000);
	observation ss; obsv_id sid;
	sid.uniqId = -1;
	filter_info inf;

	while(cat.next())
	{

		tick.tick();
		if(s.id % 1000000 == 0) { sel.sync(); selidx.sync(); runidx.sync(); }

		bool accept = filter(s, f1, f2, inf);

		if(!accept) { selidx[s.id] = -1; continue; }

		// record indices
		ss.fitsId = sid.fitsId = s.id;
		sid.sloanId = packed_id(s.run, s.col, s.field, s.objid);

		// store the info we're interested in
		ss.ra = s.ra;
		ss.dec = s.dec;
		// a hack - add u & z magnitudes to the end of the array, to keep
		// backwards compatibility, and yet have them there too
		// so the ordering is grizu
		FOR(0, 5) { ss.mag[i] = s.mag[(i+1) % 5]; ss.magErr[i] = s.magErr[(i+1) % 5]; }		
		ss.Ar = s.Ar;

		// store
		selidx[sid.fitsId] = sel.size();
		sel.push_back(ss);
		runidx.push_back(sid);
	}
	tick.close();

	cout << "\n";

	cout << "Lookup table: " << select << " [ " << sel.size() << " entries ]\n";
	cout << "Cat->Lut intex: " << selindex << " [ " << selidx.size() << " entries ]\n";
	cout << "Lut->other index: " << " [ " << runidx.size() << " entries ]\n";

	cout << "\n" << inf << "\n";
	
	cout << "Total observations in FITS files: " << s.id+1 << "\n";
	cout << "Total observations accepted     : " << sel.size() << "\n";
	cout << "\n";
	
	ASSERT(sel.size() == runidx.size());
	ASSERT(selidx.size() == s.id+1);

	cout.flush();
//#else
//	cerr << "CCfits support not compiled\n";
//	abort();
//#endif
}

struct proc_obs_info
{
	std::map<float, zero_init<double> > lnLhist;
	int lf_failed, magerr_too_small;
	proc_obs_info() { lf_failed = magerr_too_small = 0; }
};

OSTREAM(const proc_obs_info &inf)
{
	out << "# Locus fit failed                   : " << inf.lf_failed << " objects.\n";
	out << "# Magnitude error too small          : " << inf.magerr_too_small << " measurements.\n";
	out << "# Likelihood histogram:\n";
	FOREACH(inf.lnLhist) { out << (*i).first << "\t" << (*i).second << "\n"; }
	return out;
}

mobject process_observations(int obs_offset, double ra, double dec, float Ar, std::vector<obsv_mag> &obsv, proc_obs_info &inf)
{
	static plx_gri_locus plx;

	int N = obsv.size();

	// initialize
	mobject m;
	memset(&m, 0, sizeof(m));
	m.obs_offset = obs_offset;
	m.ra = ra;
	m.dec = dec;
	m.Ar = Ar;
	m.n = obsv.size();

	m.flags = MAG_MEAN;

	// for each magnitude...
	FOR(0, 5)
	{
		// sum things up
		float f[N], w[N], mag[N];
		int n = 0;
		FORj(k, 0, N) // for each observation...
		{
			obsv_mag &s = obsv[k];

			// throw out bad measurements
			if(!finite(s.mag[i]) || !finite(s.magErr[i]))
			{
				continue;
			}

			// fix rerun 137 magnitude error bugs
			if(s.magErr[i] < 0.01)
			{
				s.magErr[i] = sqrt(1./2.*(sqr(0.01) + sqr(s.magErr[i])));
				inf.magerr_too_small++;
			}
			n++;

			// calculate the magnitude, flux and inverse variance of flux
			mag[k] = s.mag[i];
			phot::fluxfromluptitude(mag[k], s.magErr[i], i, f[k], w[k], 1e9);
			w[k] = 1./sqr(w[k]);

			m.mag[i] += f[k]*w[k];
			m.magErr[i] += w[k];
		}
		m.N[i] = n;

		if(n > 0)	// averages make sense if we have at least one valid observation in this band
		{
			// average
			m.mag[i] /= m.magErr[i];

			// convert to luptitudes and stdevs
			m.magErr[i] = 1./sqrt(m.magErr[i]);
			phot::luptitude(m.mag[i], m.magErr[i], i, m.mag[i], m.magErr[i], 1e9);

			// correct for extinction
			m.mag[i] -= extinction(Ar, (i+1) % 5); // note the hack to account for grizu ordering
		} else {
			// mark magnitude as unknown, by setting sigma to infinity
			static double zero = 0.;
			m.magErr[i] = 1./zero;
			cerr << "+\n";
		}
	}

	// calculate distance, from paralax relation
#if 1
	float lnL;
	float RI = plx.ml_r_band(m.ri(), m.gr(), m.gErr(), m.rErr(), m.iErr(), &lnL);
	if(RI != -1)
	{
		float GR = plx.gr(RI);
		plx.ml_magnitudes(m.ml_g(), m.ml_r(), m.ml_i(), m.g(), m.r(), m.i(), m.gErr(), m.rErr(), m.iErr(), RI, GR);

		double Mr = plx.Mr(RI);
		m.D = stardist::D(m.ml_r(), Mr);

		ASSERT(fabs(m.ml_ri() - RI) < 0.001);
		ASSERT(fabs(m.ml_gr() - GR) < 0.001);

		//cout << lnL << " " << 0.1*int(-lnL / 0.1) << "\n";
		inf.lnLhist[0.1*int(-lnL / 0.1)]++;
	} else {
		inf.lf_failed++;
	}
#else
	sdss_star s;
	s.ra = m.ra*ctn::d2r; s.dec = m.dec*ctn::d2r;
	coordinates::equgal(s.ra, s.dec, s.l, s.b);
	FOR(0, 3) { s.mag[i+1] = m.mag[i]; s.magErr[i+1] = m.magErr[i]; } // grizu ordering hack (sdss_star has ugriz ordering)

	if(paralax(s))	// calculate the absolute magnitude, ML colors, and distances
	{
		m.D = s.earth.D;
		m.ml_mag[0] = s.ml_g;
		m.ml_mag[1] = s.ml_r;
		m.ml_mag[2] = s.ml_i;
	} else {
		inf.lf_failed++;
	}
#endif
	//cout << m << "\t" << lnL << "\n";
	//exit(-1);
	return m;
}

// defined in binner_tng.cpp
bool calculate_cartesian(V3 &earth, const mobject &m);

//#define NOMOD 1
#define NOMOD 0
void make_object_catalog(
	const std::string &uniqObjectCat, // "dm_unique_stars.dmm" [output]
	const std::string &uniqObsvCat, // "dm_starmags.dmm" [output]

	const std::string &matchGroups, // "match_groups.lut" (produced by unique.x)

	const std::string &select,   // name of the temporary object catalog (indexed by uniq ID) (DMM file) [in]
	const std::string &selindex, // name of the fitsID -> uniq ID map (DMM file) [in]
	const std::string &runindex, // name of the sloanID -> uniqID (DMM file) [in]
	
	const std::string &stagesummary_fn // name of text file for output summary
)
{
	binary_input_or_die(in, matchGroups);

	// open lookup catalogs (input)
	DMMArray<observation> sel; sel.setmaxwindows(40); sel.setwindowsize(1024*1024 / 5);
	DMMArray<int> selidx; selidx.setmaxwindows(40); selidx.setwindowsize(1024*1024 / 5);
	DMMArray<obsv_id> runidx; runidx.setmaxwindows(40); runidx.setwindowsize(1024*1024 / 5);

	sel.open(select, "r");
	selidx.open(selindex, "r");
	runidx.open(runindex, NOMOD ? "r" : "rw");
	cout << "Lookup table: " << select << " [ " << sel.size() << " entries ]\n";
	cout << "Cat->Lut intex: " << selindex << " [ " << selidx.size() << " entries ]\n";
	cout << "Lut->other index: " << " [ " << runidx.size() << " entries ]\n";

	int N, id, size;
	vector<obsv_mag> obsv;

	// open object/observation catalogs (output)
	DMMArray<mobject> out;
	DMMArray<obsv_mag> obsv_mags;
#if !NOMOD
	out.create(uniqObjectCat);
	obsv_mags.create(uniqObsvCat);
#else
	out.open(uniqObjectCat, "r");
	obsv_mags.open(uniqObsvCat, "r");
#endif
	ticker tick(10000);
	unsigned int t = time(NULL), tstart = t;
	in >> size;
	// stream through all observations, grouped by unique objects
	proc_obs_info inf;
	FORj(k, 0, size)
	{
		int uniqId = out.size();

		// load the whole group of observations belonging to an object
		double ra, dec; float Ar;
		in >> N >> ra >> dec;
		obsv.clear();
		for(int i = 0; i != N; i++)
		{
			in >> id; // source catalog id of the observation

			// check if this observation was selected as valid
			int sid = selidx[id];
			if(sid == -1) continue;
#if !NOMOD
			// load this observation
			observation &o = sel[sid];
			obsv.push_back(o);

			if(obsv.size() == 1) { ra = o.ra; dec = o.dec; Ar = o.Ar; }
#endif

			// store ID of the unique object together with this observation
#if !NOMOD
			obsv_id &oid = runidx[sid];
			oid.uniqId = uniqId;

			ASSERT(oid.fitsId == o.fitsId);
#endif
		}

#if !NOMOD
		// process the object - e.g., average magnitudes, calculate proper motions, etc.,
		// and then store it to the object catalog
		if(obsv.size() != 0)		// ignore groups with zero selected observations
		{
			// process and store
			mobject m = process_observations(obsv_mags.size(), ra, dec, Ar, obsv, inf);

			// store observation lookup table (holds magnitudes only, for now)
			FOREACH(obsv) { obsv_mags.push_back((obsv_mag&)*i); }

			out.push_back(m);
		}
#endif
		if((out.size() % 500000) == 0) { out.sync(); }

		// user interface stuff
		tick.tick();
		#if 0
		if((k % 10000) == 0)
		{
			unsigned int tnow = time(NULL);
			cout << sel.winopenstat << " " << tnow - t << "s [ " << (tnow - tstart)/60 << "min ] ";
			if(k > 20000) { cout << "[ " << io::format("% 5.2f") << ((double(size) / double(k) -1.) * (tnow - tstart)/3600.) << "hrs left]";}
			cout << "\n";
			t = tnow;
			cout.flush();
		}
		#endif
	}
	
	// convert histogram to cumulative
	double total = 0;
	FOREACH(inf.lnLhist) { total += (*i).second; }
	double cur = 0;
	FOREACH(inf.lnLhist) { cur = ((*i).second += cur); (*i).second /= total; }

	std::ofstream of(stagesummary_fn.c_str());
	of << inf;
}


void recalculate_ml_colors(mobject &m)
{
	// calculate distance, from paralax relation
	sdss_star s;
	s.ra = m.ra*ctn::d2r; s.dec = m.dec*ctn::d2r;
	coordinates::equgal(s.ra, s.dec, s.l, s.b);
	FOR(0, 3) { s.mag[i+1] = m.mag[i]; s.magErr[i+1] = m.magErr[i]; }

	if(paralax(s))	// calculate the absolute magnitude and distances
	{
//		if(abs(m.D/s.earth.D - 1) > 0.5 && m.ml_mag[1]-m.ml_mag[2] > 0.1) {
//		cerr << m.D << " " << s.earth.D << "\n";
//		cerr << m.ml_mag[0] << " " << s.ml_g << "\n";
//		cerr << m.ml_mag[1] << " " << s.ml_r << "\n";
//		cerr << m.ml_mag[2] << " " << s.ml_i << "\n";
/*		cout << m.ml_mag[1]-m.ml_mag[2] << " " << s.ml_r-s.ml_i << " ";
		cout << m.ml_mag[0]-m.ml_mag[1] << " " << s.ml_g-s.ml_r << "\n";*/
//		}
		m.D = s.earth.D;
		m.ml_mag[0] = s.ml_g;
		m.ml_mag[1] = s.ml_r;
		m.ml_mag[2] = s.ml_i;
	} else {
		m.ml_mag[0] = m.ml_mag[1] = m.ml_mag[2] = 0;
	}
}

void recalculate_absolute_magnitude(mobject &m)
{
	if(m.ml_mag[1] == 0 && m.ml_mag[2] == 0) return;

	const float ml_ri = m.ml_mag[1] - m.ml_mag[2];
	float Mr = paralax.Mr(ml_ri);
	double D = stardist::D(m.ml_mag[1], Mr);
	if(abs(D/m.D - 1) > 0.0001)
	{
		m.D = D;
	}
}

void reprocess_driver()
{
	DMMArray<mobject> uniq("dm_unique_stars.dmm", "rw");
	DMMArray<obsv_mag> obsv_mags("dm_starmags.dmm");

	ticker tick(10000);
	vector<obsv_mag> obsv;
	FOR(0, uniq.size())
	{
		mobject &m = uniq[i];

#if 0
		obsv.clear();
		obsv.reserve(m.n);

		FORj(m.obs_offset, m.obs_offset + m.n)
		{
			obsv.push_back(obsv_mags[j]);
		}
#endif
		// do stuff
		recalculate_ml_colors(m);
		//recalculate_absolute_magnitude(m);

		tick.tick();
	}
}

#define DUMP(x) cerr << #x << " = " << x << "\n";

/*
	Find the ID of the first observation for each run, and 
	store them to a simple text file. This is mainly to help
	debugging and mapping observation IDs to catalog entries.
*/
struct columns_t { int col[6]; };
inline OSTREAM(const columns_t &c) { out << c.col[0]; FOR(1, 6) { out << "\t" << c.col[i]; } return out; }
void make_run_index_offset_map(std::ostream &out, const std::string &runidxfn)
{
	DMMArray<obsv_id> runidx;
	runidx.open(runidxfn, "r");

	map<int, columns_t> runoffs;

	int oldrun = -1, oldcol = -1;
	FOR(0, runidx.size())
	{
		int run = runidx[i].sloanId.run();
		int col = runidx[i].sloanId.col();
		if(run != oldrun)
		{
			ASSERT(runoffs.count(run) == 0);
			oldrun = run;
			col = -1;
			FORj(j,0,6) { runoffs[run].col[j] = -1; }
		}
		if(col != oldcol)
		{
			ASSERT(runoffs[run].col[col] == -1) { DUMP(runoffs[run].col[col]); }
			runoffs[run].col[col] = i;

			cout << run << "\t" << col << "\t" << i << "\n";
			oldcol = col;
		}
	}

	FOREACH(runoffs) { out << (*i).first << "\t" << (*i).second << "\n"; }
}

void loadRuns(set<int> &runs, const std::string &runfile)
{
	string fn(runfile.size() ? runfile : "catalogs/runs.txt");
	cerr << "Loading runs from: " << runfile << "\n";
	text_input_or_die (in, runfile);
	load(in, runs, 0);
	in_stream.close();
}

void obs_info(
	int fitsID,
	const std::string &uniqObjectCat, // "dm_unique_stars.dmm" [output]
	const std::string &uniqObsvCat, // "dm_starmags.dmm" [output]

	const std::string &selindex, // name of the fitsID -> obsv ID map (DMM file) [in]
	const std::string &runindex // name of the sloanID -> uniqID (DMM file) [in]
)
{
	// open catalogs
	DMMArray<mobject> unique(uniqObjectCat);
	DMMArray<int> selidx(selindex);
	DMMArray<obsv_id> runidx(runindex);
	
	if(fitsID < 0 || fitsID > selidx.size())
	{
		cerr << "fitsID out of range (0, " << selidx.size() << ")\n";
		return;
	}

	// map fitsId to selId
	int selid = selidx[fitsID];
	if(selid == -1)
	{
		cerr << "This observation [fitsId = " << fitsID << "] was not included in our sample\n";
		return;
	}

	// observation information
	obsv_id sid = runidx[selid];
	ASSERT(sid.fitsId == fitsID);
	cout << "fitsId = " << sid.fitsId << "\n";
	cout << "uniqId = " << sid.uniqId << "\n";
	cout << "sloanID = " << sid.sloanId << "\n";

	// object information
	mobject m = unique[sid.uniqId];
	print_mobject(cout, m);
}

void print_mobject(std::ostream &out, const mobject &m)
{
	out << "ra = " << m.ra << ", dec = " << m.dec << "\n";
	out << "Magnitudes (ugriz):\n";
	out << "mag    = "; FOR(0, 5) { cout << m.mag[(i+4) % 5] << " "; } cout << "\n";
	out << "magErr = "; FOR(0, 5) { cout << m.magErr[(i+4) % 5] << " "; } cout << "\n";
	out << "N      = "; FOR(0, 5) { cout << m.N[(i+4) % 5] << " "; } cout << "\n";
	out << "ml_mag = "; FOR(0, 3) { cout << m.ml_mag[i] << " "; } cout << "\n";
	out << "Distances:\n";
	//out << "Mr = " << m.Mr << "\n";
	V3 p; calculate_cartesian(p, m);
	out << "\tD = " << abs(p) << "\n";
	out << "\t(x, y, z) = " << p.x << " " << p.y << " " << p.z << "\n";
}
