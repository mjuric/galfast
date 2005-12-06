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
#if 1

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

bool filter(sdss_star &s, const valarray<int> &flags, const valarray<int> &flags2)
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
	if(!finite(s.mag[2])) { return false; }		// must have r magnitude
	if(nValidBands < 2) { return false; }		// must have at least one color
	if(s.r > 22) { return false; }			// r magnitude cutoff

	return true;
}

template<typename Streamer>
void makelookup(
	std::set<int> &runs, 		// list of runs to import
	const std::string &select, 	// name of the object catalog (indexed by uniq ID) (DMM file)
	const std::string &selindex,// name of the fitsID -> uniq ID map (DMM file)
	const std::string &runindex // name of the sloanID -> uniqID (DMM file)
	)
{
#ifdef HAVE_LIBCCFITS
	// setup SDSS FITS file catalog streamer
	sdss_star s;
	valarray<int> f1(5), f2(5);	// photometry flags (FLAGS and FLAGS2 FITS fields)
	Streamer cat(runs, s, f1, f2);

	// open arrays for catalog and index
	DMMArray<star> sel;	// catalog of selected stars
	DMMArray<int> selidx;	// fitsID -> selIndex map. May be sparse, because not all FITS records will be accepted
	DMMArray<starid> runidx;// sloanID -> uniqID map - just creating a placeholder array here

	sel.create(select);
	selidx.create(selindex);
	runidx.create(runindex);

	ticker tick(10000);
	star ss; starid sid;
	sid.uniqId = -1;

	while(cat.next())
	{

		tick.tick();
		if(cat.fitsId() % 1000000 == 0) { sel.sync(); selidx.sync(); runidx.sync(); }

		bool accept = filter(s, f1, f2);

		if(!accept) { selidx[cat.fitsId()] = -1; continue; }

		// record indices
		ss.fitsId = sid.fitsId = cat.fitsId();
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
	cout << "Total stars in FITS files: " << cat.fitsId()+1 << "\n";
	cout << "Total stars accepted     : " << sel.size() << "\n";
	cout << "\n";
	cout.flush();
#else
	cerr << "CCfits support not compiled\n";
	abort();
#endif
}

mobject process_observations(int obs_offset, double ra, double dec, float Ar, std::vector<starmag> &obsv)
{
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

	// magnitudes
	FOR(0, 5)
	{
		// sum things up
		float f[N], w[N], mag[N];
		int n = 0;
		FORj(k, 0, N)
		{
			starmag &s = obsv[k];

			// throw out bad measurements
			if(!finite(s.mag[i]) || !finite(s.magErr[i]))
			{
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
	sdss_star s;
	s.ra = m.ra*ctn::d2r; s.dec = m.dec*ctn::d2r;
	coordinates::equgal(s.ra, s.dec, s.l, s.b);
	FOR(0, 3) { s.mag[i+1] = m.mag[i]; s.magErr[i+1] = m.magErr[i]; }

	if(paralax(s))	// calculate the absolute magnitude and distances
	{
		m.D = s.earth.D;
		m.ml_mag[0] = s.ml_g;
		m.ml_mag[1] = s.ml_r;
		m.ml_mag[2] = s.ml_i;
	}

	return m;
}

// defined in binner_tng.cpp
bool calculate_cartesian(V3 &earth, const mobject &m);

//#define NOMOD
void average_magnitudes()
{
	binary_input_or_die(in, "match_groups.lut");

	// open catalogs
	DMMArray<star> sel; sel.setmaxwindows(40); sel.setwindowsize(1024*1024 / 5);
	DMMArray<int> selidx; selidx.setmaxwindows(40); selidx.setwindowsize(1024*1024 / 5);
	DMMArray<starid> runidx; runidx.setmaxwindows(40); runidx.setwindowsize(1024*1024 / 5);

	sel.open("dm_tmpcat.dmm", "r");
	selidx.open("dm_tmpcat_index.dmm", "r");
#ifndef NOMOD
	runidx.open("dm_run_index.dmm", "rw");
#else
	runidx.open("dm_run_index.dmm", "r");
#endif
	int N, id, size;
	vector<starmag> obsv;

	in >> size;
	DMMArray<mobject> out;
	DMMArray<starmag> starmags;
#ifndef NOMOD
	out.create("dm_unique_stars.dmm");
	starmags.create("dm_starmags.dmm");
#else
	out.open("dm_unique_stars.dmm", "r");
	starmags.open("dm_starmags.dmm", "r");
#endif
	ticker tick(10000);
	unsigned int t = time(NULL), tstart = t;
	FORj(k, 0, size)
	{
		int uniqId = out.size();

		// load the whole group
		in >> N;
		obsv.clear();
		double ra, dec; float Ar;
		for(int i = 0; i != N; i++)
		{
			in >> id;

			// check if this observation was selected as valid
			int sid = selidx[id];
			if(sid == -1) continue;
#ifndef NOMOD
			// load this observation
			star &s = sel[sid];
			obsv.push_back(s);
			
			if(obsv.size() == 1) { ra = s.ra; dec = s.dec; Ar = s.Ar; }
#endif
			
			// store ID of the unique object together with this observation
#ifndef NOMOD
			starid &stid = runidx[sid];
			stid.uniqId = uniqId;
			
			ASSERT(stid.fitsId == s.fitsId);
#endif
		}

#ifndef NOMOD
		if(obsv.size() != 0)		// ignore groups with zero selected observations
		{
			// process and store
			mobject m = process_observations(starmags.size(), ra, dec, Ar, obsv);

			FOREACH(obsv) { starmags.push_back((starmag&)*i); }

			out.push_back(m);
		}
#endif
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
	DMMArray<starmag> starmags("dm_starmags.dmm");

	ticker tick(10000);
	vector<starmag> obsv;
	FOR(0, uniq.size())
	{
		mobject &m = uniq[i];

#if 0
		obsv.clear();
		obsv.reserve(m.n);

		FORj(m.obs_offset, m.obs_offset + m.n)
		{
			obsv.push_back(starmags[j]);
		}
#endif
		// do stuff
		recalculate_ml_colors(m);
		//recalculate_absolute_magnitude(m);

		tick.tick();
	}
}


/*
	Find the ID of the first observation for each run, and 
	store them to a simple text file. This is mainly to help
	debugging and mapping observation IDs to catalog entries.
*/
void make_run_index_offset_map(std::ostream &out, const std::string &runidxfn)
{
	DMMArray<starid> runidx;
	runidx.open(runidxfn, "r");

	map<int, int> runoffs;

	int oldrun = -1;
	FOR(0, runidx.size())
	{
		int run = runidx[i].sloanId.run();
		if(run != oldrun)
		{
			ASSERT(runoffs.count(run) == 0);
			cout << run << "\t" << i << "\n";
			runoffs[run] = i;
			oldrun = run;
		}
	}

	FOREACH(runoffs) { out << (*i).first << "\t" << (*i).second << "\n"; }
}

void loadRuns(set<int> &runs, const std::string &runfile)
{
	string fn(runfile.size() ? runfile : "catalogs/runs.txt");
	cout << "Loading runs from: " << runfile << "\n";
	text_input_or_die (in, runfile);
	load(in, runs, 0);
	in_stream.close();
}

void starinfo(int fitsID)
{
	// open catalogs
	DMMArray<star> sel("dm_tmpcat.dmm");
	DMMArray<int> selidx("dm_tmpcat_index.dmm");
	DMMArray<starid> runidx("dm_run_index.dmm");
	DMMArray<mobject> unique("dm_unique_stars.dmm");
	
	if(fitsID < 0 || fitsID > selidx.size())
	{
		cerr << "fitsID out of range (0, " << selidx.size() << ")\n";
		return;
	}

	int selid = selidx[fitsID];
	if(selid == -1)
	{
		cerr << "This observation was not included in our sample\n";
		return;
	}

	starid sid = runidx[selid];
	ASSERT(sid.fitsId == fitsID);
	cout << "uniqId = " << sid.uniqId << "\n";
	cout << "sloanID = " << sid.sloanId << "\n";

	mobject m = unique[sid.uniqId];
	print_mobject(cout, m);
}

void print_mobject(std::ostream &out, const mobject &m)
{
	out << "ra = " << m.ra << ", dec = " << m.dec << "\n";
	out << "Magnitudes:\n";
	out << "mag    = "; FOR(0, 3) { cout << m.mag[i] << " "; } cout << "\n";
	out << "magErr = "; FOR(0, 3) { cout << m.magErr[i] << " "; } cout << "\n";
	out << "N      = "; FOR(0, 3) { cout << m.N[i] << " "; } cout << "\n";
	out << "ml_mag = "; FOR(0, 3) { cout << m.ml_mag[i] << " "; } cout << "\n";
	out << "Distances:\n";
	//out << "Mr = " << m.Mr << "\n";
	V3 p; calculate_cartesian(p, m);
	out << "D = " << abs(p) << "\n";
	out << "(x, y, z) = " << p.x << " " << p.y << " " << p.z << "\n";
}
#endif