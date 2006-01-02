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

#include "dm.h"
#include "paralax.h"

#include <set>
#include <iomanip>
#include <stdlib.h>

#include "ximage.h"
#include <astro/io/fits.h>
#include <astro/system/shell.h>
#include <astro/useall.h>
using namespace std;

//////////////////////

class nformat
{
protected:
	int w, p;
public:
	nformat(int width = 10, int precision = -2)
	{
		w = width;
		if(precision == -1)
		{
			precision = width - 5;
		}
		p = precision;
	}

	friend OSTREAM(const nformat &);
};
OSTREAM(const nformat &nf)
{
	out << setw(nf.w);
	if(nf.p >= 0) out << setprecision(nf.p);
	return out;
}
typedef nformat nfmt;

OSTREAM(const mobject &m)
{
	char buf[1000];
	sprintf(buf, "%9.5f %9.5f %9.2f  %8.4f %6.4f %3d  %8.4f %6.4f %3d  %8.4f %6.4f %3d  %8.4f %6.4f %3d  %8.4f %6.4f %3d   %8.4f %8.4f %8.4f",
		m.ra, m.dec, m.D,
		m.mag[0], m.magErr[0], (int)m.N[0],
		m.mag[1], m.magErr[1], (int)m.N[1],
		m.mag[2], m.magErr[2], (int)m.N[2],
		m.mag[3], m.magErr[3], (int)m.N[3],
		m.mag[4], m.magErr[4], (int)m.N[4],
		m.ml_mag[0], m.ml_mag[1], m.ml_mag[2]);
	out << buf;

	return out;
}

OSTREAM(const obsv_mag &sm)
{
	char buf[1000];
	sprintf(buf, "%12d %8.4f %6.4f  %8.4f %6.4f  %8.4f %6.4f  %8.4f %6.4f  %8.4f %6.4f",
		sm.fitsId,
		sm.mag[0], sm.magErr[0],
		sm.mag[1], sm.magErr[1],
		sm.mag[2], sm.magErr[2],
		sm.mag[3], sm.magErr[3],
		sm.mag[4], sm.magErr[4]
		);
	return out << buf;
}

OSTREAM(const binned_runset &brs)
{
	int nstars = 0, nobs = 0;
	FOREACH(brs.pixels)
	{
		const binned_runset::pixel &p = (*i).second;
		nstars += p.uniqueN;
		nobs += p.N;
	}
	
	out << "# dx = " << brs.dx << "\n";
	out << "# runs = {"; out << brs.runs.size() << " : "; FOREACH(brs.runs) { out << " " << *i; }; out << " }\n";
	out << "# colorbins = {"; FOREACH(brs.colorbins) { out << " [" << (*i).first << ", " << (*i).second << ")"; }; out << " }\n";
	out << "# pixels = " << brs.pixels.size() << "\n";
	out << "# stars = " << nstars << "\n";
	out << "# observations = " << nobs << "\n";
	out << "#\n";
	out << "#   x          y       z    obs      volume    unq   unqVolume Nrun runs[Nrun]\n";

	FOREACH(brs.pixels)
	{
		const S3 &idx = (*i).first;
		const binned_runset::pixel &p = (*i).second;

		out << setw(8) << brs.dx*idx.x << setw(8) << brs.dx*idx.y << setw(8) << brs.dx*idx.z;
		out << setw(7) << p.N << " " << setw(11) << p.volume << setw(7) << p.uniqueN << " " << setw(11) << p.uniqueVolume;
		out << setw(4) << p.runs.size();
		FOREACHj(j, p.runs) { out << setw(6) << *j; }

		out << "\n";
	}

	return out;
}

template<typename V>
BOSTREAM(const std::set<V> &m)
{
	out << m.size(); FOREACH(m) { out << *i; }
	return out;
}

template<typename V>
BISTREAM(std::set<V> &m)
{
	size_t size;
	V val;
	in >> size; FOR(0, size) { in >> val; m.insert(val); }
	return in;
}

BOSTREAM(const binned_runset::pixel &p)
{
	return out << p.N << p.uniqueN << p.volume << p.uniqueVolume << p.runs;
}

BISTREAM(binned_runset::pixel &p)
{
	return in >> p.N >> p.uniqueN >> p.volume >> p.uniqueVolume >> p.runs;
}

BOSTREAM(const binned_runset &brs)
{
	return out << brs.dx << brs.runs << brs.colorbins << brs.pixels;
}

BISTREAM(binned_runset &brs)
{
	return in >> brs.dx >> brs.runs >> brs.colorbins >> brs.pixels;
}

//////////////////////

#include <xcat.h>

sdss::RunGeometryDB db;
std::map<int, sdss::Mask> masks;
std::map<int, sdss::RunGeometry> geoms;

bool boundsCheck(mobject &m, int run)
{
	if(!masks.count(run))
	{
		masks[run] = geoms[run] = db.getGeometry(run);
	}
	sdss::RunGeometry &geom = geoms[run];
	Radians mu, nu;
	Coordinates::equgcs(geom.node, geom.inc, rad(m.ra), rad(m.dec), mu, nu);
	return masks[run].contains(mu, nu);
}

#define DIST_RECALC 0

#if DIST_RECALC

static plx_gri_locus paralax;
bool calculate_distance(mobject &m)
{
	// calculate distance, from paralax relation
	sdss_star s;
	s.ra = m.ra*ctn::d2r; s.dec = m.dec*ctn::d2r;
	coordinates::equgal(s.ra, s.dec, s.l, s.b);
	FOR(0, 3) { s.mag[i+1] = m.mag[i]; s.magErr[i+1] = m.magErr[i]; }
	s.earth.D = sqrt(sqr(m.x) + sqr(m.y) + sqr(m.z));
	if(s.calculate_distances(paralax, false))
	{
		m.x = s.earth.x; m.y = s.earth.y; m.z = s.earth.z;
/*		m.Mr = s.Mr;
		m.ml_mag[0] = s.ml_g;
		m.ml_mag[1] = s.ml_r;
		m.ml_mag[2] = s.ml_i;*/
		return true;
	} else {
		m.x = m.y = m.z = 0;
		return false;
	}
}
#endif

bool calculate_cartesian(V3 &earth, const mobject &m)
{
	if(m.D == 0) return false;

	Radians l, b;
	coordinates::equgal(m.ra*ctn::d2r, m.dec*ctn::d2r, l, b);

	float D = m.D;

	// switching to rectangular geocentric coordinates
	float Re = D*cos(b);			// geocentric cylindrical coord. sys. distance
	earth.x = Re*cos(l);			// x distance from Earth
	earth.y = Re*sin(l);			// y distance from Earth
	earth.z = D*sin(b);			// z distance from Earth
	earth.x *= -1; earth.y *= -1;		// rotate the coordinate system to point _away_ from the galactic center (180deg rotation)

	return true;
}

static bool bin3d = true;

inline
binned_run::pixel &find_bin(binned_run &br, const mobject &m, double dx)
{
	V3 p;
	calculate_cartesian(p, m);
	
	if(bin3d)
	{
		int x = (int)floor(p.x / dx + 0.5);
		int y = (int)floor(p.y / dx + 0.5);
		int z = (int)floor(p.z / dx + 0.5);
		
		S3 idx(x, y, z);
		return br.pixels[idx];
	}
	else
	{
		double	x = (p.x + Rg) / dx,
			y = (p.y) / dx;
		double rho = sqrt(sqr(x) + sqr(y));
		S3 rz((short)floor(rho + 0.5), (short)floor(p.z / dx + 0.5), 0);
		return br.pixels[rz];
	}

	ASSERT(0);
}

int bin_run(binned_run &br, int run, double dx, ribin colorbin, double Dmin, double Dmax,
	DMMArray<mobject> &unique,
	DMMArray<obsv_id> &runindex,
	map<int, int> &runoffsets,
	int *boundsrejected_ = NULL,
	int *novolumerejected_ = NULL)
{
	int boundsrejected = 0;
	int novolumerejected = 0;

	// setup binned_run structure
	br.run = run;
	br.dx = dx;
	br.colorbin = colorbin;
	br.Dmin = Dmin;
	br.Dmax = Dmax;

	// find offset from which to start loading observations in a run
	int at = 0; obsv_id sid = runindex[0];
	if(run > 0)
	{
		at = runoffsets[run];
		sid = runindex[at];
		ASSERT(sid.sloanId.run() == run);
	}

	int n = 0; at--;
	V3 p;
	ticker tick(10000);
	for(;;)
	{
		tick.tick();

		if(++at == runindex.size()) break;
		sid = runindex[at];
		if(run > 0 && sid.sloanId.run() != run) break;
		mobject m = unique[sid.uniqId];

#if DIST_RECALC
		calculate_distance(m);
#endif

		// check if this observation matches a set of criteria
		// color selection
		float ri = m.ml_mag[1] - m.ml_mag[2]; // r-i color
		if(!(colorbin.first <= ri && ri < colorbin.second)) { continue; }
		// volume limited selection
		if(m.D < Dmin || m.D >= Dmax) { continue; }
		// check if this observation is within formal run bounds
		switch(run)
		{
		case -1: // Saggitarius dwarf cut
		//	if(!toSagDw(m)) { continue; }
			break;
		default:
			if(!boundsCheck(m, run)) { boundsrejected++; continue; }
			break;
		}

		// find the pixel
		binned_run::pixel &p = find_bin(br, m, dx);

		// check if this pixel has any volume assigned to it
		if(p.volume == 0)
		{
			novolumerejected++;
		}

		// add this observation to pixel
		p.stars.insert(sid.uniqId);
		n++;
	}

	br.loaded = true;
	if(boundsrejected_ != NULL) { *boundsrejected_ = boundsrejected; }
	if(novolumerejected_ != NULL) { *novolumerejected_ = novolumerejected; }
	
	return n;
}

void load_run_index_offset_map(map<int, int> &runoffs)
{
	input_or_die(in, "dm_run_index.map");
	int run, offs;
	
	in >> run >> offs;
	while(!in.eof())
	{
		runoffs[run] = offs;
		in >> run >> offs;
	}
}

bool color_less(const binned_run *a, const binned_run *b)
{
	if(a->colorbin.first < b->colorbin.first) return true;
	if(a->colorbin.first > b->colorbin.first) return false;
	if(a->run < b->run) return true;
	return false;
}

void merge(binned_runset &brs, vector<binned_run *> &runs)
{
	ASSERT(runs.size());
	
	// sort list by colorbin, run
	sort(runs.begin(), runs.end(), color_less);

	binned_run &first = *runs.front();
	bool shouldUnload = first.load();
	brs.dx = first.dx;

	binned_run::pixelmap tmp;

	vector<binned_run *>::iterator i = runs.begin();
	while(i != runs.end())
	{
		binned_run &br = *(*i);
		ASSERT(br.dx == brs.dx);

		FOREACHj(j, br.pixels)
		{
			const S3 &idx = (*j).first;
			binned_run::pixel &p = (*j).second;
			
			binned_run::pixel &dest = tmp[idx];
			dest.stars.insert(p.stars.begin(), p.stars.end());

			binned_runset::pixel &merged = brs.pixels[idx];
			merged.N += p.stars.size();
			merged.volume += p.volume;
			merged.runs.insert(br.run);
		}
		brs.runs.insert(br.run);
		brs.colorbins.insert(br.colorbin);

		if(shouldUnload) br.unload();

		float ri0 = br.colorbin.first;
		++i;

		if(i != runs.end()) { shouldUnload = (*i)->load(); }

		if(i == runs.end() || (*i)->colorbin.first != ri0)
		{
			// add up unique stars for each pixel
			FOREACHj(j, tmp)
			{
				const S3 &idx = (*j).first;
				binned_run::pixel &dest = (*j).second;
				binned_runset::pixel &merged = brs.pixels[idx];

				merged.uniqueN += dest.stars.size();
			}

			// clear tmp map
			tmp.clear();
		}

		cout << br.run << "\n";
	}
}

int bin(const set<int> &runs, pair<float, float> r, pair<float, float> ri)
{
	gsl_set_error_handler_off ();

	DMMArray<obsv_id> runidx("dm_run_index.dmm");
	DMMArray<mobject> unique("dm_unique_stars.dmm");
	map<int, int> runoffs;

	load_run_index_offset_map(runoffs);

	double dx;
	double Dmin, Dmax;
	int boundsr, novolr;

	plx_gri_locus paralax;
	paralax.distance_limits(Dmin, Dmax, ri.first, ri.second, r.first, r.second);
	cout << "Distance limits: " << Dmin << " " << Dmax << "\n";

	FOREACH(runs)
	{
		int run = *i;
		binned_run br;
		string fn, fnout;

		switch(run)
		{
		case -1:	// Saggitarius dwarf cut
			fn = io::format("maps/map.sagdw.%5.3f-%5.3f.%5.3f-%5.3f.bin")
				<< r.first << r.second << ri.first << ri.second;
			fnout = io::format("bins/density.sagdw.%5.3f-%5.3f.%5.3f-%5.3f.bin")
				<< r.first << r.second << ri.first << ri.second;
			break;
		default:
			fn = io::format("maps/map.%s%05d.%5.3f-%5.3f.%5.3f-%5.3f.bin")
				<< (bin3d ? "" : "rz.") << run << r.first << r.second << ri.first << ri.second;
			fnout = io::format("bins/density.%s%05d.%5.3f-%5.3f.%5.3f-%5.3f.bin")
				<< (bin3d ? "" : "rz.") << run << r.first << r.second << ri.first << ri.second;
			break;
		}
//		fn = "maps/merged.1740.bin";
//		cerr << "FN:NNNNN: " << fn << "\n";

		binary_input_or_die(in, fn);
		in >> br;
		
		if(i == runs.begin()) { dx = br.dx; }
//		cerr << setprecision(12) << br.dx << " " << setprecision(12) << br.Dmin << "=?=" << setprecision(12) << Dmin << " " << setprecision(12) << br.Dmax << "=?=" << setprecision(12) << Dmax << "\n";
		ASSERT(br.dx == dx);
		ASSERT(abs(br.Dmin - Dmin) < 0.01);
		ASSERT(abs(br.Dmax - Dmax) < 0.01);

		int n = bin_run(br, run, dx, ri, Dmin, Dmax, unique, runidx, runoffs, &boundsr, &novolr);

		binary_output_or_die(out, fnout);
		out << br;

		// debugging
		output_or_die(tout, "bin.txt");
		tout << br;

		cout << run << ": binned = " << n << " [" << boundsr << " " << novolr << "]\n";
	}

	return 0;
}

int make_run_plots(const set<int> &runs)
{
	gsl_set_error_handler_off ();

	DMMArray<obsv_id> runindex("dm_run_index.dmm");
	DMMArray<mobject> unique("dm_unique_stars.dmm");
	map<int, int> runoffs;

	load_run_index_offset_map(runoffs);

	double dx;
	double Dmin, Dmax;
	int boundsr, novolr;

	FOREACHj(j, runs)
	{
		int run = *j;
		cerr << "Run " << run << "...";

		{
			text_output_or_die(out, "run.tmp");
	
			int at = runoffs[run];
			for(;;)
			{
				if(++at == runindex.size()) break;
				obsv_id sid = runindex[at];
				if(sid.sloanId.run() != run) break;
	
				mobject &m = unique[sid.uniqId];
				out << m.ra << m.dec << nl();
			}
		}

		cerr << " Running SM...";

		::system("/usr/local/bin/sm -m mkplots.sm");
		cerr << string(io::format("./mvplots.csh %05d ") << run) << "\n";
		::system(string(io::format("./mvplots.csh %05d ") << run).c_str());
	}

	return 0;
}

// defined in raytrace.cpp
void mergeUniqVolume(binned_runset &brs, const std::string &uniqMapFn, pair<float, float> r, pair<float, float> ri);

int merge_maps(const std::string &outputPrefix, const set<int> &runs, pair<float, float> r, pair<float, float> ri, const std::string &uniqMapFn)
{
	binned_runset brs;

	vector<binned_run *> bruns;
	FOREACH(runs)
	{
		int run = *i;
		string fn = io::format("bins/density.%s%05d.%5.3f-%5.3f.%5.3f-%5.3f.bin")
			<< (bin3d ? "" : "rz.") << run << r.first << r.second << ri.first << ri.second;
		bruns.push_back(new binned_run(fn));
	}

	merge(brs, bruns);
	if(uniqMapFn.size()) { mergeUniqVolume(brs, uniqMapFn, r, ri); }

	FOREACH(bruns) { delete *i; }

	// store merged 3D volumes
	string fn = io::format("merged/txt/%s.%5.3f-%5.3f.%5.3f-%5.3f.txt")
		<< outputPrefix << r.first << r.second << ri.first << ri.second;
	output_or_die(out, fn);
	out << brs;

	fn = io::format("merged/bin/%s.%5.3f-%5.3f.%5.3f-%5.3f.bin")
		<< outputPrefix << r.first << r.second << ri.first << ri.second;
	binary_output_or_die(bout, fn);
	bout << brs;

	return 0;
}

#if 0
int make_z_slice_fits(const std::string &prefix, const binned_runset &brs, float z)
{
	// convert a slice into a FITS file
	cout << "Making FITS slice\n";

	// autodetect image size
	FOREACH(brs.pixels)
	{
		const S3 idx = (*i).first;
		if(idx.z != zz) continue;
	}

	XImage img(48, 48);
	int xc = img.x() / 2;
	int yc = img.y() / 2;

	int zz = (int)floor(z / brs.dx + 0.5);

	FOREACH(brs.pixels)
	{
		const S3 idx = (*i).first;
		if(idx.z != zz) continue;

		const binned_runset::pixel &p = (*i).second;
		img(idx.x + xc, idx.y + yc) = float(p.N) / p.volume;
	}

	fits::write(img, "density.fits");

	return zz;
}
#endif

int bin()
{
	DMMArray<obsv_id> runidx("dm_run_index.dmm");
	DMMArray<mobject> unique("dm_unique_stars.dmm");
	typedef map<int, int> romap;
	romap runoffs;

	load_run_index_offset_map(runoffs);
#if 0

	float r0 = 15, r1 = 20.5;
	float ri0 = 1, ri1 = 1.1;
	double Dmin, Dmax;
	plx_gri_locus::distance_limits(Dmin, Dmax, ri0, ri1, r0, r1);
	cout << "Distance limits: " << Dmin << " " << Dmax << "\n";
	double dx = 10;
	int boundsr, novolr;
#if 0
	// test of volume/counts compatibility
	binned_run br;
	int run = 745;
	binary_input_or_die(in, "br_vol.745.bin");
	in >> br;
	ASSERT(br.dx == dx);
	int n = bin_run(br, run, dx, make_pair(ri0, ri1), Dmin, Dmax, unique, runidx, runoffs, &boundsr, &novolr);
	cout << run << ": binned = " << n << " [" << boundsr << " " << novolr << "]\n";
	std::string fn = io::format("bins/runtest.%05d.bin") << run;
	binary_output_or_die(out, fn);
	out << br;

	return 0;
#else
	FOREACH(runoffs)
	{
		int run = (*i).first;
		binned_run br;
		std::string fn = io::format("maps/map.%05d.bin") << run;
		binary_input_or_die(in, fn);
		in >> br;
		if(i == runoffs.begin()) { dx = br.dx; }
		ASSERT(br.dx == dx);
		int n = bin_run(br, run, dx, make_pair(ri0, ri1), Dmin, Dmax, unique, runidx, runoffs, &boundsr, &novolr);
		cout << run << ": binned = " << n << " [" << boundsr << " " << novolr << "]\n";

		fn = io::format("bins/density.%05d.bin") << run;
		binary_output_or_die(out, fn);
		out << br;
	}
#endif
#else
	binned_runset brs;

	vector<binned_run *> runs;
	FOREACH(runoffs)
	{
		int run = (*i).first;
		std::string fn = io::format("bins/density.%05d.bin") << run;
		runs.push_back(new binned_run(fn));
	}

	merge(brs, runs);
	
	output_or_die(out, "merged.txt");
	out << brs;
	binary_output_or_die(bout, "merged.bin");
	bout << brs;

#ifndef XX
	// convert a slice into a FITS file
	cout << "Making FITS slice\n";
	XImage img(48, 48);
	int xc = img.x() / 2;
	int yc = img.y() / 2;
	int z = 400 / 50;
	
	FOREACH(brs.pixels)
	{
		const S3 idx = (*i).first;
		if(idx.z != z) continue;
		
		const binned_runset::pixel &p = (*i).second;
		img(idx.x + xc, idx.y + yc) = float(p.N) / p.volume;
	}

#ifdef HAVE_PKG_CCfits
	fits::write(img, "density.fits");
#else
	ASSERT(0) { cerr << "libCCfits support not compiled in.\n"; }
#endif
#endif
#endif

}
