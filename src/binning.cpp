#include "config.h"

#ifdef HAVE_LIBCCFITS

#include <astro/constants.h>
#include <astro/coordinates.h>
#include <astro/system/log.h>
#include <astro/util.h>
#include <astro/sdss/rungeometry.h>
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

#include "textstream.h"
#include "analysis.h"

#include "vrml.h"
#include "gslcc.h"
#include "ximage.h"
#include "xcat.h"
#include "integerbin.h"
#include "sparsevolume.h"

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

using namespace std;

/////////////

#if 0

template<typename T>
class mmopen_struct
{
public:
	int size;
	MemoryMapVector<T> &v;
	int mflag;

public:
	mmopen_struct(MemoryMapVector<T> &v_, int size_ = -1, const std::string &mode = "r");

	void open(std::string &fn, int offset);
};

template<typename T>
mmopen_struct<T>::mmopen_struct(MemoryMapVector<T> &v_, int size_, const std::string &mode)
	: v(v_), size(size_)
{
	     if(mode == "r" ) mflag = MemoryMap::ro;
	else if(mode == "rw") mflag = MemoryMap::rw;
	else if(mode == "w" ) mflag = MemoryMap::wo;
	else { THROW(EIOException, "Unknown mode flag (can be one of 'r', 'w' or 'rw')"); }
}

template<typename T>
void mmopen_struct<T>::open(std::string &fn, int offset)
{
	ASSERT(size > -1);
	v.open(fn, size, offset, mflag);
}

template<typename T>
mmopen_struct<T> mmopen(MemoryMapVector<T> &v, int size = -1, std::string mode = "")
{
	return mmopen_struct<T>(v, size, mode);
}

template<typename T>
ibinarystream &operator >>(ibinarystream &in, mmopen_struct<T> &mm)
{
	if(!in.filename) { THROW(EIOException, "Input binary stream is not tied to a file (cannot memory map something that is not a file)"); }
	if(mm.size == -1) { in >> mm.size; }
	int offset = in.f.tellg();
	mm.open(in.filename, offset);
	in.seekg(ios_base::cur, sizeof(T) * mm.size);
	return in;
}

#endif

/////////////

inline double frac(double x) { return x - trunc(x); }

void resample(valarray<float> &out, valarray<float> &in,
	int w0, int h0, double xo, double yo,
	int w1, int h1, double xc, double yc,
	double b)
{
	ASSERT(b >= 1);

	ind2<float>  iin(in,  w0, h0),
	            iout(out, w1, h1);
	double b2 = sqr(b);

	ASSERT(w0*h0 == in.size());
	ASSERT(w1*h1 == out.size());

	FORj(i, 0, w0) FORj(j, 0, h0)
	{
		float v = iin(i, j);
		
		double x0 = xc + (i - xo) / b;
		double y0 = yc + (j - yo) / b;
		double x1 = xc + (i + 1 - xo) / b;
		double y1 = yc + (j + 1 - yo) / b;

		int i0((int)x0), i1((int)x1), j0((int)y0), j1((int)y1);
		double wx0, wx1, wy0, wy1;

		if(i0 != i1)
		{
			wx0 = (1 - frac(x0)) * b;
			wx1 = 1 - wx0;

//			cout << wx0 << " " << wx1 << "\n";
		}
		else
		{
			wx0 = 1; wx1 = 0;
		}
		if(j0 != j1)
		{
			wy0 = (1 - frac(y0)) * b;
			wy1 = 1 - wy0;
		}
		else
		{
			wy0 = 1; wy1 = 0;
		}

		ASSERT(i0 >= 0 && i1 <= w1);
		ASSERT(j0 >= 0 && j1 <= h1);
		ASSERT(abs(wx0*wy0 + wx1*wy0 + wx0*wy1 + wx1*wy1 - 1) < 0.000001);
		ASSERT(j1 != h1 || (j1 == h1 && wx0*wy1 < 0.00001));
		ASSERT(j1 != h1 || (j1 == h1 && wx1*wy1 < 0.00001));
		ASSERT(i1 != w1 || (i1 == w1 && wx1*wy0 < 0.00001));
		ASSERT(i1 != w1 || (i1 == w1 && wx1*wy1 < 0.00001));

		iout(i0, j0) += v*wx0*wy0;

		if(i1 != w1) 		 iout(i1, j0) += v*wx1*wy0;
		if(j1 != h1) 		 iout(i0, j1) += v*wx0*wy1;
		if(i1 != w1 && j1 != h1) iout(i1, j1) += v*wx1*wy1;
	}

	double s1 = std::accumulate(&in[0], &in[in.size()], 0.);
	double s2 = std::accumulate(&out[0], &out[out.size()], 0.);
	ASSERT(abs(s1/s2 - 1) < 0.0001);
//	cout << abs(s1/s2 - 1) << "\n";
}

//
// Shrink image in by factor b and store the result in out
// while keeping the center of the central pixel of in as the center
// of the central pixel of out
//
void resample(XImage &out, XImage &in, double b)
{
	double xo = in.x() / 2 + .5;
	double yo = in.y() / 2 + .5;

	double w = ceil(in.x() / b) + 1;
	double h = ceil(in.y() / b) + 1;
	out.resize((int)w, (int)h);

	cout << w << " " << h << "\n";
	cout << out.x() << " " << out.y() << "\n";

	resample(out.array, in.array, in.x(), in.y(), in.x() / 2 + .5, in.y() / 2 + .5,
		out.x(), out.y(), out.x() / 2 + .5, out.y() / 2 + .5,
		b);
}

void resample(XImage &img, double b)
{
	XImage tmp;
	resample(tmp, img, b);
	img = tmp;
}

io::magick::IdentityFilter io::magick::IdentityFilter::I;

#include "paralax.h"

ticker tickk(10000);

struct filter_function
{
	virtual bool operator()(const sdss_star &s) = 0;
};

struct mstar_record
{
	sdss::RunGeometryDB geomDB;
	Config conf;

	int binfactor;
	double dx;		///< linear pixel size _after_ shrinking the volume binfactor times

	auto_ptr<sdss_star_cat> cat;

	void loadCatalog(const std::string &inputfile)
	{
		cout << "Loading catalog [" << inputfile << "] ...";
		cat.reset(new sdss_star_cat(inputfile));
		cout << " done [" << cat->runs.size() << " runs, " << cat->size() << " stars]\n";
	}

	/**
		\brief Extract the list of runs which appear in mdwarfs.dat.gz file. Store it to runs.txt file.
	*/
	void getUsedRuns()
	{
		ASSERT(0); // TODO: reimplement this because now data is in separate files
//		text_output_or_die (out, "runs.txt");

//		FOREACH(cat->index) { out << (*i).first << nl(); }
//		cout << cat->stars.size() << " stars in cat, " << cat->index.size() << " runs.\n";
	}

	/**
		\brief Sanity check to see if all the observed objects fall into given run
		masks and bounds.
	*/
	void checkMasksAndBounds()
	{
		set<int> runs;
		loadRuns(runs);

		double mu, nu;

		int cnt = 0;
		FOREACH(runs) {
			int rcnt = 0, rtotal = 0;
			const int run = *i;

			sdss::RunGeometry geom;
			geomDB.getGeometry(run, geom);
			sdss::Mask mask(geom);


			cout << "Checking run " << geom.run << "\n";
			sdss_star_cat::subset s(cat->getrun(run));
			FOREACH(s) {
				rtotal++;
				Radians ra = (*i).ra, dec = (*i).dec;

				Coordinates::equgcs(geom.node, geom.inc, ra, dec, mu, nu);

				if(!mask.contains(mu, nu)) {
					cnt++; rcnt++;
					cout << run << "    " << deg(ra) << " " << deg(dec) 
						<< "    " << deg(mu) << " " << deg(nu)
						<< "    " << deg(geom.muStart) << " " << deg(geom.muEnd) << "\n";
				}
			}
			cout << "# run " << run << " out of bounds: " << rcnt << " [out of " << rtotal << "]\n";
		}
		cout << "# total out of bounds: " << cnt << "\n";
	}

	/**
		\brief Load a list of runs specified in runs.txt file (just a list of runs, every run on its own line)

		runs.txt file is generated by getUsedRuns() and mdwarfs.dat file.
	*/
	void loadRuns(set<int> &runs, const std::string &runfile = "")
	{
		string fn(runfile.size() ? runfile : "catalogs/runs.txt");
		cout << "Loading runs from: " << runfile << "\n";
		text_input_or_die (in, runfile);
		load(in, runs, 0);
		in_stream.close();
	}

	/*
		Bin stars for a given set of runs and with a given dx

		counts image has to be preinitialized and null-ed.
	*/
	int bin(set<int> &runs, XImage &counts, double dx, filter_function &filter)
	{
		const int x = counts.x();
		const int y = counts.y();
		int yc = y / 2;
		int xc = x / 2;
		ticker tick(100);
		text_output_or_die(text, "binned.txt");

//		int k = 0;
		FOREACHj(run, runs)
		{
			sdss::RunGeometry geom;
			geomDB.getGeometry(*run, geom);
			sdss::Mask mask(geom);
			
			// cerr << "\nRun " << *run << ":\n";
			sdss_star_cat::subset ss(cat->getrun(*run));
			FOREACH(ss) {
				sdss_star &s = *i;
//				if(s.ri() > 1.0 && s.ri() < 1.1 && s.r < 20.5 && s.r > 15) k++;
//				if(s.ri() > 1.0 && s.ri() < 1.1 && s.r < 20.5 && s.r > 15 && !filter(s) && k==5) {cout << s << "\n"; exit(-1);}
				if(!filter(s)) continue;

				Radians mu, nu, l, b;
				coordinates::equgcs(geom.node, geom.inc, s.ra, s.dec, mu, nu);
				coordinates::equgal(s.ra, s.dec, l, b);
//				text << s.id << s.run << s.col << s.field << s.objid << deg(s.ra) << deg(s.dec) << s.r << s.i << s.rErr << s.iErr << s.ri() << s.ml_ri << s.Mr << nl();
				int px = xc + (int)floor((s.gc.rho-Rg) / dx + 0.5);
				int py = yc + (int)floor(s.gc.z / dx + 0.5);

				text << s.id << s.run << s.col << s.field << s.objid << deg(mu) << deg(nu) << deg(l) << deg(b) << 
					s.gc.z << nl();

#if 1
				if(!(0 <= px && px < x)) { cout << s << "\n: px = " << px << "\n"; }
				if(!(0 <= py && py < y)) { cout << s << "\n: py = " << py << "\n"; }
				ASSERT(0 <= px && px < x);
				ASSERT(0 <= py && py < y);
#endif
				counts(px, py) += 1.;
				tick.tick();
			}
		}
//		cerr << "k = " << k << "\n";
//		cerr << "Done\n";
//		cerr << (int)std::accumulate(&counts.array[0], &counts.array[counts.array.size()], 0.) << " objects\n";

		return (int)std::accumulate(&counts.array[0], &counts.array[counts.array.size()], 0.);
	}

	/*
		Bin stars to 3D volume for a given set of runs and with a given dx

		counts image has to be preinitialized and null-ed.
	*/
	int bin3D(const set<int> &runs, int b, sparse_volume &merged, sparse_volume_info &vi, double dx, filter_function &filter)
	{
		ticker tick(100);

		int n = 0;
		FOREACHj(run, runs)
		{
			sdss_star_cat::subset ss(cat->getrun(*run));
			sparse_volume vm;
			FOREACH(ss) {
				sdss_star &s = *i;
				if(!filter(s)) continue;

				// bin it
				V3 v(-s.earth.x, -s.earth.y, s.earth.z);
				I3 pix(floor(v / double(b)));
				set_or_add(vm, pix, 1);
				n++;

				tick.tick();
			}

			FOREACH(vm.volume)
			{
				set_or_add(merged, (*i).first, (*i).second);
				vi.runs[(*i).first].push_back(*run);
			}
		}

		return n;
	}

	/*
		Bin stars for a given set of runs and with a given dx, in x-y plane.
		NOTE: Pixel coordinates are centers

		counts image has to be preinitialized and null-ed.
	*/
	int binxy(set<int> &runs, XImage &counts, double dx, filter_function &filter)
	{
		const int x = counts.x();
		const int y = counts.y();
		int yc = y / 2;
		int xc = x / 2;
		ticker tick(100);

		FOREACHj(run, runs)
		{
			cerr << "\nRun " << *run << ":\n";
			sdss_star_cat::subset ss(cat->getrun(*run));
			FOREACH(ss) {
				sdss_star &s = *i;
				if(!filter(s)) continue;

				int px = xc + (int)floor(s.earth.x / dx + .5);
				int py = yc + (int)floor(s.earth.y / dx + .5);

//				cout << px << " " << py << "\n";

#if 1
				if(!(0 <= px && px < x)) { cout << s << "\n: px = " << px << "\n"; }
				if(!(0 <= py && py < y)) { cout << s << "\n: py = " << py << "\n"; }
				ASSERT(0 <= px && px < x);
				ASSERT(0 <= py && py < y);
#endif
				counts(px, py) += 1.;
				tick.tick();
			}
		}
//		cerr << "Done\n";
//		cerr << (int)std::accumulate(&counts.array[0], &counts.array[counts.array.size()], 0.) << " objects\n";

		return (int)std::accumulate(&counts.array[0], &counts.array[counts.array.size()], 0.);
	}

	void starinfo(int run, int id, int col = -1, int field = -1)
	{

		if(run != -1)
		{
			if(!cat->runs.count(run)) { cout << "Run " << run << " not in catalog\n"; return; }
			sdss_star_cat::subset curslice(cat->getrun(run));

			pair<int, float> Dmin(-1, 1E10), Dmax(-1, 0);
			FOR(0, curslice.size())
			{
				sdss_star &s = curslice[i];
				if(s.earth.D > Dmax.second) { Dmax = make_pair(i, s.earth.D); }
				if(s.earth.D < Dmin.second) { Dmin = make_pair(i, s.earth.D); }
				
				if(col != -1)
				{
					if(s.col == col && s.field == field)
					{
						cout << i << " -> " << s.objid << " ";
						if(s.objid == id) { id = i; col = 0; }
					}
				}
			}
			cout << "\n";

			cout << "Run statistics:\n";
			cout << "\tmin(earth.D) = " << Dmin.second << " [" << Dmin.first << "]\n";
			cout << "\tmax(earth.D) = " << Dmax.second << " [" << Dmax.first << "]\n";
			cout << curslice.size() << " of stars in run.\n\n";

			cout << cat->size() << " of stars in catalog.\n";

			if(field != -1 && col == -1)
			{
				cout << "Object not found.\n";
				return;
			}

			if(id >= curslice.size()) { cout << "Index out of range.\n"; }
			else			   { cout << curslice[id] << "\n"; }
		} else {
			cout << cat->size() << " of stars in catalog.\n";

			if(id >= cat->size()) { cout << "Index out of range.\n"; }
			else			   { cout << (*cat)[id] << "\n"; }
		}
	}











	void setbinfactor(int bf)
	{
		// in the future, dx should somehow be stored in volume-map fits files
		binfactor = bf;
		dx = binfactor*double(conf["volume.dx"]);		// in parsecs/pixel of the shrunk image
		// cout << "dx = " << dx << "\n";
	}

	mstar_record()
	: conf("mstars.conf")
	{
		setbinfactor(conf["binfactor"]);
	}

	/**
		\brief Load map of the volume covered (produced using volume.x) for specified \a run and
			\a ri bin

		\returns Total volume, in cubic parsecs
	*/
	double load_volume_map(int run, float ri, XImage &volume)
	{
		string volfile("column_maps/column_map." + str(run) + ".ri=" + str(ri, "%4.2f") + ".fits.gz");

		//
		// load and shrink the volume map
		//
		valarray<float> volume_big;
		io::fits::read(volume_big, volume.x(), volume.y(), volfile);
		integer_bin(volume.array, volume_big, volume.x(), volume.y(), binfactor, binfactor);

		// convert to cubic parsecs
		// caveat: this is now HARDCODED to match the output of volume.x, because volume.x
		//         stores the volume column maps in cubic pixels covered (so it depends on
		//         the value of dx, the linear size of the pixel (dx), as set in volume.x)
		// -- update: this is now a configuration file parameter. However, volume.x still
		// 	needs to be changed to read it's parameters from the same configuration file
		volume.array *= pow(conf["volume.dx"], 3.);

//		cerr << run << " " << volume(volume.x()/2 + 1, volume.y()/2 + 1) << " ";
//		cerr << volume(volume.x()/2 + 1, volume.y()/2 + 2) << "\n";

		if(conf["load_volume_map.save_fits"])
		{
			string volumefile("vol/vol." + str(run) + ".ri=" + str(ri, "%4.2f") + ".fits");
			io::fits::write(volume, volumefile);
			io::compress::gzip(volumefile);
		}

		// return the total volume (in cubic parsecs)
		return std::accumulate(&volume.array[0], &volume.array[volume.size()], 0.);
	}
#if 0
	/**
		\brief Load (r,z) coordinates for stars in run \a run, which are in \a ri_lo < \a ri < \a ri_hi range.
			Then bin the stars into array \a countsa, whose dimensions are \a x x \a y.

		\remark The counts image must allocated and size must be set.

		\returns Total number of stars binned
	*/
	double bin_stars_in_run(int run, float ri_lo, float ri_hi, XImage &counts)
	{
//		cerr << "Run " << run << "\n";
		//
		// load and bin the (r,z) plane data
		//
		string runfile("galactocentric/galactocentric." + str(run) + ".dat");
		gz_text_input_or_die (in, runfile);

		float rho, z, ri, D;
//		load(in, rho,0, z,2, D,3, ri,4);
		bind(in, rho,0, z,2, D,6, ri,3);

		counts.array = 0.;
		int &x = counts.x();
		int &y = counts.y();

		double Dmin, Dmax;
		mstars::distance_limits(Dmin, Dmax, ri_lo, ri_hi, conf["r_min"], conf["r_max"]);

		int yc = y / 2;
		int xc = x / 2;
		ticker tick(100);
		int i = 0;
//-		text_output_or_die(out, "out.txt");
		while(in.next())
		{
			if(ri_lo > ri || ri > ri_hi) { continue; }	// throw out everything not in our r-i range
			if(Dmin  >  D || D  > Dmax)  { continue; }	// throw out everything not in our volume limited sample

			// IMPORTANT NOTE:
			// Earth is in the central pixel of the image, bottom left corner of the pixel
			int px = xc + (int)floor((rho-Rg) / dx);
			int py = yc + (int)floor(z / dx);

#if 0
			if(!(0 <= px && px < x)) { cout << "[" << i << "] px = " << px << ", rho = " << rho << "\n"; continue; }
			if(!(0 <= py && py < x)) { cout << "[" << i << "] py = " << py << ", z = " << z << "\n"; continue; }
			ASSERT(0 <= px && px < x);
			ASSERT(0 <= py && py < y);
#endif
//			if(px + 1 == 236 && py + 1 == 856)
//			{
//				cerr << io::format(" ", true) << run << i << D << rho << z << ri << "\n";
//			}
//			if(counts(px,py) == 0) {
//				out << i << rho << z << D << ri << px << py << nl();
//			}
			counts(px, py) += 1.;
			tick.tick(); i++;
		}

		if(conf["bin_stars_in_run.save_fits"])
		{
			string countsfile("counts/rz." + str(run) + ".ri=" + str((ri_hi+ri_lo)/2, "%4.2f") + ".fits");
			io::fits::write(counts, countsfile);
			io::compress::gzip(countsfile);
		}

		return std::accumulate(&counts.array[0], &counts.array[counts.array.size()], 0.);
	}
#endif
	double load_counts_map(int run, float ri_lo, float ri_hi, XImage &counts)
	{
		string countsfile("counts/rz." + str(run) + ".ri=" + str((ri_hi+ri_lo)/2, "%4.2f") + ".fits");

		valarray<float> counts_big;
		io::fits::read(counts_big, counts.x(), counts.y(), countsfile);
		integer_bin(counts.array, counts_big, counts.x(), counts.y(), binfactor, binfactor);

		return std::accumulate(&counts.array[0], &counts.array[counts.array.size()], 0.);
	}

	/**
		\brief Calculates density for a given run and ri bin. Returns the result in \a density.
	*/
	void calculate_density_for_run(int run, float ri_lo, float ri_hi, XImage &density)
	{
		XImage volume, sigma, &counts = density;

		float ri = (ri_lo + ri_hi) / 2.;

		// load volume covered by the run
		double nvolume = load_volume_map(run, ri, volume);

		// load star counts
		counts.be_compatible_with(volume);
		sigma.be_compatible_with(volume); sigma.array = 0;
		double ncounts = load_counts_map(run, ri_lo, ri_hi, counts);

		// calculate density just for this run
		// note/warning: &density == &counts
		double dd = 0, cc = 0;
		FOR(0, counts.size()) {
			if(counts[i] > 0 && volume[i] > 0)
			{
				sigma[i] = 1. / sqrt(counts[i]);	// relative Poisson sigma
				density[i] = counts[i] / volume[i];
				dd += density[i];
				cc += 1;
			} else {
				// not measured
//				density[i] = -.1; // BUG: <-- this should be -1, 0 is just while testing
				density[i] = 0;
			}
		}
		
		cout << "mean density : " << ncounts/nvolume << "\n";
		cout << "mean img     : " << dd/cc << "\n";

		if(conf["calculate_density_for_run.save_density_fits"])
		{
			string denfile("den/den." + str(run) + ".ri=" + str(ri, "%4.2f") + ".fits");
			io::fits::write(density, denfile);
			io::compress::gzip(denfile);

			string sigmafile("den/sigma." + str(run) + ".ri=" + str(ri, "%4.2f") + ".fits");
			io::fits::write(sigma, sigmafile);
			io::compress::gzip(sigmafile);
		}
	}

	/**
		\brief Sum up all volumes covered by individual \a runs to produce total \a volume covered for bin \ri
	*/
	void total_volume(float ri, const set<int> &runs, XImage &volume)
	{
		// load volume covered by the run
		XImage run_volume;
		FOREACH(runs)
		{
			load_volume_map(*i, ri, run_volume);

			if(!volume.is_compatible_with(run_volume))
			{
				volume.be_compatible_with(run_volume);
				volume.array = 0.;
			}

			volume.array += run_volume.array;
		}
	}

	/**
		\brief Sum up all counts covered by individual \a runs to produce total \a volume covered for bin \ri
	*/
	void total_counts(float ri_lo, float ri_hi, const set<int> &runs, XImage &counts)
	{
		// load volume covered by the run
		XImage run_counts;

		run_counts.be_compatible_with(counts);
		counts.array = 0.;

		FOREACH(runs)
		{
			load_counts_map(*i, ri_lo, ri_hi, run_counts);

#if 0
			// eliminate edge pixels
			XImage cnt2;
			cnt2.be_compatible_with(run_counts);
			const int border = 1;
			FOREACH(run_counts)
			{
				cnt2(i) = *i;
				FORj(k, -border+i.x, border+1+i.x) FORj(l, -border+i.y, border+1+i.y)
				{
					if(0 > k || k >= run_counts.x()) continue;
					if(0 > l || l >= run_counts.y()) continue;
					if(run_counts(k, l) == 0) { cnt2(i) = 0; std::cerr << "#"; break; }
				}
			}
			run_counts.array = cnt2.array;
#endif

			counts.array += run_counts.array;
		}
	}

	/**
		\brief Measures density for a given \c ri bin, including measurements from runs given in \a runs.
	*/
	void calculate_density_for_ri_bin(float ri_lo, float ri_hi, const set<int> &runs, XImage &density, XImage &sigma)
	{
		double ri = (ri_lo + ri_hi) / 2;

		XImage volume, counts;

		total_volume(ri, runs, volume);					// calculate total volume covered

		counts.be_compatible_with(volume);					// counts
		total_counts(ri_lo, ri_hi, runs, counts);

		sigma.be_compatible_with(volume); sigma.array = 0;		// prepare sigma and density images
		density.be_compatible_with(volume); density.array = 0;

		// calculate cumulative maps and uncertancies, while binning at the same time
		// at this point, variables contain:
		// 	volume  - map of total volume covered in (r,z) space
		//	counts  - counts in each (r,z) bin
		//	density - will contain the density
		//	sigma   - will contain Poisson noise uncertancies
		//
		// Value of -1 in density will mean that pixel has not been measured
		// Note that density=0 is a valid value - at that pixel, no stars have been
		// detected (ergo, density is zero)

		//
		// our limit for "significant" volume sampling will be that the volume sampled was greater
		// than dV_limit cubic pc (TODO: make this a little less arbitrary)
		//
		double dV_limit = conf["volume.limit"];

		float zero = 0.;
		const float inf = 1./zero;

		// calculating density
		FOR(0, density.size())
		{
//			if(volume[i] > dV_limit)
			if(counts[i] > 0 && volume[i] > 0)
			{
				density[i] = counts[i] / volume[i];
				sigma[i] = 1. / sqrt(counts[i]);
			} else {
				// not measured
//				density[i] = -1;
				density[i] = 0;
				sigma[i] = inf;
			}
		}

		//
		// store results
		//
		if(conf["calculate_density_for_ri_bin.save_counts_fits"])
		{
			string countsfile("result/counts.ri=" + str(ri, "%4.2f") + ".fits");
			io::fits::write(counts.array, counts.x(), counts.y(), countsfile);
			io::compress::gzip(countsfile);
		}

		if(conf["calculate_density_for_ri_bin.save_volume_fits"])
		{
			string volfile("result/volume_map.ri=" + str(ri, "%4.2f") + ".fits");
			io::fits::write(volume.array, volume.x(), volume.y(), volfile);
			io::compress::gzip(volfile);
		}

		if(conf["calculate_density_for_ri_bin.save_density_fits"])
		{
			string denfile("result/den.ri=" + str(ri, "%4.2f") + ".fits");
			io::fits::write(density.array, density.x(), density.y(), denfile);
			io::compress::gzip(denfile);
		}

		if(conf["calculate_density_for_ri_bin.save_density_fits"])
		{
			string sigmafile("result/sigma.ri=" + str(ri, "%4.2f") + ".fits");
			io::fits::write(sigma.array, density.x(), density.y(), sigmafile);
			io::compress::gzip(sigmafile);
		}
	}

	void calculate_density(float ri_lo)
	{
		set<int> runs;
		loadRuns(runs);

		// bin bounds
		double ri_min = conf["ri_min"], ri_max = conf["ri_max"];
		int bins = conf["ri_bins"];
		double step = (ri_max - ri_min)/bins;
		double ri_hi = ri_lo + step;

		cout << "Processing bin " << ri_lo << " <= ri < " << ri_hi << "\n";

		XImage density, sigma;
		calculate_density_for_ri_bin(ri_lo, ri_hi, runs, density, sigma);
	}

	string filename(const string &dir, const string &prefix, float ri, int run = 0)
	{
		if(run == 0) {
			return string("./" + dir + "/" + prefix + ".ri=" + str(ri, "%4.2f") + ".fits.gz");
		} else {
			return string("./" + dir + "/" + prefix + "." + str(run) + ".ri=" + str(ri, "%4.2f") + ".fits.gz");
		}
	}

	void analyze(float ri)
	{
		XImage den, sigma, volume, counts;

		io::fits::read(den, filename("result", "den", ri) );
		io::fits::read(sigma, filename("result", "sigma", ri) );
		io::fits::read(volume, filename("result", "volume_map", ri) );
		io::fits::read(counts, filename("result", "counts", ri) );

		int &x = den.x();
		int &y = den.y();

		//
		// Fit data to a single exponential
		//

		// number of measurements (samples)
		int n = 0;
		FOREACH(den) { n += *i > 0; }

		text_output_or_die(out, io::format("analysis/density.ri=%4.2f.txt") << ri);

		gsl::vector Y(n), w(n);
		gsl::matrix X(n, 3);
		int k = 0;
		FOREACH(den)
		{
			double rho = Rg + dx*((i.x-x/2) + .5);
			double z   = dx*((i.y-y/2) + .5);
			
			if(*i <= 0) {
				out << Rg + dx*((i.x-x/2) + .5) <<  dx*((i.y-y/2) + .5) << *i << 0 << volume(i) << counts(i) << nl();
				continue;
			}

			Y(k) = log(*i);					// log density
			w(k) = *i / sqr(sigma(i));			// sigma of log density

			X(k, 0) = 1;					// constant a
			X(k, 1) = rho;					// galactocentric distance of the center of that pixel

			X(k, 2) = abs(z);					// abs z of that pixel

			out << X(k, 1) << z << *i << sigma(i) << volume(i) << counts(i) << nl();

			k++;
			ASSERT(k <= n);
		}
		ASSERT(k == n);

		gsl::vector c(3);
		gsl::matrix covariance(3, 3);
		double chisq;

		gsl::multifit (n, 3) . wlinear(X, w, Y, c, covariance, chisq);

		cout << "r-i bin:    " << ri << "\n";
		cout << "Samples:    " << k << "\n";
		cout << "Results:    l=" << -1./c(1) << "pc h=" << -1./c(2) << "pc\n";
		cout << "Chisq:      " << chisq << "\n";
		cout << "Chisq/dof:  " << chisq/n << "\n";
		cout << "Covariance:\n " << covariance << "\n";
	}

	void analyze2(float ri)
	{
		XImage den, sigma, volume, counts;

		io::fits::read(den, filename("result", "den", ri) );
		io::fits::read(sigma, filename("result", "sigma", ri) );
		io::fits::read(volume, filename("result", "volume_map", ri) );
		io::fits::read(counts, filename("result", "counts", ri) );

		int &x = den.x();
		int &y = den.y();

		//
		// Fit data to a single exponential
		//

		// number of measurements (samples)
		int n = 0;
		FOREACH(den) { n += *i > 0; }

		text_output_or_die(out, io::format("analysis/density.ri=%4.2f.txt") << ri);

		gsl::vector Y(n), w(n);
		gsl::matrix X(n, 3);
		int k = 0;
		FOREACH(den)
		{
			double rho = Rg + dx*((i.x-x/2) + .5);
			double z   = dx*((i.y-y/2) + .5);
			
			if(*i <= 0) {
				out << Rg + dx*((i.x-x/2) + .5) <<  dx*((i.y-y/2) + .5) << *i << 0 << volume(i) << counts(i) << nl();
				continue;
			}

			Y(k) = log(*i);					// log density
			w(k) = *i / sqr(sigma(i));			// sigma of log density

			X(k, 0) = 1;					// constant a
			X(k, 1) = rho;					// galactocentric distance of the center of that pixel

			X(k, 2) = abs(z);					// abs z of that pixel

			out << X(k, 1) << z << *i << sigma(i) << volume(i) << counts(i) << nl();

			k++;
			ASSERT(k <= n);
		}
		ASSERT(k == n);

		gsl::vector c(3);
		gsl::matrix covariance(3, 3);
		double chisq;

		gsl::multifit (n, 3) . wlinear(X, w, Y, c, covariance, chisq);

		cout << "r-i bin:    " << ri << "\n";
		cout << "Samples:    " << k << "\n";
		cout << "Results:    l=" << -1./c(1) << "pc h=" << -1./c(2) << "pc\n";
		cout << "Chisq:      " << chisq << "\n";
		cout << "Chisq/dof:  " << chisq/n << "\n";
		cout << "Covariance:\n " << covariance << "\n";
	}

	#include "tests.cpp"
	#include "vrml.cpp"
};

string shift(vector<string> &s, bool retempty = false)
{
	if(!s.size())
	{
		if(retempty) { return ""; }
		THROW(EAny, "Not enough arguments specified");
	}
	string x = s.front(); s.erase(s.begin());
	return x;
}


struct plx_gri_locus_filter : public filter_function
{
	double Dmin, Dmax;
	float ri0, ri1;
	plx_gri_locus paralax;

	plx_gri_locus_filter(float ri0_, float ri1_, float rmin, float rmax)
		: ri0(ri0_), ri1(ri1_)
	{
		paralax.distance_limits(Dmin, Dmax, ri0, ri1, rmin, rmax);
		cout << "Distance limits: " << Dmin << " " << Dmax << "\n";
	}

	virtual bool operator()(const sdss_star &s)
	{
#if 0
		FOR(1, 4) // this is now in schlegel.x, throw it out later
		{
			ASSERT(finite(s.mag[i]) && finite(s.magErr[i]));
			ASSERT(abs(s.magErr[i]) <= 0.2);
		}
#endif
		if(ri0  > s.ml_ri() || s.ml_ri() > ri1)   { return false; }	// throw out everything not in our r-i range
		if(Dmin > s.earth.D || s.earth.D > Dmax)  { return false; }	// throw out everything not in our volume limited sample
		return true;
	}
};

struct plx_gri_locus_filter_z : plx_gri_locus_filter
{
	double zmin, zmax;
	
	plx_gri_locus_filter_z(float ri0_, float ri1_, float rmin, float rmax, float zmin_, float zmax_)
		: plx_gri_locus_filter(ri0_, ri1_, rmin, rmax), zmin(zmin_), zmax(zmax_)
	{}

	virtual bool operator()(const sdss_star &s)
	{
		if(!plx_gri_locus_filter::operator()(s)) return false;

		if(zmin > s.earth.z || zmax < s.earth.z) return false;
		
		return true;
	}
};

#if 0

struct brightstars_z_filter : public brightstars_filter
{
	float z0, z1;

	brightstars_z_filter(float gr0_, float gr1_, float rmin, float rmax, float z0_, float z1_)
		: brightstars_filter(gr0_, gr1_, rmin, rmax), z0(z0_), z1(z1_)
	{
	}

	virtual bool operator()(const sdss_star &s)
	{
		if(gr0  > s.gr() || s.gr() > gr1) 		{ return false; }	// throw out everything not in our g-r range
		if(Dmin > s.earth.D || s.earth.D > Dmax)  { return false; }	// throw out everything not in our volume limited sample
		if(z0 > s.earth.z || s.earth.z > z1)	{ return false; } // throw out everything not in our z slice
//		if(s.earth.z < 0)	{ return false; } // throw out everything not in our z slice
		return true;
	}
};

#endif

void safe_divide(XImage &density, const XImage &counts, const XImage &volume)
{
	float zero = 0.;
	const float inf = 1./zero;

	ASSERT(volume.is_compatible_with(counts));
	density.ensure_compat(counts);

	// calculating density
	FOR(0, density.size())
	{
		if(counts[i] > 0 && volume[i] > 0)
		{
			density[i] = counts[i] / volume[i];
		} else {
			// not measured
			density[i] = 0; // -1;
		}
	}
}

void safe_divide(sparse_volume &density, const sparse_volume &counts, const sparse_volume &volume)
{
	float zero = 0.;
	const float inf = 1./zero;

	ASSERT(counts.dx == volume.dx);
	density.dx = counts.dx;
	density.r0 = counts.r0;
	density.r1 = counts.r1;

	// calculating density
	FOREACH(counts.volume)
	{
		sparse_volume::K p = (*i).first;
		float c = (*i).second;
		
		sparse_volume::volume_map::const_iterator v = volume.volume.find(p);
		if(v == volume.volume.end())
		{
			cerr << "Error: volume = 0 at " << p << " and N(p) = " << c << "\n";
			continue;
		}
		
		density.volume[p] = c / (*v).second;
	}
}


void write_sparse_volume(const sparse_volume &vm, const std::string &fn, const std::string &header_text = "")
{
	header h(header_text);
	h["dx"] = str(vm.dx);

	gz_binary_output_or_die(binf, fn);
	binf << h << vm;
}

void read_sparse_volume(sparse_volume &vm, const std::string &fn)
{
	header h;
	gz_binary_input_or_die(binf, fn);
	binf >> h >> vm;
}

int main(int argc, char **argv)
{
try
{
	mstar_record r;

	if(argc < 1) { THROW(EAny, "Insufficient number of arguments"); }
	string cmd = argv[1];

	vector<string> args;
	FORj(k, 2, argc) { args.push_back(argv[k]); }

	// ---------------------------

	if(cmd == "getUsedRuns") 	{ r.getUsedRuns(); }			/// 0.) Extract the list of runs from mdwarfs.dat file
	else if(cmd == "checkMasks")  	{ r.checkMasksAndBounds(); }		/// Sanity check (*)
	else if(cmd == "vrml")        	{ r.vrml_model(); }			/// Make a VRML model (*)
	else if(cmd == "starinfo")								/// Info about a particular star (interactive mode, type in star IDs)
	{
		int id, run;
		r.loadCatalog("catalogs");
		
		// make some catalog statistics
		while(true) {
			cout << "> ";
			std::string line;
			int a = -1, b = -1, c = -1, d = -1;
			cin >> line;
			int k = sscanf(line.c_str(), "%d %d %d %d", &a, &b, &c, &d);
			cout << k << " params (" << a << " " << b << " " << c << " " << d << ")\n";
			switch(k)
			{
			case 1: r.starinfo(-1, a); break; // global catalog ID (id)
			case 2: r.starinfo(a, b); break; // ID within a run (run, idx)
			case 3: r.starinfo(a, -1, b, c); break; // run, camcol, field
			case 4: r.starinfo(a, d, b, c); break; // run, camcol, field, objid
			}
		}
	}
	else if(cmd == "findstar")								/// Info about a particular star (interactive mode, type in star IDs)
	{
		int id, run, col, frame;
		r.loadCatalog("catalogs");
		
		// make some catalog statistics
		while(true) {
			cout << "> ";
			cin >> run >> col >> frame >> id;
			r.starinfo(run, id, col, frame);
		}
	}
	else if(cmd == "distlimits")								/// Info about a particular star (interactive mode, type in star IDs)
	{
		plx_gri_locus paralax;

		float ri0 = atof(shift(args).c_str());
		float ri1 = atof(shift(args).c_str());
		float rmin = atof(shift(args, true).c_str());
		float rmax = atof(shift(args, true).c_str());
		if(rmin == 0) rmin = 15.0;
		if(rmax == 0) rmax = 20.5;

		double Dmin, Dmax;
		paralax.distance_limits(Dmin, Dmax, ri0, ri1, rmin, rmax);
		
		cout << "# r = [" << rmin << ", " << rmax << "], r-i = [" << ri0 << ", " << ri1 << "]\n";
		cout << Dmin << " " << Dmax << "\n";
	}
	else if(cmd == "resample")
	{
		string inf = shift(args);
		string outf = shift(args);

		XImage in, out;
		io::fits::read(in, inf);

		resample(out, in, 50);

		io::fits::write(out, outf);
		
	}
	else if(cmd == "binRun")								/// Bin a run
	{
		set<int> runs;
		int run;
		string runfile;
		while(run = atoi((runfile = shift(args, true)).c_str()))
		{
			runs.insert(run);
		}
		if(!runs.size())
		{
			r.loadRuns(runs, runfile);
		}
		string suffix = runs.size() == 1 ? str(*runs.begin()) : "merged";

		float rmin = 15.0;
		float rmax = 20.5;
		float ri0 = 1.0;
		float ri1 = 1.1;
		float dx = 1;
		int x = 2*1200;

		int y = x;
		plx_gri_locus_filter volumelimit(ri0, ri1, rmin, rmax);
		cout << io::format("Limits: D_min = %.5f, D_max = %.5f\n") << volumelimit.Dmin << volumelimit.Dmax;

		XImage counts(x, y);
		int n = r.bin(runs, counts, dx, volumelimit);
		cout << "Total stars binned : " << n << "\n";

		std::string filename = io::format("counts/rz.%.3f-%.3f.%s.fits") << ri0 << ri1 << suffix;
		io::fits::write(counts, filename);
		io::compress::gzip(filename);

	}
	else if(cmd == "binxy")
	{
		float rmin = atof(shift(args).c_str());
		float rmax = atof(shift(args).c_str());
		float zmin = atof(shift(args).c_str());
		float zmax = atof(shift(args).c_str());
		float ri0 = atof(shift(args).c_str());
		float ri1 = atof(shift(args).c_str());
		float dx = atof(shift(args).c_str());
		int x = atoi(shift(args).c_str());

		set<int> runs;
		int run;
		string runfile;
		while(run = atoi((runfile = shift(args, true)).c_str()))
		{
			runs.insert(run);
		}
		if(!runs.size())
		{
			r.loadRuns(runs, runfile);
		}
		string suffix = runs.size() == 1 ? str(*runs.begin()) : "merged";

		r.loadCatalog("catalogs");

		int y = x;
		plx_gri_locus_filter_z volumelimit(ri0, ri1, rmin, rmax, zmin, zmax);
		cout << io::format("Limits: D_min = %.5f, D_max = %.5f\n") << volumelimit.Dmin << volumelimit.Dmax;

		XImage counts(x, y);
		int n = r.binxy(runs, counts, dx, volumelimit);
		cout << "Total stars binned : " << n << "\n";

		string filename = io::format("counts/xy.%.3f-%.3f.%05.0f-%05.0f.%s.fits") << ri0 << ri1 << zmin << zmax << suffix;
		io::fits::write(counts, filename);
		io::compress::gzip(filename);

	}
	else if(cmd == "denxy")
	{
		/*
		float rmin = 15.0;
		float rmax = 20.5;
		float zmin = 3500; //200;
		float zmax = 4000; //300;
		float ri0 = 0.1; //1.0;
		float ri1 = 0.15; //1.1;
		float dx = 10; //1;
		int b = 50;
		int x = 2*1200;
		*/

		float rmin = atof(shift(args).c_str());
		float rmax = atof(shift(args).c_str());
		float zmin = atof(shift(args).c_str());
		float zmax = atof(shift(args).c_str());
		float ri0 = atof(shift(args).c_str());
		float ri1 = atof(shift(args).c_str());
		float dx = atof(shift(args).c_str());
		int x = atoi(shift(args).c_str());
		int b = atoi(shift(args).c_str());

		XImage volume, counts, density;
		string suffix = shift(args, true);
		if(!suffix.size()) { suffix = "merged"; }

		string volumefn = io::format("column_maps/xy.%.3f-%.3f.%05.0f-%05.0f.%s.fits.gz") << ri0 << ri1 << zmin << zmax << suffix;
		io::fits::read(volume, volumefn);

		string countsfn = io::format("counts/xy.%.3f-%.3f.%05.0f-%05.0f.%s.fits.gz") << ri0 << ri1 << zmin << zmax << suffix;
		io::fits::read(counts, countsfn);

		resample(counts, b);
		resample(volume, b);

		safe_divide(density, counts, volume);

		// write a bunch of slice files for SM
		FOR(0, density.x())	// z slices
		{
			text_output_or_die(out, io::format("density/den.z_slice.%d.txt") << i);
			FORj(j, 0, density.y()) { out << (j - density.y()/2)*dx*b << " " << density(i, j) << " " << counts(i, j) << "\n"; }
		}

		FOR(0, density.y())	// r slices
		{
			text_output_or_die(out, io::format("density/den.r_slice.%d.txt") << i);
			FORj(j, 0, density.x()) { out << (j - density.x()/2)*dx*b << " " << density(j, i) << " " << counts(j, i) << "\n"; }
		}

		string densityfn = io::format("density/xy.density.%.3f-%.3f.%05.0f-%05.0f.%s.fits.gz") << ri0 << ri1 << zmin << zmax << suffix;
		io::fits::write(density, densityfn);
		string countsbfn = io::format("density/xy.counts.%.3f-%.3f.%05.0f-%05.0f.%s.fits.gz") << ri0 << ri1 << zmin << zmax << suffix;
		io::fits::write(counts, countsbfn);
		string volumebfn = io::format("density/xy.volume.%.3f-%.3f.%05.0f-%05.0f.%s.fits.gz") << ri0 << ri1 << zmin << zmax << suffix;
		io::fits::write(volume, volumebfn);
	}
	else if(cmd == "bin3D")								/// Bin a run
	{
		set<int> runs;
		int run;
		while(run = atoi(shift(args, true).c_str()))
		{
			runs.insert(run);
		}
		if(!runs.size())
		{
			r.loadRuns(runs);
		}
		string suffix = runs.size() == 1 ? str(*runs.begin()) : "merged";

		float rmin = 15.0;
		float rmax = 20.5;
		float ri0 = 1.0;
		float ri1 = 1.1;
		float dx = 1;
		int b = 50;
		int x = 2*1200;

		int y = x;
		plx_gri_locus_filter volumelimit(ri0, ri1, rmin, rmax);
		cout << io::format("Limits: D_min = %.5f, D_max = %.5f\n") << volumelimit.Dmin << volumelimit.Dmax;

		sparse_volume vm; sparse_volume_info vi;
		vm.dx = dx * b;
		vm.r0 = volumelimit.Dmin / vm.dx;
		vm.r1 = volumelimit.Dmax / vm.dx;
		int n = r.bin3D(runs, b, vm, vi, dx, volumelimit);
		cout << "Total stars binned : " << n << "\n";

		std::string filename = io::format("counts3d/counts.%.3f-%.3f.%s.bin") << ri0 << ri1 << suffix;
		header h("Binned star counts in geocentric galactic cartesian coordinates, in dx units");
		string sruns;
		FOREACH(runs) { sruns += str(*i) + " "; }
		h["runs"] = sruns;
		h["dx"] = str(vm.dx);
		h["hasinfo"] = 1;

		gz_binary_output_or_die(binf, filename);
		binf << h << vm << vi;
	}
	else if(cmd == "filterPlane")
	{
		set<int> runs;
		r.loadRuns(runs);
		FOREACH(runs)
		{
			int run = *i;
			
			sdss::RunGeometry geom;
			r.geomDB.getGeometry(run, geom);

			Radians len = geom.length();
			bool reject = false;
			for(double mu = 0; mu < len; mu += rad(.5))
			{
				Radians l, b;
				coordinates::gcsgal(geom.node, geom.inc, geom.muStart + mu, 0, l, b);
				if(abs(b) < rad(15)) { reject = true; break; }
			}
			if(!reject) { cout << run << "\n"; }
		}
	}
	else if(cmd == "dumpLocus")
	{
#if 0
		plx_gri_locus paralax;
		float ri = 1.1;
		float gi = 1.5;
		float riErr = .1, giErr = .1;
		float riMl = paralax.mlri(ri, gi, riErr, giErr);
		cerr << paralax.gr(ri) << " " << GSL_FN_EVAL(&paralax.mlri,ri) << "\n";
		cerr << riMl << "\n";

		exit(-1);
		/*
			Make some canonical diagrams.
		*/
	
		text_output_or_die(out, "mags.txt");
		ticker tick(1);
		FOREACHj(run, r.cat->runs)
		{
			auto_ptr<star_cat_slice> stars(r.cat->slice(*run));
			FOREACH(*stars)
			{
				sdss_star &s = *i;
				out << s.g << s.r << s.i << s.gErr << s.rErr << s.iErr << nl();
//				out << s.ri() << s.gi() << s.riErr() << s.giErr() << nl();
			}
			tick.tick();
		}

		text_output_or_die(out2, "locusfit.txt");
		for(float gi=-1; gi <= 3.3; gi += 0.1)
		{
			out2 << gi << paralax.ri_from_gi(gi) << nl();
		}

		text_output_or_die(out3, "locusscatter.txt");
		valarray<double> rms(0., 40), w(0., 40), xbar(0., 40);
		valarray<int> n(0, 40);
		FOREACHj(run, r.cat->runs)
		{
			auto_ptr<star_cat_slice> stars(r.cat->slice(*run));
			FOREACH(*stars)
			{
				sdss_star &s = *i;
				if(s.riErr() > 0.4 || s.giErr() > 0.4) { continue; }
				if(s.gi() < 0 || s.gi() > 3.5) { continue; }
				
				int bin = (int)(s.gi() / 0.1);
				double wt = 1./sqr(s.riErr()) + 1./sqr(s.giErr());
				rms[bin] += wt*sqr(s.ri() - paralax.ri_from_gi(s.gi()));
				xbar[bin] += wt*(s.ri() - paralax.ri_from_gi(s.gi()));
				w[bin] += wt;
				n[bin]++;
			}
		}
		rms = sqrt(rms / w);
		xbar /= w;
		out3 << "# g-r    rms    xbar    n" << nl();
		FOR(0, rms.size()) { out3 << i*0.1 << rms[i] << xbar[i] << n[i] << nl(); }
#endif
	}
	else if(cmd == "density")								/// 2.) Calculate density maps
	{
		float gr0 = 0.4;
		float gr1 = 0.5;
		float rmin = 15;
		float rmax = 20.5;
		float dx = 9.042877;
		int b = 20;

		XImage volume, counts, density;
		string volumefn = io::format("column_maps/cm-bright-%3.1f-%3.1f.merged.fits.gz") << gr0 << gr1;
		io::fits::read(volume, volumefn);
		string countsfn = io::format("counts/rz.%3.1f-%3.1f.fits.gz") << gr0 << gr1;
		io::fits::read(counts, countsfn);

		integer_bin(counts, b);
		integer_bin(volume, b);

		safe_divide(density, counts, volume);

		// write a bunch of slice files for SM
		FOR(0, density.x())	// z slices
		{
			text_output_or_die(out, io::format("density/den.z_slice.%d.txt") << i);
			FORj(j, 0, density.y()) { out << (j - density.y()/2)*dx*b << " " << density(i, j) << " " << counts(i, j) << "\n"; }
		}

		FOR(0, density.y())	// r slices
		{
			text_output_or_die(out, io::format("density/den.r_slice.%d.txt") << i);
			FORj(j, 0, density.x()) { out << (j - density.x()/2)*dx*b << " " << density(j, i) << " " << counts(j, i) << "\n"; }
		}

		// axis drawing
		int cx = density.x() / 2, cy = density.y() / 2;
		cx -= (int)ceil(8000. / dx / b);
		FOR(0, density.x()) { density(i, cy) = .2; }
		float step = 0;
		FOR(0, density.y()) {
			density(cx, i) = .2;
			if((i - cy) * dx * b >= step) { density(cx-1, i) = .2; step += 1000; }
		}

		string densityfn = io::format("density/density.%3.1f-%3.1f.fits.gz") << gr0 << gr1;
		io::fits::write(density, densityfn);
		string countsbfn = io::format("density/counts.%3.1f-%3.1f.fits.gz") << gr0 << gr1;
		io::fits::write(counts, countsbfn);
		string volumebfn = io::format("density/volume.%3.1f-%3.1f.fits.gz") << gr0 << gr1;
		io::fits::write(volume, volumebfn);
	}
	else if(cmd == "denMdw")								/// 2.) Calculate density maps
	{
		float rmin = 15.0;
		float rmax = 20.5;
		float ri0 = 1.0;
		float ri1 = 1.1;
		float dx = 1;
		int b = 50;
		int x = 2*1200;

		XImage volume, counts, density;
//		string volumefn = io::format("column_maps/rz.%.3f-%.3f.merged.fits.gz") << ri0 << ri1;
//		string volumefn = io::format("column_maps/rz.%.3f-%.3f.merged-orig.fits.gz") << ri0 << ri1;
//		string volumefn = io::format("column_maps/rz.1.0-1.1.745.fits.gz");

		string suffix = shift(args, true);
		if(!suffix.size()) { suffix = "merged"; }

		string volumefn = io::format("column_maps/rz.%.3f-%.3f.%s.fits.gz") << ri0 << ri1 << suffix;
		io::fits::read(volume, volumefn);

		string countsfn = io::format("counts/rz.%.3f-%.3f.%s.fits.gz") << ri0 << ri1 << suffix;
		io::fits::read(counts, countsfn);

		integer_bin(counts, b);
		integer_bin(volume, b);

		safe_divide(density, counts, volume);

		// write a bunch of slice files for SM
		FOR(0, density.x())	// z slices
		{
			text_output_or_die(out, io::format("density/den.z_slice.%d.txt") << i);
			FORj(j, 0, density.y()) { out << (j - density.y()/2)*dx*b << " " << density(i, j) << " " << counts(i, j) << "\n"; }
		}

		FOR(0, density.y())	// r slices
		{
			text_output_or_die(out, io::format("density/den.r_slice.%d.txt") << i);
			FORj(j, 0, density.x()) { out << (j - density.x()/2)*dx*b << " " << density(j, i) << " " << counts(j, i) << "\n"; }
		}

		string densityfn = io::format("density/rz.density.%.3f-%.3f.%s.fits.gz") << ri0 << ri1 << suffix;
		io::fits::write(density, densityfn);
		string countsbfn = io::format("density/rz.counts.%.3f-%.3f.%s.fits.gz") << ri0 << ri1 << suffix;
		io::fits::write(counts, countsbfn);
		string volumebfn = io::format("density/rz.volume.%.3f-%.3f.%s.fits.gz") << ri0 << ri1 << suffix;
		io::fits::write(volume, volumebfn);
	}
	else if(cmd == "den3D")								/// 2.) Calculate density maps
	{
		float rmin = 15.0;
		float rmax = 20.5;
		float ri0 = 1.0;
		float ri1 = 1.1;
		float dx = 1;
		int b = 50;
		int x = 2*1200;

		sparse_volume volume, counts, density;

		string suffix = shift(args, true);
		if(!suffix.size()) { suffix = "merged"; }

		string volumefn = io::format("binned_volumes/ri-%.1f-%.1f.%s.bin") << ri0 << ri1 << suffix;
		read_sparse_volume(volume, volumefn);
		string countsfn = io::format("counts3d/counts.%.3f-%.3f.%s.bin") << ri0 << ri1 << suffix;
		read_sparse_volume(counts, countsfn);

		safe_divide(density, counts, volume);

		string densityfn = io::format("density3d/density.%.3f-%.3f.%s.bin") << ri0 << ri1 << suffix;
		write_sparse_volume(density, densityfn);
		string countsbfn = io::format("density3d/counts.%.3f-%.3f.%s.bin") << ri0 << ri1 << suffix;
		write_sparse_volume(counts, countsbfn);
		string volumebfn = io::format("density3d/volume.%.3f-%.3f.%s.bin") << ri0 << ri1 << suffix;
		write_sparse_volume(volume, volumebfn);
	}
	else if(cmd == "analyze")								/// 3.) Analyze results
	{
		FOR(0, 10) {
			r.analyze(1.05 + double(i)/10.);
		}
	}
	else if(cmd == "runDensity")
	{
		XImage density;
		r.calculate_density_for_run(atoi(shift(args).c_str()), 1., 1.1, density);
	}
	else
	{
		THROW(EAny, "No command specified");
	}
} catch (peyton::exceptions::EAny &e)
{
	e.print();
} catch (...)
{
	cout << "Uncaught exception!\n";
}

}
#else
#include <iostream>

int main(int argc, char **argv)
{
	std::cerr << "This exe has not been compiled because of the lack of CCfits library.\n";
	return -1;
}
#endif // HAVE_LIBCCFITS
