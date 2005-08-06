
//--- libpeyton includes. Add any libpeyton includes _before_ including
//--- astro/useall.h

#include "config.h"

#ifdef HAVE_LIBCCFITS

#include "xcat.h"
#include "paralax.h"
#include "analysis.h"
#include "fitsloader.h"

#include <astro/system/options.h>
#include <astro/io/format.h>
#include <astro/coordinates.h>
#include <astro/sdss/rungeometry.h>
#include <astro/sdss/photometry.h>
#include <astro/system/fs.h>
#include <astro/useall.h>
#include <sstream>

#include <CCfits>

#include <iostream>
#include <fstream>
#include <valarray>

using namespace std;
using namespace CCfits;

#ifndef DEBUGMODE
	#define DEBUGMODE 0
#endif

/////////////

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

#include "binarystream.h"

#if 0
template<typename T>
ibinarystream & operator >>(ibinarystream &in, const mmopen_struct<T> &mmc)
{
	mmopen_struct<T> &mm = const_cast<mmopen_struct<T> &>(mmc);
	if(!in.filename.size()) { THROW(EIOException, "Input binary stream is not tied to a file (cannot memory map something that is not a file)"); }
	if(mm.size == -1) { in >> mm.size; }

	in >> mmalign();

	int offset = in.f.tellg();
	mm.open(in.filename, offset);
	in.f.seekg(sizeof(T) * mm.size, ios::cur);
	return in;
}

/////////////

struct mmalign
{
};

obinarystream &operator <<(obinarystream &out, const mmalign &mmp)
{
	int at = out.f.tellp();
	at += sizeof(int);

	int offs = at % MemoryMap::pagesize;
	if(offs) { offs = MemoryMap::pagesize - offs; }

	out << offs;
	out.f.seekp(offs, ios::cur);

	return out;
};

ibinarystream &operator >>(ibinarystream &in, const mmalign &mmp)
{
	int offs;
	in >> offs;
	in.f.ignore(offs);
	return in;
};

#endif

/////////////

sdss::RunGeometryDB db;
std::map<int, sdss::Mask> masks;
std::map<int, sdss::RunGeometry> geoms;

bool boundsCheck(sdss_star &s)
{
	if(!masks.count(s.run))
	{
		masks[s.run] = geoms[s.run] = db.getGeometry(s.run);
	}
	sdss::RunGeometry &geom = geoms[s.run];
	Radians mu, nu;
	Coordinates::equgcs(geom.node, geom.inc, s.ra, s.dec, mu, nu);
	return masks[s.run].contains(mu, nu);
}


#if 0
enum
{
	PEAKCENTER = 1 << 5,
	NOTCHECKED = 1 << 19,
	CR = 1 << 12,
	BINNED1 = 1 << 28,
	BRIGHT = 1 << 1,
	SATURATED = 1 << 18,
	EDGE = 1 << 2,
	BLENDED = 1 << 3,
	NODEBLEND = 1 << 6,
	NOPROFILE = 1 << 7,
};

enum
{
	DEBLENDED_AS_MOVING = 1 << 0,
	DEBLEND_NOPEAK = 1 << 14,
	PSF_FLUX_INTERP = 1 << 15,
	BAD_COUNTS_ERROR = 1 << 8,
	INTERP_CENTER = 1 << 12,
};
#endif

#include "phObjc.h"

bool filter(sdss_star &s, const valarray<int> &flags, const valarray<int> &flags2)
{
	//
	// Flags cuts on magnitudes
	//
	const int f1req = OBJECT1_BRIGHT | OBJECT1_SATUR;
	const int f2req = OBJECT2_DEBLENDED_AS_MOVING;

	FOR(1, 4)	// check flags for g, r, i magnitudes
	{
		if(flags[i]  & f1req) return false;
		if(flags2[i] & f2req) return false;

		if(!finite(s.mag[i]) || !finite(s.magErr[i])) { return false; }	// magnitude data and error validity cut
		if(abs(s.magErr[i]) > 0.2) return false;
	}


	//
	// Apply magnitude cuts
	//
	if(s.r > 22) { return false; }

	return true;
}

int importCatalog(int run, const string &output, int id)
{
	sdss_star s;
	valarray<int> flags1(5), flags2(5);
	s.Nobservations = -1;
	s.primaryobservation = true;

	vector<sdss_star> stars;
	stars.reserve(1400000);

	FITS::setVerboseMode(true);

	//
	// Loading
	//
	ticker tick(io::format("Importing run %d") << run, 1);
	plx_gri_locus paralax;
	int rejected = 0, distrejected = 0;

	string pattern = io::format("data/schlegel/calibObj-%06d-?-star.fits.gz") << run;
	dir inputfiles(pattern);
	int n = 0;
	FOREACH(dir::iterator, inputfiles)
	{
		cerr << "\t" << *i << "\n";
		fits_loader in(*i, s, flags1, flags2);
		while(in.next())
		{
			n++;
			preprocess_sdss_star(s);

			if(!filter(s, flags1, flags2)) { continue; }
			if(!boundsCheck(s)) { rejected++; continue; }
			if(!s.calculate_distances(paralax)) { distrejected++; continue; }

			stars.push_back(s);
		}
//		tick.tick();
	}
//	tick.close();
	cout << "Stars out of bounds: " << rejected << "\n";
	cout << "Stars rejected because of bad locus fit: " << distrejected << "\n";
	cout << "Stars accepted: " << stars.size() << "/" << n << io::format(" (%.2f\%)\n") << 100.0*stars.size()/n;

	//
	// Sorting
	//
	sort(stars.begin(), stars.end(), sdss_star_less());

	// find out how many runs are there, and where each run begins (i.e., create an index _of this run_ (zero is the beginning of the run)
	valarray<int> index(0, 7);
	FOR(0, stars.size())
	{
		stars[i].id = id; id++;

		int &ind = index[stars[i].col + 1];
		ind = max(ind, i);
	}
	index += 1;
	index[0] -= 1;

	ostringstream columns;
	FOR(0, 7) { columns << io::format(" %7d") << index[i]; }
	cerr << "Run " << run << " import complete [" << columns.str() << " ]\n\n";
	columns << "\n";

	// store
	header h("SDSS stars catalog in binary format");
	h.data["run"] = str(run);
	h.data["size"] = str(stars.size());
	h.data["columns"] = columns.str();
	binary_output_or_die(out, io::format("%s/run-%d.bin") << output << run);
	out << h << index;

	// store page aligned array of stars
	out << stars.size();
	MemoryMap::pagesizealign(out.f);
	FOREACH(vector<sdss_star>::iterator, stars) { out << *i; }

	return stars.size();
}

int main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version);

	Options opts(
		"This program has not been described",
		version,
		Authorship::majuric
	);

	//# add any arguments your program needs. eg:
	opts.argument("run", "Run number or text file with list of runs");
	opts.argument("output", "Output directory");

	// add any options your program might need. eg:
	// opts.option("meshFactor", "meshFactor", 0, "--", Option::required, "4", "Resolution decrease between radial steps");

	try {
#if !DEBUGMODE
		opts.parse(argc, argv);
#endif
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}

	/////// Start your application code here
	gsl_set_error_handler_off ();

	string input = DEBUGMODE ? (string)"data/schlegel/runs.txt" : (string)opts["run"];
	string output = DEBUGMODE ? (string)"catalogs" : (string)opts["output"];

	set<int> runs;
	int run = atoi(input.c_str());
	if(run == 0)
	{
		text_input_or_die (in, input);
		load(in, runs, 0);
	} else {
		runs.insert(run);
	}

	int id = 0;
	FOREACH(set<int>::iterator, runs)
	{
		id += importCatalog(*i, output, id);
	}
	cout << "Total stars in catalog: " << id << "\n";
}
catch(EAny &e)
{
	e.print();
}
catch(FitsException &e)
{
	cerr << "FITS exception!" << "\n";
}
}

#else
#include <iostream>

int main(int argc, char **argv)
{
	std::cerr << "This exe has not been compiled because of the lack of CCfits library.\n";
	return -1;
}
#endif
