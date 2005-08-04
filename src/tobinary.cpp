
//--- libpeyton includes. Add any libpeyton includes _before_ including
//--- astro/useall.h

#include "xcat.h"
#include "paralax.h"
#include "analysis.h"

#include <astro/system/options.h>
#include <astro/io/format.h>
#include <astro/io/gzstream/fstream.h>
#include <astro/coordinates.h>
#include <astro/sdss/rungeometry.h>
#include <astro/system/fs.h>
#include <astro/useall.h>
#include <sstream>

#include <CCfits>

#include <iostream>
#include <fstream>
#include <valarray>

using namespace std;
using namespace CCfits;

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

#include "phObjc.h"

bool filter(sdss_star &s, const valarray<int> &flags, const valarray<int> &flags2)
{
	//
	// Flags cuts on magnitudes
	//
	const int f1req = OBJECT1_BRIGHT | OBJECT1_SATUR;
	const int f2req = OBJECT2_DEBLENDED_AS_MOVING;

#if 1
	FOR(1, 4)	// check flags for g, r, i magnitudes
	{
		if(flags[i]  & f1req) return false;
		if(flags2[i] & f2req) return false;

		if(!finite(s.mag[i]) || !finite(s.magErr[i])) { return false; }	// magnitude data and error validity cut
		if(abs(s.magErr[i]) > 0.2) return false;
	}
#endif

	//
	// Apply magnitude cuts
	//
	if(s.r > 22) { return false; }

	return true;
}

void importCatalog(const string &input, const string &output)
{
	sdss_star s;
	gz_text_input_or_die(in, input);

	bind(in, s.run, 0, s.col, 1, s.field, 2);
	bind(in, s.ra, 5, s.dec, 6, s.l, 7, s.b, 8);
	bind(in, s.Ar, 9);
	bind(in, s.u, 10, s.g, 11, s.r, 12, s.i, 13, s.z, 14);
	bind(in, s.uErr, 15, s.gErr, 16, s.rErr, 17, s.iErr, 18, s.zErr, 19);

	s.objid = -1;
	s.Nobservations = -1;
	s.primaryobservation = true;
	
	valarray<int> flags(0, 5), flags2(0, 5);

	while(in.next())
	{
		vector<sdss_star> stars;
		stars.reserve(1400000);

		int run = s.run;

		//
		// Loading
		//
		ticker tick(io::format("Importing run %d") << run, 1000);
		plx_gri_locus paralax;
		int rejected = 0, distrejected = 0, n = 0;
		
		do
		{
			n++;
			preprocess_sdss_star(s);

			if(!filter(s, flags, flags2)) { continue; }
			if(!boundsCheck(s)) { rejected++; continue; }
			if(!s.calculate_distances(paralax)) { distrejected++; continue; }

			stars.push_back(s);
			tick.tick();
		} while(in.next() && run == s.run);

		tick.close();
		cout << "Stars out of bounds: " << rejected << "\n";
		cout << "Stars rejected because of bad locus fit: " << distrejected << "\n";
		cout << "Stars accepted: " << stars.size() << "/" << n << io::format(" (%.2f\%)\n") << 100.0*stars.size()/n;

		//
		// Sorting
		//
		sort(stars.begin(), stars.end(), sdss_star_less());

		// find out how many runs are there, and where each run begins (i.e., create an index)
		valarray<int> index(0, 7);
		FOR(0, stars.size())
		{
			stars[i].id = i;
			index[stars[i].col + 1] = max(index[stars[i].col + 1], i);
		}
		index += 1;
		index[0] -= 1;

		cerr << "Run " << run << " import complete.\n\n";
		ostringstream columns;
		FOR(0, 7) { columns << io::format(" %7d") << index[i]; }
		columns << "\n";

		// store
		header h("SDSS stars catalog in binary format");
		h.data["run"] = str(run);
		h.data["size"] = str(stars.size());
		h.data["columns"] = columns.str();
		binary_output_or_die(out, io::format("%s/run-%d.bin") << output << run);
		out << h << index << stars;
	}
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
	opts.argument("input", "Input mdwarf.dat file");
	opts.argument("output", "Output directory");

	// add any options your program might need. eg:
	// opts.option("meshFactor", "meshFactor", 0, "--", Option::required, "4", "Resolution decrease between radial steps");

	try {
		opts.parse(argc, argv);
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}

	/////// Start your application code here
	gsl_set_error_handler_off ();

	string input = opts["input"];
	string output = opts["output"];

	importCatalog(input, output);
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
