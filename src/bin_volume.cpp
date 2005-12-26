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

#include "analysis.h"

#include <astro/system/options.h>
#include <astro/system/log.h>
#include <astro/math.h>

#include <set>
#include <utility>
#include <fstream>

#include <astro/useall.h>
using namespace std;

// in raytrace.cpp
int bin_volumes(const std::set<int> &runs, double dx, int ndx, pair<float, float> r, pair<float, float> ri, bool justinfo);
int bin_volumes2(const std::string &outfile, const std::string &volfn, double dx, int ndx, pair<float, float> r, pair<float, float> ri);
int bin_gc_cut(
	const std::string &outfile, const std::string &volfn, 
	double dx, int ndx,
	pair<float, float> r, pair<float, float> ri,
	Radians node, Radians inc, pair<Radians, Radians> nu);
int bin_plane_cut(
	const std::string &outfile, const std::string &volfn, 
	double dx, int ndx,
	pair<float, float> r, pair<float, float> ri,
	const std::string &coordsys,
	double d1, pair<double, double> p1,
	double d2, pair<double, double> p2,
	double d3, pair<double, double> p3,
	double d0, pair<double, double> p0,
	double delta, bool earth_on_x_axis);
int bin_cylindrical(
	const std::string &outfile, const std::string &volfn, 
	double dx, int ndx,
	pair<float, float> r, pair<float, float> ri,
	Radians phi0
);

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
	opts.argument("r0");
	opts.argument("r1");
	opts.argument("ri0");
	opts.argument("ri1");
	opts.argument("run", "Run number or text file with list of runs, or binned unique volume file, if <uniqMapFn> is specified");
	opts.argument("dx", "Scale of a pixel of sampled volume [pc]");
	opts.argument("ndx", "How many sampled volume pixels to bin into one binned volume pixel [even integer]");
	opts.argument("uniqMapFn", "Filename which contains map of unique volume", Option::optional);

	// add any options your program might need. eg:
	// opts.option("meshFactor", "meshFactor", 0, "--", Option::required, "4", "Resolution decrease between radial steps");
	bool justinfo = false;
	{ // setup arguments and options
		using namespace peyton::system::opt;

		opts.option("just-info", binding(justinfo), arg_none, desc("Just print the information about what will be done, do no actual work. Works only for some subcommands."));
	}

	opts.option_arg("type", "Type of binning to perform ('gc', 'plane' and 'cyl' are defined)");
	opts.option_arg("coordsys", "Coordinate system used for x0...xn ('equ' is the default, 'gal' and 'galcart' are available)");

	opts.option_arg("node", "Great circle node [deg] (for use with --type=gc)");
	opts.option_arg("inc", "Great circle inclination [deg] (for use with --type=gc)");
	opts.option_arg("nu", "Great circle nu range [deg] - everything in [-nu, nu) will be selected (for use with --type=gc)");
	
	opts.option_arg("x1", "Coordinates of the 1st point in the plane (the object), in 'dist,ra,dec' format - eg. --x1=\"5443,342.4221,11.3223\" [pc,deg,deg] (for use with --type=plane)");
	opts.option_arg("x2", "Coordinates of the 2nd point in the plane (for use with --type=plane)");
	opts.option_arg("x3", "Coordinates of the 3rd point in the plane (for use with --type=plane)");
	opts.option_arg("x0", "Coordinates of the origin point (for use with --type=plane). This point will become (0,0,0) point in plane coordinates.");
	opts.option_arg("earthonaxis", "If nonzero, the x axis will be constructed by projecting the x1->x2 direction to vector pointing from the galactic center to Earth. (for use with --type=plane)", "0");
	opts.option_arg("delta", "Halfthickness of the plane to bin (for use with --type=plane)");

	opts.option_arg("phi0", "Zero point of galactocentric phi angle, in degrees (for use with --type=cyl)", "0");

	try {
		opts.parse(argc, argv);
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}

	/////// Start your application code here
	pair<double, double> x1, x2, x3, x0;
	double d1, d2, d3, d0;
	if(opts.found("type"))
	{
		if(opts["type"] == "gc")
		{
			ASSERT(opts.found("uniqMapFn"));
			ASSERT(opts.found("node"));
			ASSERT(opts.found("inc"));
			ASSERT(opts.found("nu"));
		}
		else if(opts["type"] == "plane")
		{
			ASSERT(opts.found("x1"));
			ASSERT(sscanf(opts["x1"].c_str(), "%lf,%lf,%lf", &d1, &x1.first, &x1.second) == 3);
			ASSERT(opts.found("x2"));
			ASSERT(sscanf(opts["x2"].c_str(), "%lf,%lf,%lf", &d2, &x2.first, &x2.second) == 3);
			ASSERT(opts.found("x3"));
			ASSERT(sscanf(opts["x3"].c_str(), "%lf,%lf,%lf", &d3, &x3.first, &x3.second) == 3);
			ASSERT(opts.found("x0"));
			ASSERT(sscanf(opts["x0"].c_str(), "%lf,%lf,%lf", &d0, &x0.first, &x0.second) == 3);
			ASSERT(opts.found("delta"));
		}
		else if(opts["type"] == "cyl")
		{
		}
		else
		{
			die("Unknown option " + opts["type"] + "\n\n" + opts.usage());
		}
	}

	if(!opts.found("uniqMapFn"))
	{
		cerr << "Binning volume\n";

		set<int> runs;
		int run = opts["run"];
		if(run == 0)
		{
			string fn = opts["run"];
			text_input_or_die (in, fn);
			load(in, runs, 0);
		} else {
			runs.insert(run);
		}

		bin_volumes(runs, opts["dx"], opts["ndx"],
			make_pair((float)opts["r0"], (float)opts["r1"]), 
			make_pair((float)opts["ri0"], (float)opts["ri1"]),
			justinfo
		);
	}
	else
	{
		if(opts["type"] == "gc")
		{
			cerr << "Binning GC volume to [" << opts["uniqMapFn"] << "]\n";
			Radians nu = rad(opts["nu"]);
			bin_gc_cut(opts["uniqMapFn"], opts["run"], opts["dx"], opts["ndx"],
				make_pair((float)opts["r0"], (float)opts["r1"]), 
				make_pair((float)opts["ri0"], (float)opts["ri1"]),
				rad(opts["node"]), rad(opts["inc"]), make_pair(-nu, nu)
			);
		}
		else if(opts["type"] == "plane")
		{
			cerr << "Binning planar volume to [" << opts["uniqMapFn"] << "]\n";
			bin_plane_cut(opts["uniqMapFn"], opts["run"], opts["dx"], opts["ndx"],
				make_pair((float)opts["r0"], (float)opts["r1"]), 
				make_pair((float)opts["ri0"], (float)opts["ri1"]),
				opts["coordsys"], d1, x1, d2, x2, d3, x3, d0, x0, opts["delta"],
				opts["earthonaxis"]
			);
		}
		else if(opts["type"] == "cyl")
		{
			cerr << "Binning cylindrical volume to [" << opts["uniqMapFn"] << "]\n";
			bin_cylindrical(opts["uniqMapFn"], opts["run"], opts["dx"], opts["ndx"],
				make_pair((float)opts["r0"], (float)opts["r1"]), 
				make_pair((float)opts["ri0"], (float)opts["ri1"]),
				rad(opts["phi0"])
			);
		}
		else
		{
			cerr << "Binning unique volume to [" << opts["uniqMapFn"] << "]\n";
			bin_volumes2(opts["uniqMapFn"], opts["run"], opts["dx"], opts["ndx"],
				make_pair((float)opts["r0"], (float)opts["r1"]), 
				make_pair((float)opts["ri0"], (float)opts["ri1"])
			);
		}
	}
}
catch(EAny &e)
{
	e.print();
}
}
