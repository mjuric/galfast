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

#include "gpu.h"

#include <astro/math.h>
#include <astro/system/options.h>

#include <fstream>
#include <boost/shared_ptr.hpp>

#include <astro/useall.h>
using namespace std;

///////////////////////////////////

void resample_texture(const std::string &outfn, const std::string &texfn, float2 crange[3], int npix[3], bool deproject, Radians l0, Radians b0);
void postprocess_catalog(const std::string &conffn, const std::string &input, const std::string &output, std::set<std::string> modules);

int main(int argc, char **argv)
{
try
{
	std::string argv0 = argv[0];
	VERSION_DATETIME(version, "$Id: simulate.cpp,v 1.19 2007/04/15 12:09:52 mjuric Exp $");
	std::string progdesc = "galfast.x, a mock star catalog simulator.";

	std::string cmd, input, output;
	std::map<std::string, boost::shared_ptr<Options> > sopts;
	Options opts(argv[0], progdesc, version, Authorship::majuric);
	opts.argument("cmd").bind(cmd).desc(
		"What to make. Can be one of:\n"
		"postprocess - \tpostprocess the mock catalog (e.g., derive photometry, add instrumental errors, etc.)\n"
		"       util - \tvarious utilities\n"
					   );
	opts.stop_after_final_arg = true;
	opts.prolog = "For detailed help on a particular subcommand, do `galfast.x <cmd> -h'";
	opts.add_standard_options();

	Radians dx = 4.;

	std::vector<std::string> modules;
	std::string infile(""), outfile("");
	sopts["postprocess"].reset(new Options(argv0 + " postprocess", progdesc + " Apply postprocessing steps to catalog sources.", version, Authorship::majuric));
	sopts["postprocess"]->argument("conf").bind(input).desc("Postprocessing (\"postprocess.conf\") configuration file");
	sopts["postprocess"]->argument("modules").bind(modules).optional().gobble().desc("List of postprocessing module configuration files (or module names)");
	sopts["postprocess"]->option("i").addname("input").bind(infile).param_required().desc("Input catalog file");
	sopts["postprocess"]->option("o").addname("output").bind(outfile).param_required().desc("Output catalog file");
	sopts["postprocess"]->add_standard_options();

	std::string util_cmd;
	sopts["util"].reset(new Options(argv0 + " util", progdesc + " Utilities subcommand.", version, Authorship::majuric));
	sopts["util"]->argument("utility").bind(util_cmd).desc(
		"Run a utility. Can be one of:\n"
		"  cudaquery - \tquery available cuda devices\n"
		" resample3d - \tresample a 3D FITS file, store output to text table\n"
//		"   footplot - \tmake a PostScript plot of the footprints\n"
	);
	sopts["util"]->stop_after_final_arg = true;
	sopts["util"]->prolog = "For detailed help on a particular subcommand, do `galfast.x util <cmd> -h'";
	sopts["util"]->add_standard_options();
	std::map<std::string, boost::shared_ptr<Options> > uopts;
	
	//
	// util sub-subcommands
	//
	uopts["cudaquery"].reset(new Options(argv0 + " util cudaquery", progdesc + " Query available CUDA devices.", version, Authorship::majuric));

	float2 crange[3] = { make_float2(0, 0), make_float2(0, 0), make_float2(0, 0) };
	int npix[3] = { 0 };
	bool deproject = false;
	double l0 = -100., b0 = -100.;
	uopts["resample3d"].reset(new Options(argv0 + " util resample3d", progdesc + " Resample 3D data cube.", version, Authorship::majuric));
	uopts["resample3d"]->argument("input").bind(input).desc("Input FITS file");
	uopts["resample3d"]->argument("output").bind(output).desc("Output text file");
	uopts["resample3d"]->option("x0").bind(crange[0].x).param_required().desc("Minimum x coordinate");
	uopts["resample3d"]->option("x1").bind(crange[0].y).param_required().desc("Maximum x coordinate");
	uopts["resample3d"]->option("y0").bind(crange[1].x).param_required().desc("Minimum y coordinate");
	uopts["resample3d"]->option("y1").bind(crange[1].y).param_required().desc("Maximum y coordinate");
	uopts["resample3d"]->option("z0").bind(crange[2].x).param_required().desc("Minimum z coordinate");
	uopts["resample3d"]->option("z1").bind(crange[2].y).param_required().desc("Maximum z coordinate");
	uopts["resample3d"]->option("nx").bind(npix[0]).param_required().desc("Number of output pixels in x dimension");
	uopts["resample3d"]->option("ny").bind(npix[1]).param_required().desc("Number of output pixels in y dimension");
	uopts["resample3d"]->option("nz").bind(npix[2]).param_required().desc("Number of output pixels in z dimension");
	uopts["resample3d"]->option("d").bind(deproject).value("true").desc("Deproject to (l,b). If this option is active, the x coordinate is l, and y is b. Projection pole can be changed using --l0 and --b0 options.");
	uopts["resample3d"]->option("l0").bind(l0).param_required().desc("Longitude of projection pole, degrees");
	uopts["resample3d"]->option("b0").bind(b0).param_required().desc("Latitude of projection pole, degrees");

#if 0 // TODO: Implement this
	std::vector<std::string> footprint_confs;
	std::string outputps = "foot.ps";
	uopts["footplot"].reset(new Options(argv0 + " util footplot", progdesc + " Make a PostScript plot of covered footprint.", version, Authorship::majuric));
	uopts["footplot"]->argument("footprints").bind(footprint_confs).gobble().desc("Footprint configuration files, or a single skygen.conf configuration file.");
	uopts["footplot"]->option("o").addname("output").bind(outputps).param_required().desc("Output PostScript file");
#endif

	//
	// Parse
	//
	Options::option_list optlist;
	parse_options(optlist, opts, argc, argv);
	if(sopts.count(cmd))
	{
		parse_options(*sopts[cmd], optlist);
		if(cmd == "util")
		{
			cmd = cmd + " " + util_cmd;
			parse_options(*uopts[util_cmd], optlist);
		}
	}
	else
	{
		ostringstream ss;
		ss << "Unrecognized subcommand `" << cmd << "'";
		print_options_error(ss.str(), opts);
		return -1;
	}

	/////// Application starts here

	if(!cuda_init())
	{
		MLOG(verb1) << "Error initializing GPU acceleration. Aborting.";
		return -1;
	}

	if(cmd == "util cudaquery")
	{
		// The real test was that we successfully passed the cuda_init() step above.
		return 0;
	}
	if(cmd == "util resample3d")
	{
		resample_texture(output, input, crange, npix, deproject, rad(l0), rad(b0));
		return 0;
	}

	ifstream in(input.c_str());
	if(!in) { THROW(EFile, "Error accessing " + input + "."); }

	if(cmd == "postprocess")
	{
		// ./galfast.x postprocess postprocess.conf [--infile=sky.cat.txt] [--outfile=sky.obs.txt] [module1.conf [module2.conf]....]
		std::set<std::string> mset;
		mset.insert(modules.begin(), modules.end());
		postprocess_catalog(input, infile, outfile, mset);
	}
	else
	{
		THROW(ENotImplemented, "Should not get to here. I'm really confused. Aborting.");
	}
	return 0;
}
catch(EAny &e)
{
	e.print();
	return -1;
}
}
