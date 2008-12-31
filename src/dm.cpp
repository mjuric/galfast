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

#include "container.h"
#include "textloader.h"

#include <iostream>

#include <astro/system/options.h>
#include <astro/useall.h>

using namespace std;

// defined in dm_common.cpp
void make_object_catalog(
	const std::string &uniqObjectCat, // "dm_unique_stars.dmm" [output]
	const std::string &uniqObsvCat, // "dm_starmags.dmm" [output]

	const std::string &matchGroups, // "match_groups.lut" (produced by unique.x)

	const std::string &select,  // name of the temporary object catalog (indexed by uniq ID) (DMM file) [in]
	const std::string &selindex,// name of the fitsID -> uniq ID map (DMM file) [in]
	const std::string &runindex,// name of the sloanID -> uniqID (DMM file) [in]

	const std::string &stagesummary_fn // name of text file for output summary
);

// defined in dm_common.cpp
void obs_info(
	int fitsID,
	const std::string &uniqObjectCat, // "dm_unique_stars.dmm" [output]
	const std::string &uniqObsvCat, // "dm_starmags.dmm" [output]

	const std::string &selindex, // name of the fitsID -> obsv ID map (DMM file) [in]
	const std::string &runindex // name of the sloanID -> uniqID (DMM file) [in]
);

void reprocess_driver();

main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version, "$Id: dm.cpp,v 1.11 2006/07/13 23:27:30 mjuric Exp $");
	Options opts(
		argv[0],
		"Create the unique object catalog.",
		version, Authorship::majuric
	);

	std::string
		inputCatalog = "/home/mjuric/projects/galaxy/workspace/catalogs/test.txt",
		stage = "help",
		stagesummary_fn("dm-summary.txt"),
		group_dir("."),
		lutfn("match_groups.lut"),
		indexfn("match_index.dmm"),
		cattype("fits"),
		tmp_dir("."),
		tmp_prefix("uniq_lut"),
		object_dir("."),
		object_prefix("uniq"),
		run_offset_map_fn("uniq_run_index.map"),
		fitsfilepat("fits/calibObj-%06d-?-star.fits.gz");

	std::vector<std::string> stages;
	container(stages) = "mklookup", "mkobject", "mkrunidx", "asciidump", "info";
	std::string allstages = join(", ", stages);
	stages.push_back("all");
	stages.push_back("help");
	std::set<std::string> cattypes;
	container(cattypes) = "fits", "text";
	int info_fitsid = 0;

	opts.argument("inputCatalog", "Input catalog file"); //, Option::optional);

	opts.option("stage").bind(stage).desc("Execute a stage of matching algorithm. Can be one of " + allstages + ". Special value 'all' causes all stages to be executed in sequence, value 'none' executes nothing and prints a summary of parameters.");
	opts.option("type").bind(cattype).desc("Type of the input catalog (valid choices are 'fits' and 'text'");

	opts.option("group-dir").bind(group_dir).desc("Directory where the group table and index will be stored");
	opts.option("group-table-file").bind(lutfn).desc("Name of the group table file");
	opts.option("group-index-file").bind(indexfn).desc("Name of the group index file");

	opts.option("tmp-dir").bind(tmp_dir).desc("Directory where the temporary lookup file, used while creating the object catalog, will be stored.");
	opts.option("tmp-lookup-prefix").bind(tmp_prefix).desc("Prefix of temporary lookup files.");

	opts.option("fits-files-pattern").bind(fitsfilepat).desc("sscanf-formatted shell pattern to resolve path to files of a run (used with --type=fits)");

	opts.option("object-dir").bind(object_dir).desc("Directory where the output object catalog files will be stored.");
	opts.option("object-catalog-prefix").bind(object_prefix).desc("Prefix of output object catalog files.");

	opts.option("info-fitsid").bind(info_fitsid).desc("fitsID of the observation for which to get the info (use with --stage=info).");

	
	try {
		parse_options(opts, argc, argv);
		if(opts.found("inputCatalog")) { inputCatalog = opts["inputCatalog"]; }

		if(find(stages.begin(), stages.end(), stage) == stages.end())
			{ THROW(EOptions, "Argument to option --stage must be one of " + join(", ", stages) + "."); }
		if(!cattypes.count(cattype))
			{ THROW(EOptions, "Argument to option --cattype must be one of " + join(", ", cattypes) + "."); }
		if(inputCatalog != "-" && !file_exists(inputCatalog))
			{ THROW(EOptions, "Input catalog file '" + inputCatalog + "' does not exist or is inaccessible"); }
	} catch(EOptions &e) {
		opts.usage(cerr);
		e.print();
		exit(-1);
	}

#if 0
	cout << "stage = " << stage << "\n";
	cout << "inputCatalog = " << inputCatalog << "\n";
	cout << "group_dir = " << group_dir << "\n";
	cout << "lutfn = " << lutfn  << "\n";
	cout << "indexfn = " << indexfn << "\n";
	cout << "cattype = " << cattype << "\n";
	cout << "tmp_dir = " << tmp_dir << "\n";
	cout << "tmp_prefix = " << tmp_prefix << "\n";
#endif


	gsl_set_error_handler_off ();

	std::string select_fn   = tmp_dir + '/' + tmp_prefix + "_tmp_obsv_cat.dmm";
	std::string selindex_fn = tmp_dir + '/' + tmp_prefix + "_tmp_catid2selidx.dmm";
	std::string selindices_fn = tmp_dir + '/' + tmp_prefix + "obsv_indices.dmm";

	std::string object_fn      = object_dir + '/' + object_prefix + "_objects.dmm";
	std::string object_obsv_fn = object_dir + '/' + object_prefix + "_observations.dmm";

	std::string match_lut_path = group_dir + '/' + lutfn;

	if(stage == "mklookup" || stage == "all")
	{
		if(cattype == "fits")
		{
		#ifdef HAVE_PKG_CCfits
			cerr << "\n[1/3] Creating matching lookup table\n";
			std::set<int> runs;
			loadRuns(runs, inputCatalog);
			fits_set_streamer fitscat(fitsfilepat, runs);
			makelookup(fitscat, select_fn, selindex_fn, selindices_fn);
		#else
			cerr << "FITS support not compiled in!\n";
			abort();
		#endif
		} else {
			ifstream inf(inputCatalog.c_str());
			text_file_streamer cat(inputCatalog == "-" ? cin : inf, false);
			makelookup(cat, select_fn, selindex_fn, selindices_fn);
		}
	}

	if(stage == "mkobject" || stage == "all")
	{
		make_object_catalog(object_fn, object_obsv_fn,
			match_lut_path,
			select_fn, selindex_fn, selindices_fn,
			stagesummary_fn);
	}

	if(stage == "mkrunidx" || stage == "all")
	{
		cout << "Run offset map file: " << run_offset_map_fn << "\n";
		std::ofstream out(run_offset_map_fn.c_str());
		make_run_index_offset_map(out, selindices_fn);
	}

	if(stage == "info")
	{
		obs_info(info_fitsid,
			object_fn,
			object_obsv_fn,
			selindex_fn,
			selindices_fn
		);
	}

	if(stage == "help")
	{
		opts.help(cerr);
		return -1;
	}
#if 0
	reprocess_driver();
#endif

} catch (peyton::exceptions::EAny &e)
{
	e.print();
} catch (...)
{
	std::cout << "Uncaught exception!\n";
}

}
