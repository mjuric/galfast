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

#include <set>
#include <utility>
#include <fstream>

#include <astro/useall.h>
using namespace std;

// in binner_tng.cpp
int merge_maps(const std::string &outputPrefix, const set<int> &runs, pair<float, float> r, pair<float, float> ri, const std::string &uniqMapFn = "");

int main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version);

	Options opts(
		"\
	Bla bla",
		version,
		Authorship::majuric
	);

	opts.argument("r0");
	opts.argument("r1");
	opts.argument("ri0");
	opts.argument("ri1");
	opts.argument("run", "Run number or a file which contains runs to flatten");
	opts.argument("prefix");
	opts.argument("uniqMapFn", "Filename which contains map of unique volume", Option::optional);

	try {
		opts.parse(argc, argv);
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}

	/////// Now start with the real busyness //////
	
	// parse the command line
	set<int> runs;
	int run = opts["run"];
	if(run == 0)
	{
		text_input_or_die (in, (std::string)opts["run"]);
		load(in, runs, 0);
	} else {
		runs.insert(run);
	}

	pair<float, float>
		ri (opts["ri0"], opts["ri1"]),
		r  (opts["r0"],  opts["r1"]);

	merge_maps(opts["prefix"], runs,
	    make_pair((float)opts["r0"], (float)opts["r1"]), 
	    make_pair((float)opts["ri0"], (float)opts["ri1"]),
	    opts["uniqMapFn"]
	    );
}
catch(EAny &e)
{
	e.print();
}
}
