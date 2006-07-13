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

// in raytrace.cpp
int merge_volumes(const std::string &outputPrefix, const std::set<int> &runs);

int main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version, "$Id: merge_volumes.cpp,v 1.2 2006/07/13 23:27:30 mjuric Exp $");

	Options opts(
		argv[0],
		"\
	Bla bla",
		version,
		Authorship::majuric
	);

	opts.argument("run", "Run number or a file which contains runs whose volumes to merge");
	opts.argument("prefix", "Output merged volume file");

	parse_options(opts, argc, argv);

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

	merge_volumes(opts["prefix"], runs);
}
catch(EAny &e)
{
	e.print();
}
}
