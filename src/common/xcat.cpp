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

#include "xcat.h"
#include "binarystream.h"
#include "analysis.h"

#include <astro/util.h>

#include <fstream>
#include <valarray>
#include <iostream>

#include <astro/useall.h>
using namespace std;

/////////////////

sdss_star_cat::sdss_star_cat(std::string dir_)
	: dir(dir_)
{
	text_input_or_die (in, dir + "/runs.txt");
	load(in, runs, 0);

	// add all files to memory map
	valarray<int> cols; header h; int size; int n = 0;
	FOREACH(runs)
	{
		cerr << *i << " ";

		// load headers
		std::string fn(dir + "/run-" + str(*i) + ".bin");
		binary_input_or_die(in, fn);
		in >> h >> cols;

		// add memory map file
		in >> size;
		MemoryMap::pagesizealign(in.f);
		dmm.addblock(size, in.f.tellg(), fn);

		// store subset
		runsubsets[*i] = subset(begin() + n, begin() + n + size);
		n += size;
	}

	// set the memory map parameters
	setmaxwindows(40);
	setwindowsize(1024*1024 / 5);

	cerr << "\n";
}

sdss_star_cat::subset sdss_star_cat::getrun(int run)
{
	return runsubsets[run];
}

