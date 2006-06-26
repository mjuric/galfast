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

#include "paralax.h"

#include <astro/useall.h>
using namespace std;

//plx_gri_locus plx;

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

	opts.argument("toColor", "Which color is the output (then, the other one must be specified on input). Valid values are 'ri' and 'gr'");
	opts.argument("color");

	try {
		opts.parse(argc, argv);
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}
#if 0
	/////// Now start with the real busyness //////
	float color = opts["color"];
	if(opts["toColor"].vstring() == "gr")
	{
		cout << "For gr=" << color << ", ri=" << gr_to_ri(color) << "\n";
	}
	if(opts["toColor"].vstring() == "ri")
	{
		cout << "For ri=" << color << ", gr=" << plx::gr(color) << "\n";
	}
#else
	ASSERT(0); // this code is unfinished
#endif
}
catch(EAny &e)
{
	e.print();
}
}
