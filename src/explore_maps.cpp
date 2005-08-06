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
#include <astro/io/fits.h>

#include "dm.h"
#include "ximage.h"
#include "projections.h"

#include <set>
#include <utility>
#include <fstream>

#include <astro/useall.h>
using namespace std;

void lambertw(double dx, Radians l0, Radians b0, Radians l1, Radians b1)
{
	DMMArray<mobject> arr("dm_unique_stars.dmm");
	int w = (int)(4. / dx)+1;
	int h = (int)(4. / dx)+1;
	cerr << "Image dimension: " << w << " x " << h << "\n";
	XImage img(w,h);

	lambert map(rad(0), rad(90));
#if 1
	ticker tick(10000);
	FOREACH(DMMArray<mobject>::iterator, arr)
	{
		tick.tick();
		mobject &m = *i;

		Radians l, b;
		coordinates::equgal(rad(m.ra), rad(m.dec), l, b);

		if(!coordinates::inBox(l, b, l0, b0, l1, b1)) continue;
		if(m.mag[1] > 21.0) continue;
		if(!between(m.gr(), 0.2, 0.3)) continue;
		if(!between(m.mag[1], 20., 21.)) continue;

		double x, y;
		map.convert(l, b, x, y);
		x = (int)((x + 2) / dx);
		y = (int)((y + 2) / dx);

		img((int)x, (int)y)++;
	}
#endif	
	text_output_or_die(out, "counts.txt");
	FOREACH(XImage::iterator, img)
	{
		if(*i == 0) continue;
		out << (0.5 + i.x)*dx - 2 << (0.5 + i.y)*dx - 2 << *i << nl();
	}
	out_stream.close();

#ifdef HAVE_LIBCCFITS
	fits::write(img, "counts.fits");
#else
	ASSERT(0) { cerr << "libCCfits support not compiled in.\n"; }
#endif
}

void work(double dx, double dy, Radians l0, Radians b0, Radians l1, Radians b1)
{
	DMMArray<mobject> arr("dm_unique_stars.dmm");
	int w = (int)((l1 - l0) / dx)+1;
	int h = (int)((b1 - b0) / dy)+1;
	cerr << "Image dimension: " << w << " x " << h << "\n";
	XImage img(w,h);

	ticker tick(10000);
	FOREACH(DMMArray<mobject>::iterator, arr)
	{
		tick.tick();
		mobject &m = *i;

		Radians l, b;
		coordinates::equgal(rad(m.ra), rad(m.dec), l, b);

		if(!coordinates::inBox(l, b, l0, b0, l1, b1)) continue;
		if(m.mag[1] > 21.0) continue;
		if(!between(m.gr(), 0.2, 0.3)) continue;
		if(!between(m.mag[1], 20., 21.)) continue;

		int x = (int)((l - l0) / dx);
		int y = (int)((b - b0) / dy);
		img(x, y)++;
	}
#ifdef HAVE_LIBCCFITS
	fits::write(img, "counts.fits");
#else
	ASSERT(0) { cerr << "libCCfits support not compiled in.\n"; }
#endif
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
#if 0
	opts.argument("r0");
	opts.argument("r1");
	opts.argument("ri0");
	opts.argument("ri1");
	opts.argument("run", "Run number or text file with list of runs");
	opts.argument("dx", "Scale of a pixel of sampled volume [pc]");
	opts.argument("ndx", "How many sampled volume pixels to bin into one binned volume pixel [even integer]");
#endif
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
#if 0
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
#endif
//	work(rad(1), rad(1), rad(0), rad(0), rad(359.999), rad(80));
	lambertw(.0333333, rad(0), rad(0), rad(359.999), rad(80));
}
catch(EAny &e)
{
	e.print();
}
}
