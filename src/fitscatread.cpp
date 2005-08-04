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

#define NO_SDSS_STAR_CAT

#include "dm.h" 
#include "fitsloader.h"
#include "analysis.h"
#include "xcat.h"
#include "paralax.h"

#include <astro/constants.h>
#include <astro/coordinates.h>
#include <astro/sdss/rungeometry.h>
#include <astro/sdss/photometry.h>
#include <astro/math.h>
#include <astro/math/vector.h>
#include <astro/io/format.h>
#include <astro/system/memorymap.h>
#include <astro/system/fs.h>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics_double.h>

#include <sys/types.h>
#include <sys/stat.h>

#include <cmath>
#include <map>
#include <fstream>
#include <valarray>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <time.h>

#include <astro/useall.h>
using namespace std;

int main(int argc, char **argv)
{
	int run = 3919;
	
	sdss_star s; valarray<int> f1(5), f2(5);
	FOR(1, 7)
	{
		string file = io::format("data/schlegel/calibObj-%06d-%d-star.fits.gz") << run << i;
		fits_loader fl(file, s, f1, f2);
		while(fl.next())
		{
			preprocess_sdss_star(s);
			if(0.1 < s.gr() && s.gr() < 1.2 &&
			   18.5 < s.r && s.r < 19)
			{
				cout << deg(s.ra) << " " << deg(s.dec) << " " << s.col + 1;
				FORj(j,0,5)
				{
					cout << " " << s.mag[j] << " " << s.magErr[j];
				}
				cout << "\n";
			}
		}
	}
}
