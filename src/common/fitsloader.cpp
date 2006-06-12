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

#ifdef HAVE_PKG_CCfits

#include "fitsloader.h"
#include <astro/sdss/photometry.h>
#include <astro/io/format.h>

#include <astro/useall.h>
using namespace std;
using namespace CCfits;
 
void fits_loader::read_vector(ExtHDU &table, const std::string &column, valarray<float> &array)
{
	int status = 0; int index;
	fits_get_colnum(table.fitsPointer(), CASEINSEN, (char *)column.c_str(), &index, &status);
	int veclen = 5;
	array.resize(rows*veclen);
	fits_read_col(table.fitsPointer(), TFLOAT, index, 1, 1, rows*veclen, NULL, &array[0], NULL, &status);
	fits_report_error(stdout, status);
	//cerr << array[veclen*(rows-1) + 4] << "\n";
	//exit(-1);
}
void fits_loader::read_vector_int(ExtHDU &table, const std::string &column, valarray<int> &array)
{
	int status = 0; int index;
	fits_get_colnum(table.fitsPointer(), CASEINSEN, (char *)column.c_str(), &index, &status);
	int veclen = 5;
	array.resize(rows*veclen);
	fits_read_col(table.fitsPointer(), TLONG, index, 1, 1, rows*veclen, NULL, &array[0], NULL, &status);
	fits_report_error(stdout, status);
	//cerr << array[veclen*(rows-1) + 4] << "\n";
	//exit(-1);
}
	
fits_loader::fits_loader(const string &fitsfile, sdss_star &s_, valarray<int> &f1_, valarray<int> &f2_, bool justSniff_)
	: s(s_), at(0), extinction(5), psfflux(5), psfflux_ivar(5), f1(f1_), f2(f2_), justSniff(justSniff_),
	in(new FITS(fitsfile, Read, 1, false))
{
	ExtHDU &table = in->extension(1);		// get the first extension
	rows = table.rows();

	if(justSniff) return;
	
	// load the data
	table.column("RUN").read(run, 0, rows);
	table.column("CAMCOL").read(camcol, 0, rows);
	table.column("FIELD").read(field, 0, rows);
	table.column("ID").read(objid, 0, rows);
	table.column("RA").read(ra, 0, rows);
	table.column("DEC").read(dec, 0, rows);

	read_vector(table, "EXTINCTION", extinction);
	read_vector(table, "PSFFLUX", psfflux);
	read_vector(table, "PSFFLUX_IVAR", psfflux_ivar);
	read_vector(table, "COLC", colc);
	
	read_vector_int(table, "FLAGS", flags);
	read_vector_int(table, "FLAGS2", flags2);
//	table.column("NCHILD").read(nchild, 0, rows);
}

bool fits_loader::next()
{
	if(at >= rows) return false;

	if(!justSniff)
	{
		s.run = run[at]; s.col = camcol[at]; s.field = field[at]; s.objid = objid[at];
		s.ra = ra[at]; s.dec = dec[at];
		s.colc = colc[5*at+2];
		s.Ar = extinction[5*at+2];
	
		FOR(0, 5)
		{
			float f = psfflux[5*at+i];
			float fErr = 1./sqrt(psfflux_ivar[5*at+i]);
	
			phot::luptitude(f, fErr, i, s.mag[i], s.magErr[i], 1e+9);
	
			f1[i] = flags[5*at+i];
			f2[i] = flags2[5*at+i];
		}
	}

	++at;
	return true;
}


///////////

fits_set_streamer::fits_set_streamer(const std::string &filepattern, const std::set<int> &runs_)
	: filepat(filepattern), runs(runs_), at(NULL), first(true), fitsid(-1)
{
	FITS::setVerboseMode(true);
}

void fits_set_streamer::nextfile()
{
	// advance to next run
	string pattern = io::format(filepat) << *at;
	inputfiles.reset(new dir(pattern));
	curfile = inputfiles->begin();
	ASSERT(curfile != inputfiles->end());

	cur.reset(new fits_loader(*curfile, s, f1, f2));
}

bool fits_set_streamer::next()
{
	if(at == runs.end()) return false;

	if(first) { first = false; at = runs.begin(); nextfile(); }

	bool ret = cur->next();
	if(ret) { ++fitsid; s.id = fitsid; return true; }

	curfile++;
	if(curfile == inputfiles->end())
	{
		++at;
		if(at == runs.end()) return false;
		nextfile();
	} else {
		cur.reset(new fits_loader(*curfile, s, f1, f2));
	}

	return next();
}

#endif // HAVE_LIBCCFITS
