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

#ifndef common_fitsloader_h__
#define common_fitsloader_h__

#include "config.h"

#ifdef HAVE_PKG_CCfits
 
#include <valarray>
#include <memory>
#include <CCfits>
#include <set>
#include <astro/system/fs.h>

#include "xcat.h"

struct fits_loader
{
	sdss_star &s;
	std::valarray<int> &f1, &f2;

	int rows;
	bool justSniff;
protected:
	std::auto_ptr<CCfits::FITS> in;

	int at;

	std::valarray<int> run, camcol, field, objid;
	std::valarray<double> ra, dec;
	std::valarray<float> rowc, colc;
	std::valarray<float> extinction, psfflux, psfflux_ivar;

	std::valarray<int> flags, flags2;

	// a special function to read in the arrays, because CCfits readArrays method is
	// *EXTREMELY* slow
	void read_vector(CCfits::ExtHDU &table, const std::string &column, std::valarray<float> &array);
	void read_vector_int(CCfits::ExtHDU &table, const std::string &column, std::valarray<int> &array);
public:
	fits_loader(const string &fitsfile, sdss_star &s_, std::valarray<int> &f1_, std::valarray<int> &f2_, bool justSniff = false);

	bool next();
};

class fits_set_streamer : public catalog_streamer
{
protected:
	std::auto_ptr<fits_loader> cur;
	std::set<int> runs;
	std::set<int>::iterator at;

	std::auto_ptr<peyton::system::dir> inputfiles;
	peyton::system::dir::iterator curfile;
	bool first;

	int fitsid;

	void nextfile();
public:
	fits_set_streamer(const std::set<int> &runs_);
	bool next();
	//int catalogId() const { return fitsid; }	/// fits ID of the last record loaded
};

#endif // HAVE_LIBCCFITS

#endif
