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

#ifndef _textloader_h__
#define _textloader_h__

#include "analysis.h"

class text_file_streamer : public catalog_streamer
{
protected:
	itextstream in;
	int id;
	bool astrometry_only;
public:
	text_file_streamer(std::istream &datafile, bool astrometry_only = false);
	bool next()
	{
		if(in.next()) {
//			cerr << "R: " << s.id << " " << s.run << " " << s.ra << " " << s.dec << "\n";
//			s.print(std::cerr); exit(-1);
			s.id = id;
			++id; return true;
		}
		return false;
	}
};

text_file_streamer::text_file_streamer(std::istream &datafile, bool astrometry_only_)
	: in(datafile), id(0), astrometry_only(astrometry_only_)
{
	s.r = 0.;
	s.objid = 0;
	s.field = 0;

	bind(in, s.ra, 1-1);
	bind(in, s.dec, 2-1);
	bind(in, s.run, 3-1);

	if(!astrometry_only)
	{
		bind(in, s.col, 5-1);
		bind(in, s.field, 6-1);
		bind(in, s.objid, 7-1);

		bind(in, s.Ar, 10-1);

		bind(in, s.u, 11-1);
		bind(in, s.g, 12-1);
		bind(in, s.r, 13-1);
		bind(in, s.i, 14-1);
		bind(in, s.z, 15-1);

		bind(in, s.uErr, 16-1);
		bind(in, s.gErr, 17-1);
		bind(in, s.rErr, 18-1);
		bind(in, s.iErr, 19-1);
		bind(in, s.zErr, 20-1);
        }
}

#endif
