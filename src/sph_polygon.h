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

#ifndef sph_polygon_h__
#define sph_polygon_h__

#include <list>

/**
	sph_contour -- A simple (nonintersecting) contour on a sphere

	You should not construct this class directly, but obtain a reference
	to it using sph_polygon::add_contour().
*/
class sph_contour
{
public:
	std::vector<gpc_vertex> c;	// contour vertices
	gpc_vertex p;			// a point inside the contour (NOTE: this used to be here to determine the what is the "inside"; soon switching to CCW ordering of vertices for this purpose)
	peyton::Radians dx;		// resolution at which this contour should be polygonized when projected

public:

	sph_contour(bool allsky = false)
	{
		reset(allsky);
	}

	sph_contour &operator =(const sph_contour &a)
	{
		c = a.c;
		p = a.p;
		dx = a.dx;
		return *this;
	}

	void reset(bool allsky = false)
	{
		dx = 0.;		// try to "guess" the right resolution
		p.x = p.y = 0;
		if(allsky)
		{
			c.push_back(p);	// add a single point to the contour (doesn't matter what this point actually is)
		}
	}

	bool allsky() const { return c.size() == 1; }

public:	// projected gpc_polygon conversion routines

	// Project and resample the contour onto projection proj, sampling the edges along geodesics with resolution dx
	gpc_polygon project(const peyton::math::lambert &proj) const;
};

/**
	sph_polygon -- An abstraction of a spherical polygon

	The polygon can consist of many contours, which will be assembled
	(in order of their appearance in the 'contours' list) using gpc_op
	operator. I.e., some of the contours may be holes, etc...
*/
class sph_polygon
{
public:
	std::list<std::pair<gpc_op, sph_contour> > contours;
public:
	// Project the polygon onto the requested projection
	gpc_polygon project(const peyton::math::lambert &proj) const;

	sph_contour &add_contour(gpc_op op = GPC_UNION)
	{
		contours.push_back(std::make_pair(op, sph_contour()));
		return contours.back().second;
	}
};

#endif // sph_polygon_h__
