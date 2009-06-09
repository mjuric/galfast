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

// Interface to GPC library, with binary input/output goodies

#ifndef gpc_cpp__h__
#define gpc_cpp__h__

#include <map>

extern "C"
{
	#include "../gpc/gpc.h"
}

#include <astro/io/binarystream.h>
#include <astro/math.h>
#include "projections.h"

BLESS_POD(gpc_vertex);

// binary serialization routines for gpc_polygons and aux. data structures
BOSTREAM2(const gpc_vertex_list& p);
BISTREAM2(gpc_vertex_list& p);
BOSTREAM2(const gpc_polygon& p);
BISTREAM2(gpc_polygon& p);

void poly_bounding_box(double &x0, double &x1, double &y0, double &y1, const gpc_polygon &p);
gpc_polygon poly_rect(double x0, double x1, double y0, double y1);
double contour_area(const gpc_vertex_list &c);
double polygon_area(const gpc_polygon &p);
bool in_contour(const gpc_vertex &t, const gpc_vertex_list &vl);
bool in_polygon(const gpc_vertex &t, const gpc_polygon &p);

// polygon split into little rectangular pieces
struct partitioned_skymap
{
	peyton::Radians dx; 			// lambert grid angular resolution
	peyton::Radians x0, x1, y0, y1;		// lambert survey footprint bounding box

	struct pixel_t {
		gpc_polygon poly;
		double coveredArea;	// Area within the pixel covered by the polygon
		double pixelArea;	// Area of the entire rectangular lambert pixel
	};
	
	std::map<std::pair<int, int>, pixel_t> skymap;	// a map of rectangular sections of the sky, for fast is-point-in-survey-area lookup

	bool in(peyton::Radians x, peyton::Radians y);
};
BOSTREAM2(const partitioned_skymap::pixel_t &m);
BISTREAM2(partitioned_skymap::pixel_t &m);

BOSTREAM2(const partitioned_skymap &m);
BISTREAM2(partitioned_skymap &m);

void xgpc_write(const std::string &fn, const gpc_polygon &sky, const peyton::math::lambert &proj);
gpc_polygon xgpc_read(const std::string &fn, peyton::math::lambert &proj, bool *hadProjection = NULL);

#endif
