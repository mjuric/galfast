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

#include "gpc_cpp.h"
#include <cmath>

using namespace std;

// binary serialization routines for gpc_polygons and aux. data structures
BOSTREAM2(const gpc_vertex_list& p)
{
	out << p.num_vertices;
	FOR(0, p.num_vertices)
	{
		out << p.vertex[i];
	}
	return out;
}

BISTREAM2(gpc_vertex_list& p)
{
	in >> p.num_vertices;
	p.vertex = (gpc_vertex*)malloc(p.num_vertices*sizeof(gpc_vertex));
	FOR(0, p.num_vertices)
	{
		in >> p.vertex[i];
	}
	return in;
}

BOSTREAM2(const gpc_polygon& p)
{
	out << p.num_contours;
	FOR(0, p.num_contours)
	{
		out << p.hole[i];
		out << p.contour[i];
	}
	return out;
}

BISTREAM2(gpc_polygon& p)
{
	in >> p.num_contours;
	p.contour = (gpc_vertex_list*)malloc(p.num_contours*sizeof(gpc_vertex_list));
	p.hole = (int*)malloc(p.num_contours*sizeof(int));
	FOR(0, p.num_contours)
	{
		in >> p.hole[i];
		in >> p.contour[i];
	}
	return in;
}

void poly_bounding_box(double &x0, double &x1, double &y0, double &y1,
	const gpc_polygon &p)
{
	// find minimum and maximum vertex
	x0 = x1 = p.contour[0].vertex[0].x;
	y0 = y1 = p.contour[0].vertex[0].y;
	FOR(0, p.num_contours)
	{
		const int n = p.contour[i].num_vertices;
		FORj(j, 0, n)
		{
			const double x = p.contour[i].vertex[j].x;
			const double y = p.contour[i].vertex[j].y;
			x0 = std::min(x0, x);
			x1 = std::max(x1, x);
			y0 = std::min(y0, y);
			y1 = std::max(y1, y);
		}
	}
}

// create a gpc_polygon rectangle
gpc_polygon poly_rect(double x0, double x1, double y0, double y1)
{
	static gpc_vertex v[4];
	static gpc_vertex_list vl = {4, v};
	static gpc_polygon rect = {1, NULL, &vl};

	v[0].x = x0; v[0].y = y0;
	v[1].x = x1; v[1].y = y0;
	v[2].x = x1; v[2].y = y1;
	v[3].x = x0; v[3].y = y1;

	return rect;
}

// calculate and return the area of a gpc_polygon
double polygon_area(const gpc_polygon &p)
{
	double A = 0;
	FOR(0, p.num_contours)
	{
		int n = p.contour[i].num_vertices;
		gpc_vertex *v = p.contour[i].vertex;
		double cA = 0;
		FORj(j, 0, n)
		{
			gpc_vertex &a = v[j];
			gpc_vertex &b = (j + 1 == n) ? v[0] : v[j+1];

			cA += a.x*b.y - b.x*a.y;
		}
		cA = abs(cA) * (p.hole[i] ? -1 : 1);
		A += cA;
	}
	A *= 0.5;
	return abs(A);
}


BOSTREAM2(const partitioned_skymap &m)
{
	out << m.dx;
	out << m.x0 << m.x1 << m.y0 << m.y1;	// sky footprint bounding rectangle
	out << m.skymap;			// skymap lookup table

	return out;
}

BISTREAM2(partitioned_skymap &m)
{
	in >> m.dx;
	in >> m.x0 >> m.x1 >> m.y0 >> m.y1;
	in >> m.skymap;

	return in;
}

BOSTREAM2(const partitioned_skymap::pixel_t &m)
{
	return out << m.poly << m.area;
}

BISTREAM2(partitioned_skymap::pixel_t &m)
{
	return in >> m.poly >> m.area;
}
