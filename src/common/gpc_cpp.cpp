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

#define MJURIC_IMPORT_GPC_IO

#include "gpc_cpp.h"
#include "projections.h"
#include <cmath>

#include <astro/exceptions.h>
#include <astro/useall.h>

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
	if(p.num_contours == 0)
	{
		x0 = x1 = y0 = y1 = 0;
		return;
	}

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
#if 1
double polygon_area(const gpc_polygon &p)
{
	double A = 0;
	FOR(0, p.num_contours)
	{
		double cA = contour_area(p.contour[i]);
		if(p.hole[i]) cA *= -1;
		A += cA;
	}
	return fabs(A);
}
#else
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
#endif

// calculate and return the area of a gpc_vertex_list
double contour_area(const gpc_vertex_list &c)
{
	gpc_vertex *v = c.vertex;
	double cA = 0;
	FORj(j, 0, c.num_vertices)
	{
		gpc_vertex &a = v[j];
		gpc_vertex &b = (j + 1 == c.num_vertices) ? v[0] : v[j+1];

		cA += a.x*b.y - b.x*a.y;
	}
	return 0.5*fabs(cA);
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
	return out << m.poly << m.coveredArea << m.pixelArea;
}

BISTREAM2(partitioned_skymap::pixel_t &m)
{
	return in >> m.poly >> m.coveredArea >> m.pixelArea;
}

// ASCII file serialization with projection information
void xgpc_write(const std::string &fn, const gpc_polygon &sky, const peyton::math::lambert &proj)
{
	FILE *ofp = fopen(fn.c_str(), "w");
	fprintf(ofp, "# lambert %.8f %.8f\n", deg(proj.l0), deg(proj.phi1));

	gpc_write_polygon(ofp, 1, const_cast<gpc_polygon*>(&sky));

	fclose(ofp);
}

gpc_polygon xgpc_read(const std::string &fn, peyton::math::lambert &proj, bool *hadProjection)
{
	FILE *fp = fopen(fn.c_str(), "r");
	if(fp == NULL) { THROW(EIOException, "Could not open footprint file [" + fn + "]"); }

	// check for projection information embedded in a first line comment
	char c = getc(fp);
	bool hp;
	if(c != '#')
	{
		ungetc(c, fp);
		hp = false;
	}
	else
	{
		double l0, b0;
		if(fscanf(fp, " lambert %lf %lf", &l0, &b0) == 2)
		{
			proj = peyton::math::lambert(rad(l0), rad(b0));
			hp = true;
		}
		else
		{
			hp = false;
		}
	}

	if(hadProjection)
		*hadProjection = hp;
	else if(!hp) {
		fclose(fp);
		THROW(EAny, "The .xgpc file " + fn + "had no projection information. File corrupted?");
	}

	gpc_polygon sky;
	gpc_read_polygon(fp, 1, &sky);

	fclose(fp);
	return sky;
}

//
// Simple crossing number algorithm, valid for non-self intersecting polygons.
// See http://softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm for details.
//
bool in_contour(const gpc_vertex &t, const gpc_vertex_list &vl)
{
	int cn = 0;
	const int n = vl.num_vertices;
	for(int i = 0; i != n; i++)
	{
		gpc_vertex &a = vl.vertex[i];
		gpc_vertex &b = (i + 1 == n) ? vl.vertex[0] : vl.vertex[i+1];

		if   (((a.y <= t.y) && (b.y > t.y))    		// an upward crossing
		   || ((a.y > t.y)  && (b.y <= t.y)))  		// a downward crossing
		{
			// compute the actual edge-ray intersect x-coordinate
			double vt = (float)(t.y - a.y) / (b.y - a.y);
			if (t.x < a.x + vt * (b.x - a.x)) 	// t.x < intersect
				++cn;				// a valid crossing of y=t.y right of t.x
		}
	}

	return (cn&1);    // 0 if even (out), and 1 if odd (in)	
}

// test if a given gpc_vertex in inside of a given gpc_polygon
bool in_polygon(const gpc_vertex &t, const gpc_polygon &p)
{
	bool in = false;
	FOR(0, p.num_contours)
	{
		bool inc = in_contour(t, p.contour[i]);
		in = (in != inc);
	}
	return in;
}

// test if a given point is inside of the partitioned_skymap
bool partitioned_skymap::in(peyton::Radians x, peyton::Radians y)
{
	int X = (int)((x - x0)/dx);
	int Y = (int)((y - y0)/dx);

	std::map<std::pair<int, int>, pixel_t>::iterator it = skymap.find(std::make_pair(X, Y));
	if(it == skymap.end()) { return false; }

	gpc_polygon &poly = it->second.poly;
	gpc_vertex vtmp = { x, y };
	if(!in_polygon(vtmp, poly))
	{
		return false;
	}
}
