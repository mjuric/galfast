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

#ifndef __raytrace_h
#define __raytrace_h

#include <astro/math/vector.h>
#include <astro/types.h>
 
#include <algorithm>

//
// Simple 3 x 3 matrix
//
class M3
{
protected:
	peyton::math::V3 v[3]; // columns
public:
	M3(const double d = 0) { v[0] = d; v[1] = d; v[2] = d;}
	M3(double a, double b, double c) { v[0] = v[1] = v[2] = 0.; v[0][0] = a; v[1][1] = b; v[2][2] = c; }
	M3(const peyton::math::V3 &x, const peyton::math::V3 &y, const peyton::math::V3 &z) { v[0] = x; v[1] = y; v[2] = z; }
	double &operator ()(int r, int c) { return v[c][r]; }
	double operator ()(int r, int c) const { return v[c][r]; }
	peyton::math::V3 &operator()(int c) { return v[c]; }

	M3 &rotate(const peyton::math::V3 &axis, peyton::Radians angle);
};

struct Plane
{
	double a, b, c, w;
	Plane(double a_ = 0, double b_ = 0, double c_ = 0, double w_ = 0)
		: a(a_), b(b_), c(c_), w(w_)
	{}

	Plane(const peyton::math::V3 &v, double w_ = 0)
		: a(v.x), b(v.y), c(v.z), w(w_)
	{}

	double operator()(double x, double y, double z) { return a*x + b*y + c*z - w; }
};

class plane_transformer
{
public:
	Plane p;
	peyton::math::V3 x[3];	// original points used to define the plane p
	double delta;	// stuff +- delta from central plane will be accepted

	M3 M;	// earthcentric -> plane c.s.
	peyton::math::V3 t0;	// origin of plane c.s. in earthcentric coords
public:
	void setup(const peyton::math::V3 &x1, const peyton::math::V3 &x2, const peyton::math::V3 &x3, const peyton::math::V3 &origin, double delta_, bool earth_on_x_axis);
	peyton::math::V3 toPlane(const peyton::math::V3 &v);
};

void setupPlaneTransformer(
	plane_transformer &pt,
	const std::string &coordsys,
	double d1, std::pair<double, double> p1,
	double d2, std::pair<double, double> p2,
	double d3, std::pair<double, double> p3,
	double d0, std::pair<double, double> p0,
	double delta, bool earth_on_x_axis);

void ecenequ(const peyton::math::V3 &v, double &d, peyton::Radians ra, peyton::Radians dec);
peyton::math::V3 equecen(const double d, const peyton::Radians ra, const peyton::Radians dec);

#endif
