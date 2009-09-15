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
#include "simulate_base.h"
#include "sph_polygon.h"

#include <astro/exceptions.h>
#include <astro/system/log.h>
#include <astro/system/config.h>
#include <astro/coordinates.h>
#include <astro/useall.h>

template<typename T> inline void memzero(T &x) { memset(&x, 0, sizeof(x)); }

gpc_polygon make_circle(double x0, double y0, double r, double dx)
{
	gpc_polygon p = {0, 0, NULL};
	int n = (int)(2*ctn::pi*r / dx);
	gpc_vertex v[n];
	gpc_vertex_list c = { n, v };
	FOR(0, n)
	{
		v[i].x = x0+r*cos(i*dx/r);
		v[i].y = y0+r*sin(i*dx/r);
	}
	gpc_add_contour(&p, &c, false);
	return p;
}

// Adds a spherical coordinates rectangle contour(s) to the polygon
int add_lonlat_rect(sph_polygon &poly, Radians l0, Radians b0, Radians l1, Radians b1, Radians dl, gpc_op op = GPC_UNION)
{
	// Reduce longitudes to [0, 2pi) domain
	l0 = peyton::math::modulo(l0, ctn::twopi);
	l1 = peyton::math::modulo(l1, ctn::twopi);
	if(fabs(l0 - l1) < 100*EPSILON_OF(l0))
	{
		// full band -- generate as two >half-sky contours (the 0.1*pi term is to ensure the contours
		// overlap and will be properly merged by gpc_* routines later on).
		int k = 0;
		k += add_lonlat_rect(poly, l0, b0, l0 + ctn::pi + 0.1*ctn::pi, b1, dl, op);
		k += add_lonlat_rect(poly, l0 + ctn::pi, b0, l0 + 0.1*ctn::pi, b1, dl, op);
		return k;
	}
	if(l0 > l1)
	{
		l0 -= ctn::twopi;
	}

	sph_contour &c = poly.add_contour(op);
	c.dx = dl;

	// Construct the rectangle (in spherical coordinates)
	int nl = (int)ceil((l1 - l0)/dl);
	dl = (l1 - l0) / nl;
	FOR(0, nl+1) { gpc_vertex vv = { l0 + dl*i, b0 }; c.c.push_back(vv); }	// Lower edge
	FOR(0, nl+1) { gpc_vertex vv = { l1 - dl*i, b1 }; c.c.push_back(vv); }	// Upper edge
	c.p.x = l0 + 0.5*(l1 - l0);	// A point _inside_ the polygon
	c.p.y = b0 + 0.5*(b1 - b0);

	return 1;
}

/**
	Return the polygon area in degrees (sugar).
*/
inline double polygon_area_deg(const gpc_polygon &sky)
{
	return polygon_area(sky)*sqr(180./ctn::pi);
}

/**
	great_circle -- Unit sphere great circle abstraction

	Representation of a great circle on a celestial sphere, bound by
	two endpoints. Undefined if the endpoints are the opposing poles
	of the sphere

	Based on formulae from http://williams.best.vwh.net/avform.htm
*/
struct great_circle
{
	double clon1, slon1, clat1, slat1, clon2, slon2, clat2, slat2, d;

	great_circle(Radians lon1, Radians lat1, Radians lon2, Radians lat2)
	{
		clon1 = cos(lon1);
		slon1 = sin(lon1);
		clat1 = cos(lat1);
		slat1 = sin(lat1);

		clon2 = cos(lon2);
		slon2 = sin(lon2);
		clat2 = cos(lat2);
		slat2 = sin(lat2);

		d = 2*asin(sqrt(sqr(sin((lat1-lat2)/2)) + clat1*clat2*sqr(sin((lon1-lon2)/2))));
	}

	double dist() const
	{
		/*
			Return the distance along the great circle between
			(lon1, lat1) and (lon2, lat2)
		*/
		return d;
	}

	/**
		Return an intermediate point along the great circle.

		f = 0 maps to (lon1, lat1)
		f = 1 maps to (lon2, lat2)
	*/
	const great_circle &intermediate(Radians &lon, Radians &lat, const double f) const
	{
		double A = sin((1.-f)*d)/sin(d);
		double B = sin(f*d)/sin(d);
		double x = A*clat1*clon1 +  B*clat2*clon2;
		double y = A*clat1*slon1 +  B*clat2*slon2;
		double z = A*slat1       +  B*slat2;
		lat = atan2(z,sqrt(x*x+y*y));
		lon = atan2(y,x);
		return *this;
	}
};

/**
	project_geodesic_aux -- Support for project_geodesic(). See project_geodesic_aux() documentation.
*/
void project_geodesic_aux(std::map<Radians, gpc_vertex> &out, const great_circle &gc, const lambert &proj, gpc_vertex a, gpc_vertex b, Radians fa, Radians fb, Radians dxmin)
{
	// this can happen if both ends of the interval have hit the pole
	if(a.x == b.x && a.y == b.y) { return; }

	// check if we've reached the desired precision in lambert plane
	Radians dlam = sqrt(sqr(a.x - b.x) + sqr(a.y - b.y));
	bool lam_done = dlam < dxmin;
	
	// check if we reached the desired precision in sky coordinates
	Radians dsky = (fb - fa)*gc.dist();
	bool sky_done = dsky < dxmin;

	// exit if precision reached, or we've reached the maximum resolution
	if(lam_done && sky_done) { return; }
	if(dsky < rad(1./3600.))
	{
		// we pass through here when subdividing an interval with one
		// end being the pole
		if(a.x != -10 && b.x != -10 || !sky_done)
		{
			MLOG(verb1) << "WARNING: Could not reach the desired precision when projecting sky polygon.";
		}
		return;
	}

	// compute the geodesic midpoint
	double fm = 0.5*(fa + fb);
	gpc_vertex p;
	gc.intermediate(p.x, p.y, fm);
	bool polem = !proj.project(p.x, p.y, p.x, p.y);
	if(polem)
	{
		 // signal we've hit the pole here
		p.x = p.y = -10;
	}
	out[fm] = p;

	// sample subintervals
	project_geodesic_aux(out, gc, proj, a, p, fa, fm, dxmin);
	project_geodesic_aux(out, gc, proj, p, b, fm, fb, dxmin);
}

/**
	project_geodesic -- project and sample a geodesic onto a Lambert equal area map

	This routine should be able to handle the corner cases where the south pole of the lambert
	projection lies on the geodesic.

	Adjacent sampled points on the geodesic are guaranteed to be not more than dxmin radians
	apart, both in projection plane and on the celestial sphere. If for any reason the routine
	fails to converge to this accuracy, a warning will be issued.
*/
void project_geodesic(std::vector<gpc_vertex> &vv, size_t &pole_at, gpc_vertex a, gpc_vertex b, const lambert proj, double dxmin)
{
	great_circle gc(a.x, a.y, b.x, b.y);

	// if this is a long geodesic, subdivide to 1deg scale before recursive refinment
	int nsamp = (int)ceil(gc.dist() / rad(1.));
	ASSERT(nsamp > 0);
	Radians df = 1. / nsamp;

	static const gpc_vertex polemarker = { -10, -10 };

	std::map<Radians, gpc_vertex> out;
	gpc_vertex pprev = a, pa, pb;
	FOR(1, nsamp+1)
	{
		out.clear();

		double f = i*df;
		gc.intermediate(pb.x, pb.y, f);
		great_circle gc2(pprev.x, pprev.y, pb.x, pb.y);

		// project and refine the geodesic
		bool polea = !proj.project(pa.x, pa.y, pprev.x, pprev.y);
		if(polea) { pa = polemarker; }
		pprev = pb;
		bool poleb = !proj.project(pb.x, pb.y, pb.x, pb.y);
		if(poleb) { pb = polemarker; }

		out[0.] = pa;
		if(!polea || !poleb)
		{
			project_geodesic_aux(out, gc2, proj, pa, pb, 0., 1., dxmin);
		}

		// add sampled vertices, and set a flag if we detected a pole
		vv.reserve(vv.size() + out.size());
		FOREACH(out)
		{
			if(i->second.x != -10)
			{
				vv.push_back(i->second);
			}
			else
			{
				if(pole_at == -1)
				{
					pole_at = vv.size();
				}
				else
				{
					if(pole_at != 0 && pole_at != vv.size())
					{
						THROW(EAny, "Multiple south pole crossings not allowed (self-intersecting polygon!).");
					}
				}
			}
		}
	}
}

/**
	sph_distance -- Compute the distance between two points on a unit sphere.
*/
Radians sph_distance(const gpc_vertex &a, const gpc_vertex &b)
{
	Radians lon1 = a.x, lat1 = a.y;
	Radians lon2 = b.x, lat2 = b.y;

	return 2*asin(sqrt(sqr(sin((lat1-lat2)/2)) + cos(lat1)*cos(lat2)*sqr(sin((lon1-lon2)/2))));
}

/**
	sph_triangle_area -- Compute the area of a spherical triangle.
*/
Radians sph_triangle_area(const gpc_vertex &A, const gpc_vertex &B, const gpc_vertex &C)
{
	Radians
		b = sph_distance(A, C),
		a = sph_distance(C, B),
		c = sph_distance(B, A),
		s = (a+b+c)/2;
	
	Radians E = 4*atan(sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2)));

	return E;
}

/**
	gpc_polygon -- Convert a contour in spherical coordinates to a (gpc_)polygon in projected coordinates

	This function should correctly handle the corner cases when the south pole is enclosed by the contour,
	or when a vertex of the contour is at the south pole.

	WARNING: The input must be well-behaved (no self-intersecting contours, or vertices that are used
	         more than once).
*/
gpc_polygon spherical_to_lambert(int nvert, const gpc_vertex *vsph, const gpc_vertex &inptsph, const lambert &proj, Radians dx)
{
	// Algorithm:
	// 0. Detect special cases:
	//	0a. All sky (polygon with one vertex)
	//	0b. Empty (polygon with no vertices)
	//	0c. Two-angle (error)
	// 1. Convert the vertices to output projection
	// 2. Figure out if inpt is inside the projected polygon
	//	2a. Yes: The polygon does not enclose the south pole
	//	2b. No: Poly topology inverted => encloses the south pole. Clip it.

	// 0.
	if(nvert == 1)	// All sky
	{
		// TODO: should I add a small margin here, so that the polygon fully encloses the r=2 circle?
		gpc_polygon allsky = make_circle(0, 0, 2., dx);
		return allsky;
	}
	if(nvert == 0)	// Empty footprint
	{
		gpc_polygon foot = { 0, NULL, NULL };
		return foot;
	}
	if(nvert == 2) // Two-angle
	{
		ASSERT(!(nvert == 2));
	}

	// 1.
	size_t pole_hit = (size_t)-1;
	std::vector<gpc_vertex> vv;
	FOR(0, nvert)
	{
		int j = (i+1) % nvert;
		const gpc_vertex &a = vsph[i], &b = vsph[j];

		project_geodesic(vv, pole_hit, a, b, proj, dx);
	}

	if(pole_hit != -1)
	{
		std::vector<gpc_vertex> v2;
		size_t p = pole_hit;
		std::copy(&vv[0], &vv[p], back_inserter(v2));

		// insert an arc with r=2, connecting a to b
		gpc_vertex a = vv[p ? p-1 : vv.size()-1];
		gpc_vertex b = vv[p % vv.size()];
		Radians phia = atan2(a.y, a.x);
		Radians phib = atan2(b.y, b.x);
		if(phia > phib)
		{
			std::swap(phia, phib);
		}
		int nsamp = (int)ceil(2.*(phib-phia) / dx);
		Radians dphi = (phib-phia) / nsamp;
		FOR(0, nsamp+1)
		{
			Radians phi = phia + dphi*i;
			gpc_vertex vv = { 2.*cos(phi), 2.*sin(phi) };
			v2.push_back(vv);
		}

		std::copy(vv.begin()+p, vv.end(), back_inserter(v2));

		vv.swap(v2);
	}

	gpc_polygon foot = { 0, NULL, NULL };
	gpc_vertex_list vl = { vv.size(), &vv[0] };
	gpc_add_contour(&foot, &vl, 0);

	// 2.
	gpc_vertex inpt;
	if(!proj.project(inpt.x, inpt.y, inptsph.x, inptsph.y))
	{
		inpt.x = inpt.y = -10;
	}
	Radians A = polygon_area(foot, true);
	if(A < 0)
	{
		// TODO: Once I'm convinced that the sign of area always changes if the south is
		// included in the footprint, inpt mechanism will be removed
		ASSERT(!in_polygon(inpt, foot));

		// 2b
		Radians dx = rad(0.05);
		gpc_polygon southpole = make_circle(0, 0, 2., dx);
		gpc_polygon_clip(GPC_DIFF, &southpole, &foot, &foot);

		gpc_free_polygon(&southpole);
	}

	size_t nv = 0;
	for(int i=0; i != foot.num_contours; i++)
	{
		nv += foot.contour[i].num_vertices;
	}

	DLOG(verb2) << "Area of projected contour: " << polygon_area_deg(foot) << "\n";
	DLOG(verb2) << "     [Contours, vertices]: " << foot.num_contours << ", " << nv << "\n";

	return foot;
}

gpc_polygon sph_contour::project(const peyton::math::lambert &proj) const
{
	return spherical_to_lambert(c.size(), &c[0], p, proj, dx);
}

gpc_polygon sph_polygon::project(const peyton::math::lambert &proj) const
{
	gpc_polygon sky;
	memzero(sky);

	//
	// Project the spherical polygons onto north/south hemispheres
	//
	FOREACH(contours)
	{
		gpc_polygon tmp = i->second.project(proj);
		gpc_polygon_clip(i->first, &sky, &tmp, &sky);
		gpc_free_polygon(&tmp);
	}

	return sky;
}

// Project and resample the list of contours onto two hemispheres defined by the projection proj (north), using
// project() to do the sampling, and possibly adding a margin along the equator to each hemisphere
std::pair<gpc_polygon, gpc_polygon> project_to_hemispheres(const std::list<sph_polygon> &foot, const peyton::math::lambert &proj, Radians dx)
{
	std::pair<gpc_polygon, gpc_polygon> sky, ns;
	memzero(sky);

	lambert sproj(modulo(proj.l0 + ctn::pi, ctn::twopi), -proj.phi1);

	FOREACH(foot)
	{
		ns.first  = i->project(proj);
		ns.second = i->project(sproj);

		gpc_polygon_clip(GPC_UNION, &sky.first,  &ns.first,  &sky.first);
		gpc_polygon_clip(GPC_UNION, &sky.second, &ns.second, &sky.second);

		gpc_free_polygon(&ns.first);
		gpc_free_polygon(&ns.second);
	}
	DLOG(verb2) << "Hemisphere areas, preclip (north, south) = (" << polygon_area_deg(sky.first) << ", " << polygon_area_deg(sky.second) << ")\n";

	// DEBUGGING
	// writeGPCPolygon("north.foot.txt", sky.first,  proj,  true);
	// writeGPCPolygon("south.foot.txt", sky.second, sproj, true);

	//
	// Clip north/south along the equator
	//
	Radians margin = 1./cos(0.5*dx) - 1.;	// Ensure the polygon makeHemisphereMaps will create
						// will be superscribed, not inscribed, in the unit circle
	gpc_polygon boundary = make_circle(0, 0, sqrt(2.) + margin, dx);

	gpc_polygon_clip(GPC_INT, &sky.first, &boundary, &sky.first);
	gpc_polygon_clip(GPC_INT, &sky.second, &boundary, &sky.second);

	gpc_free_polygon(&boundary);

	DLOG(verb2) << "Hemisphere map margin: " << margin << " (== " << deg(margin) << " deg)";
	DLOG(verb2) << "Hemisphere areas (north, south) = (" << polygon_area_deg(sky.first) << ", " << polygon_area_deg(sky.second) << ")\n";

	return sky;
}

int makeBeamMap(std::list<sph_polygon> &out, Radians l, Radians b, Radians r, Radians rhole, Radians dx);

int load_footprint_rect(std::list<sph_polygon> &out, const peyton::system::Config &cfg)
{
	// coordinate system (default: galactic)
	std::string coordsys;
	cfg.get(coordsys, "coordsys", std::string("gal"));

	// transformer from the chosen to galactic coordinate system
	peyton::coordinates::transform xxx2gal= peyton::coordinates::get_transform(coordsys, "gal");
	if(xxx2gal == NULL) { THROW(EAny, "Don't know how to convert from coordinate system '" + coordsys + "' to galactic coordinates."); }

	// load rectangle bounds and sampling scale dx
	std::string rect; double l0, l1, b0, b1;
	cfg.get(rect, "rect", std::string("0 1 0 1"));
	std::istringstream ss(rect);
	ss >> l0 >> l1 >> b0 >> b1;
	l0 = rad(l0); b0 = rad(b0);
	l1 = rad(l1); b1 = rad(b1);

	// sampling resolution
	Radians dx;
	cfg.get(dx, "dx", 0.);
	dx = rad(dx);

	// Reduce longitudes to [0, 360) domain
	l0 = peyton::math::modulo(l0, ctn::twopi);
	l1 = peyton::math::modulo(l1, ctn::twopi);
	// The complicated if is there to ensure that the common case
	// of l0=l1=0 (all sky) gets printed out as l0=0, l1=360
	if(l0 >= l1) { if(l0 == 0) { l1 = ctn::twopi; } else { l0 -= ctn::twopi; } }

	if(fabs((l1 - l0)/ctn::twopi - 1.) < 1e-5)
	{
		// ring -- delegate to beam
		Radians lon0, lat0;
		xxx2gal(0, ctn::halfpi, lon0, lat0);
		return makeBeamMap(out, lon0, lat0, ctn::halfpi - b0, ctn::halfpi - b1, dx);
	}
	MLOG(verb1) << "Footprint: lon=(" << deg(l0) << " to " << deg(l1) << ") lat=(" << deg(b0) << " to " << deg(b1) << "), coordsys=" << coordsys;

	// autosetting of dx
	if(dx == 0.)
	{
		dx = std::min(l1 - l0, b1 - b0) / 10.;	// 1/10th of the smaller dimension ...
		dx = std::min(dx, rad(0.1));			// ... but not bigger than 0.1deg
		MLOG(verb2) << "Autoset dx (rect): " << deg(dx) << "\n";
	}

	out.push_back(sph_polygon());
	sph_contour &c = out.back().add_contour();
	c.dx = dx;

	MLOG(verb2) << "Sampling resolution (dx, degrees)             = " << deg(dx);

	// Construct the rectangle (in spherical coordinates)
	int nl = (int)ceil((l1 - l0)/dx);
	double dl = (l1 - l0) / nl;

	bool southpole = fabs(-b0/ctn::halfpi - 1.) < 1e-5;
	bool northpole = fabs( b1/ctn::halfpi - 1.) < 1e-5;

	// Lower edge
	if(!southpole)
	{
		FOR(0, nl+1) { gpc_vertex vv = { l0 + dl*i, b0 }; c.c.push_back(vv); }
	}
	else
	{
		gpc_vertex vv = { l0, b0 }; c.c.push_back(vv);
//		gpc_vertex vv = { rad(180+15), rad(-20) }; c.c.push_back(vv);
	}

	// Upper edge
	if(!northpole)
	{
		FOR(0, nl+1) { gpc_vertex vv = { l1 - dl*i, b1 }; c.c.push_back(vv); }	// Upper edge
	}
	else
	{
		if(southpole)
		{
			// a case such as l0=30,b0=-90 l1=60,b1=90
			{ gpc_vertex vv = { l1, 0. }; c.c.push_back(vv); }
			{ gpc_vertex vv = { l1, b1 }; c.c.push_back(vv); }
			{ gpc_vertex vv = { l0, 0. }; c.c.push_back(vv); }
		}
		else
		{
			gpc_vertex vv = { l1, b1 }; c.c.push_back(vv);
		}
	}

	// A point _inside_ the polygon
	c.p.x = l0 + 0.5*(l1 - l0);
	c.p.y = b0 + 0.5*(b1 - b0);

	// transform to radians, Galactic coordinate system
	FOREACH(c.c)
	{
		xxx2gal(i->x, i->y, i->x, i->y);
	}
	xxx2gal(c.p.x, c.p.y, c.p.x, c.p.y);

	return 1;
}

int makeBeamMap(std::list<sph_polygon> &out, Radians l, Radians b, Radians r, Radians rhole, Radians dx)
{
	std::ostringstream ss;
	if(rhole) { ss << ", hole = " << deg(rhole) << "deg"; }
	MLOG(verb1) << "Footprint: Beam towards (l, b) = " << deg(l) << " " << deg(b) << ", radius = " << deg(r) << "deg" << ss.str();

	// autosetting of dx
	int npts;
	if(dx == 0.)
	{
		// Set polygon sampling as 1/1000th of 2*pi*min(r, hole)
		npts = 3600; // == 0.1 deg longitudinal resolution
		dx = ctn::twopi / npts;
		MLOG(verb2) << "Autosetting dx (beam): " << deg(dx) << "\n";
	}
	else
	{
		npts = (int)ceil(ctn::twopi / dx);
		dx = ctn::twopi / npts;
	}

	out.push_back(sph_polygon());
	sph_contour &c = out.back().add_contour();
	c.dx = dx;

	// test for full-sky
	if(fabs(r/ctn::pi - 1.) < 1e-5)
	{
		// Do we have a hole?
		if(rhole > 1e-5)
		{
			// this handles the corner case where the user sets r=180, rhole=nonzero
			// Convert it to r=180-rhole beam around the antipode
			l = modulo(l + ctn::pi, ctn::twopi);
			b = -b;
			r = ctn::pi - rhole;
			rhole = 0.;
		}
		else
		{
			// full sky
			c.reset(true);
			c.dx = dx;
			return 1;
		}
	}

	// construct the footprint in lambert coordinates centered on
	// the direction of the pencil beam
	lambert pproj(rad(90.), rad(90.));
	Radians dummy;
	pproj.project(r,     dummy, ctn::pi, ctn::pi/2. - r);
	if(rhole) pproj.project(rhole, dummy, ctn::pi, ctn::pi/2. - rhole);

	lambert proj(l, b);
	gpc_vertex v;
	FOR(0, npts)
	{
		proj.deproject(v.x, v.y, r*cos(dx*i), r*sin(dx*i));
		c.c.push_back(v);
	}

	// set the inside point
	c.p.x = l; c.p.y = b;

	if(rhole != 0.)
	{
		// add the hole contour
		sph_contour &c = out.back().add_contour(GPC_DIFF);
		c.dx = dx;
		FOR(0, npts)
		{
			proj.deproject(v.x, v.y, rhole*cos(dx*i), rhole*sin(dx*i));
			c.c.push_back(v);
		}

		// set the inside point
		c.p.x = l; c.p.y = b;
	}
	return 1;	// the number of added polygons
}

#define CFG_THROW(str) THROW(EAny, str)

int load_footprint_beam(std::list<sph_polygon> &foot, const peyton::system::Config &cfg)
{
	// coordinate system (default: galactic)
	std::string coordsys;
	cfg.get(coordsys, "coordsys", std::string("gal"));

	// transform the chosen direction to Galactic coordinate system
	peyton::coordinates::transform xxx2gal= peyton::coordinates::get_transform(coordsys, "gal");
	if(xxx2gal == NULL) { THROW(EAny, "Don't know how to convert from coordinate system '" + coordsys + "' to galactic coordinates."); }

	// beam direction and radius
	double l, b, r, rhole = 0;
	if(!cfg.count("footprint_beam")) { CFG_THROW("Configuration key footprint_beam must be set"); }
	std::istringstream ss2(cfg["footprint_beam"]);
	ss2 >> l >> b >> r >> rhole;

	DLOG(verb1) << "Radius, hole radius    = " << r << " " << rhole;

	// sampling resolution
	Radians dx;
	cfg.get(dx, "dx", 0.);
	dx = rad(dx);

	// transform to galactic coordinates (and radians)
	xxx2gal(rad(l), rad(b), l, b);

	return makeBeamMap(foot, l, b, rad(r), rad(rhole), dx);
}

