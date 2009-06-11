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

#include <astro/coordinates.h>
#include <astro/sdss/rungeometry.h>
#include <astro/system/options.h>

#include <fstream>
#include <sys/time.h>

#include "gpc_cpp.h"
#include "io.h"
#include "analysis.h"
#include "simulate.h"

#include <astro/useall.h>
using namespace std;

///////////////////////////////////


gpc_polygon make_polygon(const RunGeometry &geom, const lambert &proj, double dx)
{
	Mask mask(geom);
	Radians mu0 = geom.muStart;
	Radians mu1 = geom.muStart + mask.length();

	vector<double> ex, ey;

	gpc_polygon p = {0, 0, NULL};

	FORj(col, 0, 6)
	{
		Radians nu0 = mask.lo(col);
		Radians nu1 = mask.hi(col);
		Radians l, b, mu, nu;

		// bottom edge
		for(mu = mu0; mu < mu1; mu += dx)
		{
			coordinates::gcsgal(geom.node, geom.inc, mu, nu0, l, b);
			ex.push_back(l);
			ey.push_back(b);
		}
		// right edge
		for(nu = nu0; nu < nu1; nu += dx)
		{
			coordinates::gcsgal(geom.node, geom.inc, mu1, nu, l, b);
			ex.push_back(l);
			ey.push_back(b);
		}

		// top edge
		for(mu = mu1; mu > mu0; mu -= dx)
		{
			coordinates::gcsgal(geom.node, geom.inc, mu, nu1, l, b);
			ex.push_back(l);
			ey.push_back(b);
		}

		// left edge
		for(nu = nu1; nu > nu0; nu -= dx)
		{
			coordinates::gcsgal(geom.node, geom.inc, mu0, nu, l, b);
			ex.push_back(l);
			ey.push_back(b);
		}

		// convert to contour in lambert coordinates
		gpc_vertex v[ex.size()];
		FOR(0, ex.size())
		{
			proj.convert(ex[i], ey[i], v[i].x, v[i].y);
		}

		// TODO: check if the projection antipode is inside of the polygon
		
		// add contour to polygon
		gpc_vertex_list c = { ex.size(), v };
		gpc_add_contour(&p, &c, false);

		ex.clear();
		ey.clear();
	}
	return p;
}

class striped_polygon
{
public:
	int n;
	double dx;
	vector<gpc_polygon> stripes;
	double x0, x1;
public:
	striped_polygon(gpc_polygon &p, int n);
	~striped_polygon();
};

void sm_write(const std::string &fn, const gpc_polygon &p);

striped_polygon::striped_polygon(gpc_polygon &p, int n_)
	: n(n_)
{
	// split the polygon p into n stripes, for faster lookup
	// when doing poin-in-polygon searches.
	double y0, y1;
	poly_bounding_box(x0, x1, y0, y1, p);
	y0 += .01*y0; y1 += .01*y1;

	// split into n polygons
	dx = (x1-x0) / n;
	stripes.resize(n);
	FOR(0, n)
	{
		double xa = x0 +     i * dx;
		double xb = x0 + (i+1) * dx;
		gpc_vertex v[] = {{xa, y0}, {xb, y0}, {xb, y1}, {xa, y1}};
		gpc_vertex_list vl = {4, v};
		gpc_polygon rect = {1, NULL, &vl};

		gpc_polygon_clip(GPC_INT, &p, &rect, &stripes[i]);
//		cerr << i << " " << xa << " " << xb << " " << i * (x1-x0) / n << " " << (i+1) * (x1-x0) / n << "\n";

/*		sm_write("stripe.gpc.txt", stripes[i]);
		exit(0);*/
	}
}

striped_polygon::~striped_polygon()
{
	FOR(0, stripes.size()) { gpc_free_polygon(&stripes[i]); }
}

bool in_polygon(const gpc_vertex& t, striped_polygon &p)
{
	// find the right poly stripe
	int s = (int)floor((t.x - p.x0) / p.dx );
	ASSERT(s < p.n)
	{
		cerr << t.x << " " << t.y << " " << s << " " << p.n << "\n";
	}
	return in_polygon(t, p.stripes[s]);
}

double seconds()
{
	timeval tv;
	int ret = gettimeofday(&tv, NULL);
	assert(ret == 0);
	return double(tv.tv_sec) + double(tv.tv_usec)/1e6;
}

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

gpc_tristrip triangulatePoly(gpc_polygon sky)
{
	double dx = rad(5);
	double x0, x1, y0, y1, xa, xb, ya, yb;
//	std::vector<gpc_polygon> skymap;	// a map of rectangular sections of the sky, for fast is-point-in-survey-area lookup
	std::vector<gpc_vertex_list> tristrips;
 	poly_bounding_box(x0, x1, y0, y1, sky);
	for(double x = x0; x < x1; x += dx) // loop over all x values in the bounding rectangle
	{
		double xa = x, xb = x+dx;
		double xysum = 0.;
		for(double y = y0; y < y1; y += dx) // loop over all y values for a given x
		{
			double ya = y, yb = y+dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon poly;
			gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
			if(poly.num_contours == 0) continue; // if there are no observations in this direction

			gpc_tristrip tri = { 0, NULL };
			gpc_polygon_to_tristrip(&poly, &tri);
			tristrips.insert(tristrips.begin(), tri.strip, tri.strip + tri.num_strips);
//			skymap.push_back(poly); // store the polygon into a fast lookup map
		}
		cerr << "#";
	}

	gpc_tristrip tri;
	tri.num_strips = tristrips.size();
	tri.strip = (gpc_vertex_list*)malloc(sizeof(gpc_vertex_list) * tri.num_strips);
	FOR(0, tristrips.size())
	{
		tri.strip[i] = tristrips[i];
	}

	return tri;
}

void sm_write(const std::string &fn, const gpc_polygon &p)
{
	// dump polygons in SM compatible format
	ofstream out(fn.c_str());
	FOR(0, p.num_contours)
	{
		gpc_vertex *v = p.contour[i].vertex;
		FORj(j, 0, p.contour[i].num_vertices)
		{
			out << v[j].x << " " << v[j].y << " " << i << "\n";
		}
		out << "#\n";
	}
}

gpc_polygon make_polygon(const std::vector<double> &x, const std::vector<double> &y)
{
	ASSERT(x.size() == y.size());

	gpc_polygon p = {0, 0, NULL};
	int n = x.size();
	gpc_vertex v[n];
	gpc_vertex_list c = { n, v };
	FOR(0, n)
	{
		v[i].x = x[i];
		v[i].y = y[i];
	}
	gpc_add_contour(&p, &c, false);
	return p;
}

void writeGPCPolygon(const std::string &output, const gpc_polygon &sky, const lambert &proj, bool smOutput)
{
	int nvert = 0;
	FOR(0, sky.num_contours) { nvert += sky.contour[i].num_vertices; }
	MLOG(verb1) << "Statistics: " << polygon_area(sky)*sqr(deg(1)) << "deg2 area, "
			<< sky.num_contours << " contours, " << nvert << " vertices";

	// store the footprint polygon
	if(smOutput)
	{
		sm_write(output, sky);
		MLOG(verb1) << "Output stored in SM-readable format in " << output << " file.";
	}
	else
	{
		xgpc_write(output, sky, proj);
		MLOG(verb1) << "Output stored in .xgpc format in " << output << " file.";
	}

}

void poly_print_info(const gpc_polygon &poly, const char *name = "poly")
{
	DLOG(verb1) << name << ": ncontours=" << poly.num_contours;
	for(int i = 0; i != poly.num_contours; i++)
	{
		gpc_vertex_list &contour = poly.contour[i];
		DLOG(verb1) << name << ": contour " << i << ": nvertices=" << contour.num_vertices
			<< " area= " << contour_area(poly.contour[i])*sqr(deg(1.))
			<< (poly.hole[i] ? " (hole)" : "");
	}
}

// Splits an all-sky map into two nort/south hemisphere maps, possibly including an overlapping margin 
void makeHemisphereMaps(gpc_polygon &nsky, gpc_polygon &ssky, lambert &sproj, const lambert &nproj, gpc_polygon allsky, Radians dx = rad(.1), Radians margin = rad(0.))
{
	const Radians south_pole_epsilon2 = sqr(rad(0.03));

	gpc_polygon northBoundary = make_circle(0, 0, sqrt(2.) + margin, dx);
	gpc_polygon southBoundary = make_circle(0, 0, sqrt(2.) - margin, dx);

	// North sky (easy)
	gpc_polygon_clip(GPC_INT, &allsky, &northBoundary, &nsky);

	// South sky (hard)
	// Algorithm:
	//  1) Clip with equator
	//  2) Transform all vertices to south projection
	//  3) Reconstruct the polygon
	//  4) Try to detect & fix pileup of vertices near the south pole
	gpc_polygon_clip(GPC_DIFF, &allsky, &southBoundary, &ssky);

	sproj = lambert(modulo(nproj.l0 + ctn::pi, ctn::twopi), -nproj.phi1);
	DLOG(verb1) << "South projection pole: " << deg(sproj.l0) << " " << deg(sproj.phi1);
	gpc_vertex pole = {0., 0.};
	int at = 0;
	for(int i = 0; i != ssky.num_contours; i++)
	{
		gpc_vertex_list &contour = ssky.contour[i];
		bool south_pole = true;
//		int atv = 0; gpc_vertex prev = pole;
		for(int k = 0; k != contour.num_vertices; k++)
		{
			Radians l, b;
			double &x = contour.vertex[k].x;
			double &y = contour.vertex[k].y;

			nproj.inverse(x, y, l, b);
			sproj.convert(l, b, x, y);

#if 0
			// merge adjacent points
			Radians d2 = sqr(x - prev.x) + sqr(y - prev.y);
			const Radians epsilon2 = sqr(rad(0.005));
			if(d2 > epsilon2)
			{
				contour.vertex[atv] = contour.vertex[k];
				atv++;
				prev = contour.vertex[k];
			}
			else
			{
				std::cerr << "atv=" << atv << " k=" << k << " merging " << prev.x << "," << prev.y << " w. " << x << "," << y << " d=" << sqrt(d2) << "\n";
			}
#endif
			south_pole &= sqr(x) + sqr(y) < south_pole_epsilon2;
		}
//		contour.num_vertices = atv;

		// if the contour contains the projection pole,
		// invert it's hole significance
		if(in_contour(pole, contour))
		{
			ssky.hole[i] = !ssky.hole[i];
			south_pole &= ssky.hole[i];
		}
		else
		{
			south_pole = false;
		}

		if(south_pole)
		{
			DLOG(verb1) << "South pole detected, removing contour " << i;
			free(ssky.contour[i].vertex);
		}
		else
		{
			ssky.hole[at]    = ssky.hole[i];
			ssky.contour[at] = ssky.contour[i];
			at++;
		}
	}
	ssky.num_contours = at;

	double  area = polygon_area(allsky) * sqr(deg(1.));
	double narea = polygon_area(nsky) * sqr(deg(1.));
	double sarea = polygon_area(ssky) * sqr(deg(1.));
	double tarea = narea + sarea;
	double relerr   = fabs(tarea / area - 1.);
	double abserr   = fabs(tarea - area);

	#define failed  (relerr > 1e-4 && abserr > sqr(0.01))
	if(1 || margin == 0 && failed)
	{
		poly_print_info(nsky, "north");
		poly_print_info(ssky, "south");
		DLOG(verb1) << "(area, narea+sarea, relerr, abserr, narea, sarea) = "
			<< area << " " << tarea << " " << relerr << " " << abserr
			<< "     " << narea << " " << sarea;
		assert(margin != 0 || !failed);
	}
	#undef failed
}

gpc_polygon makeBeamMap(Radians l, Radians b, Radians r, Radians rhole, const lambert &proj);
gpc_polygon makeRectMap(const peyton::system::Config &cfg, lambert &proj)
{
	// coordinate system (default: galactic)
	std::string coordsys;
	cfg.get(coordsys, "coordsys", std::string("gal"));

	// transformer from the chosen to galactic coordinate system
	peyton::coordinates::transform xxx2gal= peyton::coordinates::get_transform(coordsys, "gal");
	if(xxx2gal == NULL) { THROW(EAny, "Don't know how to convert from coordinate system '" + coordsys + "' to galactic coordinates."); }

	// set up the lambert projection (in galactic coordinates), with the pole being the
	// pole of the chosen coordinate system
	Radians lp, bp;
	xxx2gal(rad(90), rad(90), lp, bp);
	proj = lambert(lp, bp);

	// load rectangle bounds and sampling scale dx
	std::string rect; double l0, l1, b0, b1;
	cfg.get(rect, "rect", std::string("0 1 0 1"));
	std::istringstream ss(rect);
	ss >> l0 >> l1 >> b0 >> b1;
	double dx;
	cfg.get(dx, "polydx", .1);

	// Rationalize l
	l0 = peyton::math::modulo(l0, 360);
	l1 = peyton::math::modulo(l1, 360);
	// The complicated if is there to ensure that the common case
	// of l0=l1=0 (all sky) gets printed out as l0=0, l1=360
	if(l0 >= l1) { if(l0 == 0) { l1 = 360; } else { l0 -= 360; } }
	if(fabs((l1 - l0) - 360) < 1e-5)
	{
		DLOG(verb1) << "Generating a ring -- delegating to beam footprint generation routine.";
		Radians r = rad(90. - b0), rhole = rad(90. - b1);
		return makeBeamMap(lp, bp, r, rhole, proj);
	}

	MLOG(verb1) << "Projection pole (l, b)               = " << deg(lp) << " " << deg(bp);
	MLOG(verb1) << "Observed rectangle (l0, l1, b0, b1)  = " << l0 << " " << l1 << " " << b0 << " " << b1;
	MLOG(verb1) << "Sampling resolution (dx)             = " << dx;

	// South pole causes a divergence in lambert routines
	if(b0 <= -89.99)
	{
		b0 = -89.99;
		MLOG(verb1) << "Footprint includes or approaches the south pole. Setting b0=" << b0 << " to avoid numerical issues.";
	}

	// Upper/lower edge
	std::vector<gpc_vertex> vup, vdown;
	gpc_vertex v;
	for(double l = l0; true; l += dx)
	{
		if(l > l1) { l = l1; }

		//std::cerr << l << " " << l1 << " " << dx << "\n";
		v.x = rad(l); v.y = rad(b1);
		xxx2gal(v.x, v.y, v.x, v.y); proj.convert(v.x, v.y, v.x, v.y);
		vup.push_back(v);

		v.x = rad(l); v.y = rad(b0);
		xxx2gal(v.x, v.y, v.x, v.y); proj.convert(v.x, v.y, v.x, v.y);
		vdown.push_back(v);

		if(l >= l1) { break; }
	}

	// Note: There's no reason to sample the sides, as they're going to be
	// straight lines in lambert projection because of the choice of coordinate
	// system pole as the projection pole.

	// Put together the rectangle
	std::vector<gpc_vertex> vv;
	vv.reserve(vup.size() + vdown.size());
	//std::cerr << "b1 = " << b1 << "\n";
	if(b1 <= 89.99999999)
	{
		vv.insert(vv.end(), vup.begin(), vup.end());
	}
	else
	{
		// insert a single vertex for the pole
		MLOG(verb1) << "Footprint includes the north pole. Reducing upper edge to a single vertex.";
		vv.push_back((gpc_vertex){0, 0});
	}
	vv.insert(vv.end(), vdown.rbegin(), vdown.rend());

	// construct the GPC polygon and write it to file
	gpc_polygon sky = {0, 0, NULL};
	gpc_vertex_list c = { vv.size(), &vv[0] };
	gpc_add_contour(&sky, &c, false);

	return sky;
}

// Create rectangular footprint based on cfg, store to output, possibly with sm if smOutput=true
void makeRectMap(const std::string &output, peyton::system::Config &cfg, bool smOutput)
{
	lambert proj;
	gpc_polygon sky = makeRectMap(cfg, proj);

//	gpc_polygon nsky, ssky; makeHemisphereMaps(nsky, ssky, proj, sky); exit();

	writeGPCPolygon(output, sky, proj, smOutput);
	gpc_free_polygon(&sky);
}

gpc_polygon makeBeamMap(Radians l, Radians b, Radians r, Radians rhole, const lambert &proj)
{

	ostringstream ss;
	if(rhole) { ss << ", hole = " << deg(rhole) << "deg"; }
	MLOG(verb1) << "Footprint: Beam towards (l, b) = " << deg(l) << " " << deg(b) << ", radius = " << deg(r) << "deg" << ss.str();

	double x, y;
	std::vector<double> lx, ly, hx, hy;
	lambert bproj(l, b);

	// Convert spherical (angular) radii to lambert projection radii
	lambert pproj(rad(0), rad(90));
	if(fabs(r/ctn::pi - 1) < 1e-5) // Full sky ?
	{
		//MLOG(verb1) << "DANGER: Detected attempt to compute full sky. There may be a small region around the south pole with no stars present. Aborting preventively.";
		//abort();
		r = 1.99999999;
	}
	else if(r > ctn::pi)
	{
		MLOG(verb1) << "ERROR: Beam radius (r = " << deg(r) << ") bigger than 180deg (full sky). Aborting.";
		abort();
	}
	else
	{
		pproj.convert(ctn::pi, ctn::pi/2. - r, x, r);
	}
	if(rhole)
	{
		pproj.convert(ctn::pi, ctn::pi/2. - rhole, x, rhole);
	}

	/* determine polygon sampling resolution */
	Radians rmax = rhole > 0 ? std::max(r, rhole) : r;
	Radians dx = ctn::twopi * rmax / 1000.;
	dx = std::min(dx, rad(.1)); // keep minimum .1 deg resolution along the circumference
	Radians dphi = dx / rmax;

	for(double phi = 0; phi < 2. * ctn::pi; phi += dphi)
	{
		Radians lp, bp;

		x = r * cos(phi);
		y = r * sin(phi);
		bproj.inverse(x, y, lp, bp);
		proj.convert(lp, bp, x, y);
		lx.push_back(x);
		ly.push_back(y);

		x = rhole * cos(phi);
		y = rhole * sin(phi);
		bproj.inverse(x, y, lp, bp);
		proj.convert(lp, bp, x, y);
		hx.push_back(x);
		hy.push_back(y);
	}

	gpc_polygon sky = make_polygon(lx, ly);

	if(rhole)
	{
		//std::cerr << "Clipping the hole.\n";
		gpc_polygon poly, hole = make_polygon(hx, hy);
		gpc_polygon_clip(GPC_DIFF, &sky, &hole, &poly);
		sky = poly;
	}

	// Sanity check -- the areas must be what we expect.
	double area_expected = ctn::pi * (sqr(r) - sqr(rhole));
	double area_poly = polygon_area(sky);
	double err = fabs(area_poly/area_expected - 1);
	double aerr = area_poly - area_expected;
	if(err > 1e-5)
	{
		MLOG(verb1) << "Polygonized footprint area differs from the expected by more than prescribed tolerance. Aborting preventively.";
		MLOG(verb1) << "Details: expected=" << area_expected*sqr(deg(1)) << " polygonized=" << area_poly*sqr(deg(1)) << " relerr=" << err << " abserr=" << aerr*sqr(deg(1));
		abort();
	}

//	gpc_polygon nsky, ssky; makeHemisphereMaps(nsky, ssky, proj, sky); exit();

	return sky;
}

#define CFG_THROW(str) THROW(EAny, str)

gpc_polygon makeBeamMap(const peyton::system::Config &cfg, lambert &proj)
{
	// coordinate system (default: galactic)
	std::string coordsys;
	cfg.get(coordsys, "coordsys", std::string("gal"));

	// transform the chosen direction to Galactic coordinate system
	peyton::coordinates::transform xxx2gal= peyton::coordinates::get_transform(coordsys, "gal");
	if(xxx2gal == NULL) { THROW(EAny, "Don't know how to convert from coordinate system '" + coordsys + "' to galactic coordinates."); }

	// projection
	std::string pole; double l0, b0;
	cfg.get(pole,	"projection",	std::string("90 90"));
	std::istringstream ss(pole);
	ss >> l0 >> b0;

	// beam direction and radius
	double l, b, r, rhole = 0;
	if(!cfg.count("footprint_beam")) { CFG_THROW("Configuration key footprint_beam must be set"); }
	std::istringstream ss2(cfg["footprint_beam"]);
	ss2 >> l >> b >> r >> rhole;

	DLOG(verb1) << "Projection pole (l, b) = " << l0 << " " << b0;
	DLOG(verb1) << "Radius, hole radius    = " << r << " " << rhole;

	// transform to galactic coordinates (and radians)
	xxx2gal(rad(l), rad(b), l, b);
	xxx2gal(rad(l0), rad(b0), l0, b0);

	proj = lambert(l0, b0);
	return makeBeamMap(l, b, rad(r), rad(rhole), proj);

}

inline double polygon_area_deg(const gpc_polygon &sky)
{
	return polygon_area(sky)*sqr(180./ctn::pi);
}

gpc_polygon clip_zone_of_avoidance(gpc_polygon &sky, double bmin, const peyton::math::lambert &proj)
{
	static const int npts = 360;
	Radians dx = ctn::twopi / npts;

	if(abs(proj.phi1) == bmin) { bmin *= 0.999; }	// slightly reduce ZOA to ensure the pole is out of it

	// construct ZoA contours
	std::vector<double> x(npts), y(npts);
	gpc_polygon mask, contour1, contour2;
	FOR(0, npts) { proj.convert(i*dx,  bmin, x[i], y[i]); }
	contour1 = make_polygon(x, y);
	FOR(0, npts) { proj.convert(i*dx, -bmin, x[i], y[i]); }
	contour2 = make_polygon(x, y);

	DLOG(verb2) << "A(contour1) = " << polygon_area_deg(contour1) << "\n";
	DLOG(verb2) << "A(contour2) = " << polygon_area_deg(contour2) << "\n";
	DLOG(verb2) << "    A(band) = " << polygon_area_deg(contour2)-polygon_area_deg(contour1) << "\n";

	// determine the topology of ZoA curve and accordingly construct
	// the clipping mask
	if(abs(proj.phi1) < bmin)
	{
		// pole is within the zone
		gpc_polygon circle = make_circle(0., 0., 2.001, dx);		// whole sky

		// construct the mask
		gpc_polygon tmp;
		gpc_polygon_clip(GPC_DIFF, &circle, &contour1, &tmp);
		gpc_polygon_clip(GPC_DIFF, &tmp, &contour2, &mask);

		gpc_free_polygon(&tmp);
		gpc_free_polygon(&circle);
	}
	else
	{
		// pole is out of the zone. Topology is a simple band.
		if(proj.phi1 < -bmin) { std::swap(contour1, contour2); } // ensure contour1 is the inner contour

		// construct the mask
		gpc_polygon_clip(GPC_DIFF, &contour2, &contour1, &mask);
	}
	DLOG(verb2) << "A(mask) = " << polygon_area_deg(mask) << "\n";

	gpc_polygon sky2;
	gpc_polygon_clip(GPC_DIFF, &sky, &mask, &sky2);
	DLOG(verb2) << "A(sky) = " << polygon_area_deg(sky) << "\n";
	DLOG(verb2) << "A(sky2) = " << polygon_area_deg(sky2) << "\n";

	gpc_free_polygon(&mask);
	gpc_free_polygon(&contour1);
	gpc_free_polygon(&contour2);

	return sky2;
}

void makeSkyMap(std::set<int> &runs, const std::string &output, const lambert &proj, Radians b0 = rad(0.))
{
	Radians dx = rad(.25); /* polygon sampling resolution in radians */
	cerr << "footprint = " << output << ", dx = " << dx << " radians\n";

	RunGeometryDB db;
	double x, y;
	gpc_polygon sky = {0, 0, NULL};
	proj.convert(rad(0.), b0, x, y);
	double r = sqrt(x*x+y*y);
	cerr << "Excluding r > " << r << " from the origin of lambert projection.\n";
	gpc_polygon circle = make_circle(0., 0., r, dx);

	cerr << "Processing " << runs.size() << " runs.\n";

	int k = 0;
	FOREACH(runs)
	{
		const RunGeometry &geom = db.getGeometry(*i);
 		cerr << ".";

		gpc_polygon rpoly = make_polygon(geom, proj, dx);
		gpc_polygon_clip(GPC_INT, &rpoly, &circle, &rpoly);
		gpc_polygon_clip(GPC_UNION, &sky, &rpoly, &sky);

	}
	cerr << "\n";

	int nvert = 0;
	FOR(0, sky.num_contours) { nvert += sky.contour[i].num_vertices; }
	cerr << "total [" << polygon_area(sky)*sqr(deg(1)) << "deg2 area, "
	     << sky.num_contours << " contours, " << nvert << " vertices]\n";

	// store the footprint polygon
	sm_write(output + ".foot.txt", sky);
	xgpc_write(output, sky, proj);
//	FILE *ofp = fopen(output.c_str(), "w");
//	gpc_write_polygon(ofp, 1, &sky);
//	fclose(ofp);

	// free memory
	gpc_free_polygon(&sky);
	gpc_free_polygon(&circle);
}

namespace footprints
{
	std::string output;
	bool smOutput = false;

	int pencilBeam(peyton::system::Config &cfg)
	{
		lambert proj;
		gpc_polygon sky = makeBeamMap(cfg, proj);
		writeGPCPolygon(output, sky, proj, smOutput);
		gpc_free_polygon(&sky);

		return 0;
#if 0
		// projection
		std::string pole; double l0, b0;
		cfg.get(pole,	"projection",	std::string("90 90"));
		std::istringstream ss(pole);
		ss >> l0 >> b0;
		lambert proj(rad(l0), rad(b0));

		// beam direction and radius
		double l, b, r, rhole = 0;
		if(!cfg.count("footprint_beam")) { CFG_THROW("Configuration key footprint_beam must be set"); }
		std::istringstream ss2(cfg["footprint_beam"]);
		ss2 >> l >> b >> r >> rhole;

		MLOG(verb1) << "Projection pole (l, b) = " << l0 << " " << b0;
		DLOG(verb1) << "Radius, hole radius    = " << r << " " << rhole;
		makeBeamMap(output, rad(l), rad(b), rad(r), rad(rhole), proj, smOutput);

		return 0;
#endif
	}

	int rectangle(peyton::system::Config &cfg)
	{
		makeRectMap(output, cfg, smOutput);
		return 0;
#if 0
		// coordinate system (default: galactic)
		std::string coordsys;
		cfg.get(coordsys, "coordsys", std::string("gal"));

		// transformer from the chosen to galactic coordinate system
		peyton::coordinates::transform xxx2gal= peyton::coordinates::get_transform(coordsys, "gal");
		if(xxx2gal == NULL) { THROW(EAny, "Don't know how to convert from coordinate system '" + coordsys + "' to galactic coordinates."); }

		// set up the lambert projection (in galactic coordinates), with the pole being the
		// pole of the chosen coordinate system
		Radians lp, bp;
		xxx2gal(rad(90), rad(90), lp, bp);
		lambert proj(lp, bp);

		// load rectangle bounds and sampling scale dx
		std::string rect; double l0, l1, b0, b1;
		cfg.get(rect, "rect", std::string("0 1 0 1"));
		std::istringstream ss(rect);
		ss >> l0 >> l1 >> b0 >> b1;
		double dx;
		cfg.get(dx, "dx", .1);

		// Rationalize l
		l0 = peyton::math::modulo(l0, 360);
		l1 = peyton::math::modulo(l1, 360);
		// The complicated if is there to ensure that the common case
		// of l0=l1=0 (all sky) gets printed out as l0=0, l1=360
		if(l0 >= l1) { if(l0 == 0) { l1 = 360; } else { l0 -= 360; } }

		MLOG(verb1) << "Projection pole (l, b)               = " << deg(lp) << " " << deg(bp);
		MLOG(verb1) << "Observed rectangle (l0, l1, b0, b1)  = " << l0 << " " << l1 << " " << b0 << " " << b1;
		MLOG(verb1) << "Sampling resolution (dx)             = " << dx;

		// South pole causes a divergence in lambert routines
		if(b0 <= -89.99999999)
		{
			b0 = -89.9999;
			MLOG(verb1) << "Footprint includes the south pole. Setting b0=" << b0 << " to avoid numerical issues.";
		}

		// Upper/lower edge
		std::vector<gpc_vertex> vup, vdown;
		gpc_vertex v;
		for(double l = l0; true; l += dx)
		{
			if(l > l1) { l = l1; }

			//std::cerr << l << " " << l1 << " " << dx << "\n";
			v.x = rad(l); v.y = rad(b1);
			xxx2gal(v.x, v.y, v.x, v.y); proj.convert(v.x, v.y, v.x, v.y);
			vup.push_back(v);

			v.x = rad(l); v.y = rad(b0);
			xxx2gal(v.x, v.y, v.x, v.y); proj.convert(v.x, v.y, v.x, v.y);
			vdown.push_back(v);

			if(l >= l1) { break; }
		}

		// Note: There's no reason to sample the sides, as they're going to be
		// straight lines in lambert projection because of the choice of coordinate
		// system pole as the projection pole.

		// Put together the rectangle
		std::vector<gpc_vertex> vv;
		vv.reserve(vup.size() + vdown.size());
		//std::cerr << "b1 = " << b1 << "\n";
		if(b1 <= 89.99999999)
		{
			vv.insert(vv.end(), vup.begin(), vup.end());
		}
		else
		{
			// insert a single vertex for the pole
			MLOG(verb1) << "Footprint includes the north pole. Reducing upper edge to a single vertex.";
			vv.push_back((gpc_vertex){0, 0});
		}
		vv.insert(vv.end(), vdown.rbegin(), vdown.rend());

		// construct the GPC polygon and write it to file
		gpc_polygon sky = {0, 0, NULL};
		gpc_vertex_list c = { vv.size(), &vv[0] };
		gpc_add_contour(&sky, &c, false);
		writeGPCPolygon(output, sky, proj, smOutput);
		gpc_free_polygon(&sky);

		return 0;
#endif
	}
}

//void make_skymap(partitioned_skymap &m, Radians dx, const std::string &skypolyfn);
void pdfinfo(std::ostream &out, const std::string &pdffile);

inline uint32_t rng_mwc(uint32_t *xc)
{
	#define c (xc[0])
	#define x (xc[1])
	#define a (xc[2])

	uint64_t xnew = (uint64_t)a*x + c;
//	printf("%016llx\n", xnew);
	c = xnew >> 32;
	x = (xnew << 32) >> 32;
	return x;

	#undef c
	#undef x
	#undef a
}

void test_mwc_rng()
{
	return;

	uint32_t xc[3] = { 42, 332, 8834973 };
	FOR(0, 50)
	{
		printf("%08x\n", rng_mwc(xc));
	}
	exit(-1);
}

#if 0
void test_pm_conversions2()
{
//	return;

	float  x = 100, y = 200, z = 300;
	float vx = 30, vy = 20, vz = 10;
/*	float  x = 0, y = 100, z = 0;
	float vx = 0, vy = 30, vz = 0;*/
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = x, y, z, |r|\n", x, y, z, sqrt(x*x + y*y + z*z));
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = vx, vy, vz, |v|\n", vx, vy, vz, sqrt(vx*vx + vy*vy + vz*vz));

	// convert to l,b,pm_lb
	Radians l, b; float r, vl, vb, vr;
	    xyz2lbr(l, b, r, x, y, z);
	vel_xyz2lbr(vl, vb, vr, vx, vy, vz, l, b);
	printf("(% 12.7f, % 12.7f, % 12.7f) = l, b, r\n", l, b, r);
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = vl, vb, vr, mu\n", vl, vb, vr, sqrt(vl*vl + vb*vb));

	if(0) {
		// test if cartesian<->celestial conversions work
		float ex, ey, ez, evx, evy, evz;
		    lbr2xyz(ex, ey, ez, l, b, r);
		vel_lbr2xyz(evx, evy, evz, vl, vb, vr, l, b);
		printf("(% 12.7f, % 12.7f, % 12.7f) = ex, ey, ez\n", ex, ey, ez);
		printf("(% 12.7f, % 12.7f, % 12.7f) = evx, evy, evz\n", evx, evy, evz);
		exit(0);
	}

	// convert to ra,dec,pm_ra,pm_dec
	Radians ra, dec; float vra, vdec;
	peyton::coordinates::galequ(l, b, ra, dec);
	pm_galequ(vra, vdec, l, b, vl, vb);
	printf("(% 12.7f, % 12.7f, % 12.7f) = ra, dec, r\n", ra, dec, r);
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = vra, vdec, vr, mu\n", vra, vdec, vr, sqrt(vra*vra + vdec*vdec));

	// convert to equatorial cartesian
	float ex, ey, ez, evx, evy, evz;
	    lbr2xyz(ex, ey, ez, ra, dec, r);
	vel_lbr2xyz(evx, evy, evz, vra, vdec, vr, ra, dec);
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = ex, ey, ez, |er|\n", ex, ey, ez, sqrt(ex*ex + ey*ey + ez*ez));
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = evx, evy, evz, |ev|\n", evx, evy, evz, sqrt(evx*evx + evy*evy + evz*evz));

	// advance the trajectory in equatorial cartesian space
	float dt = 1.;
	ex += dt*evx;
	ey += dt*evy;
	ez += dt*evz;
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = moved ex, ey, ez, |er|\n", ex, ey, ez, sqrt(ex*ex + ey*ey + ez*ez));

	// advance the trajectory in galactic cartesian space
	x += dt*vx;
	y += dt*vy;
	z += dt*vz;
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = moved x, y, z, |r|\n", x, y, z, sqrt(x*x + y*y + z*z));

	// convert back to equatorial coordinates
	Radians ra2, dec2, l2, b2; float r2, x2, y2, z2;
	xyz2lbr(ra2, dec2, r2, ex, ey, ez);
	printf("(% 12.7f, % 12.7f, % 12.7f) = ra2, dec2, r2\n", ra2, dec2, r2);
	// convert back to galactic coordinates
	peyton::coordinates::equgal(ra2, dec2, l2, b2);
	printf("(% 12.7f, % 12.7f) = l2, b2, r2\n", l2, b2, r2);
	// convert to galactic cartesian
	lbr2xyz(x2, y2, z2, l2, b2, r2);
	printf("\n");
	printf("(% 12.7f, % 12.7f, % 12.7f) = x, y, z  from equatorial\n", x2, y2, z2);

	// pray that they're the same
	printf("(% 12.7f, % 12.7f, % 12.7f) = x, y, z  from galactic\n", x, y, z);

#if 0
	double x, y;
	peyton::coordinates::equgal(0., rad(90), x, y);
	printf("pole       = %.16lf %.16lf\n", deg(x), deg(y));

/*	double l = rad(139.09170669), b = rad(-89.95);
	float vl = 108.6, vb = -75.1;*/
	double l = rad(139.09170669), b = rad(-61.75582115);
	float vl = 61.095611, vb = -46.828220;

	double ra, dec;
	peyton::coordinates::galequ(l, b, ra, dec);
	printf("     radec = %.16lf %.16lf\n", deg(ra), deg(dec));

	printf("   vlb, mu = %.16f %.16f %.16f\n", vl, vb, sqrt(sqr(vl*cos(b)) + sqr(vb)));
	float vra, vdec;
	pm_galequ(vra, vdec, l, b, vl, vb);
	printf("vradec, mu = %.16f %.16f %.16f\n", vra, vdec, sqrt(sqr(vra*cos(dec)) + sqr(vdec)));
	pm_equgal(vl, vb, ra, dec, vra, vdec);
	printf("       vlb = %.16f %.16f\n", vl, vb);	
#endif
	exit(0);
}
#endif

#include "simulate_base.h"
void pm_galequ(float &vra, float &vdec, float l, float b, float vl, float vb);
void pm_equgal(float &vl, float &vb, float ra, float dec, float vra, float vdec);
void vel_xyz2lbr(float &vl, float &vb, float &vr, const float vx, const float vy, const float vz, const float l, const float b);
void vel_lbr2xyz(float &vx, float &vy, float &vz, const float vl, const float vb, const float vr, const float l, const float b);
void xyz2lbr(Radians &l, Radians &b, float &r, float x, float y, float z);
void lbr2xyz(float &x, float &y, float &z, Radians l, Radians b, float r);

static const double AU = 149597870691.; /* 1AU in meters; JPL DE405 value; http://en.wikipedia.org/wiki/Astronomical_unit */
static const Radians as = ctn::twopi/(360*3600);
static const double pc = AU / as;
static const double yr = 31557600; /* 1 Julian year in seconds; see http://en.wikipedia.org/wiki/Julian_year_(astronomy) */
static const float kms_per_masyrpc = as/yr * pc / 1000.; /* ~4.74 km/s @ 1kpc is 1mas/yr */
static const double kms_to_pcyr = 1000.*yr/pc;

/* convert cartesian velocities (in km/s) to proper motions (in mas/yr) + radial velocity (km/s) */
void vel_xyz2pmr(float &pml, float &pmb, float &vr, const float vx, const float vy, const float vz, const float l, const float b, const float r)
{
	vel_xyz2lbr(pml, pmb, vr, vx, vy, vz, l, b);

	// proper motion in mas/yr
	pml /= kms_per_masyrpc * r*1e-3;
	pmb /= kms_per_masyrpc * r*1e-3;
}

/* convert proper motions (in mas/yr) + radial velocity (km/s) to cartesian velocities (in km/s) */
void vel_pmr2xyz(float &vx, float &vy, float &vz, float pml, float pmb, const float vr, const float l, const float b, const float r)
{
	pml *= kms_per_masyrpc * r*1e-3;
	pmb *= kms_per_masyrpc * r*1e-3;

	vel_lbr2xyz(vx, vy, vz, pml, pmb, vr, l, b);
}

/* convert velocities in km/s to pc/yr */
void kms2pcyr(float &vx2, float &vy2, float &vz2, float vx, float vy, float vz)
{
	vx2 = vx * kms_to_pcyr;
	vy2 = vy * kms_to_pcyr;
	vz2 = vz * kms_to_pcyr;
}

// advance the trajectory a given number of years
void advance_trajectory(Radians &l, Radians &b, float &r, float &vl, float &vb, float &vr, float dt)
{
	// convert to cartesian
	float x, y, z, vx, vy, vz, vxpc, vypc, vzpc;
	    lbr2xyz( x,  y,  z,  l,  b,  r);		// convert position
	vel_pmr2xyz(vx, vy, vz, vl, vb, vr, l, b, r);	// convert velocities
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = x, y, z, |r|\n", x, y, z, sqrt(x*x + y*y + z*z));
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = vx, vy, vz, |v| (km/s)\n", vx, vy, vz, sqrt(vx*vx + vy*vy + vz*vz));

	// advance the trajectory
	kms2pcyr(vxpc, vypc, vzpc, vx, vy, vz);	// convert velocities to parsecs/yr
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = vx, vy, vz, |v| (pc/Myr)\n", 1e6*vxpc, 1e6*vypc, 1e6*vzpc, 1e6*sqrt(vxpc*vxpc + vypc*vypc + vzpc*vzpc));
	x += dt*vxpc;
	y += dt*vypc;
	z += dt*vzpc;
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f) = new x, y, z, |r|\n", x, y, z, sqrt(x*x + y*y + z*z));

	// convert back to celestial
	    xyz2lbr(l, b, r, x, y, z);
	vel_xyz2pmr(vl, vb, vr, vx, vy, vz, l, b, r);
	printf("(% 12.7f, % 12.7f, % 12.7f) (% 12.7f, % 12.7f, % 12.7f) = new (astrom, r) (pm, vr) (|pm|)\n",   deg(l), deg(b), r, vl, vb, vr, sqrt(vl*vl + vb*vb));
	printf("\n");
}

void test_pm_conversions()
{
	return;

	Radians l, b, ra, dec, l2, b2;
	float r, vl, vb, vr, vra, vdec, vl2, vb2, re, vre;
	float dt;

	/*************** Setup **************/
/*	l = rad(11); b = rad(77); r = 100.;     vl = 20.; vb = 20.; vr = 30.;	dt = 1e6;*/
/*	l = rad(139.09170669); b = rad(-61.75582115); r = 100.;     vl = 61.095611; vb = -46.828220; vr = 25.557892;	dt = 1e6;*/
/*	l = rad(0); b = rad(0); r = 100.;     vl = 0; vb = 10.; vr = 0.;	dt = 1e6;*/
	l = rad(50.92620493); b = rad(89.84075570); r = 654.;     vl = -5.8; vb = 4.2; vr = -20.6;	dt = 1e6;
	l = rad(55.26824585); b = rad(89.43553778); r = 654.;     vl = -4.05; vb = 5.00; vr = -51.54;	dt = 1e6;
	l = rad(139.09170669);  b=rad(-61.75582115); r = 222.;    vl = 167.18; vb = 96.21; vr = -0.28;  dt=1e6;

	peyton::coordinates::galequ(l, b, ra, dec);	// convert coordinates to equatorial
	pm_galequ(vra, vdec, l, b, vl, vb);		// convert proper motion to equatorial
	re = r; vre = vr;

	printf("Input values:\n");
	printf("Gal (pos, pm&vr, |pm|): (% 12.7f, % 12.7f, % 12.7f) (% 12.7f, % 12.7f, % 12.7f) (% 12.7f)\n",  deg(l),   deg(b),  r,  vl,   vb,  vr, sqrt( vl*vl  +   vb*vb));
	printf("Equ (pos, pm&vr, |pm|): (% 12.7f, % 12.7f, % 12.7f) (% 12.7f, % 12.7f, % 12.7f) (% 12.7f)\n", deg(ra), deg(dec), re, vra, vdec, vre, sqrt(vra*vra + vdec*vdec));
	printf("\n");

	/*********** Work in galactic coordinate system *************/
	printf("Moving the star for %.0fMyr, using Gal. position and pm:\n", dt/1e6);
	advance_trajectory(l, b, r, vl, vb, vr, dt);

	/*********** Work in equatorial coordinate system *************/
	printf("Moving the star for %.0fMyr, using Equ. position and pm:\n", dt/1e6);
	advance_trajectory(ra, dec, re, vra, vdec, vre, dt);

	/*********** Compare the results *************/
	printf("Final coordinates and pm:\n");
	peyton::coordinates::equgal(ra, dec, l2, b2);
	pm_equgal(vl2, vb2, ra, dec, vra, vdec);
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f, % 12.7f, % 12.7f) =  l, b, r, vl, vb, vr  (advanced in galactic)\n",   deg(l), deg(b), r, vl, vb, vr);
	printf("(% 12.7f, % 12.7f, % 12.7f, % 12.7f, % 12.7f, % 12.7f) =  l, b, r, vl, vb, vr  (advanced in equatorial)\n", deg(l2), deg(b2), re, vl2, vb2, vre);

	exit(0);
}

void test_kin();
void test_tags();
void test_otable();
void test_pm_conversions();
int main(int argc, char **argv)
{
try
{
#define ASCII_ESC 27
//	printf( "%c[2J", ASCII_ESC );
//	printf( "%c[1m Bla bla %c[m", ASCII_ESC, ASCII_ESC );
//	abort();

	test_pm_conversions();
//	test_kin();
	test_otable();
	test_mwc_rng();
//	test_tags(); return 0;

	std::string argv0 = argv[0];
	VERSION_DATETIME(version, "$Id: simulate.cpp,v 1.19 2007/04/15 12:09:52 mjuric Exp $");
	std::string progdesc = "simulate.x, a mock star catalog simulator.";

	std::string cmd, input, output;
	std::map<std::string, boost::shared_ptr<Options> > sopts;
	Options opts(argv[0], progdesc, version, Authorship::majuric);
	opts.argument("cmd").bind(cmd).desc(
		"What to make. Can be one of:\n"
		"       foot - \tcalculate footprint given a config file\n"
		"  footprint - \tcalculate footprint of a set of runs on the sky (deprecated)\n"
		"       beam - \tcalculate footprint of a single conical beam (deprecated)\n"
		"        pdf - \tcalculate cumulative probability density functions (CPDF) for a given model and footprint\n"
		"    pdfinfo - \tget information about the contents of a .pdf.bin file\n"
		"    catalog - \tcreate a mock catalog given a set of CPDFs\n"
		"postprocess - \tpostprocess the mock catalog (e.g., derive photometry, add instrumental errors, etc.)\n"
					   );
	opts.stop_after_final_arg = true;
	opts.prolog = "For detailed help on a particular subcommand, do `simulate.x <cmd> -h'";
	opts.add_standard_options();

	Radians dx = 4.;

	sopts["footprint"].reset(new Options(argv0 + " footprint", progdesc + " Footprint polygon generation subcommand.", version, Authorship::majuric));
	sopts["footprint"]->argument("conf").bind(input).desc("Footprint configuration file");
	sopts["footprint"]->prolog = 
		"Input run list filename is read from $footprint_runs configuration variable. "
		"Output filename is read from $footprint confvar.";
	sopts["footprint"]->add_standard_options();

	sopts["foot"].reset(new Options(argv0 + " foot", progdesc + " Footprint polygon generation subcommand.", version, Authorship::majuric));
	sopts["foot"]->argument("conf").bind(input).desc("Footprint configuration file (input)");
	sopts["foot"]->argument("output").bind(footprints::output).desc(".xgpc.txt footprint output file (output)");
	sopts["foot"]->option("s").bind(footprints::smOutput).addname("sm_out").value("true").desc("Generate SM-readable output instead of .xgpc format");
	sopts["foot"]->add_standard_options();

	sopts["beam"].reset(new Options(argv0 + " beam", progdesc + " Conical beam footprint polygon generation subcommand.", version, Authorship::majuric));
	sopts["beam"]->argument("conf").bind(input).desc("Footprint configuration file (input)");
	sopts["beam"]->argument("output").bind(footprints::output).desc(".xgpc.txt footprint output file (output)");
	sopts["beam"]->option("s").bind(footprints::smOutput).addname("sm_out").value("true").desc("Generate SM-readable output instead of .xgpc format");
	sopts["beam"]->prolog = "Deprecated (will be removed in future versions). Use '" + argv0 + " foot' instead";
	sopts["beam"]->add_standard_options();

	std::string footfn, modelfn;
	sopts["pdf"].reset(new Options(argv0 + " pdf", progdesc + " Cumulative probability density function (CPDF) generation subcommand.", version, Authorship::majuric));
	sopts["pdf"]->argument("conf").bind(input).desc("CPDF (\"sky\") configuration file (input)");
	sopts["pdf"]->argument("footprint").bind(footfn).desc(".xgpc.txt footprint file (input) ");
	sopts["pdf"]->argument("model").bind(modelfn).desc("model.XXX.conf model configuration file (input) ");
	sopts["pdf"]->argument("output").bind(output).desc("CPDF file (output) ");
	sopts["pdf"]->add_standard_options();

	std::string pdffile;
	sopts["pdfinfo"].reset(new Options(argv0 + " pdfinfo", progdesc + " Print information about a .pdf.bin file.", version, Authorship::majuric));
	sopts["pdfinfo"]->argument("pdf").bind(pdffile).desc(".pdf.bin file (input)");
	sopts["pdfinfo"]->add_standard_options();

	bool simpleOutput = true;
	sopts["catalog"].reset(new Options(argv0 + " catalog", progdesc + " Star catalog generation subcommand.", version, Authorship::majuric));
	sopts["catalog"]->argument("conf").bind(input).desc("Catalog (\"sim\") configuration file (input)");
	sopts["catalog"]->argument("pdf").bind(pdffile).desc(".pdf.bin file (input)");
	sopts["catalog"]->argument("output").bind(output).desc("Generated catalog prefix (output)");
//	sopts["catalog"]->option("s").bind(simpleOutput).addname("simple").value("true").desc("Generate simple .txt output");
	sopts["catalog"]->option("d").bind(simpleOutput).addname("dmm").value("false").desc("Generate DMM output (deprecated)");
	sopts["catalog"]->prolog =
		"If -d is specified, the output is by default stored to $(output)/uniq_objects.dmm and $(output)/uniq_observations.dmm DMM files.\n"
		"Otherwise, textual output, possibly compressed depending on file extension, is stored $(output).";
	sopts["catalog"]->add_standard_options();

/*	std::string catalog, in_module = "textin", out_module = "textout";*/
	std::vector<std::string> modules;
	std::string infile("sky.cat.txt"), outfile("sky.obs.txt");
	sopts["postprocess"].reset(new Options(argv0 + " postprocess", progdesc + " Apply postprocessing steps to catalog sources.", version, Authorship::majuric));
	sopts["postprocess"]->argument("conf").bind(input).desc("Postprocessing (\"postprocess.conf\") configuration file");
	sopts["postprocess"]->argument("modules").bind(modules).optional().gobble().desc("List of postprocessing module configuration files (or module names)");
	sopts["postprocess"]->option("i").addname("input").bind(infile).param_required().desc("Input catalog file");
	sopts["postprocess"]->option("o").addname("output").bind(outfile).param_required().desc("Output catalog file");
/*	sopts["postprocess"]->option("i").bind(in_module).addname("inmodule").param_required().desc("Input file reading module");
	sopts["postprocess"]->option("o").bind(out_module).addname("outmodule").param_required().desc("Output file writing module");*/
	sopts["postprocess"]->add_standard_options();

	//
	// Parse
	//
	Options::option_list optlist;
	parse_options(optlist, opts, argc, argv);
	if(sopts.count(cmd))
	{
		parse_options(*sopts[cmd], optlist);
	}
	else
	{
		ostringstream ss;
		ss << "Unrecognized subcommand `" << cmd << "'";
		print_options_error(ss.str(), opts);
		return -1;
	}

	/////// Start your application code here

	if(!cuda_init())
	{
		MLOG(verb1) << "Error initializing GPU acceleration. Aborting.";
		return -1;
	}

	if(cmd == "pdfinfo")
	{
		// ./simulate.x pdfinfo sky.bin.pdf
		pdfinfo(cout, pdffile);
		return 0;
	}

	ifstream in(input.c_str());
	if(!in) { THROW(EFile, "Error accessing " + input + "."); }
	if(cmd == "footprint")
	{
		Config cfg; cfg.load(in);

		// load run list
		if(!cfg.count("footprint_runs")) { CFG_THROW("Configuration key footprint_runs must be set"); }
		std::string runsFn = cfg["footprint_runs"];
		text_input_or_die(in, runsFn);
		std::set<int> runs;
		load(in, runs, 0);

		// output filename
		if(!cfg.count("footprint")) { CFG_THROW("Configuration key footprint_runs must be set"); }
		output = cfg["footprint"];

		// projection
		std::string pole; double l0, b0;
		cfg.get(pole,	"projection",	std::string("90 90"));
		std::istringstream ss(pole);
		ss >> l0 >> b0;
		lambert proj(rad(l0), rad(b0));

		// lower limit in latitude
		cfg.get(b0, "footprint_min_lat", 0.);

		makeSkyMap(runs, output, proj, rad(b0));

		return 0;
	}
	if(cmd == "beam")
	{
		Config cfg; cfg.load(in);
		return footprints::pencilBeam(cfg);
	}
	if(cmd == "foot")
	{
		Config cfg; cfg.load(in);

		if(!cfg.count("type")) { THROW(EAny, "The footprint configuration file must contain the 'type' keyword."); }

		if(cfg["type"] == "beam") { return footprints::pencilBeam(cfg); }
		if(cfg["type"] == "rect") { return footprints::rectangle(cfg); }

		THROW(EAny, "The requested footprint type " + cfg["type"] + " is unknown.");
	}
	if(cmd == "pdf")
	{
		// ./simulate.x pdf north.conf north.pdf.bin
		model_pdf pdf(in, input);
		std::ofstream oout(output.c_str());
		if(!oout) { THROW(EFile, "Cannot access " + output + " for output."); }

		pdf.construct_mpdf(footfn, modelfn);

		io::obstream out(oout);
		out << pdf;

		DLOG(verb2) << io::binary::manifest;
	}
	else if(cmd == "catalog")
	{
		// turn off GSL's error handler or else locusfitting routines
		// may barf
		gsl_set_error_handler_off();
	
		// ./simulate.x catalog sim.conf dmmwriter.conf
		sky_generator skygen(in, pdffile);

		if(!simpleOutput)
		{
			THROW(EAny, "DMM output support has been discontinued");
		}
		else
		{
			MLOG(verb1) << "Simple text file output to " << output << ".";
			flex_output out(output.c_str());
			star_output_to_textstream cat_out(out.out());
			skygen.montecarlo(cat_out);
		}
	}
	else if(cmd == "postprocess")
	{
		// turn off GSL's error handler or else locusfitting routines
		// may barf
		gsl_set_error_handler_off();

		// ./simulate.x postprocess postprocess.conf [--infile=sky.cat.txt] [--outfile=sky.obs.txt] [module1.conf [module2.conf]....]
		std::set<std::string> mset;
		mset.insert(modules.begin(), modules.end());
		postprocess_catalog(input, infile, outfile, mset);
	}
	else
	{
		THROW(ENotImplemented, "Should not get to here. I'm really confused and aborting.");
	}
	return 0;
}
catch(EAny &e)
{
	e.print();
	return -1;
}
}
