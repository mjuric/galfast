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
#define COMPILING_SIMULATE

#include <astro/io/binarystream.h>
#include <astro/math/vector.h>
#include <astro/math.h>
#include <astro/util.h>
#include <astro/sdss/rungeometry.h>
#include <astro/system/options.h>

#include <fstream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "projections.h"
#include "model.h"
#include "gpc_cpp.h"
#include "io.h"

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
//		if(inc) { cerr << "+"; }
	}
//	cerr << "\n";
	return in;
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

#include <sys/time.h>
double seconds()
{
	timeval tv;
	int ret = gettimeofday(&tv, NULL);
	assert(ret == 0);
	return double(tv.tv_sec) + double(tv.tv_usec)/1e6;
}

#if 0
void testIntersection(const std::string &prefix)
{
	lambert proj(rad(90), rad(90));
	double l = rad(300);
	double b = rad(60);

	gpc_polygon allsky;
	FILE *fp = fopen((prefix + ".gpc.txt").c_str(), "r");
	gpc_read_polygon(fp, 1, &allsky);
	fclose(fp);
	striped_polygon sky(allsky, 100);
//	gpc_polygon &sky = allsky;

	// specific point test
	double x, y;
	proj.convert(l, b, x, y);
	x = 0.272865; y = -0.387265;
	//x = 0.280607; y = -0.387462;
	x = 1.05179; y = 0.807653;

	pair<double, double> tmp(x, y);
	cerr << "Inside: " << in_polygon((gpc_vertex const&)tmp, sky) << "\n";
//exit(0);
	// monte-carlo tests
	double begin = seconds();

	srand(23);
	RunGeometryDB db;
	for(int k = 0; k != 500; k++) {
	FOREACH(db.db)
	{
		const RunGeometry &geom = (*i).second;
		Mask mask(geom);

		Radians mu0 = geom.muStart;
		Radians mu1 = geom.muStart + mask.length();

		FORj(col, 0, 6)
		{
			Radians nu0 = mask.lo(col);
			Radians nu1 = mask.hi(col);
			Radians l, b, mu, nu;
			mu = math::rnd(mu0, mu1);
			nu = math::rnd(nu0, nu1);
			coordinates::gcsgal(geom.node, geom.inc, mu, nu, l, b);
			if(b < 0) { continue; }

			proj.convert(l, b, x, y);
			pair<double, double> tmp(x, y);
			if(!in_polygon((gpc_vertex const&)tmp, sky))
			{
				cerr << "-";
// 				cerr << "Test failed for:\n";
// 				cerr << " l, b = " << deg(l) << " " << deg(b) << "\n";
// 				cerr << " mu, nu = " << deg(mu) << " " << deg(nu) << "\n";
// 				cerr << " x, y = " << x << " " << y << "\n";
// 				cerr << " run, col = " << geom.run << " " << col << "\n";
// 				cerr << " mu0, mu1 = " << deg(mu0) << " " << deg(mu1) << "\n";
// 				cerr << " nu0, nu1 = " << deg(nu0) << " " << deg(nu1) << "\n";
			}
		}
	}
		cerr << ".";
	}
	cerr << "\n";
	cerr << "time: " << seconds() - begin << "\n";
	gpc_free_polygon(&allsky);
}
#endif

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
#if 0	
	gpc_tristrip tri = { 0, NULL };
	gpc_polygon_to_tristrip(&sky, &tri);
#endif

	return tri;
}

#if 0
gpc_polygon loadSky(const std::string &prefix)
{
	gpc_polygon sky;
	FILE *fp = fopen((prefix + ".gpc.txt").c_str(), "r");
	gpc_read_polygon(fp, 1, &sky);
	fclose(fp);

	return sky;
}
#endif

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

void makeBeamMap(std::string &output, Radians l, Radians b, Radians r, Radians rhole, const lambert &proj, bool smOutput)
{
	Radians dx = rad(.1); /* polygon sampling resolution in radians */

	double x, y;
	MLOG(verb1) << "Beam towards (l, b) = " << deg(l) << " " << deg(b) << ", radius = " << deg(r) << "deg, hole = " << deg(rhole);

	std::vector<double> lx, ly, hx, hy;
	Radians dphi = dx / r;
	lambert bproj(l, b);
	Radians rprev(r);
	lambert pproj(rad(0), rad(90));
	pproj.convert(ctn::pi, ctn::pi/2. - r, x, y); r = y;
	DLOG(verb1) << rprev << " " << r << " " << x << " " << y;
	if(rhole) {
		cerr << rhole << " -> ";
		pproj.convert(ctn::pi, ctn::pi/2. - rhole, x, y); rhole = y;
		cerr << rhole << " " << x << " " << y << "\n";
	}
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
		std::cerr << "Clipping the hole.\n";
		gpc_polygon poly, hole = make_polygon(hx, hy);
		gpc_polygon_clip(GPC_DIFF, &sky, &hole, &poly);
		sky = poly;
	}

	writeGPCPolygon(output, sky, proj, smOutput);

	// free memory
	gpc_free_polygon(&sky);
}

void makeSkyMap(std::set<int> &runs, const std::string &output, const lambert &proj, Radians b0 = rad(0.))
{
	Radians dx = rad(.25); /* polygon sampling resolution in radians */
//	Radians dx = rad(1); /* polygon sampling resolution - good for footprint plots */
	cerr << "footprint = " << output << ", dx = " << dx << " radians\n";

	RunGeometryDB db;
	double x, y;
	gpc_polygon sky = {0, 0, NULL};
/*	proj.convert(rad(0.), rad(-30.), x, y);
	cerr << sqrt(x*x+y*y) << "\n";*/
	proj.convert(rad(0.), b0, x, y);
	double r = sqrt(x*x+y*y);
	cerr << "Excluding r > " << r << " from the origin of lambert projection.\n";
	gpc_polygon circle = make_circle(0., 0., r, dx);

	cerr << "Processing " << runs.size() << " runs.\n";

	int k = 0;
	FOREACH(runs)
	{
		const RunGeometry &geom = db.getGeometry(*i);
// 		cerr << geom.run << " ";
 		cerr << ".";

		//if(geom.run != 752) continue;
		//if(geom.run != 752 && geom.run != 756) continue;

		gpc_polygon rpoly = make_polygon(geom, proj, dx);
		gpc_polygon_clip(GPC_INT, &rpoly, &circle, &rpoly);
		gpc_polygon_clip(GPC_UNION, &sky, &rpoly, &sky);

// 		double A = polygon_area(rpoly);
// 		gpc_free_polygon(&rpoly);
// 
// 		int nvert = 0;
// 		FOR(0, sky.num_contours) { nvert += sky.contour[i].num_vertices; }
// 		cerr << " [" << A*sqr(deg(1)) << "] [" << sky.num_contours << " contours, " << nvert << " vertices]\n";
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

#ifdef COMPILE_SIMULATE_X

#define CFG_THROW(str) THROW(EAny, str)

#include "simulate.h"

namespace footprints
{
	std::string output;
	bool smOutput = false;

	int pencilBeam(peyton::system::Config &cfg)
	{
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
	}

	int rectangle(peyton::system::Config &cfg)
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
	}
}

//void make_skymap(partitioned_skymap &m, Radians dx, const std::string &skypolyfn);
void pdfinfo(std::ostream &out, const std::string &pdffile);

void test_kin();
void test_tags();
void test_otable();
int main(int argc, char **argv)
{
try
{
	test_kin();
	test_otable();
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

	std::string catalog, in_module = "textin", out_module = "textout";
	std::vector<std::string> modules;
	sopts["postprocess"].reset(new Options(argv0 + " postprocess", progdesc + " Apply postprocessing steps to catalog sources.", version, Authorship::majuric));
	sopts["postprocess"]->argument("conf").bind(input).desc("Postprocessing (\"postprocess.conf\") configuration file");
	sopts["postprocess"]->argument("input").bind(catalog).desc("Input catalog file");
	sopts["postprocess"]->argument("output").bind(output).desc("Output catalog file");
	sopts["postprocess"]->argument("modules").bind(modules).optional().gobble().desc("List of postprocessing module configuration files");
	sopts["postprocess"]->option("i").bind(in_module).addname("inmodule").param_required().desc("Input file reading module");
	sopts["postprocess"]->option("o").bind(out_module).addname("outmodule").param_required().desc("Output file writing module");
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
#if 0
	if(cmd == "pskymap")
	{
		partitioned_skymap sky;
		dx = rad(dx);

		make_skymap(sky, dx, input);
		{ std::ofstream f(output.c_str()); io::obstream out(f); out << sky; }

		return 0;
	}
	else
#endif
	if(cmd == "pdfinfo")
	{
		// ./simulate.x pdfinfo sky.bin.pdf
		pdfinfo(cout, pdffile);
		return 0;
	}

//	std::cerr << "cmd=" << cmd << "\n";
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
//		peyton::system::print_trace();
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
			std::cerr << "DMM output.\n";
			star_output_to_dmm cat_out(output + "/uniq_objects.dmm", output + "/uniq_observations.dmm", true);
			skygen.montecarlo(cat_out);
		}
		else
		{
			MLOG(verb1) << "Simple text file output to " << output << ".";
//			std::ofstream out(output.c_str());
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

		// ./simulate.x postprocess postprocess.conf file.in.txt file.out.txt [module1 [module2]....]
		//observe_catalog(input, catalog, output);
		modules.push_back(in_module);
		modules.push_back(out_module);
		postprocess_catalog(input, catalog, output, modules);
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
#else
int main(int argc, char **argv)
{
	cerr << "simulate.x not compiled. Do ./configure with --enable-simulate first\n";
	return -1;
}
#endif
