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
 
#ifdef HAVE_BOOST

#include <sstream>
#include <astro/system/options.h>
#include <astro/system/shell.h>
#include <astro/system/fs.h>
#include <astro/io/binarystream.h>
#include <boost/shared_ptr.hpp>

#include "dm.h"
#include "projections.h"
#include "raytrace.h"
#include "model.h"
#include "conic_volume_map.h"
#include "container.h"

#include "gpc_cpp.h"

#include <astro/useall.h>
using namespace std;

class RunMap
{
protected:
	std::map<int, int> rm;
	std::map<int, int> colmap;
public:
	RunMap(const std::string &runmap = "schlegel.cat.txt")
	{
		if(Filename(runmap).exists())
		{
			text_input_or_die(in, runmap);
			in.skip(1); // skip over the header line

			int run, end_id, c[7];
			bind(in, run, 0, c[0], 1, c[1], 2, c[2], 3, c[3], 4, c[4], 5, c[5], 6, c[6], 7);
	
			while(in.next())
			{
				rm[c[6]] = run;
				FOR(1, 7)
				{
					colmap[c[i]] = i-1;
				}
			}
		} else {
			std::cerr << "WARNING: RunMap file [" << runmap << "] was not found.\n";
		}
	}

	int operator[](int fitsId)
	{
		if(rm.size() == 0) return -1; // runmap not loaded

		ASSERT(fitsId >= 0);
		std::map<int, int>::iterator k = rm.upper_bound(fitsId);
		ASSERT(k != rm.end());
		return (*k).second;
	}
	
	int column(int fitsId, int *reloffs = NULL)
	{
		if(rm.size() == 0) // runmap not loaded
		{
			if(reloffs == NULL) { *reloffs = -1; }
			return -1;
		}

		ASSERT(fitsId >= 0);
		std::map<int, int>::iterator k = colmap.upper_bound(fitsId);
		ASSERT(k != rm.end());
		if(reloffs)
		{
			*reloffs = fitsId;

			if(colmap.begin() != k)
			{
				std::map<int, int>::iterator j = k; j--;
				*reloffs -= (*j).first;
			}
		}
		return (*k).second;
	}
};

class selector;

class driver
{
public:
	DMMArray<mobject> arr;
	DMMArray<obsv_mag> obsv_mags;

protected:
	struct selector_d
	{
		ofstream *out;
		selector *s;
		
		selector_d() : out(NULL), s(NULL) {}
	};

	list<selector_d> selectors;
	
	ticker tick;
	
	int att;
public:
	driver(const std::string &objectCat, const std::string &obsvCat) : tick(10000)
	{
/*		arr.open("dm_unique_stars.dmm", "r");
		obsv_mags.open("dm_starmags.dmm", "r");*/
		arr.open(objectCat, "r");
		obsv_mags.open(obsvCat, "r");

		arr.setmaxwindows(2);
		obsv_mags.setmaxwindows(2);
	}

	void run();
	int at() const { return att; }
	selector *newSelector(const std::string &filename);
	
	~driver();
};

RunMap rm;

void increment(std::map<int, int> &h, int b)
{
	if(h.count(b)) { h[b]++; }
	else { h[b] = 1; }
}

double pixel_area(int i, int j, int npix, int side, const gnomonic &proj);

std::string &filterName(int filter)
{
	static std::map<int, std::string> *names = new std::map<int, std::string>();
	return (*names)[filter];
}
struct nameSetter { nameSetter(int val, const char *name) { filterName(val) = name; } };

#define DEFINE_FILTER(filt, val) \
	static const int filt = val; \
	nameSetter ns_##filt(val, #filt);

DEFINE_FILTER(F_BEAM,		0x01);
DEFINE_FILTER(F_RECT,		0x02);
DEFINE_FILTER(F_RI,		0x04);
DEFINE_FILTER(F_R,		0x08);
DEFINE_FILTER(F_ML_RI,		0x10);
DEFINE_FILTER(F_ML_R,		0x20);
DEFINE_FILTER(F_GAL_RECT,	0x40);
DEFINE_FILTER(F_SPARSE,		0x80);
DEFINE_FILTER(F_RUN,		0x100);
DEFINE_FILTER(F_GR,		0x200);
DEFINE_FILTER(F_D,		0x400);
DEFINE_FILTER(T_COORD,		0x800);
DEFINE_FILTER(F_ML_GR,		0x1000);
DEFINE_FILTER(F_COLOR,		0x2000);
DEFINE_FILTER(T_DISTANCE,	0x4000);
DEFINE_FILTER(F_CART,		0x8000);
DEFINE_FILTER(F_SIGMA,		0x10000);
DEFINE_FILTER(F_LIMIT,		0x20000);
DEFINE_FILTER(T_GMIRROR,	0x40000);
DEFINE_FILTER(T_MALMQUIST_CORRECTION,	0x80000);

class selector
{
public:
/*	#define F_BEAM		0x01
	#define F_RECT		0x02
	#define F_RI		0x04
	#define F_R		0x08
	#define F_ML_RI		0x10
	#define F_ML_R		0x20
	#define F_GAL_RECT	0x40
	#define F_SPARSE	0x80
	#define F_RUN		0x100
	#define F_GR		0x200
	#define F_D		0x400
	#define T_COORD		0x800
	#define F_ML_GR		0x1000
	#define F_COLOR		0x2000
	#define T_DISTANCE	0x4000
	#define F_CART		0x8000
	#define F_SIGMA		0x10000
	#define F_LIMIT		0x20000
	#define T_GMIRROR	0x40000*/
	int filters;
	
	// supporting per-mobject fields
	Radians lon, lat;	// transformed longitude, latitude (see T_COORD)
	
	Radians ra, dec, radius;	// for select_beam
	Radians ra0, dec0, ra1, dec1;	// for select_rect
	pair<float, float> m_ri;
	pair<float, float> m_r;
	pair<float, float> m_ml_ri;
	pair<float, float> m_ml_r;
	Radians l0, b0, l1, b1;
	int nskip;
	int run;
	pair<float, float> m_gr;
	pair<float, float> m_D;
	Radians t_node, t_inc;
	enum {Equatorial, Galactic} t_gc_coordsys;
	pair<float, float> m_ml_gr;
	struct m_color_t
	{
		pair<float, float> range;
		sdss_color color;
	};
	vector<m_color_t> m_color_vec;
	struct m_sigma_t
	{
		pair<float, float> range;
		sdss_color color;
	} m_sigma;

	struct m_cart_t
	{
		bool use[3];
		V3 v0, v1;
		
		m_cart_t() { use[0] = use[1] = use[2] = false; }
	} m_cart;
	int nread, nlimit;

	driver &db;
public:
	#define OM_STARS	1
	#define OM_OBSV		2
	#define OM_MAP		3
	#define OM_CMD		4
	#define OM_VOLMAP	5
	#define OM_PLANECUT	6
	#define OM_CYLINDRICAL	7
	#define OM_CATSTATS	8
	#define OM_HISTOGRAM	9
	#define OM_CATSTATS2	10
	#define OM_PROJECTION	11
	#define OM_CONICMAP	12

	int outputMode;

	struct om_gnomonic_t
	{
		// gnomonic projectors for the sides, top and bottom
		std::auto_ptr<peyton::math::gnomonic> proj[6];

		int npix;
		// key == (x coordinate, y coordinate, side). side == 1,2,3,4,5,6 (5 and 6 are top/bottom)
		std::map<S3, zero_init<int>, less_S3> maps;
		om_gnomonic_t()
			: npix(256)
			{
				proj[0].reset(new gnomonic(0, 0)); // front
				proj[1].reset(new gnomonic(rad(90), 0)); // left
				proj[2].reset(new gnomonic(rad(180), 0)); // back
				proj[3].reset(new gnomonic(rad(270), 0)); // right
				proj[4].reset(new gnomonic(0, rad(90))); // top
				proj[5].reset(new gnomonic(0, rad(-90))); // bottom

				int side = 0;
				//double f = pixel_area(128, 128, npix, side, *proj[side]);
				double f = pixel_area(64, 64, npix, side, *proj[side]);
			}

		void bin_side(double &x, double &y, Radians l, Radians b, int side)
		{
			ASSERT(side >= 0 && side <= 3);
			proj[side]->convert(l, b, x, y);
		}

		bool bin_top_bottom(double &x, double &y, Radians l, Radians b, int &side)
		{
			int xside = b > 0 ? 4 : 5;
			proj[xside]->convert(l, b, x, y);
			if(std::abs(x) <= 1 && std::abs(y) <= 1) {
				side = xside;
				return true;
			}
			return false;
		}

		void bin(Radians l, Radians b)
		{
			// this is planetarium view from the inside, so switch
			// the direction of l
			l = 2*ctn::pi - l;

			int side = (int)(modulo(l + rad(45), 2*ctn::pi) / (ctn::pi/2.));
			double x, y;
			if(abs(b) < rad(35) || !bin_top_bottom(x, y, l, b, side))
			{
				bin_side(x, y, l, b, side);
			}

			ASSERT(std::abs(x) <= 1 && std::abs(y) <= 1)
			{
				cerr << setw(15) << setprecision(10);
				cerr << "side = " << side << "\n";
				cerr << x << " " << y << "\n";
				cerr << deg(l) << " " << deg(b) << "\n";
				
				proj[4]->convert(l, b, x, y);
				cerr << x << " " << y << "\n";
				cerr << deg(l) << " " << deg(b) << "\n";
			}

			// bin
			short i = (short)(((x + 1.) / 2.) * npix);
			short j = (short)(((y + 1.) / 2.) * npix);
			if(i == npix) { i--; }
			if(j == npix) { j--; }
			if(i == -1) { i++; }
			if(j == -1) { j++; }
			ASSERT(i >= 0 && i < npix);
			ASSERT(j >= 0 && j < npix) {
				cerr << setw(15) << setprecision(10);
				cerr << "side = " << side << "\n";
				cerr << x << " " << y << "\n";
				cerr << deg(l) << " " << deg(b) << "\n";
				
				proj[4]->convert(l, b, x, y);
				cerr << x << " " << y << "\n";
				cerr << deg(l) << " " << deg(b) << "\n";
			}
			
			// store
			maps[S3(i, j, side)]++;
		}
	} om_gnomonic;

	struct om_lambert_t
	{
		double dx;
		double x0, y0, x1, y1;
		double l0, phi1;

		double w() { return x1 - x0; }
		double h() { return y1 - y0; }

		om_lambert_t()
			: dx(.0333333),
			  x0(-2), y0(-2), x1(2), y1(2),
			  l0(90), phi1(90)
			{}
	} om_lambert;

	conic_volume_map om_conic;
	bool om_conic_doserialize;
	std::ofstream om_conic_bout;
	
	auto_ptr<galactic_model> mq;

	struct om_cyl_t
	{
	public:
		Radians phi0;

		pair<float, float> r, ri;
		double dx;
		int ndx;
	public:
		om_cyl_t() : phi0(0) {}
	} om_cyl;
	
	struct om_cmd_t
	{
	public:
//		enum { MAGNITUDE, COLOR } xtype, ytype;

		double dx, dy;
		double x0, x1, y0, y1;

		int mx, my;
		sdss_color cx, cy;
	public:
		void set(const std::string &x, const std::string &y)
		{
			cx.set(x);
			cy.set(y);
/*			//cerr << x << "|" << y << "\n";
			if(x.size() >= 2) { cx.set(x); xtype = COLOR; }
			else { mx = sdss_color::bandIdx(x[0]); xtype = MAGNITUDE; }
			
			if(y.size() >= 2) { cy.set(y); ytype = COLOR; }
			else { my = sdss_color::bandIdx(y[0]); ytype = MAGNITUDE; }*/
		}
		om_cmd_t()
		{
			set("gr", "r");
			//	0.033333, 0.1,
			//	0.066666, 0.2,
			dx = 0.0333333; dy = 0.1;
			//x0 = -0.2; x1 = 2; y0 = 14; y1 = 22;
			x0 = -100; x1 = 100; y0 = -100; y1 = 100;
		};
	} om_cmd;

	// for OM_VOLMAP and OM_PLANECUT
	binned_runset brs;
	struct om_planecut_t {
		pair<float, float> ri, r;
		plane_transformer pt;
		std::string coordsys, mfn;
		pair<double, double> x1, x2, x3, x0;
		double d1, d2, d3, d0;
		double dx;
		int ndx;
		bool earth_on_x_axis;
	} om_planecut;

	// for OM_CATSTATS2
	struct stat_elem {
		double result;
		stat_elem() : result(0.) {}
		
		virtual void finalize() {}
		virtual void print(std::ostream &out) { out << result; }
		
		virtual void add(const mobject &s) = 0;
		virtual ~stat_elem() {};
		
		static bool
		create(boost::shared_ptr<selector::stat_elem> &e,
			std::string cmd, std::istream &in);
	};
	std::vector<boost::shared_ptr<stat_elem> > om_stats2;

	// for OM_CATSTATS
	struct {
		typedef map<int, int> nappmap;
		nappmap napp;
	} om_stats;

	// for OM_HISTOGRAM
	struct om_hist_t {
		map<int, int> h;
		sdss_color c;
		float dx;
		float x0, x1;
		
		om_hist_t() : x0(-1e10), x1(1e10), dx(.02), c("ri") {}
	} om_hist;
protected:
	lambert lambert_map;
	typedef map<S2, int, less_S2> binnedmap;
	binnedmap binned;
	typedef map<S3, int, less_S3> volmap;
	volmap volume;

	int nsparsereject;
public:
	ostream &out;
	std::string filename;	// filename of the output stream (or the description, if it's not a file)
	int nselected;	// number of records selected
	std::map<int, std::map<int, int> > nrejected;

	bool running;
public:
	selector(ostream &out_, driver &db_, std::string fn_)
		: running(true), db(db_), nselected(0), filters(0), outputMode(OM_STARS), out(out_), filename(fn_), nsparsereject(0),
		  t_gc_coordsys(Equatorial)
		{
			ASSERT(!out.fail());
		}

	int select_beam(Radians ra, Radians dec, Radians radius);
	int select_rect(Radians ra0, Radians dec0, Radians ra1, Radians dec1);
	int select_gal_rect(Radians l0_, Radians b0_, Radians l1_, Radians b1_);

	bool parse(const std::string &text);
	bool parse(istream &ss);
	bool parse(const std::string &cmd, istream &ss);
	void start();
	bool select(mobject m);
	void action(mobject m);
	void finish();
};

class stat_mean : public selector::stat_elem
{
public:
	sdss_color c;
	int n;
	stat_mean(istream &in)
		: n(0)
	{
		in >> c;
		std::cerr << "Color: " << c.name << "\n";
	}
	virtual void add(const mobject &m)
	{
		double v = m.field(c);
		result += v;
		++n;
	}
	virtual void finalize()
	{
		result /= n;
	}
};

bool selector::stat_elem::create(
	boost::shared_ptr<selector::stat_elem> &e, std::string stat, std::istream &in)
{
	if(stat == "mean") { e.reset(new stat_mean(in)); }
	else
		{ die("Unknown statistical function '" + stat + "'\n"); }

	return true;
}

driver::~driver()
{
	FOREACH(selectors)
	{
		delete (*i).out;
		delete (*i).s;
	}
}

OSTREAM(const selector::om_gnomonic_t &g)
{
	out << "# npix = " << g.npix << "\n";
	out << "#";
	return out;
}

OSTREAM(const selector::om_lambert_t &l)
{
	out << "# dx = " << deg(l.dx) << " deg [" << l.dx << " rad]\n";
	out << "# (lambda0, phi1) = (" << l.l0 << ", " << l.phi1 << ")\n";
	out << "# x = [" << l.x0 << ", " << l.x1 << ")\n";
	out << "# y = [" << l.y0 << ", " << l.y1 << ")";
	return out;
}

OSTREAM(const selector::m_color_t &l)
{
	out << "# " << l.color << " = " << '[' << l.range.first << ", " << l.range.second << ")\n";
	return out;
}

OSTREAM(const selector::m_sigma_t &l)
{
	out << "# sigma(" << l.color << ") = " << '[' << l.range.first << ", " << l.range.second << ")\n";
	return out;
}

OSTREAM(const selector::om_cmd_t &l)
{
	out << "# ";
//	if(l.xtype == selector::om_cmd_t::COLOR) { out << l.cx; } else { out << sdss_color::bandName(l.mx); };
	out << l.cx;
	out << " vs. ";
//	if(l.ytype == selector::om_cmd_t::COLOR) { out << l.cy; } else { out << sdss_color::bandName(l.my); };
	out << l.cy;
	out << "\n";
	out << "# dx = " << l.dx << ", dy = " << l.dy << "\n";
	out << "# x = [" << l.x0 << ", " << l.x1 << ")\n";
	out << "# y = [" << l.y0 << ", " << l.y1 << ")";
	return out;
}

OSTREAM(const selector::om_hist_t &l)
{
	out << "# Histogram of " << l.c << " in dx = " << l.dx << " bins.\n";
	out << "# x = [" << l.x0 << ", " << l.x1 << ")";
	return out;
}

/*
http://www.truevision3d.com/phpBB2/viewtopic.php?p=41760#41760
    
        0        1       2
    |  CE      -CF      -D  | 0
M = | -BDE+AF   BDF+AE  -BC | 1
    |  ADE+BF  -ADF+BE   AC | 2

where A,B are the cosine and sine of the X-axis rotation axis,
C,D are the cosine and sine of the Y-axis rotation axis,
E,F are the cosine and sine of the Z-axis rotation axis. 
*/

void eulerAngles(Radians &phi, Radians &theta, Radians &psi, M3 &mat)
{
	double A, B, C, D, angle_x, angle_y, angle_z;
	
	angle_y = D = -asin( mat(0, 2)); /* Calculate Y-axis angle */
	C = cos( angle_y );

	if ( fabs( C ) > 0.005 ) /* Gimball lock? */
	{
		double trxa = mat(2, 2) / C; /* No, so get X-axis angle */
		double trya = -mat(1, 2) / C;

		angle_x = atan2( trya, trxa );

		trxa = mat(0, 0) / C; /* Get Z-axis angle */
		trya = -mat(0, 1) / C;

		angle_z = atan2( trya, trxa );
	}
	else /* Gimball lock has occurred */
	{
		angle_x = 0; /* Set X-axis angle to zero */

		double trxa = mat(1, 1); /* And calculate Z-axis angle */
		double trya = mat(1, 0);

		angle_z = atan2( trya, trxa );
	}

	phi = modulo(angle_x, ctn::pi2);
	theta = modulo(angle_y, ctn::pi2);
	psi = modulo(angle_z, ctn::pi2);
}

void selector::start()
{
	lambert_map = lambert(rad(om_lambert.l0), rad(om_lambert.phi1));
	
	switch(outputMode)
	{
	case OM_CONICMAP: {
		break;
		}
	case OM_CYLINDRICAL: {
		double dx = om_cyl.dx;
		int ndx = om_cyl.ndx;
		pair<float, float> ri = om_cyl.ri, r = om_cyl.r;
		float phi0 = deg(om_cyl.phi0);

		std::string fn = io::format("cache/cylindrical.%6.3f-%6.3f.%5.3f-%5.3f.dx=%6.3f.phi0=%03.0f.bin")
			<< r.first << r.second << ri.first << ri.second << dx*ndx << phi0;

		if(!Filename(fn).exists())
		{
			std::string cmd = io::format(
				"$BIN/bin_volume.x --type=cyl --phi0=%f "
				"%6.3f %6.3f %5.3f %5.3f merged.bin %f %d %s"
			) << phi0 << r.first << r.second << ri.first << ri.second << dx << ndx << fn;

			out << "# Generating volume map file using: " << cmd << "\n";
			cerr << "Calling external command:\n  " << cmd << "\n- - - - - - - -\n";
			shell(cmd);
			cerr << "- - - - - - - -\n";
		}
		ASSERT(Filename(fn).exists());

		// load volume map and copy it into a binned runset
		binned_run br;
		binary_input_or_die(in, fn);
		in >> br;
		brs.dx = br.dx;
		FOREACH(br.pixels)
		{
			const S3 &k = (*i).first;
			binned_run::pixel &p = (*i).second;
			brs.pixels[k].uniqueVolume = p.volume;
		}

		out << "# volume map loaded from " << fn << "\n";
		out << "# outputing cylindrical volume map\n";
		out << "#    phi0 = " << deg(om_cyl.phi0) << " deg\n";
		out << "#\n";
		out.flush();
		} break;
	
	case OM_PLANECUT: {
		om_planecut_t &o = om_planecut;
		Radians phi, theta, psi;
		eulerAngles(phi, theta, psi, o.pt.M);
		std::string fn = io::format("cache/plane.r=%6.3f-%6.3f,ri=%5.3f-%5.3f,dx=%6.3f,origin=%.5f,%.5f,%.5f,euler_rot=%.5f,%.5f,%.5f.bin")
			<< o.r.first << o.r.second << o.ri.first << o.ri.second << o.dx*o.ndx
			<< o.pt.t0[0] << o.pt.t0[1] << o.pt.t0[2] // earthcentric origin
			<< deg(phi) << deg(theta) << deg(psi);

		if(!Filename(fn).exists())
		{
			std::string cmd = io::format(
				"$BIN/bin_volume.x --type=plane "
				"%6.3f %6.3f %5.3f %5.3f merged.bin %f %d \"%s\" "
				"--coordsys=%s --x1=%f,%f,%f --x2=%f,%f,%f --x3=%f,%f,%f --x0=%f,%f,%f --delta=%f --earthonaxis=%d"
			)
			<< o.r.first << o.r.second << o.ri.first << o.ri.second << o.dx << o.ndx << fn
			<< o.coordsys
			<< o.d1 << o.x1.first << o.x1.second
			<< o.d2 << o.x2.first << o.x2.second
			<< o.d3 << o.x3.first << o.x3.second
			<< o.d0 << o.x0.first << o.x0.second
			<< 1e+10 << (int)o.earth_on_x_axis;

			out << "# Generating volume map file using: " << cmd << "\n";
			cerr << "Calling external command:\n  " << cmd << "\n- - - - - - - -\n";
			shell(cmd);
			cerr << "- - - - - - - -\n";
		}
		ASSERT(Filename(fn).exists());

		// load volume map created with 'bin_volume --type=plane' and
		// copy it into a binned runset
		binned_run br;
		binary_input_or_die(in, fn);
		in >> br;
		ASSERT(std::abs(br.dx/(o.dx*o.ndx) - 1) < 1e-5)
		{
			std::cerr << "br.dx = " << br.dx << "\n";
			std::cerr << "o.dx = " << o.dx << "\n";
			std::cerr << "o.ndx = " << o.ndx << "\n";
		}

		// copy binned_run structure to binned_runset structure (a backwards compat. quirk)
		brs.dx = br.dx;
		FOREACH(br.pixels)
		{
			const S3 &k = (*i).first;
			binned_run::pixel &p = (*i).second;
			brs.pixels[k].uniqueVolume = p.volume;
		}

		out << "# volume map loaded from " << fn << "\n";
		out << "# outputing plane cut\n";
		out << "#    coordsys = " << o.coordsys << "\n";
		out << "#    x1       = (" << o.d1 << ", " << o.x1.first << ", " << o.x1.second << ") == " << o.pt.x[0] << " earthcentric\n";
		out << "#    x2       = (" << o.d2 << ", " << o.x2.first << ", " << o.x2.second << ") == " << o.pt.x[1] << " earthcentric\n";
		out << "#    x3       = (" << o.d3 << ", " << o.x3.first << ", " << o.x3.second << ") == " << o.pt.x[2] << " earthcentric\n";
		out << "#    origin   = (" << o.d0 << ", " << o.x0.first << ", " << o.x0.second << ") == " << o.pt.t0 << " earthcentric\n";
		out << "#    plane    = " << o.pt.p.a << " " << o.pt.p.b << " " << o.pt.p.c << " " << o.pt.p.w << " (ax+by+cz+w=0 parameters)\n";
		out << "#    earth_on_x_axis    = " << o.earth_on_x_axis << "\n";
		if(o.mfn.size()) { out << "#    transformation matrix file = " << o.mfn << "\n"; }
		out << "#\n";
		} break;
	}
	
	out.flush();

}

int selector::select_beam(Radians ra_, Radians dec_, Radians radius_)
{
	ra = ra_; dec = dec_; radius = radius_;
	filters |= F_BEAM;
}

int selector::select_rect(Radians ra0_, Radians dec0_, Radians ra1_, Radians dec1_)
{
	ra0 = ra0_; dec0 = dec0_; ra1 = ra1_; dec1 = dec1_;
	filters |= F_RECT;
}

int selector::select_gal_rect(Radians l0_, Radians b0_, Radians l1_, Radians b1_)
{
	l0 = l0_; b0 = b0_; l1 = l1_; b1 = b1_;
	filters |= F_GAL_RECT;
}

#define FILTER(f) if(selected && (filters & (lastfilter = f)))
#define TRANSFORM(f) if(filters & f)
bool selector::select(mobject m)
{
	bool selected = true;
	int lastfilter = 0;
	int lastfilterinstance = 0;

	//print_mobject(cout, m);

	lon = rad(m.ra);
	lat = rad(m.dec);

	// rectangular coordinates
	V3 v;
	Radians l, b;
	coordinates::equgal(lon, lat, l, b);
	v.celestial(m.D, l, b);
	v.x = -v.x;
	v.y = -v.y; // convert to earthcentric galactic

	// transformation filters
	TRANSFORM(T_GMIRROR)
	{
		// mirror around the l=0 meridian
		l =  2*ctn::pi-l;
		//if(l < 0) { l += 2*ctn::pi; }
		// update lon/lat, XYZ
		coordinates::galequ(l, b, lon, lat);
		v.y = -v.y;
	}
	TRANSFORM(T_COORD)
	{
		// rotate coordinate system in coordinate system requested
		switch(t_gc_coordsys)
		{
		case Galactic:
			// node and inclination are given wrt. galactic coordinate system
			coordinates::galgcs(t_node, t_inc, lon, lat, lon, lat);
			break;
		case Equatorial:
			// node and inclination are given wrt. equatorial coordinate system
			Radians ra = lon, dec = lat;
			//lon = rad(14.99);
			//lat = rad(78.36);
			//coordinates::galequ(lon, lat, ra, dec);
			//std::cerr << deg(lon) << " " << deg(lat) << " " << deg(ra) << " " << deg(dec) << "\n";
			coordinates::equgcs(t_node, t_inc, ra, dec, lon, lat);
			//std::cerr << deg(lon) << " " << deg(lat) << "\n";
			//exit(0);
			break;
		default:
			ASSERT(0) { std::cerr << "Unknown coordinate system for 'transform' (" << t_gc_coordsys << ")\n"; }
		}
	}
	TRANSFORM(T_DISTANCE)
	{
		// recalculate distance
		const float ml_ri = m.ml_mag[1] - m.ml_mag[2];
		float Mr = paralax.Mr(ml_ri);
		double D = stardist::D(m.ml_mag[1], Mr);
		m.D = D;
	}
#if 1
	TRANSFORM(T_MALMQUIST_CORRECTION)
	{
		// modify the distance to correct for Malmquist correction
		double del2 = sqr(0.46*0.3);
		double d1 = m.D * exp(del2);
		V3 v1; v1.celestial(d1, l, b);

		double frac = mq->rho(v1.x, v1.y, v1.z, 0.) / mq->rho(v.x, v.y, v.z, 0.);
		double corr = exp(3.5 * del2) * frac;
		m.D *= corr;
	}
#endif
	// selection filters
	FILTER(F_LIMIT)
	{
		++nread;
		selected = nread <= nlimit;
		if(nread == nlimit)
		{
			running = false;
		}
	}
	FILTER(F_BEAM)
	{
		selected = coordinates::distance(lon, lat, ra, dec) < radius;
	}
	FILTER(F_RECT)
	{
		selected = coordinates::inBox(lon, lat, ra0, dec0, ra1, dec1);
	}
	FILTER(F_GAL_RECT)
	{
		double l, b;
		coordinates::equgal(rad(m.ra), rad(m.dec), l, b);
		selected = coordinates::inBox(l, b, l0, b0, l1, b1);
	}
	FILTER(F_RI)
	{
		float ri = m.ri();
		selected = m_ri.first <= ri && ri < m_ri.second;
	}
	FILTER(F_COLOR)
	{
		FOREACH(m_color_vec)
		{
			m_color_t &m_color = *i;
			float c = m.field(m_color.color);
			selected = between(c, (pair<float,float>)m_color.range);
/*			if(c > 1 && selected)
				std::cerr << m_color.range.first << " " << m_color.range.second << " " << c << "\n";*/
			if(!selected) { lastfilterinstance = i - m_color_vec.begin(); break; }
		}
	}
	FILTER(F_SIGMA)
	{
		float s = m.sigma(m_sigma.color);
		selected = between(s, (pair<float,float>)m_sigma.range);
	}
	FILTER(F_CART)
	{
		FOR(0, 3)
		{
			if(!m_cart.use[i]) continue;

			selected = between(v[i], m_cart.v0[i], m_cart.v1[i]);
			if(!selected) { break; }
		}
	}
	FILTER(F_GR)
	{
		float gr = m.gr();
		selected = m_gr.first <= gr && gr < m_gr.second;
	}
	FILTER(F_R)
	{
		float r = m.mag[1];
		selected = m_r.first <= r && r < m_r.second;
	}
	FILTER(F_D)
	{
		selected = m_D.first <= m.D && m.D < m_D.second;
	}
	FILTER(F_ML_RI)
	{
		float ri = m.ml_mag[1] - m.ml_mag[2];
		selected = m_ml_ri.first <= ri && ri < m_ml_ri.second;
	}
	FILTER(F_ML_GR)
	{
		float gr = m.ml_mag[0] - m.ml_mag[1];
		selected = m_ml_gr.first <= gr && gr < m_ml_gr.second;
	}
	FILTER(F_ML_R)
	{
		float r = m.ml_mag[1];
		selected = m_ml_r.first <= r && r < m_ml_r.second;
	}
	FILTER(F_RUN)
	{
		selected = false;
		FOR(m.obs_offset, m.obs_offset + m.n)
		{
			obsv_mag &sm = db.obsv_mags[i];
			if(rm[sm.fitsId] == run) { selected = true; break; }
		}
	}
	
	// Note: thisone has to be _last_
	FILTER(F_SPARSE)
	{
		selected = ((nselected + nsparsereject) % nskip) == 0;
		if(!selected) nsparsereject++;
	}

	if(!selected)
	{
		if(lastfilter)
		{
			// record which filter deselected this object
			if(nrejected.count(lastfilter)) { nrejected[lastfilter][lastfilterinstance]++; }
			else nrejected[lastfilter][lastfilterinstance] = 0;
		}
		return false;
	}

	action(m);
	nselected++;

	return true;
}
#undef FILTER
#undef TRANSFORM

void selector::action(mobject m)
{
//	std::cerr << "In action\n";
	
	// SM doesn't like infinities
	FOR(0, 3) { if(!isfinite(m.magErr[i])) { m.magErr[i] = 0; } }

	double l, b; V3 v;
	coordinates::equgal(rad(m.ra), rad(m.dec), l, b);
	v.celestial(m.D, l, b);
	v.x = -v.x;
	v.y = -v.y; // convert to earthcentric galactic

	switch(outputMode)
	{
	case OM_STARS:
		out << setw(12) << db.at() << " " << m;// << "\n";
		out << setw(12) << setprecision(8) << deg(l) << " ";
		out << setw(12) << setprecision(8) << deg(b) << " ";
	//	FOR(0, 5) { out << setw(8) << m.mag[(i+1) % 5] << " "; }
	//	out << setw(10) << m.color(sdss_color("p1s")) << " ";
	//	out << setw(10) << m.color(sdss_color("p2s"));
	//	out << setw(12) << " " << v;
 		out << "\n";
		break;
	case OM_OBSV:
		FOR(m.obs_offset, m.obs_offset + m.n)
		{
			obsv_mag &sm = db.obsv_mags[i];
			char buf[30]; sprintf(buf, "%6.3f", m.Ar);
			int reloffs;
			int col = rm.column(sm.fitsId, &reloffs);
			out 	<< setw(12) << db.at()
				<< setw(6) << rm[sm.fitsId]
				<< setw(3) << col
				<< setw(10) << reloffs
				<< setw(12) << setprecision(8) << m.ra
				<< setw(12) << setprecision(8) << m.dec
				<< " " << buf
				<< " " << sm << "\n";
		}
		break;
	case OM_PROJECTION: {
		om_gnomonic.bin(l, b);
		} break;
	case OM_CONICMAP: {
		om_conic.bin(l, b, m);
		break;
		};
	case OM_MAP: {
		double x, y;
		lambert_map.convert(l, b, x, y);
		if(!between(x, om_lambert.x0, om_lambert.x1)) break;
		if(!between(y, om_lambert.y0, om_lambert.y1)) break;

		S2 k(
			(int)(floor(x / om_lambert.dx)),
			(int)(floor(y / om_lambert.dx))
		);

		if(binned.count(k) == 0) { binned[k] = 1; }
		else { binned[k]++; }
		
		} break;
	case OM_VOLMAP: {
		V3 v;

		v.celestial(m.D, lon, lat);
		S3 k = floor(v / brs.dx + 0.5);

		if(brs.pixels.count(k) == 0)
		{
			cerr << k << " " << v << " not there\n";
		}
		else { brs.pixels[k].uniqueN++; }

		} break;
	case OM_PLANECUT: {
		V3 v;

		v.celestial(m.D, l, b);
		v.x = -v.x;
		v.y = -v.y; // convert to earthcentric galactic
		v = om_planecut.pt.toPlane(v); // convert to planar

		if(abs(v.z) > om_planecut.pt.delta) break; // out of plane thinkness

		S3 k = floor(v / brs.dx + 0.5);
//		std::cerr << "Binning to " << k << "\n";
		if(brs.pixels.count(k) == 0)
		{
			cerr << m.D << " " << m.ml_r() << " " << m.ml_ri() << "\n";
			cerr << k << " " << v << " not there\n";
		}
		else { brs.pixels[k].uniqueN++; }

		} break;
	case OM_CYLINDRICAL: {
		V3 v;

		v.celestial(m.D, l, b);
		v.x = Rg - v.x;
		v.y = -v.y; // convert to galactocentric galactic

		double rho = brs.dx * floor(v.rho()/brs.dx + 0.5);
		Radians phi = modulo(v.phi() - om_cyl.phi0, 2*ctn::pi);
		double rphi = phi*rho;

		S3 k = floor(V3(rho, rphi, v.z) / brs.dx + 0.5);
		if(brs.pixels.count(k) == 0)
		{
			cerr << k << " " << v << " not there\n";
		}
		else { brs.pixels[k].uniqueN++; }

		} break;
	case OM_CMD: {
		double x, y;
		x = m.field(om_cmd.cx);
		y = m.field(om_cmd.cy);
		//cerr << om_cmd.cx << " " << x << "\n";
		if(!between(x, om_cmd.x0, om_cmd.x1)) break;
		if(!between(y, om_cmd.y0, om_cmd.y1)) break;

		S2 k(
			(int)(floor(x / om_cmd.dx)),
			(int)(floor(y / om_cmd.dy))
		);

//		if(abs(x) > 40) { cerr << x << " " << y << " " << m << "\n"; }

		if(binned.count(k) == 0) { binned[k] = 1; }
		else { binned[k]++; }
		
		} break;
	case OM_CATSTATS: {
		if(om_stats.napp.count(m.n))
		{
			om_stats.napp[m.n]++;
		} else {
			om_stats.napp[m.n] = 1;
		}
		} break;
	case OM_CATSTATS2: {
		FOREACH(om_stats2)
		{
			(*i)->add(m);
		}
		}
	case OM_HISTOGRAM: {
		float x = m.field(om_hist.c);
		if(!between(x, om_hist.x0, om_hist.x1)) break;

		int b = (int)floor(x / om_hist.dx + 0.5);
		increment(om_hist.h, b);

		} break;
	}
}

extern "C"
{
	#include "gpc/gpc.h"
}

double poly_area(const gpc_vertex *v, int n)
{
	double A = 0;
	FORj(j, 0, n)
	{
		const gpc_vertex &a = v[j];
		const gpc_vertex &b = (j + 1 == n) ? v[0] : v[j+1];

		A += a.x*b.y - b.x*a.y;
	}

	A *= 0.5;
	return abs(A);
}

// calculate pixel area in rad^2 at position i,j of the gnomonic image
double pixel_area(int i, int j, int npix, int side, const gnomonic &proj)
{
	static lambert lnorth(90, 90), lsouth(90, -90);
	lambert &eqa = side == 5 ? lsouth : lnorth;

	double x = ((double)i) / npix * 2 - 1;
	double y = ((double)j) / npix * 2 - 1;
	ASSERT(-1 <= x && x <= 1);
	ASSERT(-1 <= y && y <= 1);
	double dx = 1. / npix;

	// take a small triangle around x, y in gnomonic coords.
	// and map it to lambert coords
	gpc_vertex c[3] = { {x, y}, {x+dx, y}, {x, y+dx} };
	double area0 = poly_area(c, 3);
	double l, b;
	FOR(0, 3)
	{
		proj.inverse(c[i].x, c[i].y, l, b);
		eqa.convert(l, b, c[i].x, c[i].y);
	}
	
	// return the conversion factor from unit pixel
	// in gnomonic coordinates to square degree
	double area1 = poly_area(c, 3);
	double f = area0 / area1;
	return area1;
}

void selector::finish()
{
	if(outputMode == OM_MAP)
	{
		FOREACH(binned)
		{
			int v = (*i).second;
			const S2 &k = (*i).first;

			out << (0.5 + k.x)*om_lambert.dx << " " << (0.5 + k.y)*om_lambert.dx << " " << v << "\n";
		}
	}
	else if(outputMode == OM_CONICMAP)
	{
		om_conic.text_serialize(out);
		if(om_conic_doserialize)
		{
			{
				io::obstream out(om_conic_bout);
				out << om_conic;
			}
			om_conic_bout.flush();
		}
#if 1 // debugging
		conic_volume_map_interpolator magerrs;
		magerrs.initialize(om_conic);
#endif
	}
	else if(outputMode == OM_PROJECTION)
	{
		FOREACH(om_gnomonic.maps)
		{
			double n = (*i).second;
			const S3 &k = (*i).first;
			int side = k.z;

			// convert to stars per unit square degree
			n /= pixel_area(k.x, k.y, om_gnomonic.npix, side, *om_gnomonic.proj[side]);
			n /= sqr(180/ctn::pi);

			out << k.x << " " << k.y << " " << k.z << " " << n << "\n";
		}
	}
	else if(outputMode == OM_CMD)
	{
		FOREACH(binned)
		{
			int v = (*i).second;
			const S2 &k = (*i).first;

			out << (0.5 + k.x)*om_cmd.dx << " " << (0.5 + k.y)*om_cmd.dy << " " << v << "\n";
		}
	}
	else if(outputMode == OM_VOLMAP || outputMode == OM_PLANECUT || outputMode == OM_CYLINDRICAL)
	{
		binned_runset bout;
		bout.dx = brs.dx;
		FOREACH(brs.pixels)
		{
			S3 k = (*i).first;
			binned_runset::pixel &p = (*i).second;
//			k.z = 0;
			bout.pixels[k].uniqueVolume += p.uniqueVolume;
			bout.pixels[k].uniqueN += p.uniqueN;
		}
		out << bout;
	}
	else if(outputMode == OM_CATSTATS)
	{
		out << "# Statistics:\n";
		int total = 0, stars = 0;
		typedef map<int,int> nappmap;
		FOREACH(om_stats.napp)
		{
			out << (*i).first << " " << (*i).second << "\n";
			total += (*i).first * (*i).second;
			stars += (*i).second;
		}
		out << "# total observations: " << total << "\n";
		out << "# total stars:        " << stars << "\n";
	}
	else if(outputMode == OM_CATSTATS2)
	{
		out << "# Statistics:\n";
		bool first = true;
		FOREACH(om_stats2)
		{
			if(!first) { out << "\t"; }
			else { first = false; }
			(*i)->finalize();
			(*i)->print(out);
		}
		out << "\n";
	}
	else if(outputMode == OM_HISTOGRAM)
	{
		out << "#\n";
		int first = (int)floor(om_hist.x0 / om_hist.dx + 0.5);
		int last = (int)floor(om_hist.x1 / om_hist.dx + 0.5);
		FOR(first, last)
		{
			out << i*om_hist.dx << " ";
			out << (om_hist.h.count(i) ? om_hist.h[i] : 0) << "\n";
		}

/*		FOREACH(om_hist.h)
		{
			out << (*i).first*om_hist.dx << " " << (*i).second << "\n";
		}*/
	}
	out.flush();
}

template<typename T>
bool set_param(istream &ss, ostream &out, const string &cmd, const string &param_name, T& param)
{
	if(cmd != param_name) return false;

	ss >> param; ASSERT(!ss.fail());
	out << "# " << param_name << " = " << param << "\n";
	out << "#\n";

	return true;
}

template<typename T, int N>
std::string join(const char *separator, const T array[N])
{
	std::string res = str(array[0]);
	FOR(1, N)
	{
		res += separator;
		res += str(array[i]);
	}
	return res;
}

// "quoted string" wrapper for std::string IO
class qstring
{
public:
	std::string &s;
	qstring(std::string &s_) : s(s_) {}
};

ISTREAM(const qstring &cqs)
{
	string &s = const_cast<qstring &>(cqs).s;
	in >> s;
	if(!in) return in;
	if(s[0] != '"') return in;

	std::string tmp;
	getline(in, tmp, '"');
	if(!in) return in;

	s = s.substr(1, s.size()-1) + tmp.substr(0, s.size()-1);
	return in;
}

bool selector::parse(const std::string &text)
{
	istringstream in(text.c_str());
	return parse(in);
}

bool selector::parse(istream &in)
{
	string cmd;
	char line[1000];
	while(!in.eof())
	{
		in.getline(line, sizeof(line));
		if(in.eof()) break;

		stringstream ss(line);
		ss >> cmd;
		if(cmd[0] == '#' || ss.fail()) continue;

		if(!parse(cmd, ss))
		{
			// this selector signaled it's finished with parsing
			return false;
		}
	}
	return true;
}

bool selector::parse(const std::string &cmd, istream &ss)
{
	if(cmd == "end")
	{
		return false;
	}

	if(cmd == "beam")
	{
		double ra, dec, radius;
		ss >> ra >> dec >> radius; ASSERT(!ss.fail());
		select_beam(rad(ra), rad(dec), rad(radius));
		out << "# beam filter active\n";
		out << "# ra     = " << ra << "\n";
		out << "# dec    = " << dec << "\n";
		out << "# radius = " << radius << "\n";
		out << "#\n";
	}
	else if(cmd == "rect")
	{
		double ra0, dec0, ra1, dec1;
		ss >> ra0 >> dec0 >> ra1 >> dec1; ASSERT(!ss.fail());
		select_rect(rad(ra0), rad(dec0), rad(ra1), rad(dec1));
		out << "# rect filter active\n";
		out << "# (ra0, dec0)    = " << ra0 << " " << dec0 << "\n";
		out << "# (ra1, dec1)    = " << ra1 << " " << dec1 << "\n";
		out << "#\n";
	}
	else if(cmd == "gal_rect")
	{
		double l0, b0, l1, b1;
		ss >> l0 >> b0 >> l1 >> b1; ASSERT(!ss.fail());
		select_gal_rect(rad(l0), rad(b0), rad(l1), rad(b1));
		out << "# gal_rect filter active\n";
		out << "# (l0, b0)    = " << l0 << " " << b0 << "\n";
		out << "# (l1, b1)    = " << l1 << " " << b1 << "\n";
		out << "#\n";
	}
	else if(cmd == "gal_beam")
	{
		double l, b, radius;
		ss >> l >> b >> radius; ASSERT(!ss.fail());
		Radians ra, dec;
		coordinates::galequ(rad(l), rad(b), ra, dec);
		select_beam(ra, dec, rad(radius));
		out << "# gal_beam filter active\n";
		out << "# l     = " << l << "\n";
		out << "# b    = " << b << "\n";
		out << "# radius = " << radius << "\n";
		out << "#\n";
	}
	else if(cmd == "obsv")
	{
		outputMode = OM_OBSV;
		out << "# outputing observations\n";
		out << "#\n";
	}
	else if(cmd == "projection")
	{
		// projection gnomonic [dx=256]
		outputMode = OM_PROJECTION;
		int npix; std::string proj;

		ss >> proj;
		ASSERT(!ss.fail());
		ASSERT(proj == "gnomonic");
		ss >> npix;
		if(!ss.fail()) { om_gnomonic.npix = npix; }

		out << "# outputing binned 6-sided gnomonic map\n";
		out << om_gnomonic << "\n";
		out << "#\n";
	}
	else if(cmd == "conicmap")
	{
		// conicmap <dx> <d0> <d1> <starsperpix> <galaxy.conf> <north.pskymap.bin> <south.pskymap.bin>
		outputMode = OM_CONICMAP;
		ss >> om_conic.dx >> om_conic.d0 >> om_conic.d1 >> om_conic.starsperpix;
		ss >> om_conic.northfn >> om_conic.southfn >> om_conic.modelfn;
		ASSERT(!ss.fail());

		RAD(om_conic.dx);
		out << "# outputing binned conical map with variable radial pixel lengths\n";
		out << om_conic << "\n";
		out << "#\n";
		om_conic.initialize_for_galmodel();
		out << "# Conic support files load successfull.\n";
		out << "# " << om_conic.sky_north.size() << " beams (north).\n";
		out << "# " << om_conic.sky_south.size() << " beams (south).\n";
		out << "# " << om_conic.npix << " pixels allocated.\n";
		out << "#\n";
	}
	else if(cmd == "conicmap2")
	{
		const char *errmsg = 
		  "conicmap2 <north.pskymap.bin> <south.pskymap.bin>    <dfield> <d0> <d1> <dd>   <vfield>  [binary_output_file]";
		ss >> om_conic.northfn >> om_conic.southfn;
		ss >> om_conic.d_field;
		ss >> om_conic.d0 >> om_conic.d1 >> om_conic.dd;
		ss >> om_conic.val_field;
		ASSERT(!ss.fail());
		std::string om_conic_binaryoutputfile;
		ss >> om_conic_binaryoutputfile;

		outputMode = OM_CONICMAP;
		//om_conic.configure_binning(om_conic.d_field, om_conic.val_field, &conic_volume_map::cb_equidistant);
		om_conic.configure_binning(om_conic.d_field, om_conic.val_field,
			&conic_volume_map::cb_equidistant,
			&conic_volume_map::bb_bin_collect,
			&conic_volume_map::fin_median);
		om_conic.initialize();

		out << "# outputing binned conical map with variable radial pixel lengths\n";
		out << om_conic << "\n";
		if(om_conic_doserialize = om_conic_binaryoutputfile.size())
		{
			out << "# serializing conical map to " << om_conic_binaryoutputfile << " file\n";
			om_conic_bout.open(om_conic_binaryoutputfile.c_str());
			ASSERT(om_conic_bout.good());
		}
		out << "#\n";
	}
	else if(cmd == "malmquist")
	{
		const char *errmsg = "malmquist <galaxy.conf>";
		std::string modelfn;
		ss >> modelfn;
		ASSERT(!ss.fail());

		std::ifstream in(modelfn.c_str());
		ASSERT(!in.fail());
		mq.reset(galactic_model::load(in));
		ASSERT(mq.get() != NULL);

		filters |= T_MALMQUIST_CORRECTION;

		out << "# correcting Malmquist bias using " << modelfn << " as model.\n";
		out << "#\n";
	}
	else if(cmd == "lambert")
	{
		outputMode = OM_MAP;
		double l0, phi1, dx;
		ss >> l0 >> phi1;
		if(!ss.fail()) { om_lambert.l0 = l0; om_lambert.phi1 = phi1; }
		ss >> dx;
		if(!ss.fail()) { om_lambert.dx = rad(dx); }
		out << "# outputing binned lambert map\n";
		out << om_lambert << "\n";
		out << "#\n";
	}
	else if(set_param(ss, out, cmd, "lambert.l0", om_lambert.l0)) {}
	else if(set_param(ss, out, cmd, "lambert.phi1", om_lambert.phi1)) {}
	else if(cmd == "cmd")
	{
		// cmd [<xcolor|xmag> <ycolor|ymag> [dx dy [x0 x1 y0 y1]]]
		outputMode = OM_CMD;
		std::string x, y; // where x and y can be colors or magnitudes (eg. "gr" or "r")
		ss >> x >> y;
		if(!ss.fail())
		{
			om_cmd.set(x, y);
			float dx, dy;
			ss >> dx >> dy;
			if(!ss.fail())
			{
				om_cmd.dx = dx;
				om_cmd.dy = dy;
				float c0, c1, m0, m1;
				ss >> c0 >> c1 >> m0 >> m1;
				if(!ss.fail())
				{
					om_cmd.x0 = c0;
					om_cmd.x1 = c1;
					om_cmd.y0 = m0;
					om_cmd.y1 = m1;
				}
			}
		}
		out << "# outputing color-magnitude diagram\n";
		out << om_cmd << "\n";
		out << "#\n";
	}
	else if(set_param(ss, out, cmd, "cmd.dx", om_cmd.dx)) {}
	else if(set_param(ss, out, cmd, "cmd.dy", om_cmd.dy)) {}
	else if(set_param(ss, out, cmd, "cmd.x0", om_cmd.x0)) {}
	else if(set_param(ss, out, cmd, "cmd.x1", om_cmd.x1)) {}
	else if(set_param(ss, out, cmd, "cmd.y0", om_cmd.y0)) {}
	else if(set_param(ss, out, cmd, "cmd.y1", om_cmd.y1)) {}
	else if(cmd == "volmap")
	{
		std::string fn;
		ss >> fn; ASSERT(!ss.fail());

		// load volume map created with 'bin_volume --type=gc' and
		// copy it into a binned runset
		binned_run br;
		binary_input_or_die(in, fn);
		in >> br;
		brs.dx = br.dx;
		FOREACH(br.pixels)
		{
			const S3 &k = (*i).first;
			binned_run::pixel &p = (*i).second;
			brs.pixels[k].uniqueVolume = p.volume;
		}

		outputMode = OM_VOLMAP;
		out << "# volume map loaded from " << fn << "\n";
		out << "# outputing volume map\n";
		out << "#\n";
	}
	else if(cmd == "statistics")
	{
		outputMode = OM_CATSTATS;
		out << "# outputing statistics\n";
		out << "#\n";
	}
	else if(cmd == "stat")
	{
		outputMode = OM_CATSTATS2;
		out << "# outputing statistics\n";
		out << "#\n";
		
		std::string stat;
		ss >> stat;
		boost::shared_ptr<stat_elem> se;
		stat_elem::create(se, stat, ss);
		om_stats2.push_back(se);
	}
	else if(cmd == "histogram")
	{
		// histogram [<xcolor|xmag> [dx [x0 x1]]]
		outputMode = OM_HISTOGRAM;
		std::string x; // where x can be color or magnitude
		ss >> x;
		if(!ss.fail())
		{
			om_hist.c = x;
			float dx;
			ss >> dx;
			if(!ss.fail())
			{
				om_hist.dx = dx;
				float x0, x1;
				ss >> x0 >> x1;
				if(!ss.fail())
				{
					om_hist.x0 = x0;
					om_hist.x1 = x1;
				}
			}
		}
		out << om_hist << "\n";
		out << "#\n";
	}
	else if(cmd == "cylindrical")
	{
		// Syntax: cylindrical ri0 ri1 r0 r1 dx ndx [phi0]
		// volume map file will be autogenerated and stored to cache/ directory
		std::string fn;
		ss >> om_cyl.r.first >> om_cyl.r.second;
		ss >> om_cyl.ri.first >> om_cyl.ri.second;
		ss >> om_cyl.dx >> om_cyl.ndx; ASSERT(!ss.fail());
		ss >> om_cyl.phi0;
		cerr << om_cyl.phi0 << "\n";
		om_cyl.phi0 = rad(om_cyl.phi0);

		// push automatically volume and RI limits
		double d0, d1;
		paralax.distance_limits(d0, d1, om_cyl.ri.first, om_cyl.ri.second, om_cyl.r.first, om_cyl.r.second);
		std::string text  = io::format("color RI %.10f %.10f\n") << om_cyl.ri.first << om_cyl.ri.second;
		            text += io::format("D %.10f %.10f\n") << d0 << d1;
		parse(text);
		out.flush();

		outputMode = OM_CYLINDRICAL;
	}
	else if(cmd == "import")
	{
		// import other_file.sel [cont|stop]
		std::string fn, ifnotexist;
		ss >> qstring(fn);
		ASSERT(!ss.fail());
		ss >> ifnotexist;

		ifstream in(fn.c_str());
//		if(!Filename(fn).exists())
		if(!in)
		{
			if(ifnotexist == "cont")
			{
				std::cerr << "Error opening include file [" << fn << "].\n";
				abort();
			} else {
				out << "# File [" << fn << "] does not exist. Ignoring\n";
			}
		} else {
			out << "# Importing file [" << fn << "].\n";
			out << "# >-----------------\n#\n";
			
			parse(in);
	
			out << "# -----------------<\n";
			out << "# Import of [" << fn << "] complete.\n";
		}
		out << "#\n";
}
	else if(cmd == "planecut")
	{
		//
		// Syntax: planecut <ri[2]> <r[2]> <coordsys> <dx> <ndx> [<x0> <x1> <x2> <origin> <earth_on_x_axis> [matrixoutputfile]]
		// x0, x1, x2, origin(==x0) are (dist, ra, dec) coordinates if coordsys == equ
		// If origin info is not specified, the binning coordinate system is galactocentric
		// and coordsys has to be equal to 'galcart'
		//
		om_planecut_t &o = om_planecut;
		ss
			>> o.r.first >> o.r.second
			>> o.ri.first >> o.ri.second
			>> o.coordsys >> o.dx >> o.ndx
			>> o.d1;

		if(!ss.fail())
		{
			ss
				>> o.x1.first >> o.x1.second
				>> o.d2 >> o.x2.first >> o.x2.second
				>> o.d3 >> o.x3.first >> o.x3.second
				>> o.d0 >> o.x0.first >> o.x0.second
				>> o.earth_on_x_axis;
			ASSERT(!ss.fail());
			ss	>> o.mfn;
		} else {
			ASSERT(o.coordsys == "galcart");
			// earthcentric galactic coordinate system plane
			istringstream ss("0 0 0    1 0 0   0 1 0   8000 0 0  0  dbg.matrix.txt");
			ss
				>> o.d1 >> o.x1.first >> o.x1.second
				>> o.d2 >> o.x2.first >> o.x2.second
				>> o.d3 >> o.x3.first >> o.x3.second
				>> o.d0 >> o.x0.first >> o.x0.second
				>> o.earth_on_x_axis;
			ss	>> o.mfn;
		}

		// initialize plane_transformer
		setupPlaneTransformer(o.pt, o.coordsys, o.d1, o.x1, o.d2, o.x2, o.d3, o.x3, o.d0, o.x0, 1e10, o.earth_on_x_axis);

		if(o.mfn.size()) // store matrix and origin - for reading with SM
		{
			text_output_or_die(mout, o.mfn);
			mout << "# M_0 M_1 M_2 t0" << nl();
			FOR(0, 3) { mout << o.pt.M(0)[i] << o.pt.M(1)[i] << o.pt.M(2)[i] << o.pt.t0[i] << nl(); }
		}

		// push automatically volume and RI limits
		double d0, d1;
		paralax.distance_limits(d0, d1, o.ri.first, o.ri.second, o.r.first, o.r.second);
		std::string text  = io::format("color RI %.10f %.10f\n") << o.ri.first << o.ri.second;
		            text += io::format("D %.10f %.10f\n") << d0 << d1;
		parse(text);
		out.flush();

		outputMode = OM_PLANECUT;
	}
	else if(cmd == "sparse")
	{
		ss >> nskip; ASSERT(!ss.fail());
		filters |= F_SPARSE;
		out << "# outputing every " << nskip << " star\n";
		out << "#\n";
	}
	else if(cmd == "limit")
	{
		nread = 0;
		ss >> nlimit; ASSERT(!ss.fail());
		filters |= F_LIMIT;
		out << "# limiting scan to " << nlimit << " stars\n";
		out << "#\n";
	}
	else if(cmd == "run")
	{
		ss >> run; ASSERT(!ss.fail());
		filters |= F_RUN;
		out << "# run filter active\n";
		out << "# run   = " << run << "\n";
		out << "#\n";
	}
	else if(cmd == "color")
	{
		std::string color;
		double c1, c2;
		ss >> color >> c1 >> c2; ASSERT(!ss.fail());
		m_color_t m_color;
		m_color.color.set(color);
		m_color.range = make_pair(c1, c2);

		m_color_vec.push_back(m_color);
		filters |= F_COLOR;

		out << "# color filter active\n";
		out << m_color;
		out << "#\n";
	}
	else if(cmd == "sigma")
	{
		std::string color;
		double s1, s2;
		ss >> color >> s1 >> s2; ASSERT(!ss.fail());
		m_sigma.color.set(color);
		m_sigma.range = make_pair(s1, s2);
		filters |= F_SIGMA;
		out << "# sigma filter active\n";
		out << m_sigma;
		out << "#\n";
	}
	else if(cmd == "ri")
	{
		double ri0, ri1;
		ss >> ri0 >> ri1; ASSERT(!ss.fail());
		m_ri = make_pair(ri0, ri1);
		filters |= F_RI;
		out << "# ri filter active\n";
		out << "# [ri0, ri1)    = " << ri0 << " " << ri1 << "\n";
		out << "#\n";
	}
	else if(cmd == "gr")
	{
		double gr0, gr1;
		ss >> gr0 >> gr1; ASSERT(!ss.fail());
		m_gr = make_pair(gr0, gr1);
		filters |= F_GR;
		out << "# gr filter active\n";
		out << "# [gr0, gr1)    = " << gr0 << " " << gr1 << "\n";
		out << "#\n";
	}
	else if(cmd == "r")
	{
		double r0, r1;
		ss >> r0 >> r1; ASSERT(!ss.fail());
		m_r = make_pair(r0, r1);
		filters |= F_R;
		out << "# r filter active\n";
		out << "# [r0, r1)    = " << r0 << " " << r1 << "\n";
		out << "#\n";
	}
	else if(cmd == "D")
	{
		float D0, D1;
		ss >> D0 >> D1; ASSERT(!ss.fail());
		m_D = make_pair(D0, D1);
		filters |= F_D;
		out << "# D filter active\n";
		out << "# [D0, D1)    = " << D0 << " " << D1 << "\n";
		out << "#\n";
	}
	else if(cmd == "cartesian")
	{
		char coord; double a, b;
		ss >> coord >> a >> b; ASSERT(!ss.fail());
 		int idx = coord - 'x'; ASSERT(idx >= 0 && idx < 3);
		m_cart.use[idx] = true;
		m_cart.v0[idx] = a;
		m_cart.v1[idx] = b;
		filters |= F_CART;
		out << "# earthcentric cartesian " << coord << " filter active\n";
		out << "# [" << coord << "0, " << coord << "1)    = " << a << " " << b << "\n";
		out << "#\n";
	}
	else if(cmd == "ml_ri")
	{
		double ri0, ri1;
		ss >> ri0 >> ri1; ASSERT(!ss.fail());
		m_ml_ri = make_pair(ri0, ri1);
		filters |= F_ML_RI;
		out << "# ml_ri filter active\n";
		out << "# [ri0, ri1)    = " << ri0 << " " << ri1 << "\n";
		out << "#\n";
	}
	else if(cmd == "ml_gr")
	{
		double gr0, gr1;
		ss >> gr0 >> gr1; ASSERT(!ss.fail());
		m_ml_gr = make_pair(gr0, gr1);
		filters |= F_ML_GR;
		out << "# ml_gr filter active\n";
		out << "# [gr0, gr1)    = " << gr0 << " " << gr1 << "\n";
		out << "#\n";
	}
	else if(cmd == "ml_r")
	{
		double r0, r1;
		ss >> r0 >> r1; ASSERT(!ss.fail());
		m_ml_r = make_pair(r0, r1);
		filters |= F_ML_R;
		out << "# ml_r filter active\n";
		out << "# [r0, r1)    = " << r0 << " " << r1 << "\n";
		out << "#\n";
	}
	else if(cmd == "gal_mirror")
	{
		filters |= T_GMIRROR;
		out << "# Galactic l=0 line mirror active\n";
		out << "#\n";
	}
	else if(cmd == "transform")
	{
		double node, inc;
		ss >> node >> inc; ASSERT(!ss.fail());

		std::string coordsys;
		if(!(ss >> coordsys)) { coordsys = "equ"; }
		if(coordsys == "equ") { t_gc_coordsys = selector::Equatorial; }
		else if(coordsys == "gal") { t_gc_coordsys = selector::Galactic; }
		else { ASSERT(0) { std::cerr << "Unknown coordinate system [" << coordsys << "] for 'transform'\n"; } }

		t_node = rad(node); t_inc = rad(inc);
		filters |= T_COORD;
		out << "# GC transformation active\n";
		out << "# {node, inc, coord_system}    = " << node << ", " << inc << ", " << coordsys << "\n";
		out << "#\n";
	}
	else if(cmd == "paralax")
	{
		ASSERT(!paralax.paralaxModified())
		{
			std::cerr << "Paralax already changed but you're attempting to change it again.\n";
		}
		vector<double> Mr;
		double coef;
		while(ss >> coef) { Mr.push_back(coef); }
		ASSERT(!ss.fail());
		paralax.setParalaxCoefficients(Mr);
	}
	else
	{
		cerr << "Error - unknown command: " << cmd << "\n";
		exit(-1);
	}

	if(paralax.paralaxModified() && !(filters & T_DISTANCE))
	{
		// If paralax was changed in any way (eg., through PARALAX env.var.)
		// force distance recalculation
		//
		filters |= T_DISTANCE;
		out << "# Paralax relation changed to:\n";
		out << "#   Mr(i)  = " << join(", ", paralax.getParalaxCoefficients()) << "\n";
		out << "# Distances will be recalculated.\n";
		out << "#\n";
	}
	out.flush();
	return true;
}

selector *driver::newSelector(const std::string &filename)
{
	selector_d s;
	if(filename == "-")
	{
		s.s = new selector(cout, *this, "<STDOUT>");
	} else {
		s.out = new ofstream(filename.c_str());
		ASSERT(!s.out->fail());
		s.s = new selector(*s.out, *this, filename);
	}

	selectors.push_back(s);

	return s.s;
}

void driver::run()
{
	FOREACHj(j, selectors)
	{
		(*j).s->start();
	}

	for(att = 0; att != arr.size(); att++)
	{
		tick.tick();

		bool running = false;
		FOREACHj(j, selectors)
		{
			if((*j).s->running)
			{
				(*j).s->select(arr[att]);
				running = true;
			}
		}
		if(!running) { break; }
	}

	FOREACHj(j, selectors)
	{
		(*j).s->finish();
		std::cerr << (*j).s->filename << " : " << (*j).s->nselected << " records selected.\n";
		FOREACH((*j).s->nrejected)
		{
			FOREACHj(k, (*i).second)
			{
				std::cerr << "\t" << (*k).second << " rejected due to " << filterName((*i).first) << ", instance " << (*k).first << "\n";
			}
		}
	}
}

std::string shift(int &argc, char **argv)
{
	std::string s;
	if(argc <= 1) return s;
	
	s = argv[1];
	FOR(2, argc)
	{
		argv[i-1] = argv[i];
	}
	--argc;
	return s;
}

int make_run_plots(const set<int> &runs);

int main(int argc, char **argv)
{
#if 0
	std::string fn = "magerr.conicmap.bin";
	std::ifstream iin(fn.c_str()); ASSERT(iin) { std::cerr << "Could not open magnitude error map " << fn << "\n"; }
	io::ibstream in(iin);

	conic_volume_map magerrs_cvm;
	if(!(in >> magerrs_cvm)) { ASSERT(in >> magerrs_cvm); }

	conic_volume_map_interpolator magerrs;
	magerrs.initialize(magerrs_cvm);

	ofstream out("conic_interp_dump.txt");
	magerrs.text_dump(out);

	return 0;
#endif

#if 0
	set<int> runs;
	text_input_or_die (in, "catalogs/runs.txt");
	load(in, runs, 0);
	make_run_plots(runs);
	return -1;
#endif

#if 0
	runMap rm;
	cout << rm[atof(argv[1])] << "\n";
	return -1;
#endif

#if 0
	// some completely random code, to be deleted later
	// used to determine (ra,dec) of galactic coord.sys. ascending node
	Radians ra, dec, l, b;
	l = rad(33); b = 0;
	coordinates::galequ(l, b, ra, dec);
	cout << setprecision(15) << deg(ra) << " " << setprecision(15) << deg(dec) << "\n";
	return -1;
#endif
try
{
	gsl_set_error_handler_off ();
	
	VERSION_DATETIME(version, "$Id: selector.cpp,v 1.19 2006/12/11 05:19:55 mjuric Exp $");
	Options opts(
		argv[0],
		"Unique object & observation database query tool.",
		version, Authorship::majuric
	);

	std::string obj_cat_fn = "uniq_objects.dmm",
		obs_cat_fn = "uniq_observations.dmm",
		query_fn = "-";

	opts.argument("queryFile").bind(query_fn).optional().desc("Query specification file or - for stdin");
	opts.option("s").bind(obj_cat_fn).param_required().addname("objCat").desc("Object catalog");
	opts.option("o").bind(obs_cat_fn).param_required().addname("obsCat").desc("Observation catalog");

	parse_options(opts, argc, argv);

	// open the "database"
	driver db(obj_cat_fn, obs_cat_fn);

	// load the query
	ifstream qf;
	if(query_fn != "-") { qf.open(query_fn.c_str()); }
	if(!qf.good()) { cerr << "Cannot open [" << query_fn << "]. Aborting\n"; exit(-1); }

	istream &in = query_fn != "-" ? qf : cin;

	// some information
	cerr << "Object catalog:       " << obj_cat_fn << "\n";
	cerr << "Observations catalog: " << obs_cat_fn << "\n";
	cerr << "Query file:           " << query_fn << "\n";
	
	string cmd;
	char line[1000];
	selector *s = NULL;
	while(!in.eof())
	{
		in.getline(line, sizeof(line));
		if(in.eof()) break;

		stringstream ss(line);
		ss >> cmd;
		if(cmd[0] == '#' || ss.fail()) continue;

		if(cmd == "new")
		{
			std::string filename;
			ss >> filename; ASSERT(!ss.fail());
			s = db.newSelector(filename);
		}
		else
		{
			if(s == NULL)	// autocreate default selector
			{
				s = db.newSelector("-");
			}
			if(!s->parse(cmd, ss))
			{
				// this selector signaled it's finished with parsing
				s = NULL;
			}
		}
	}

	// run query
	db.run();
}
catch(EAny &e)
{
	e.print();
}
}

#else // HAVE_BOOST
#include <iostream>
int main()
{
	std::cerr << "Compiled without Boost support. selector.x inoperable.\n";
	return -1;
}

#endif
