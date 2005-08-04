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

#include <sstream>
#include <astro/system/options.h>
#include <astro/system/shell.h>
#include <astro/system/fs.h>

#include "dm.h"
#include "projections.h"
#include "raytrace.h"

#include <astro/useall.h>
using namespace std;

class RunMap
{
protected:
	std::map<int, int> rm;
	std::map<int, int> colmap;
public:
	RunMap()
	{
		text_input_or_die(in, "schlegel.cat.txt");
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
	}

	int operator[](int fitsId)
	{
		ASSERT(fitsId >= 0);
		std::map<int, int>::iterator k = rm.upper_bound(fitsId);
		ASSERT(k != rm.end());
		return (*k).second;
	}
	
	int column(int fitsId, int *reloffs = NULL)
	{
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
	DMMArray<starmag> starmags;

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
	driver() : tick(10000)
	{
		arr.open("dm_unique_stars.dmm", "r");
		starmags.open("dm_starmags.dmm", "r");
	}

	void run();
	int at() const { return att; }
	selector *newSelector(const std::string &filename);
	
	~driver();
};

RunMap rm;
#if 0
static char band_names[] = "grizu";

	#define MAG_G	0
	#define MAG_R	1
	#define MAG_I	2
	#define MAG_Z	3
	#define MAG_U	4
#endif

void increment(std::map<int, int> &h, int b)
{
	if(h.count(b)) { h[b]++; }
	else { h[b] = 1; }
}

class selector
{
public:
	#define F_BEAM		0x01
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
	pair<float, float> m_ml_gr;
	struct m_color_t
	{
		pair<float, float> range;
		sdss_color color;
	} m_color;
	struct m_sigma_t
	{
		pair<float, float> range;
		sdss_color color;
	} m_sigma;
	plx_gri_locus paralax;
	struct m_cart_t
	{
		bool use[3];
		V3 v0, v1;
		
		m_cart_t() { use[0] = use[1] = use[2] = false; }
	} m_cart;

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
	
	int outputMode;
	
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
	int n;	// number of records selected

	lambert lambert_map;
	typedef map<S2, int, less_S2> binnedmap;
	binnedmap binned;
	typedef map<S3, int, less_S3> volmap;
	volmap volume;

	ostream &out;
public:
	selector(ostream &out_, driver &db_)
		: db(db_), n(0), filters(0), outputMode(OM_STARS), out(out_)
		{
			ASSERT(!out.fail());
		}

	int select_beam(Radians ra, Radians dec, Radians radius);
	int select_rect(Radians ra0, Radians dec0, Radians ra1, Radians dec1);
	int select_gal_rect(Radians l0_, Radians b0_, Radians l1_, Radians b1_);

	int record_count() { return n; }	

	bool parse(const std::string &text);
	bool parse(const std::string &cmd, istream &ss);
	void start();
	bool select(mobject m);
	void action(mobject m);
	void finish();
};

driver::~driver()
{
	FOREACH(list<selector_d>::iterator, selectors)
	{
		delete (*i).out;
		delete (*i).s;
	}
}

OSTREAM(const selector::om_lambert_t &l)
{
	out << "# dx = " << l.dx << "\n";
	out << "# (lambda0, phi1) = (" << l.l0 << ", " << l.phi1 << ")\n";
	out << "# x = [" << l.x0 << ", " << l.x1 << ")\n";
	out << "# y = [" << l.y0 << ", " << l.y1 << ")";
	return out;
}

OSTREAM(const sdss_color &c)
{
/*	if(c.first < 0) return out << c.name;
	return out << sdss_color::bandName(c.first) << "-" << sdss_color::bandName(c.second);*/
	return out << c.name;
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
				"~/projects/galaxy/optimized/src/bin_volume.x --type=cyl --phi0=%f "
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
		FOREACH(binned_run::pixelmap::iterator, br.pixels)
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
		std::string fn = io::format("cache/plane.r=%6.3f-%6.3f,ri=%5.3f-%5.3f,dx=%6.3f,euler_rot=%.5f,%.5f,%.5f.bin")
			<< o.r.first << o.r.second << o.ri.first << o.ri.second << o.dx*o.ndx
			<< deg(phi) << deg(theta) << deg(psi);

		if(!Filename(fn).exists())
		{
			std::string cmd = io::format(
				"~/projects/galaxy/optimized/src/bin_volume.x --type=plane "
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
		brs.dx = br.dx;
		FOREACH(binned_run::pixelmap::iterator, br.pixels)
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

#define FILTER(f) if(selected && (filters & f))
#define TRANSFORM(f) if(filters & f)
bool selector::select(mobject m)
{
	bool selected = true;

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
	TRANSFORM(T_COORD)
	{
		// rotate coordinate system in ra/dec
		coordinates::equgcs(t_node, t_inc, lon, lat, lon, lat);
	}
	TRANSFORM(T_DISTANCE)
	{
		// recalculate distance
		const float ml_ri = m.ml_mag[1] - m.ml_mag[2];
		float Mr = paralax.Mr(ml_ri);
		double D = stardist::D(m.ml_mag[1], Mr);
		m.D = D;
	}

	// selection filters
	FILTER(F_SPARSE)
	{
		selected = (n % nskip) == 0;
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
		float c = m.color(m_color.color);
		selected = between(c, (pair<float,float>)m_color.range);
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
			starmag &sm = db.starmags[i];
			if(rm[sm.fitsId] == run) { selected = true; break; }
		}
	}

	if(!selected) return false;

	action(m);
	return true;
}
#undef FILTER
#undef TRANSFORM

void selector::action(mobject m)
{
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
			starmag &sm = db.starmags[i];
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
	case OM_HISTOGRAM: {
		float x = m.color(om_hist.c);
		if(!between(x, om_hist.x0, om_hist.x1)) break;

		int b = (int)floor(x / om_hist.dx + 0.5);
		increment(om_hist.h, b);

		} break;
	}
}

void selector::finish()
{
	if(outputMode == OM_MAP)
	{
		FOREACH(binnedmap::iterator, binned)
		{
			int v = (*i).second;
			const S2 &k = (*i).first;

			out << (0.5 + k.x)*om_lambert.dx << " " << (0.5 + k.y)*om_lambert.dx << " " << v << "\n";
		}
	}
	else if(outputMode == OM_CMD)
	{
		FOREACH(binnedmap::iterator, binned)
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
		FOREACH(binned_runset::pixelmap::iterator, brs.pixels)
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
		FOREACH(nappmap::iterator, om_stats.napp)
		{
			out << (*i).first << " " << (*i).second << "\n";
			total += (*i).first * (*i).second;
			stars += (*i).second;
		}
		out << "# total observations: " << total << "\n";
		out << "# total stars:        " << stars << "\n";
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

bool selector::parse(const std::string &text)
{
	istringstream in(text.c_str());

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
	else if(cmd == "lambert")
	{
		outputMode = OM_MAP;
		double l0, phi1;
		ss >> l0 >> phi1;
		if(!ss.fail()) { om_lambert.l0 = l0; om_lambert.phi1 = phi1; }
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
		FOREACH(binned_run::pixelmap::iterator, br.pixels)
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
	else if(cmd == "planecut")
	{
		//
		// Syntax: planecut <ri[2]> <r[2]> <coordsys> <dx> <ndx> <x0> <x1> <x2> <origin> <earth_on_x_axis> [matrixoutputfile]
		// x0, x1, x2, origin(==x0) are (dist, ra, dec) coordinates if coordsys == equ
		//
		std::string fn, mfn;
		om_planecut_t &o = om_planecut;
		ss
			>> o.r.first >> o.r.second
			>> o.ri.first >> o.ri.second
			>> o.coordsys >> o.dx >> o.ndx
			>> o.d1 >> o.x1.first >> o.x1.second
			>> o.d2 >> o.x2.first >> o.x2.second
			>> o.d3 >> o.x3.first >> o.x3.second
			>> o.d0 >> o.x0.first >> o.x0.second
			>> o.earth_on_x_axis;
		ASSERT(!ss.fail());
		ss	>> o.mfn;

		// initialize plane_transformer
		setupPlaneTransformer(o.pt, o.coordsys, o.d1, o.x1, o.d2, o.x2, o.d3, o.x3, o.d0, o.x0, 1e10, o.earth_on_x_axis);

		if(mfn.size()) // store matrix and origin - for reading with SM
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
		m_color.color.set(color);
		m_color.range = make_pair(c1, c2);
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
	else if(cmd == "transform")
	{
		double node, inc;
		ss >> node >> inc; ASSERT(!ss.fail());
		t_node = rad(node); t_inc = rad(inc);
		filters |= T_COORD;
		out << "# GC transformation active\n";
		out << "# {node, inc}    = " << node << ", " << inc << "\n";
		out << "#\n";
	}
	else if(cmd == "paralax")
	{
		ss >> paralax.Mrc[0] >> paralax.Mrc[1] >> paralax.Mrc[2] >> paralax.Mrc[3] >> paralax.Mrc[4];
		ASSERT(!ss.fail());
		filters |= T_DISTANCE;
		out << "# Paralax relation changed to:\n";
		out << "#   Mr(i)  = " << join<double, 5>(", ", &paralax.Mrc[0]) << "\n";
		out << "#\n";
	}
	else
	{
		cerr << "Error - unknown command: " << cmd << "\n";
		exit(-1);
	}
	out.flush();
	return true;
}

selector *driver::newSelector(const std::string &filename)
{
	selector_d s;
	if(filename == "-")
	{
		s.s = new selector(cout, *this);
	} else {
		s.out = new ofstream(filename.c_str());
		ASSERT(!s.out->fail());
		s.s = new selector(*s.out, *this);
	}

	selectors.push_back(s);

	return s.s;
}

void driver::run()
{
	FOREACHj(list<selector_d>::iterator, j, selectors)
	{
		(*j).s->start();
	}

	for(att = 0; att != arr.size(); att++)
	{
		tick.tick();

		FOREACHj(list<selector_d>::iterator, j, selectors)
		{
			(*j).s->select(arr[att]);
		}
	}
	
	FOREACHj(list<selector_d>::iterator, j, selectors)
	{
		(*j).s->finish();
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
	VERSION_DATETIME(version);

	driver db;

	ifstream qf;
	if(argc > 1) { qf.open(argv[1]); }
	if(!qf.good()) { cerr << "Cannot open " << argv[1] << ". Aborting\n"; exit(-1); }
	
	istream &in = argc > 1 ? qf : cin;

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

	db.run();
}
catch(EAny &e)
{
	e.print();
}
}
