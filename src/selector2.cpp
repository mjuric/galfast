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

#include "dm.h"
#include "projections.h"

#include <astro/useall.h>
using namespace std;

class RunMap
{
protected:
	std::map<int, int> rm;
public:
	RunMap()
	{
		text_input_or_die(in, "schlegel.cat.txt");
		
		int run, end_id;
		bind(in, run, 0, end_id, 7);
		
		while(in.next())
		{
			rm[end_id] = run;
		}
	}

	int operator[](int fitsId)
	{
		ASSERT(fitsId >= 0);
		std::map<int, int>::iterator k = rm.upper_bound(fitsId);
		ASSERT(k != rm.end());
		return (*k).second;
	}
};

class filter
{
	bool operator()(mobject &o) = 0;

	filter *clone() = 0;

	istream read(istream &ss) = 0;
	ostream write(ostream &ss) const = 0;
};
ISTREAM(filter &p) { return p.read(in); }
OSTREAM(const filter &p) { return p.write(in); }

class filter_factory
{
	bool register(const std::string &cmd, filter &f);
	filter *query(const std::string &cmd);
};

struct property_base
{
	istream read(istream &ss) = 0;
};

ISTREAM(property_base &p) { return p.read(in); }

template<typename T>
struct property : public property_base
{
	T value;

	istream read(istream &ss) { return ss >> value; }
	operator (typename T)() { return value; }
	T operator =(const &T v) { value = v; return v; }

	property(filter *parent, const std::string &name);
};

struct beam_filter : public filter
{
	property<Radians> ra, dec, radius;

	bool operator()(mobject &o)
	{
		return coordinates::distance(rad(m.ra), rad(m.dec), ra, dec) < radius;
	}

	filter *clone() { return new beam_filter; }
	
	beam_filter()
		: ra(this, "ra"), dec(this, "dec"), radius(this, "radius")
		{}

	istream read(istream &ss) { return ss >> value; }
};

ss >> get_filter("beam").property("ra")

class selector
{
	map<string, filter *> filters;

	filter &get_filter(std::string &cmd)
	{
		if(!filters.count())
		{
			map[cmd] = filterFactory.create(cmd);
		}
		return *map[cmd];
	}

	while(!eof())
	{
		string cmd;
		ss >> cmd;

		if(cmd == "set")
		{
			string f, p;
			cmd >> f >> p;
			ss >> get_filter(f).property(p);
		} else {
			ss >> get_filter(f);
		}
	}
	
	FOREACH(filters)
	{
		
	}
}

class selector
{
protected:
	mobject *mm;
	DMMArray<mobject> arr;
	DMMArray<starmag> starmags;
	RunMap rm;
	int att, n;
public:
	#define M_BEAM		0x01
	#define M_RECT		0x02
	#define M_RI		0x04
	#define M_R		0x08
	#define M_ML_RI		0x10
	#define M_ML_R		0x20
	#define M_GAL_RECT	0x40
	#define M_SPARSE	0x80
	#define M_RUN		0x100
	#define M_GR		0x200
	int filters;

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

	ticker tick;
protected:
	void reset();

public:
	selector(mobject *m = NULL) : mm(NULL), att(0), n(0), tick(10000), filters(0)
	{
		if(m != NULL) { init(*m, "dm_unique_stars.dmm"); }
	}

	void init(mobject &m, const std::string &dmmfile);

	int select_beam(Radians ra, Radians dec, Radians radius);
	int select_rect(Radians ra0, Radians dec0, Radians ra1, Radians dec1);
	int select_gal_rect(Radians l0_, Radians b0_, Radians l1_, Radians b1_);

	int record_count() { return n; }	
	int at() { return att; }
	bool next();
};

void selector::reset()
{
	att = 0;
	n = 0;
}

void selector::init(mobject &m, const std::string &dmmfile)
{
	arr.open(dmmfile, "r");
	starmags.open("dm_starmags.dmm", "r");
	mm = &m;
}

int selector::select_beam(Radians ra_, Radians dec_, Radians radius_)
{
	ra = ra_; dec = dec_; radius = radius_;
	filters |= M_BEAM;
	reset();
}

int selector::select_rect(Radians ra0_, Radians dec0_, Radians ra1_, Radians dec1_)
{
	ra0 = ra0_; dec0 = dec0_; ra1 = ra1_; dec1 = dec1_;
	filters |= M_RECT;
	reset();
}

int selector::select_gal_rect(Radians l0_, Radians b0_, Radians l1_, Radians b1_)
{
	l0 = l0_; b0 = b0_; l1 = l1_; b1 = b1_;
	filters |= M_GAL_RECT;
	reset();
}

//bool thruRun;

#define FILTER(f) if(selected && (filters & f))
bool selector::next()
{
	ASSERT(mm);

//	if(att == 0) { thruRun = false; }

	for(;;)
	{
		tick.tick();

		++att;
		if(att == arr.size()) return false;

		mobject &m = arr[att];
		bool selected = true;

		FILTER(M_RUN)
		{
			selected = false;
			FOR(m.obs_offset, m.obs_offset + m.n)
			{
				starmag &sm = starmags[i];
//				cerr << "\t" << sm.fitsId << " " << rm[sm.fitsId] << " " << run << "\n";
				if(rm[sm.fitsId] == run) { selected = true; break; }
			}
//			if(thruRun) { att == arr.size() - 1; continue; }
		}
		FILTER(M_SPARSE)
		{
			selected = (att % nskip) == 0;
		}
		FILTER(M_BEAM)
		{
			selected = coordinates::distance(rad(m.ra), rad(m.dec), ra, dec) < radius;
		}
		FILTER(M_RECT)
		{
			selected = coordinates::inBox(rad(m.ra), rad(m.dec), ra0, dec0, ra1, dec1);
		}
		FILTER(M_GAL_RECT)
		{
			double l, b;
			coordinates::equgal(rad(m.ra), rad(m.dec), l, b);
			selected = coordinates::inBox(l, b, l0, b0, l1, b1);
		}
		FILTER(M_RI)
		{
			float ri = m.ri();
			selected = m_ri.first <= ri && ri < m_ri.second;
		}
		FILTER(M_GR)
		{
			float gr = m.gr();
			selected = m_gr.first <= gr && gr < m_gr.second;
		}
		FILTER(M_R)
		{
			float r = m.mag[1];
			selected = m_r.first <= r && r < m_r.second;
		}
		FILTER(M_ML_RI)
		{
			float ri = m.ml_mag[1] - m.ml_mag[2];
			selected = m_ml_ri.first <= ri && ri < m_ml_ri.second;
		}
		FILTER(M_ML_R)
		{
			float r = m.ml_mag[1];
			selected = m_ml_r.first <= r && r < m_ml_r.second;
		}

		if(selected) break;
	}

	++n;
	*mm = arr[att];
	return true;
}

#define OM_STARS	1
#define OM_OBSV		2
#define OM_MAP		3
#define OM_CMD		4

int outputMode = OM_STARS;

struct om_lambert_t
{
	double dx;
	double x0, y0, x1, y1;
	double l0, phi1;
	
	double w() { return x1 - x0; }
	double h() { return y1 - y0; }
}
om_lambert = 
{
	.0333333,
	-2, -2, 2, 2,
	0, 90
};

OSTREAM(const om_lambert_t &l)
{
	out << "# dx = " << l.dx << "\n";
	out << "# (lambda0, phi1) = (" << l.l0 << ", " << l.phi1 << ")\n";
	out << "# x = [" << l.x0 << ", " << l.x1 << ")\n";
	out << "# y = [" << l.y0 << ", " << l.y1 << ")";
	return out;
}

struct om_cmd_t
{
public:
	int c1, c2;
	int m;

	double dx, dy;
	double x0, x1, y0, y1;

	#define MAG_G	0
	#define MAG_R	1
	#define MAG_I	2
public:
	int bandIdx(const char c)
	{
		switch(c)
		{
			case 'g' : return 0;
			case 'r' : return 1;
			case 'i' : return 2;
			default:
				cerr << "Unknown band (" << c << ").\n";
				ASSERT(0);
		}
	}

	float color(const mobject &m) const { return m.mag[c1] - m.mag[c2]; }
	float magnitude(const mobject &m) const { return m.mag[this->m]; }
} om_cmd =
{
	MAG_G, MAG_R,
	MAG_R,
	0.033333, 0.1,
	-0.2, 2, 14, 22
};
char band_names[] = "gri";

OSTREAM(const om_cmd_t &l)
{
	out << "# " << band_names[l.c1] << "-" << band_names[l.c2] << " vs. " << band_names[l.m] << "\n";
	out << "# dx = " << l.dx << ", dy = " << l.dy << "\n";
	out << "# x = [" << l.x0 << ", " << l.x1 << ")\n";
	out << "# y = [" << l.y0 << ", " << l.y1 << ")";
	return out;
}

void printout(mobject &m, selector &beam, DMMArray<starmag> &starmags)
{
	lambert lambert_map(rad(om_lambert.l0), rad(om_lambert.phi1));
	map<S2, int, less_S2> binned;

	while(beam.next())
	{
		// SM doesn't like infinities
		FOR(0, 3) { if(!isfinite(m.magErr[i])) { m.magErr[i] = 0; } }

		double l, b;
		coordinates::equgal(rad(m.ra), rad(m.dec), l, b);

		switch(outputMode)
		{
		case OM_STARS:
			cout << setw(12) << beam.at() << " " << m;// << "\n";
			cout << setw(12) << setprecision(8) << deg(l) << " ";
			cout << setw(12) << setprecision(8) << deg(b) << "\n";
			break;
		case OM_OBSV:
			FOR(m.obs_offset, m.obs_offset + m.n)
			{
				starmag &sm = starmags[i];
				char buf[30]; sprintf(buf, "%6.3f", m.Ar);
				cout 	<< setw(12) << beam.at()
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
		case OM_CMD: {
			double x, y;
			x = om_cmd.color(m);
			y = om_cmd.magnitude(m);
			if(!between(x, om_cmd.x0, om_cmd.x1)) break;
			if(!between(y, om_cmd.y0, om_cmd.y1)) break;

			S2 k(
				(int)(floor(x / om_cmd.dx)),
				(int)(floor(y / om_cmd.dy))
			);

			if(binned.count(k) == 0) { binned[k] = 1; }
			else { binned[k]++; }
			
			} break;
		}
	}

	if(outputMode == OM_MAP)
	{
		FOREACH(binned)
		{
			int v = (*i).second;
			const S2 &k = (*i).first;

			cout << (0.5 + k.x)*om_lambert.dx << " " << (0.5 + k.y)*om_lambert.dx << " " << v << "\n";
		}
	}
	else if(outputMode == OM_CMD)
	{
		FOREACH(binned)
		{
			int v = (*i).second;
			const S2 &k = (*i).first;

			cout << (0.5 + k.x)*om_cmd.dx << " " << (0.5 + k.y)*om_cmd.dy << " " << v << "\n";
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

template<typename T>
bool set_param(istream &ss, const string &cmd, const string &param_name, T& param)
{
	if(cmd != param_name) return false;

	ss >> param; ASSERT(!ss.fail());
	cout << "# " << param_name << " = " << param << "\n";
	cout << "#\n";

	return true;
}

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

	mobject m;
	selector beam(&m);


	ifstream qf;
	if(argc > 1) { qf.open(argv[1]); }
	if(!qf.good()) { cerr << "Cannot open " << argv[1] << ". Aborting\n"; exit(-1); }
	
	istream &in = argc > 1 ? qf : cin;

	string cmd;
	char line[1000];
	while(!in.eof())
	{
		in.getline(line, sizeof(line));
		if(in.eof()) break;

		stringstream ss(line);
		ss >> cmd;
		if(cmd[0] == '#' || ss.fail()) continue;

		if(cmd == "beam")
		{
			double ra, dec, radius;
			ss >> ra >> dec >> radius; ASSERT(!ss.fail());
			beam.select_beam(rad(ra), rad(dec), rad(radius));
			cout << "# beam filter active\n";
			cout << "# ra     = " << ra << "\n";
			cout << "# dec    = " << dec << "\n";
			cout << "# radius = " << radius << "\n";
			cout << "#\n";
		}
		else if(cmd == "rect")
		{
			double ra0, dec0, ra1, dec1;
			ss >> ra0 >> dec0 >> ra1 >> dec1; ASSERT(!ss.fail());
			beam.select_rect(rad(ra0), rad(dec0), rad(ra1), rad(dec1));
			cout << "# rect filter active\n";
			cout << "# (ra0, dec0)    = " << ra0 << " " << dec0 << "\n";
			cout << "# (ra1, dec1)    = " << ra1 << " " << dec1 << "\n";
			cout << "#\n";
		}
		else if(cmd == "gal_rect")
		{
			double l0, b0, l1, b1;
			ss >> l0 >> b0 >> l1 >> b1; ASSERT(!ss.fail());
			beam.select_gal_rect(rad(l0), rad(b0), rad(l1), rad(b1));
			cout << "# gal_rect filter active\n";
			cout << "# (l0, b0)    = " << l0 << " " << b0 << "\n";
			cout << "# (l1, b1)    = " << l1 << " " << b1 << "\n";
			cout << "#\n";
		}
		else if(cmd == "gal_beam")
		{
			double l, b, radius;
			ss >> l >> b >> radius; ASSERT(!ss.fail());
			Radians ra, dec;
			coordinates::galequ(rad(l), rad(b), ra, dec);
			beam.select_beam(ra, dec, rad(radius));
			cout << "# gal_beam filter active\n";
			cout << "# l     = " << l << "\n";
			cout << "# b    = " << b << "\n";
			cout << "# radius = " << radius << "\n";
			cout << "#\n";
		}
		else if(cmd == "obsv")
		{
			outputMode = OM_OBSV;
			cout << "# outputing observations\n";
			cout << "#\n";
		}
		else if(cmd == "lambert")
		{
			outputMode = OM_MAP;
			cout << "# outputing binned lambert map\n";
			cout << om_lambert << "\n";
			cout << "#\n";
		}
		else if(cmd == "cmd")
		{
			outputMode = OM_CMD;
			cout << "# outputing color-magnitude diagram\n";
			cout << om_cmd << "\n";
			cout << "#\n";
		}
		else if(set_param(ss, cmd, "cmd.dx", om_cmd.dx)) {}
		else if(set_param(ss, cmd, "cmd.dy", om_cmd.dy)) {}
		else if(set_param(ss, cmd, "cmd.x0", om_cmd.x0)) {}
		else if(set_param(ss, cmd, "cmd.x1", om_cmd.x1)) {}
		else if(set_param(ss, cmd, "cmd.y0", om_cmd.y0)) {}
		else if(set_param(ss, cmd, "cmd.y1", om_cmd.y1)) {}
		else if(cmd == "sparse")
		{
			ss >> beam.nskip; ASSERT(!ss.fail());
			beam.filters |= M_SPARSE;
			cout << "# outputing every " << beam.nskip << " star\n";
			cout << "#\n";
		}
		else if(cmd == "run")
		{
			ss >> beam.run; ASSERT(!ss.fail());
			beam.filters |= M_RUN;
			cout << "# run filter active\n";
			cout << "# run   = " << beam.run << "\n";
			cout << "#\n";
		}
		else if(cmd == "ri")
		{
			double ri0, ri1;
			ss >> ri0 >> ri1; ASSERT(!ss.fail());
			beam.m_ri = make_pair(ri0, ri1);
			beam.filters |= M_RI;
			cout << "# ri filter active\n";
			cout << "# [ri0, ri1)    = " << ri0 << " " << ri1 << "\n";
			cout << "#\n";
		}
		else if(cmd == "gr")
		{
			double gr0, gr1;
			ss >> gr0 >> gr1; ASSERT(!ss.fail());
			beam.m_gr = make_pair(gr0, gr1);
			beam.filters |= M_GR;
			cout << "# gr filter active\n";
			cout << "# [gr0, gr1)    = " << gr0 << " " << gr1 << "\n";
			cout << "#\n";
		}
		else if(cmd == "r")
		{
			double r0, r1;
			ss >> r0 >> r1; ASSERT(!ss.fail());
			beam.m_r = make_pair(r0, r1);
			beam.filters |= M_R;
			cout << "# r filter active\n";
			cout << "# [r0, r1)    = " << r0 << " " << r1 << "\n";
			cout << "#\n";
		}
		else if(cmd == "ml_ri")
		{
			double ri0, ri1;
			ss >> ri0 >> ri1; ASSERT(!ss.fail());
			beam.m_ml_ri = make_pair(ri0, ri1);
			beam.filters |= M_ML_RI;
			cout << "# ml_ri filter active\n";
			cout << "# [ri0, ri1)    = " << ri0 << " " << ri1 << "\n";
			cout << "#\n";
		}
		else if(cmd == "ml_r")
		{
			double r0, r1;
			ss >> r0 >> r1; ASSERT(!ss.fail());
			beam.m_ml_r = make_pair(r0, r1);
			beam.filters |= M_ML_R;
			cout << "# ml_r filter active\n";
			cout << "# [r0, r1)    = " << r0 << " " << r1 << "\n";
			cout << "#\n";
		}
		else
		{
			cerr << "Error - unknown command: " << cmd << "\n";
			exit(-1);
		}
		cout.flush();
	}

	DMMArray<starmag> starmags("dm_starmags.dmm");
	printout(m, beam, starmags);
}
catch(EAny &e)
{
	e.print();
}
}
