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

#include "conic_volume_map.h"

#include <astro/useall.h>
using namespace std;

conic_volume_map::conic_volume_map() : 
	proj_north(rad(90), rad(90)), proj_south(rad(-90), rad(-90)),
	bin_beam(NULL), construct_beam(NULL), finalize_beam(NULL)
{
}

BOSTREAM2(const conic_volume_map::pix &p)
{
	return out << p.d << p.dd << p.dOmega << p.N << p.data;
}

BISTREAM2(conic_volume_map::pix &p)
{
	return in >> p.d >> p.dd >> p.dOmega >> p.N >> p.data;
}

OSTREAM(const conic_volume_map &l)
{
	out << "# Binning radially in " << l.d_field << "\n";
	out << "# Radial binning limits and step (d0, d1, dd) = " << l.d0 << ", " << l.d1 << ", " << l.dd << "\n";
	out << "# Binning data field " << l.val_field << "\n";
	out << "# Loaded north sky from: " << l.northfn << "\n";
	out << "# Loaded south sky from: " << l.southfn;
	out << "# Angular resolution = " << deg(l.dx) << " deg [" << l.dx << " rad]\n";
	if(l.model.get() != NULL)
	{
		out << "# Loaded galactic model from: " << l.modelfn << "\n";
		out << "# Requested stars per pixel = " << l.starsperpix << "\n";
	}
	return out;
}

void conic_volume_map::initialize_for_galmodel()
{
	configure_binning("distance", "distance");

	// load the galactic model
	ifstream modelin(modelfn.c_str()); ASSERT(modelin.good());
	model.reset(galactic_model::load(modelin));
	ASSERT(model.get() != NULL);
	
	initialize();
}

void conic_volume_map::initialize()
{
	// load north and south volume maps
	{ ifstream f(northfn.c_str()); io::ibstream in(f); in >> skym_north; ASSERT(in); dx = skym_north.dx; }
	{ ifstream f(southfn.c_str()); io::ibstream in(f); in >> skym_south; ASSERT(in); ASSERT(skym_south.dx == dx); }

	// create pixels
	npix = 0;
	std::cerr << "Preparing conical pixels, north ...\n";
	create_pixels(proj_north, skym_north, sky_north);
	std::cerr << "Preparing conical pixels, south ...\n";
	create_pixels(proj_south, skym_south, sky_south);
}

io::obstream &conic_volume_map::serialize(io::obstream &out) const
{
	// north/south projections
	out << proj_north << proj_south;

	// sky
	out << dx;
	out << skym_north << skym_south;

	// gridding metadata
	out << d0 << d1 << dd;
	out << starsperpix;
	out << d_field;
	out << val_field;

	// maps
	out << sky_north << sky_south;
	out << npix;
	
	return out;
}

io::ibstream &conic_volume_map::unserialize(io::ibstream &in)
{
	// north/sinh projections
	in >> proj_north >> proj_south;

	// sky
	in >> dx;
	in >> skym_north >> skym_south;

	// gridding metadata
	in >> d0 >> d1 >> dd;
	in >> starsperpix;
	in >> d_field;
	in >> val_field;

	// maps
	in >> sky_north >> sky_south;
	in >> npix;
	
	return in;
}

ostream &conic_volume_map::text_serialize(ostream &out)
{
	out << "# Printing northern hemisphere\n";
	print_hemisphere(out, proj_north, skym_north, sky_north);

	out << "# Printing southern hemisphere\n";
	print_hemisphere(out, proj_south, skym_south, sky_south);

	return out;
}

void conic_volume_map::create_pixel(beam_t &beam, double dprev, double d, double area)
{
	double dpix = (d + dprev) / 2.;		// pixel center
	double ddpix = d - dprev;		// radial pixel length
	beam[dprev] = pix(dpix, ddpix, area);	// index by pixel _start_
}

#if 0
void conic_volume_map::construct_mag_equidistant_beam(beam_t &beam, Radians l, Radians b, double area)
{
	// create pixels equidistant in magnitude, from d0 == r0 to d1 == r1

}
#endif

void conic_volume_map::cb_equidistant(beam_t &beam, Radians l, Radians b, double area)
{
	// create pixels equidistant in magnitude, from d0 == r0 to d1 == r1, dd = dm
	double d = d0, dprev;
	while(d < d1)
	{
		dprev = d;
		d += dd;
		if(d > d1)
		{
			d -= dd; dd = d1 - d; d = d1;
		}

		create_pixel(beam, dprev, d, area);
	}
}

void conic_volume_map::cb_galmodel(beam_t &beam, Radians l, Radians b, double area)
{
	// construct beam - this basically integrates along the beam
	// and could be done very accurately using sample_integral() from simulate.cpp,
	// but such accuracy is not needed here.
	const double dOmega = sqr(dx);
	V3 unit(cos(l)*cos(b), sin(l)*cos(b), sin(b));

	double dstep = 1.;	// integration step
	double d = d0, dprev = d0;	// top edge of previous pixel
	double ninpix = 0, ntot = 0;
	while(d < d1)	// integrate from d0 to d1
	{
		d += dstep;
		if(d > d1) // adjust the final pixel to be aligned with the volume edge
		{
			d -= dstep; dstep = d1 - d; d = d1;
		}

		V3 v = unit*(d-0.5*dstep);
		double x = 8000 - v.x, y = -v.y, z = v.z;
		double rho = model->rho(x, y, z, 0.);
		double dn = rho * sqr(d-0.5*dstep)*dstep * dOmega;

		if(dn > 0.1*starsperpix)
		{
			d -= dstep;
			dstep *= .5;
			continue;
		}

		ninpix += dn; ntot += dn;
		if(ninpix >= starsperpix)
		{
			create_pixel(beam, dprev, d, area);

			ASSERT(fabs(ninpix/starsperpix - 1) < 0.3) { std::cerr << "ninpix = " << ninpix; }
			dprev = d;
			ninpix = 0;
		}

		if(dn < 0.02*starsperpix)
		{
			dstep *= 2.;
			continue;
		}
	}
	// final pixel
	if(ninpix != 0)
	{
		create_pixel(beam, dprev, d, area);
		ASSERT(ninpix < starsperpix) { std::cerr << "last ninpix = " << ninpix; }
	}
}

void conic_volume_map::create_pixels(lambert &proj, partitioned_skymap &skym, sky_t &sky)
{
	// create the pixels
	sky.clear();

	double tarea = 0;
	ticker tick(1);
	FOREACH(skym.skymap)
	{
		const std::pair<int, int> &p = (*i).first;
		double area = (*i).second.area;

		double lx = skym.x0 + (0.5 + p.first)*dx;
		double ly = skym.y0 + (0.5 + p.second)*dx;

		Radians l, b;
		proj.inverse(lx, ly, l, b);

		(this->*construct_beam)(sky[p], l, b, area);
		tick.tick();

		npix += sky[p].size();
		tarea += area;
	}
	tick.close();

	std::cerr << "Total area: " << tarea * sqr(deg(1.)) << "\n";
}

bool conic_volume_map::bin(Radians l, Radians b, const mobject &m)
{
	// convert (l,b) to map coordinates, and map to proper hemisphere
	double x, y;
	proj_north.convert(l, b, x, y);
	double r2 = sqr(x) + sqr(y);
	if(r2 < 2) { return bin(x, y, m, skym_north, sky_north); }

	proj_south.convert(l, b, x, y);
	r2 = sqr(x) + sqr(y);
	if(r2 < 2) { return bin(x, y, m, skym_south, sky_south); }
	
	ASSERT(0);
}

#define ITER(x) typeof((x).begin())
bool conic_volume_map::bb_bin_sum(beam_t &beam, Radians l, Radians b, const mobject &m)
{
	double D = m.field(d_field);

	// find the apropriate D bin in the beam
	ITER(beam) i = beam.upper_bound(D);
	if(i == beam.begin()) { return false; }
	--i;
	pix &p = (*i).second;

	p.N++;
	double val = m.field(val_field);
	if(p.data.size()) { p.data[0] += val; } else { p.data.push_back(val); }
	//std::cerr << p.data << " " << m.field(val_field) << "\n";
	ASSERT(abs((p.d - D)/p.dd) <= 0.5) {
		std::cerr << l << " " << b << " " << D << " !-> " << p;
	}
	return true;
}

void conic_volume_map::fin_average(beam_t &beam, Radians l, Radians b, double area)
{
	FOREACH(beam)
	{
		pix &p = (*i).second;
		if(p.N == 0) { continue; }

		p.data[0] /= p.N;
	}
}

bool conic_volume_map::bb_bin_collect(beam_t &beam, Radians l, Radians b, const mobject &m)
{
	double D = m.field(d_field);

	// find the apropriate D bin in the beam
	ITER(beam) i = beam.upper_bound(D);
	if(i == beam.begin()) { return false; }
	--i;
	pix &p = (*i).second;

	p.N++;
	double val = m.field(val_field);
	p.data.push_back(val);
	//std::cerr << p.data << " " << m.field(val_field) << "\n";
	ASSERT(abs((p.d - D)/p.dd) <= 0.5) {
		std::cerr << l << " " << b << " " << D << " !-> " << p;
	}
	return true;
}

void conic_volume_map::fin_median(beam_t &beam, Radians l, Radians b, double area)
{
	FOREACH(beam)
	{
		pix &p = (*i).second;
		if(p.N == 0) { continue; }

		// calculate median
		sort(p.data.begin(), p.data.end());
		int size = p.data.size(), mid = size/2;
		double median;
		if(size % 2)
		{
			median = p.data[mid];
		} else {
			median = (p.data[mid - 1] + p.data[mid]) / 2.;
		}

		// store result
		p.data = std::vector<double>(1, median);
	}
}

bool conic_volume_map::get_sky_and_index(int &skyidx, int &X, int &Y, const Radians l, const Radians b)
{
	partitioned_skymap *msky; sky_t *sky;

	// convert (l,b) to map coordinates, and map to proper hemisphere
	double x, y;
	proj_north.convert(l, b, x, y);
	double r2 = sqr(x) + sqr(y);
	if(r2 < 2)
	{
		sky = &sky_north; msky = &skym_north; skyidx = 0;
	}
	else
	{
		proj_south.convert(l, b, x, y);
		r2 = sqr(x) + sqr(y);
		if(r2 < 2)
		{
			sky = &sky_south; msky = &skym_south; skyidx = 1;
		}
		else
		{
			ASSERT(0);
		}
	}

	X = (int)((x - msky->x0) / dx);
	Y = (int)((y - msky->y0) / dx);

	return sky->count(make_pair(X, Y));
}

bool conic_volume_map::get_beam(beam_t *beam, Radians l, Radians b)
{
	int skyidx, X, Y;
	bool exists = get_sky_and_index(skyidx, X, Y, l, b);
	if(!exists) { return false; }

	sky_t &sky = skyidx == 0 ? sky_north : sky_south;
	beam = &sky[make_pair(X, Y)];
	return true;
}

bool conic_volume_map::bin(double x, double y, const mobject &m, partitioned_skymap &msky, sky_t &sky)
{
	double D = m.field(d_field);
	if(D <= d0 || D >= d1) { return false; }

	int X = (int)((x - msky.x0) / dx);
	int Y = (int)((y - msky.y0) / dx);
	bool exists = sky.count(make_pair(X, Y));
	if(!exists) { std::cerr << "XY bin not there: " << X << " " << Y << " " << D << "\n"; return false; }
	beam_t &beam = sky[make_pair(X, Y)];

	return (this->*bin_beam)(beam, x, y, m);
}

static const double radsqr2deg = sqr(deg(1.));
OSTREAM(const conic_volume_map::pix &p)
{
	return out << p.d << " " << p.N << " " << p.dOmega*radsqr2deg << " " << p.dd;
}

void conic_volume_map::print_hemisphere(ostream &out, lambert &proj, partitioned_skymap &skym, sky_t &sky)
{
	FOREACH(sky)
	{
		conic_volume_map::beam_t &beam = (*i).second;
		
		int X = (*i).first.first;
		int Y = (*i).first.second;
		const double x = skym.x0 + (0.5 + X)*dx;
		const double y = skym.y0 + (0.5 + Y)*dx;
		Radians l, b;
		proj.inverse(x, y, l, b);
		l = std::fmod(l + 2.*ctn::pi, 2.*ctn::pi);

		V3 unit(cos(l)*cos(b), sin(l)*cos(b), sin(b));
		double area = skym.skymap[(*i).first].area;

		// do any final processing
		(this->*finalize_beam)(beam, l, b, area);

		FOREACH(beam)
		{
			if((*i).second.N == 0) { continue; }
#if 0
			out << deg(l) << " " << deg(b) << " " << (*i).second << " " << X << " " << Y << " " << area*sqr(deg(1.)) << "\n";
#else
			pix &p = (*i).second;
			double d = p.d, N = p.N;
			double V = (pow(d - 0.5*p.dd, 3.) - pow(d + 0.5*p.dd, 3.)) / 3. * p.dOmega;	// exact pixel volume

			V3 v = unit*d;
			v.x = Rg - v.x;
			v.y = -v.y; // convert to galactocentric galactic

			double r = v.rho();
			const double phi0 = rad(180.);
			Radians phi = modulo(v.phi() - phi0, 2*ctn::pi);
			double rphi = phi*r;

			out << r << " " << rphi << " " << v.z << " 0 0 " << N << " " << V << "    "
				<< p.data[0] << " "
				<< setprecision(10) << deg(l) << " " << deg(b) << " " << (*i).second << " " 
				<< x << " " << y << " " << area*sqr(deg(1.)) << setprecision(6) << "\n";
#endif
		}
	}
}

////////////////////////////////////////
///  conic_volume_map_interpolator
////////////////////////////////////////

conic_volume_map_interpolator::conic_volume_map_interpolator()
: vm(NULL), use_median_beam_for_all(false)
{
}

void conic_volume_map_interpolator::initialize(conic_volume_map &vm)
{
	this->vm = &vm;
	set_aux(north, vm.sky_north);
	set_aux(south, vm.sky_south);
	
	// construct median spline which will be used in areas
	// where there is insufficient data
	std::vector<double> x, y, tmp;

	median_beam.dmin = 15; median_beam.dmax = 22;
	double dx = (median_beam.dmax - median_beam.dmin) / 100.;
	for(double d = median_beam.dmin; d <= median_beam.dmax; d += dx)
	{
		beams_t beams[2] = {north, south};

		tmp.clear();
		FORj(s, 0, 2) FOREACH(beams[s])
		{
			beam_t &beam = (*i).second;
			tmp.push_back(beam.spl(d));
		}

		sort(tmp.begin(), tmp.end());
		int size = tmp.size(), mid = size/2;
		double median = size % 2 ? tmp[mid] : (tmp[mid - 1] + tmp[mid]) / 2.;

		x.push_back(d);
		y.push_back(median);
	}
	median_beam.spl.construct(&x[0], &y[0], x.size());

	ofstream ff("conic_interp_dump.txt");
	text_dump_median_beam(ff);
}

void conic_volume_map_interpolator::text_dump_median_beam(std::ostream &out)
{
	beam_t &beam = median_beam;
	double dx = (beam.dmax - beam.dmin) / 100.;
	for(double d = beam.dmin; d <= beam.dmax; d += dx)
	{
		out << "0 " << d << " " << beam.spl(d) << "\n";
	}
}

void conic_volume_map_interpolator::text_dump(std::ostream &out)
{
	int k = 0;
	beams_t beams[2] = {north, south};
	FORj(s, 0, 2)
	{
		FOREACH(beams[s])
		{
			beam_t &beam = (*i).second;
			beam.dmin = 15; beam.dmax = 22;
			double dx = (beam.dmax - beam.dmin) / 100.;
			for(double d = beam.dmin; d <= beam.dmax; d += dx)
			{
				out << k << " " << d << " " << beam.spl(d) << "\n";
			}
			k++;
		}
		std::cerr << "Dumped " << s << "\n";
	}
}

void conic_volume_map_interpolator::set_aux(beams_t &beams, conic_volume_map::sky_t &sky)
{
	beams.clear();
	FOREACH(sky)
	{
		conic_volume_map::beam_t &beam = (*i).second;

		std::vector<double> x, y;
		double dmin, dmax, dOmega;
		
		// anchor at (0,0)
		x.push_back(0.);
		y.push_back(0.);

		FOREACHj(j, beam)
		{
			conic_volume_map::pix &p = (*j).second;
			//if(p.N == 0) { continue; }
			if(p.N < 5) { continue; }

			x.push_back(p.d);
			y.push_back(p.data[0]);
			
			if(j == beam.begin())
			{
				dOmega = p.dOmega;
				dmin = p.d - 0.5 * p.dd;
			}
			dmax = p.d + 0.5 * p.dd;
		}

		if(x.size() < 3) { continue; }

		beam_t &sbeam = beams[(*i).first];
		sbeam.dmin = dmin;
		sbeam.dmax = dmax;
		sbeam.dOmega = dOmega;
		sbeam.spl.construct(&x[0], &y[0], x.size());
	}
}

double conic_volume_map_interpolator::interpolate(Radians l, Radians b, double D)
{
	if(!use_median_beam_for_all)
	{
		// find this beam
		int skyidx, X, Y;
	
		bool exists = vm->get_sky_and_index(skyidx, X, Y, l, b);
		beams_t &sky = skyidx == 0 ? north : south;
		if(exists && sky.count(make_pair(X, Y)))
		{
			return sky[make_pair(X, Y)].spl(D);
		}
	
		// find adjacent beams and calculate the mean
		double val = 0; double wt = 0;
		FORj(dX, -1, 2)
		{
			int Xn = X + dX;
			FORj(dY, -1, 2)
			{
				int Yn = Y + dY;
				if(!sky.count(make_pair(Xn, Yn))) { continue; }
				beam_t &beam = sky[make_pair(Xn, Yn)];
				val += beam.spl(D)*beam.dOmega;
				wt += beam.dOmega;
			}
		}
		if(wt != 0)  { return val / wt; }
	}
	
	// calculate from median
	return median_beam.spl(D);
}
