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

#include "raytrace.h"
#include "interval_arithmetic.h"
#include "dm.h"
#include "projections.h"

#include <astro/math/vector.h>
#include <astro/coordinates.h>
#include <astro/math.h>
#include <astro/util.h>
#include <astro/sdss/rungeometry.h>
#include <astro/io/fits.h>
#include <astro/system/options.h>

#include <astro/useall.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>

#include "ximage.h"
#include "floodfill.h"
#include "analysis.h"
#include "binarystream.h"

using namespace std;

#if 0
#if 0
struct mobject
{
	int obs_offset;			// offset in starmags array to where the obsv. of this object begin
	int n;				// total number of observations available

	float Ar;

	float mag[5];			// corrected averaged magnitudes
	float magErr[5];		// note: magnitude ordering is grizu
	float N[5];			// number of observations of this object, which were used to calculate mag[]
	short flags;			// how exactly we calculated mag[] and magErr[]

	float ml_mag[3];		// magnitudes deduced by max.likelihood fitting to locus
	float D;			// geocentric distance

	double ra, dec;
#endif

struct galaxy_model
{
	double den0, den1, l0, l1, h0, h1, z0;	// thin and thick disk parameters
	double denhalo, r0, alpha;	// halo parameters
	float ri0, ri1;

	double dx;			// scale for this model

	double density(const V3 &p);
};

double galaxy_model::density(const V3 &p)
{
	double rho = p.rho(),
	       r = abs(p),
	       z = p.z;

	double den = den0*exp(rho/l0 + (z-z0)/h0) +
	             den1*exp(rho/l1 + (z-z0)/h1) +
	             denhalo*pow(r/r0, alpha);
	return den;
}

int galaxy_model::number(double den, double V)
{
	double n = den*V;
	int N = (int)n;

	// add Poisson noise
	// TODO: draw a number from N(mu=n, sigma=sqrt(n)) and add it to n

	// the fraction is the probability of having a particle there
	n -= N;
	float prob = float(rand())/RAND_MAX;
	if(prob < n)
	{
		N++;
	}

	return N;
}

int galaxy_model::star()
{
	// generate stellar colors
	float ri = ri0 + (ri1-ri0)*rnd();
	// find the absolute magnitude:
	float Mr = paralax.mr(ri);
	// assuming we're on the main sequence:
	float gr = paralax.gr(ri);
	// R magnitude from distance
	float r = stardist::m(D, s.Mr);
}

struct real_star
{
	float ri, gr, r;
};

class generate_model : public bin3d_pixelize_volume
{
public:
	map<float, galaxy_model> ri_models;	// disk and halo models for stars of a given ri color (spectral type)
public:
	cylindrical_pixelize_volume(const std::string &volfn, float dx, int n, pair<float, float> r, pair<float, float> ri, Radians phi0_ = 0)
		: bin3d_pixelize_volume(volfn, dx, n, r, ri), phi0(phi0_)
		{}

	double 	rho,	// galactocentric radius (set in action())
		D;	// earthcentric distance (set in action())
	mobject m;
	sdss_star ss;

	bool observe(real_star s)
	{
		// find r magnitude
		m.mag[0] = s.r;
		m.mag[1] = s.gr + s.r;
		m.mag[2] = s.r - s.ri;

		// add photometric errors, based on some braindead model of photometric errors
		// first approximation - all are drawn from some gaussian distribution, irrespective of magnitude
		m.magErr[0] = m.magErr[1] = m.magErr[2] = .02;
		FOR(0, 3)
		{
			m.mag[i] += gsl_ran_gaussian(rng, m.magErr[i]);
			ss.mag[i+1] = m.mag[i]; ss.magErr[i+1] = m.magErr[i]; // for ML color calculation
		}

		if(paralax(s))	// calculate the absolute magnitude and distances
		{
			m.D = s.earth.D;
			m.ml_mag[0] = s.ml_g;
			m.ml_mag[1] = s.ml_r;
			m.ml_mag[2] = s.ml_i;
		} else {
			m.D = 0;
		}

		// calculate position
		ecenequ(V3(x, y, z), D, m.ra, m.dec);
		DEG(m.ra); DEG(m.dec);

		return true;
	}

	virtual double action(float ddv)
	{
		// generate stars for each model, in this pixel
		rho = sqrt(sqr(x+Rg)+sqr(y));	// galactocentric radius
		D = sqrt(sqr(x)+sqr(y)+sqr(z));	// earthcentric distance
		m.ra = -1;
		FOREACH(ri_models)
		{
			const galaxy_model &model = (*i).second;
			double den = model.density(rho, z);
			int N = model.number(den, ddv);
			while(N)
			{
				real_star &s = model.star();
				if(observe(m, s))
				{
					// << write the star to DMM
				}
			}
		}
	}
};

cat_generator::generate(std::string &fn)
{
	binned_run br;
	binary_input_or_die(in, fn);
	in >> br;
	brs.dx = br.dx;
	FOREACH(br.pixels)
	{
		V3 v = (*i).first;
		v *= b
		binned_run::pixel &p = (*i).second;
		p.volume;
		
		double den = density(p);
		int N = number(den, V);
		if(N == 0) { return 0; }

		brs.pixels[k].uniqueVolume = p.volume;
	}
}

class cat_generator
{
public:
	double dx;	// volume pixel dx, used when generating star coordinates

protected:
	int N;
	double V;
	V3 p;

	int number(double den);
public:
	int init(const V3 &p, double v);
//	bool next(mobject &m);
};

int den_generator::init(const V3 &p_, double v)
{
	srand(time(NULL));

	V = v;
	p = p_;

	double den = density(p);
	N = number(den);
	
	return N;
}

#if 0
bool den_generator::next(mobject &m)
{
	if(N == 0) return false;
	--N;
	
	// assign a random coordinate within the pixel
	V3 pos;
	rand(pos, -dx/2, +dx/2);
	pos += p;

	// convert to equatorial coordinates
	ecenequ(pos, m.D, m.ra, m.dec);
	
	return true;
}
#endif

/////////////////////////////
#endif

struct BitImage
{
	const int nbits;

	int w, h;
	int x0, y0;
	double dx;

	valarray<int> data;
	
	BitImage(double dx_, int w_, int h_)
		: nbits(sizeof(int)*8), dx(dx_), w(w_), h(h_), x0(w/2), y0(h/2), data(((long long)w)*h/nbits + 1)
	{
		cerr << "w*h: " << ((long long)w)*h << "\n";
		cerr << "w*h/nbits + 1: " << ((long long)w)*h/nbits + 1 << "\n";
		cerr << "nbits: " << nbits << "\n";
		cerr << "Size: " << data.size() << "\n";
	}

	void index(int x, int y, int &byte, int &bit) const
	{
		long long ind = ((long long)(y+y0))*h + (x+x0);
		byte = ind / nbits;
		bit = ind % nbits;

		if(byte < 0 || byte >= data.size())
		{
			cerr << "Error: x = " << x << ", y = " << y << ", byte = " << byte << ", data.size() = " << data.size() << "\n";
			ASSERT(!(ind < 0 || ind >= data.size()));
		}
	}

	bool get(int x, int y) const
	{
		int byte, bit;
		index(x, y, byte, bit);

		int v = data[byte];
		return v & (1 << bit);
	}
	
	void set(int x, int y)
	{
		int byte, bit;
		index(x, y, byte, bit);

		int &v = data[byte];
		v |= (1 << bit);
	}
};

void makeSkyMap()
{
	Radians dx = rad(.125 * 1/60.); /* sampling resolution in radians */
	cerr << "dx = " << dx << "\n";
	int l = int(2*sqrt(2.)/dx + 1) + 2;
	cerr << "l = " << l << "\n";
	BitImage north(dx, l, l),
		south(dx, 1, 1);
	lambert lnorth(rad(90), rad(90)), lsouth(rad(90), rad(-90));

	RunGeometryDB db;
	double x, y;
	FOREACH(RunGeometryDB::db_t::iterator, db.db)
	{
		const RunGeometry &geom = (*i).second;
		Mask mask(geom);
		Radians muEnd = geom.muStart + mask.length() + dx;
		cerr << geom.run << " ";
		FORj(col, 0, 6)
		{
			for(Radians nu = mask.lo(col) + dx/2; nu < mask.hi(col) + dx; nu += dx)
			{
				for(Radians mu = geom.muStart + dx/2; mu < muEnd; mu += dx)
				{
					Radians ra, dec;
					coordinates::gcsequ(geom.node, geom.inc, mu, nu, ra, dec);

					dec *= -1;
					if(dec <= 0) continue;

					BitImage &img = dec > 0 ? north : south;
					lambert &proj = dec > 0 ? lnorth : lsouth;

					proj.convert(ra, dec, x, y);
//					proj.convert(0, 0, x, y);
					int I = (int)floor(x/dx + 0.5);
					int J = (int)floor(y/dx + 0.5);
					if(abs(I) > img.x0 || abs(J) > img.y0)
					{
						cerr << "Problem: " << I << ", " << J;
						cerr << " " << deg(ra) << ", " << deg(dec);
						cerr << " " << x << ", " << y << "\n";
						exit(-1);
					}
					img.set(I, J);
				}
			}
			cerr << col+1;
		}
		cerr << "\n";
	}

	int sum = 0;
	FOR(0, north.data.size()) {
		int k = north.data[i];
		int m = 1;
		FORj(j, 0, 32) {
			sum += (k & m) != 0;
			m <<= 1;
		}
	}
	cerr << "North area: " << sum/*sqr(deg(dx))*/ << " square degrees";
}

int main(int argc, char **argv)
{
	makeSkyMap();
	return 0;
}
#if 0

void binBitImage(Radians ndx, XImage &img, const BitImage &b)
{
	Radians dx = b.dx;
	int I0 = img.x() / 2;
	int J0 = img.y() / 2;
	double A = sqr(dx);
	FORj(i, 0, b.w)
	{
		FORj(j, 0, b.h)
		{
			if(!b.get(i-b.x0, j-b.y0)) continue;
			Radians x = (i - b.x0)*dx;
			Radians y = (j - b.y0)*dx;

			int I = (int)floor((x - dx/2.) / ndx + .5) + I0;
			int J = (int)floor((y - dx/2.) / ndx + .5) + J0;

			img(I, J) += A;
		}
	}
}

struct gmake
{
	XImage img;
	Radians dx;
};

struct bixel
{
	Radians x[4], y[4], l[4], b[4];	// lambert and galactic coordinates of bixel side centers
	Radians dx, dy, dm;	// bixel size
	double m[2];	// magnitude space

	void init(double d0, double d1, Radians x0, Radians y0, Radians dx, Radians dy)
	{
		m[0] = m0; m[1] = m1; dm = m[1] - m[0];

		this->dx = dx; this->dy = dy;
		x[0] = x0; x[1] = x0-dx/2; x[2] = x0; x[3] = x0+dx/2;
		y[3] = y0; y[0] = y0-dy/2; y[1] = y0; y[2] = y0+dy/2;

		inverse();
	}

	void inverse()
	{
		FOR(0, 4) { proj.inverse(x[i], y[i], l[i], b[i]); }
	}

	void split(int what, bixel &a, bixel &b)
	{
		switch(what)
		{
		case 0:
			a.init(x[0]-dx/4., y[3], dx/2, dy);	// left
			b.init(x[0]+dx/4., y[3], dx/2, dy);	// right
			break;
		case 1:
			a.init(x[0], y[3]-dx/4., dx, dy/2);	// top
			b.init(x[0], y[3]+dx/4., dx, dy/2);	// bottom
			break;
		case 2:
			a = *this; a.m[1] = m[0] + dm/2;	// closer
			b = *this; a.m[0] = m[0] + dm/2;	// farther
			break;
		}
	}
};

class galaxy_model
{
public:
	virtual bool delta(bool split[], const bixel &b) = 0;
	virtual bool generate_stars(const bixel &b) = 0;
	virtual bool next_star() = 0;
};

class disk_model : public galaxy_model
{
public:
	double den0, h, l, z0;		// disk parameters
	float ri0, ri1;
	
	double eps;
public:
	double density(const V3 &p);

	virtual bool delta(bool split[], const bixel &b);
};

double disk_model::density(const double r, const double z)
{
	return den0 * exp(r/l + (z-z0)/h);
}

plx_gri_locus paralax;
double disk_model::delta(bool split[], const bixel &b)
{
	// convert bixel side points to galactic coordinates
	double d0 = paralax.
}

int galaxy_model::number(double den, double V)
{
	double n = den*V;
	int N = (int)n;

	// add Poisson noise
	// TODO: draw a number from N(mu=n, sigma=sqrt(n)) and add it to n

	// the fraction is the probability of having a particle there
	n -= N;
	float prob = float(rand())/RAND_MAX;
	if(prob < n)
	{
		N++;
	}

	return N;
}

void ()
{
	DMMArray<mobject> out;
	out.create("simulation.dmm");

	int I0 = img.x() / 2;
	int J0 = img.y() / 2;

	lambert proj(rad(90), rad(90));

	m0 = 14; m1 = 23; dm = 1;
	epsilon = 0.01;

	// for each pixel in the image
	FOREACH(img)
	{
		Radians x = dx*(i.x - I0);
		Radians y = dx*(i.y - J0);

		// process the beam
		stack<bixel> bixels;
		bixel b;
		for(double m = m0; m <= m1; m += md)
		{
			b.init(m, m+md, x, y, dx, dx);

			// for each model
			FOREACHj(mod, models)
			{
				galaxy_model &model = *mod;
				bixels.push_back(b);
	
				// process the bixel
				bool delta[3];
				vector<bixels> tmpbix;
				while(bixels.size())
				{
					bixel b = bixels.top();
					bixels.pop();

					// check for tolerances
					if(model.delta(delta, b))
					{
						// split the bixel
						tmpbix.push_back(b);
						int at0 = 0, at = 0;
						bixel a, b;
						FORj(j, 0, 3)
						{
							if(!delta[j]) continue;

							at = tmpbix.size();
							FOR(at0, at)
							{
								tmpbix[i].split(j, a, b);
								tmpbix.push_back(a);
								tmpbix.push_back(b);
							}
							at0 = at;
						}

						// push split pixels to stack
						FOR(at0, at) { bixels.push_back(tmpbix[i]); }
						tmpbix.clear();
		
						continue;
					}

					// generate stars
					int N = model.generate_stars(b);
					mobject &obj;
					FORj(n, 0, N)
					{
						model.next_star(obj);
						out.push_back(obj);
					}
				}
			}
		}
	}
}

#endif















