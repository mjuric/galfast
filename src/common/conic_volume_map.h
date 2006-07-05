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

#ifndef conic_volume_map__h
#define conic_volume_map__h

#include <string>
#include <memory>
#include <map>
#include <vector>

#include <astro/io/binarystream.h>
#include <astro/math.h>

#include "dm.h"
#include "projections.h"
#include "model.h"
#include "gpc_cpp.h"

struct conic_volume_map
{
	// filenames for input
	std::string northfn, southfn, modelfn;
	
	// projections are _assumed_ to be north/south hemispheres
	peyton::math::lambert proj_north, proj_south;

	typedef peyton::Radians Radians;

	// input
	Radians dx;						// angular pixel scale (lambert coordinates)
	partitioned_skymap skym_north, skym_south;
	std::auto_ptr<galactic_model> model;			// galactic model used to determine radial pixel lengths
	double d0, d1, dd;					// distance limits, distance step
	int starsperpix;					// approximate number of stars per pixel to aim for

	sdss_color d_field;					// mobject field to be used as distance measure
	sdss_color val_field;					// mobject field to be used as value (for statistical worker functions)

	// output
	struct pix
	{
		double d;	/// center of the pixel (in whatever distance units were used to construct it)
		double dd;	/// radial extend of the pixel (in whatever distance units were used to construct it)

		double dOmega;	/// Angular area of the pixel. The volume is then d^2*dOmega*dd (assuming d, dd are the distance)
		int N;		/// number of stars found in the pixel

		std::vector<double> data;	/// Data being aggregated / result of aggregation. What it is depends on val_field

		pix(double d_ = 0., double dd_ = 0., double dOmega_ = 0.) : d(d_), dd(dd_), dOmega(dOmega_), N(0) {}
		pix(const pix &p) : data(p.data) { d = p.d; dd = p.dd; dOmega = p.dOmega; N = p.N; }
	};
	typedef std::map<double, pix> beam_t;
	typedef std::map<std::pair<int, int>, beam_t> sky_t;

	sky_t sky_north, sky_south;
	int npix;

protected:
	typedef bool (conic_volume_map::*bin_beam_fn)(beam_t &beam, Radians x, Radians y, const mobject &m);
	typedef void (conic_volume_map::*construct_beam_fn)(beam_t &beam, Radians l, Radians b, double area);
	typedef void (conic_volume_map::*finalize_beam_fn)(beam_t &beam, Radians l, Radians b, double area);

	bin_beam_fn bin_beam;
	construct_beam_fn construct_beam;
	finalize_beam_fn finalize_beam;

public:
	// methods
	conic_volume_map();

	void configure_binning(const sdss_color &d, const sdss_color &v, construct_beam_fn cbfn = &conic_volume_map::cb_galmodel, bin_beam_fn bbfn = &conic_volume_map::bb_bin_sum, finalize_beam_fn fbfn = &conic_volume_map::fin_average)
	{
		d_field = d; val_field = v;
		bin_beam = bbfn; construct_beam = cbfn; finalize_beam = fbfn;
	}
	void initialize();
	void initialize_for_galmodel();

	peyton::io::obstream &serialize(peyton::io::obstream &out) const;
	peyton::io::ibstream &unserialize(peyton::io::ibstream &out);
	std::ostream &text_serialize(std::ostream &out);

	bool bin(Radians l, Radians b, const mobject &m);
	bool get_beam(beam_t *beam, Radians l, Radians b);
	bool get_sky_and_index(int &skyidx, int &X, int &Y, const Radians l, const Radians b);

protected:
	void create_pixels(peyton::math::lambert &proj, partitioned_skymap &sky, sky_t &sky);
	void create_pixel(beam_t &beam, double dprev, double d, double area);

	bool bin(double x, double y, const mobject &m, partitioned_skymap &msky, sky_t &sky);

	void print_hemisphere(std::ostream &out, peyton::math::lambert &proj, partitioned_skymap &skym, sky_t &sky);

public: // various binning type functions
	// bin construction
	void cb_galmodel(beam_t &beam, Radians l, Radians b, double area);
	void cb_equidistant(beam_t &beam, Radians l, Radians b, double area);

	// binning
	bool bb_bin_sum(beam_t &beam, Radians l, Radians b, const mobject &m);
	bool bb_bin_collect(beam_t &beam, Radians l, Radians b, const mobject &m);

	// final bin processing functions
	void fin_average(beam_t &beam, Radians l, Radians b, double area);
	void fin_median(beam_t &beam, Radians l, Radians b, double area);
};

// spline interpolation of the conic volume map
struct conic_volume_map_interpolator
{
public:
	conic_volume_map *vm;
	struct beam_t {
		double dmin, dmax;
		spline spl;
		double dOmega;	/// solid angle
	};
	typedef std::map<std::pair<int, int>, beam_t> beams_t;
	beam_t median_beam;

	beams_t north, south;
public:
	conic_volume_map_interpolator();

	void initialize(conic_volume_map &vm);
	void text_dump(std::ostream &out);
	void text_dump_median_beam(std::ostream &out);

	double interpolate(peyton::Radians l, peyton::Radians b, double D);
	double operator()(peyton::Radians l, peyton::Radians b, double D) { return interpolate(l, b, D); }
	
	bool empty() { return vm == NULL; }

protected:
	void set_aux(beams_t &beams, conic_volume_map::sky_t &sky);
};

OSTREAM(const conic_volume_map::pix &p);
OSTREAM(const conic_volume_map &l);

inline BOSTREAM2(conic_volume_map &v) { return v.serialize(out); }
inline BISTREAM2(conic_volume_map &v) { return v.unserialize(in); }

BOSTREAM2(const conic_volume_map::pix &p);
BISTREAM2(conic_volume_map::pix &p);

#endif
