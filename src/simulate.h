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

#ifndef simulate_h__
#define simulate_h__

#include "projections.h"
#include "model.h"
#include "dm.h"
#include <string>

#include <astro/system/memorymap.h>
#include <boost/shared_ptr.hpp>

class cumulative_dist
{
public:
	typedef std::pair<double, double> value_type;
	std::vector<value_type> hist;
public:
	double operator()(double prob);
	void construct(const std::vector<double> &x, const std::vector<double> &y);
};

//
// some supporting datastructures
//

// cumulative marginal pdf in lambert y coodinate, for a given lambert x coordinate
// (Sy::pcum == probability that a star 
// within the observation footprint and given the model and x coordinate,
// will have the value of the lambert coordinate y *less* than Sx::Y)
struct Sy
{
	double pcum;
	int Y;
	std::vector<double> ripdf;
	Sy(double _pcum = 0, int _Y = 0) : pcum(_pcum), Y(_Y) {}
private:
	Sy(const Sy&);
	Sy& operator=(const Sy&);
};
inline bool operator<(const Sy& a, const Sy& b) { return a.pcum < b.pcum; }

// cumulative marginal pdf in lambert x coodinate (Sx::pcum == probability that a star 
// within the observation footprint and given the model will have a 
// the value of the lambert coordinate x *less* than Sx::X)
struct Sx
{
	double pcum;
	int X;
	std::vector<boost::shared_ptr<Sy> > ypdf;
	Sx(double _pcum = 0, int _X = 0) : pcum(_pcum), X(_X) {}
private:
	Sx(const Sx&);
	Sx& operator=(const Sx&);
};
inline bool operator<(const Sx& a, const Sx& b) { return a.pcum < b.pcum; }

// set of marginal pdfs in lambert x. It's the trunk of a tree
// of marginal PDFs, starting with mpdf in X, then Y, and
// finaly r-i. First branches of the tree are Sx structures,
// followed by Sy structures (stored in Sx::ypdf), followed by
// doubles of r-i mpdfs (in Sy::ripdf vector).
typedef std::vector<boost::shared_ptr<Sx> > XPDF;


struct star_output_function;
class model_pdf
{
public:
	struct star // storage structure for Monte Carlo generated stars
	{
		double x, y, ri, m;
		int X(const model_pdf *gsim) const { return (int)((x - gsim->x0)/gsim->dx); }
		int Y(const model_pdf *gsim) const { return (int)((y - gsim->y0)/gsim->dx); }
		int RI(const model_pdf *gsim) const { return (int)((ri - gsim->ri0)/gsim->dri); }
	};

public:
	std::string prefix;	// prefix for input/output datafiles
	galactic_model *model;	// models for this simulation

	double dx; 		// model grid angular resolution
	double dri;		// model CMD resolution
	double ri0, ri1;	// color limits
	double m0, m1;		// magnitude limits
	double dm;		// model CMD magnitude resolution

public: // internal variables - usually you don't need to touch these (for northern sky)
	peyton::math::lambert proj;	// lambert projector object (by default, centered at north pole)
	double x0, x1, y0, y1;		// lambert survey footprint bounding box

	std::map<std::pair<int, int>, gpc_polygon> skymap;	// a map of rectangular sections of the sky, for fast is-point-in-survey-area lookup
	XPDF xpdf; 						// marginal cumulative probability distribution function
	double N; 						// Grand total - number of stars in the whole field
public:
	model_pdf(const std::string &pfx, galactic_model *model);

	// PDF generation functions
	void precalculate_mpdf();
	void magnitude_mpdf(cumulative_dist &mspl, double x, double y, double ri);
	void store(const std::string &outfn);

	// montecarlo generation functions
	bool load(const std::string &prefix);
	bool draw_position(star &s, gsl_rng *rng);
	void draw_magnitudes(std::vector<model_pdf::star> &stars, gsl_rng *rng);

protected:
	double ri_mpdf(std::vector<double> &pdf, const double x, const double y);
};

class sky_generator
{
public:
	std::vector<model_pdf *> pdfs;
	std::vector<boost::shared_ptr<std::vector<model_pdf::star> > > stars;	// generated stars
public:
	void add_pdf(model_pdf &pdf);
	void montecarlo(unsigned int K, star_output_function &sf, gsl_rng *rng);
};

struct star_output_function
{
	virtual void output (const std::vector<model_pdf::star> &stars, peyton::math::lambert &proj) = 0;
};

struct star_output_to_dmm : star_output_function
{
public:
	proc_obs_info po_info;
	peyton::system::DMMArray<mobject> out;
	peyton::system::DMMArray<obsv_mag> starmags;
public:
	star_output_to_dmm(const std::string &objcat, const std::string &obscat, bool create);
	void close();

	virtual void output (const std::vector<model_pdf::star> &stars, peyton::math::lambert &proj);
};

struct star_output_to_textstream : star_output_function
{
protected:
	std::ostream &out;
public:
	star_output_to_textstream(std::ostream &out_) : out(out_) {}

	virtual void output (const std::vector<model_pdf::star> &stars, peyton::math::lambert &proj);
};

#endif // simulate_h__
