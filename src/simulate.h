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
#include "conic_volume_map.h"

#include <string>

#include <astro/system/memorymap.h>
#include <astro/math.h>
#include <boost/shared_ptr.hpp>
#include <astro/io/binarystream.h>

class cumulative_dist
{
public:
	typedef std::pair<double, double> value_type;
	std::vector<value_type> hist;
public:
	double operator()(double prob);
	void construct(const std::vector<double> &x, const std::vector<double> &y);
	bool empty() const { return hist.empty(); }
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
		double x, y;		// Lambert coordiates on the sky
		double ri, m;		// color and apparent magnitude
		int X(const partitioned_skymap &gsim) const { return (int)((x - gsim.x0)/gsim.dx); }
		int Y(const partitioned_skymap &gsim) const { return (int)((y - gsim.y0)/gsim.dx); }
		int RI(const model_pdf &gsim) const { return (int)((ri - gsim.ri0)/gsim.dri); }
	};

public: // serialized object state
	double ri0, ri1;	// color limits -- NOTE: this should be converted to absolute magnitude
	double dri;		// color resolution -- NOTE: this should be converted to absolute magnitude
	double m0, m1;		// apparent magnitude limits
	double dm;		// apparent magnitude resolution

	peyton::math::lambert proj;	// lambert projector object (loaded from .xgpc.txt file)
	partitioned_skymap skymap;	// a map of rectangular sections of the sky, for fast is-point-in-survey-area lookup

	XPDF xpdf; 			// marginal cumulative probability distribution function
	double N; 			// Grand total - number of stars in the whole field

protected: // temporary internal object state (NOT serialized)
	std::auto_ptr<galactic_model> model;	// model from which to generate this PDF
	std::string pdfname;	// name of this model realization (usually the .conf or .pdf.bin filename)

public:
	model_pdf(const std::string &pdfname = "UNNAMED");
	model_pdf(std::istream &in, const std::string &pdfname = "UNNAMED");

	// PDF generation functions
	void construct_mpdf(const std::string &footfn, const std::string &modelfn);
	void magnitude_mpdf(cumulative_dist &mspl, double x, double y, double ri);
	peyton::io::obstream &serialize(peyton::io::obstream &out) const;

	// montecarlo generation functions
	peyton::io::ibstream &unserialize(peyton::io::ibstream &in);
	bool draw_position(star &s, gsl_rng *rng);
	void draw_magnitudes(std::vector<model_pdf::star> &stars, gsl_rng *rng);

	// accessors
	const std::string &name() const { return pdfname; }
	galactic_model &galmodel() { return *model.get(); }
	
protected:
	double ri_mpdf(std::vector<double> &pdf, const double x, const double y);
};
inline BOSTREAM2(const model_pdf &pdf) { return pdf.serialize(out); }
inline BISTREAM2(model_pdf &pdf) { return pdf.unserialize(in); }

class sky_generator
{
protected:
	gsl_rng *rng;

	std::vector<boost::shared_ptr<model_pdf> > pdfs;
	std::vector<boost::shared_ptr<std::vector<model_pdf::star> > > stars;	// generated stars

	int nstars;	///< number of stars to generate (read from config file)
#if 0
	conic_volume_map magerrs_cvm;
	conic_volume_map_interpolator magerrs;
	
	double Ar;	///< extinction, conf["Ar"], default = 0
protected:
	double constant_photo_error;	///< photometric error to add to all magnitudes
	double paralax_dispersion;	///< Gaussian dispersion of photometric paralax relation
	int flags;
	static const int APPLY_PHOTO_ERRORS		= 0x00000001;
#endif
public:
	sky_generator(std::istream &in, const std::string &pdfs);

	void add_pdf(boost::shared_ptr<model_pdf> &pdf);
	void add_pdf(model_pdf &pdf);

	void montecarlo(star_output_function &sf);

	~sky_generator();
protected:
	int montecarlo_batch(star_output_function &out, int Ktotal, const std::vector<double> &modelCPDF, bool allowMisses);
	//void observe(const std::vector<model_pdf::star> &stars, peyton::math::lambert &proj, star_output_function &sf);
	void draw_stars(const std::vector<model_pdf::star> &stars, galactic_model &model, peyton::math::lambert &proj, star_output_function &sf);
	void draw_companion(float &gb, float &rb, float &ib, peyton::Radians l, peyton::Radians b, double dm /*distance modulus*/);
};

struct star_output_function
{
	typedef peyton::Radians Radians;

	virtual void output(Radians ra, Radians dec, double Ar, std::vector<std::pair<observation, obsv_id> > &obsvs) = 0;
	//virtual void output(Radians l, Radians b, double ri, double r, sstruct &t) = 0;
	virtual void output(sstruct &t) = 0;
	virtual void output_header(const sstruct::factory_t &factory) {};
	virtual ~star_output_function() {}
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

	virtual void output(Radians ra, Radians dec, double Ar, std::vector<std::pair<observation, obsv_id> > &obsvs);
	virtual void output(sstruct &t);
	//virtual void output(Radians l, Radians b, double ri, double r, sstruct &t);
};

struct star_output_to_textstream : star_output_function
{
protected:
	std::ostream &out;
public:
	star_output_to_textstream(std::ostream &out_) : out(out_) {}

	virtual void output(Radians ra, Radians dec, double Ar, std::vector<std::pair<observation, obsv_id> > &obsvs);
	virtual void output(sstruct &t);
	//virtual void output(Radians l, Radians b, double ri, double r, sstruct &t);
	virtual void output_header(const sstruct::factory_t &factory);
};

void observe_catalog(const std::string &conffn, const std::string &input, const std::string &output);

#endif // simulate_h__
