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

#ifdef COMPILE_SIMULATE_X
#define COMPILING_SIMULATE

#include "gsl/gsl_randist.h"

#include "gpc_cpp.h"

#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <sstream>

#include "simulate.h"
#include "projections.h"
#include "model.h"
#include "paralax.h"
#include "analysis.h"
#include "dm.h"
#include "io.h"

#include <vector>
#include <map>
#include <string>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

using namespace boost::lambda;

// Get the autoconf/automake set datadir
// TODO: This should be moved into its separate file
const std::string &datadir()
{
	static std::string dd;
	static const char *dd_hardcoded = DATADIR;
	if(dd.empty())
	{
		EnvVar ev("DATADIR");
		dd = ev ? ev.c_str() : dd_hardcoded;
		MLOG(verb2) << "datadir=" << dd << (ev ? " (initializes from $DATADIR)" : "");
	}
	return dd;
}

class osink;
class opipeline_stage
{
	protected:
		std::set<std::string> prov, req;

		osink *nextlink;
	public:
		void chain(osink *nl) { nextlink = nl; }
		virtual int run(gsl_rng *rng) = 0;

	public:
		static boost::shared_ptr<opipeline_stage> create(const std::string &name);
		virtual const std::string &name() const = 0;
		virtual const std::string &type() const { static std::string s("stage"); return s; }

	public:
		virtual bool init(const Config &cfg) = 0;
		const std::set<std::string> &provides() const { return prov; }
		const std::set<std::string> &requires() const { return req; }

		bool satisfied_with(const std::set<std::string> &haves);
		bool provides_any_of(const std::set<std::string> &needs, std::string &which);
		virtual int priority() { return 0; }

	public:
		bool inits(const std::string &cfgstring) { return inits(cfgstring.c_str()); }
		bool inits(const char *cfgstring)
		{
			std::istringstream ss(cfgstring);
			Config cfg;
			cfg.load(ss);
			return init(cfg);
		}

		opipeline_stage() : nextlink(NULL)
		{
		}
};

bool opipeline_stage::satisfied_with(const std::set<std::string> &haves)
{
	FOREACH(req)
	{
		if(!haves.count(*i)) {
			DLOG(verb2) << "Failed on: " << *i;
			return false;
		}
	}
	return true;
}

bool opipeline_stage::provides_any_of(const std::set<std::string> &needs, std::string &which)
{
	FOREACH(needs)
	{
		if(prov.count(*i)) { which = *i; return true; }
	}
	return false;
}

class osink : public opipeline_stage
{
	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng) = 0;

	public:
		virtual int run(gsl_rng *rng) { THROW(ENotImplemented, "We should have never gotten here"); } // we should never reach this place

	public:
		osink() : opipeline_stage()
		{
			req.insert("_source");
		}
};

// add Fe/H information
class os_FeH : public osink
{
protected:
	double A[2], sigma[3], offs[3];
	double Hmu, muInf, DeltaMu;

public:
	virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
	virtual bool init(const Config &cfg);
	virtual const std::string &name() const { static std::string s("FeH"); return s; }

	os_FeH() : osink()
	{
		prov.insert("FeH");
		req.insert("comp");
		req.insert("XYZ[3]");
	}
};

size_t os_FeH::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input
	//	- all stars are main sequence
	for(size_t i=0; i != count; i++)
	{
		sstruct &s = in[i];

		// fetch prerequisites
		const int component = s.component();
		float Z = s.XYZ()[2];
		float &FeH = s.FeH();

		switch(component)
		{
			case BahcallSoneira_model::THIN:
			case BahcallSoneira_model::THICK: {
				// choose the gaussian to draw from
				double p = gsl_rng_uniform(rng)*(A[0]+A[1]);
				int i = p < A[0] ? 0 : 1;

				// calculate mean
				double muD = muInf + DeltaMu*exp(-fabs(Z)/Hmu);		// Bond et al. A2
				double aZ = muD - 0.067;

				// draw
				FeH = gsl_ran_gaussian(rng, sigma[i]);
				FeH += aZ + offs[i];
/*				if(FeH > 0 && Z > 900) { std::cerr << Z << " : " << FeH << " " << muD << " " << aZ << " " << offs[i] << " " << i << "\n"; }
				if(FeH > 0 && Z > 900) { std::cerr << muInf << " " << DeltaMu << " " << abs(Z/1000) << " " << Hmu << " " << -abs(Z/1000)/Hmu << " " << i << "\n"; }*/
			} break;
			case BahcallSoneira_model::HALO:
				FeH = offs[2] + gsl_ran_gaussian(rng, sigma[2]);
				break;
			default:
				THROW(ENotImplemented, "We should have never gotten here");
		}
	}

	return nextlink->push(in, count, rng);
}

bool os_FeH::init(const Config &cfg)
{
	cfg.get(A[0],     "A0",     0.63);
	cfg.get(sigma[0], "sigma0", 0.20);
	cfg.get(offs[0],  "offs0",  0.00);
	cfg.get(A[1],     "A1",     0.37);
	cfg.get(sigma[1], "sigma1", 0.20);
	cfg.get(offs[1],  "offs1",  0.14);
	
	cfg.get(Hmu,      "Hmu",     500.);
	cfg.get(muInf,    "muInf",  -0.82);
	cfg.get(DeltaMu,  "deltaMu", 0.55);

	// renormalize disk gaussian amplitudes to sum up to 1
	double sumA = A[0] + A[1];
	A[0] /= sumA; A[1] /= sumA;

	// Halo configuration
	cfg.get(sigma[2], "sigmaH",  0.30);
	cfg.get(offs[2],  "offsH",  -1.46);

	// Output model parameters
	MLOG(verb1) << "Normalized disk amplitudes  (A[0], A[1]): "<< A[0] << " " << A[1];
	MLOG(verb1) << "Disk sigma          (sigma[0], sigma[1]): "<< sigma[0] << " " << sigma[1];
	MLOG(verb1) << "Disk offsets          (offs[0], offs[1]): "<< offs[0] << " " << offs[1];
	MLOG(verb1) << "Disk median Z dep. (muInf, deltaMu, Hmu): "<< muInf << " " << DeltaMu << " " << Hmu;
	MLOG(verb1) << "Halo gaussian              (muH, sigmaH): "<< offs[2] << " " << sigma[2];

	return true;
}



// add Fe/H information
class os_fixedFeH : public osink
{
	protected:
		float FeH;

	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual const std::string &name() const { static std::string s("fixedFeH"); return s; }

		os_fixedFeH() : osink(), FeH(0)
		{
			prov.insert("FeH");
		}
};

size_t os_fixedFeH::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input
	//	- all stars are main sequence
	for(size_t i=0; i != count; i++)
	{
		sstruct &s = in[i];

		// fetch prerequisites
		s.FeH() = FeH;
	}

	return nextlink->push(in, count, rng);
}

bool os_fixedFeH::init(const Config &cfg)
{
	if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	cfg.get(FeH, "FeH", 0.f);

	return true;
}

void print_matrix(gsl_matrix *m)
{
/* print matrix the hard way */
  printf("Matrix m\n");
  for (int i=0;i<m->size1;i++)
    {
      for (int j=0;j<m->size2;j++)
	{
	  fprintf(stderr, "%f ",gsl_matrix_get(m,i,j));
	}
      fprintf(stderr, "\n");
    }
  fprintf(stderr, "\n");
}

struct trivar_gauss
{
	gsl_matrix *A;
	gsl_vector *Z;

	trivar_gauss()
	{
		A = gsl_matrix_alloc(3, 3);
		Z = gsl_vector_alloc(3);

		gsl_matrix_set_zero(A);
	}

	void set(double s11, double s12, double s13, double s22, double s23, double s33)
	{
		// populate A (assumes the upper triang is already 0), calculate Cholesky decomp
		gsl_matrix_set(A, 0, 0, sqr(s11));
		gsl_matrix_set(A, 1, 0,     s12 ); gsl_matrix_set(A, 1, 1, sqr(s22));
		gsl_matrix_set(A, 2, 0,     s13 ); gsl_matrix_set(A, 2, 1,     s23 ); gsl_matrix_set(A, 2, 2, sqr(s33));
		//print_matrix(A);

		int status = gsl_linalg_cholesky_decomp(A);
		ASSERT(status == 0);
		//print_matrix(A); std::cerr << "status=" << status << "\n";
	}

	void draw(gsl_vector *y, gsl_rng *rng, bool zero = false)
	{
		gsl_vector_set(Z, 0, gsl_ran_ugaussian(rng));
		gsl_vector_set(Z, 1, gsl_ran_ugaussian(rng));
		gsl_vector_set(Z, 2, gsl_ran_ugaussian(rng));

		if(zero) { gsl_vector_set_zero(y); }
		//gsl_blas_dgemv(CblasNoTrans, 1., A, Z, 1., y);
		gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, A, Z);
		gsl_vector_add(y, Z);
		//std::cout << "XXXXX: " << y->data[0] << " " << y->data[1] << " " << y->data[2] << "\n";
	}

	~trivar_gauss()
	{
		gsl_vector_free(Z);
		gsl_matrix_free(A);
	}
};

// add kinematic information
class os_kinTMIII : public osink
{
	protected:
		typedef std::vector<double> dvec;

		trivar_gauss tri_rnd;

		double fk, DeltavPhi;
		dvec 	vPhi1, vPhi2, vR, vZ,
			sigmaPhiPhi1, sigmaPhiPhi2, sigmaRR, sigmaZZ, sigmaRPhi, sigmaZPhi, sigmaRZ,
			HvPhi, HvR, HvZ,
			HsigmaPhiPhi, HsigmaRR, HsigmaZZ, HsigmaRPhi, HsigmaZPhi, HsigmaRZ;
		dvec *diskEllip[6], *haloEllip[6], *diskMeans[3], *haloMeans[3];

	public:
		void add_dispersion(double v[3], double Rsquared, double Z, dvec *ellip[6], gsl_rng *rng);
		void compute_means(double v[3], double Rsquared, double Z, dvec *means[3]);

		void get_disk_kinematics(double v[3], double Rsquared, double Z, gsl_rng *rng, bool &firstGaussian);
		void get_halo_kinematics(double v[3], double Rsquared, double Z, gsl_rng *rng);

	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual const std::string &name() const { static std::string s("kinTMIII"); return s; }

		os_kinTMIII() : osink()
		{
			prov.insert("vcyl[3]");
			req.insert("comp");
			req.insert("XYZ[3]");
		}
};

void test_kin()
{
	return;
// 	Radians l, b;
// 	equgal(0, ctn::pi/2., l, b);
// 	printf("%.10f\n", deg(l));
// 	exit(-1);

	os_kinTMIII o;

	Config cfg;
	o.init(cfg);

	double v[3];
	bool firstGaussian;

	float XYZ[3] = {8000, 0, 200};
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

	// R^2 and Z in kpc
	const double Rsquared = 1e-6 * (sqr(XYZ[0]) + sqr(XYZ[1]));
	const double Z = 1e-3 * XYZ[2];

	trivar_gauss tri_rnd;
	tri_rnd.set(20, 10, 20,   20, 20,   30);
	double vvv[3] = { 0., 0., 0. };
	gsl_vector *vv = &gsl_vector_view_array(vvv, 3).vector;
	tri_rnd.draw(vv, rng);
	exit(-1);

	std::cout << "# XYZ = " << XYZ[0] << "," << XYZ[1] << "," << XYZ[2] << " (pc)\n";
	for(int i = 0; i != 1; i++)
	{
		o.get_disk_kinematics(v, Rsquared, Z, rng, firstGaussian);
		std::cout << v[0] << " " << v[1] << " " << v[2] << " " << firstGaussian << "\n";
	}
//	std::cerr << "XYZ = " << XYZ[0] << "," << XYZ[1] << "," << XYZ[2] << "\n";
//	std::cerr << "  v = " << v[0] << "," << v[1] << "," << v[2] << "\n";
	exit(-1);
}

size_t os_kinTMIII::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input
	double tmp[3]; bool firstGaussian;
	for(size_t i=0; i != count; i++)
	{
		sstruct &s = in[i];

		// fetch prerequisites
		const int component = s.component();
		float *XYZ   = s.XYZ();
		float *v     = s.vcyl();
		const double Rsquared = 1e-6 * (sqr(XYZ[0]) + sqr(XYZ[1]));
		const double Z = 1e-3 * XYZ[2];

		switch(component)
		{
			case BahcallSoneira_model::THIN:
			case BahcallSoneira_model::THICK:
				//getDiskKinematics(v[0], v[1], v[2], Z, rng);
				get_disk_kinematics(tmp, Rsquared, Z, rng, firstGaussian);
				break;
			case BahcallSoneira_model::HALO:
				//getHaloKinematics(v[0], v[1], v[2], Z, rng);
				get_halo_kinematics(tmp, Rsquared, Z, rng);
				break;
			default:
				THROW(ENotImplemented, "We should have never gotten here");
		}
		v[0] = tmp[0]; v[1] = tmp[1]; v[2] = tmp[2];
	}

	return nextlink->push(in, count, rng);
}

template<typename T>
	int split(T& arr, const std::string &text)
{
	typename T::value_type tmp;
	std::stringstream ss(text);

	arr.clear();
	while(ss >> tmp) { arr.push_back(tmp); }
	return arr.size();
}

template<typename T>
	T split(const std::string &text)
{
	T ret;
	split(ret, text);
	return ret;
}

inline double modfun(double Rsquared, double Z, double a, double b, double c, double d, double e)
{
	return a + b*pow(fabs(Z), c) + d*pow(Rsquared, 0.5*e);
}

void os_kinTMIII::add_dispersion(double v[3], double Rsquared, double Z, dvec *ellip[6], gsl_rng *rng)
{
	// compute velocity dispersions at this position, and draw from trivariate gaussian
	// NOTE: ADDS THE RESULT TO v, DOES NOT ZERO v BEFORE !!
	double sigma[6];
	FOR(0, 6)
	{
		dvec &p = *ellip[i];
		sigma[i] = modfun(Rsquared, Z, p[0], p[1], p[2], p[3], p[4]);
	}
	gsl_vector *vv = &gsl_vector_view_array(v, 3).vector;
	tri_rnd.set(sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5]);
	tri_rnd.draw(vv, rng);
}

void os_kinTMIII::compute_means(double v[3], double Rsquared, double Z, dvec *means[3])
{
	// returns means in v[3]
	FOR(0, 3)
	{
		dvec &p = *means[i];
		v[i] = modfun(Rsquared, Z, p[0], p[1], p[2], p[3], p[4]);
	}
}

void os_kinTMIII::get_disk_kinematics(double v[3], double Rsquared, double Z, gsl_rng *rng, bool &firstGaussian)
{
	// set up which gaussian are we drawing from
	double p = gsl_rng_uniform(rng);
	if(firstGaussian = (p < fk))
	{
		// first gaussian
		diskMeans[1] = &vPhi1;
		diskEllip[3] = &sigmaPhiPhi1;
	}
	else
	{
		// second gaussian
		diskMeans[1] = &vPhi2;
		diskEllip[3] = &sigmaPhiPhi2;
	}

	compute_means(v, Rsquared, Z, diskMeans);
	// truncate v_phi > 0
	if(v[1] > 0.) { v[1] = 0.; }

	add_dispersion(v, Rsquared, Z, diskEllip, rng);
}

void os_kinTMIII::get_halo_kinematics(double v[3], double Rsquared, double Z, gsl_rng *rng)
{
	compute_means(v, Rsquared, Z, haloMeans);
	add_dispersion(v, Rsquared, Z, haloEllip, rng);
}

template<typename T> inline OSTREAM(const std::vector<T> &v) { FOREACH(v) { out << *i << " "; }; return out; }

bool os_kinTMIII::init(const Config &cfg)
{
	cfg.get(fk           , "fk"           , 3.0);
	cfg.get(DeltavPhi    , "DeltavPhi"    , 34.0);
	fk = fk / (1. + fk);	// renormalize to probability of drawing from the first gaussian

	cfg.get(vR           , "vR"           , split<dvec>("0 0 0 0 0"));	diskMeans[0] = &vR;
	cfg.get(vPhi1        , "vPhi"         , split<dvec>("-194 19.2 1.25 0 0"));	diskMeans[1] = &vPhi1;
	cfg.get(vZ           , "vZ"           , split<dvec>("0 0 0 0 0"));	diskMeans[2] = &vZ;
	cfg.get(sigmaRR      , "sigmaRR"      , split<dvec>("40 5 1.5 0 0"));	diskEllip[0] = &sigmaRR;
	cfg.get(sigmaRPhi    , "sigmaRPhi"    , split<dvec>("0 0 0 0 0"));	diskEllip[1] = &sigmaRPhi;
	cfg.get(sigmaRZ      , "sigmaRZ"      , split<dvec>("0 0 0 0 0"));	diskEllip[2] = &sigmaRZ;
	cfg.get(sigmaPhiPhi1 , "sigmaPhiPhi1" , split<dvec>("12 1.8 2 0 11"));	diskEllip[3] = &sigmaPhiPhi1;	// dynamically changeable
	cfg.get(sigmaPhiPhi2 , "sigmaPhiPhi2" , split<dvec>("34 1.2 2 0 0"));
	cfg.get(sigmaZPhi    , "sigmaZPhi"    , split<dvec>("0 0 0 0 0"));	diskEllip[4] = &sigmaZPhi;
	cfg.get(sigmaZZ      , "sigmaZZ"      , split<dvec>("25 4 1.5 0 0"));	diskEllip[5] = &sigmaZZ;

	// v2 is v1 + DeltavPhi, which is what this does.
	vPhi2 = vPhi1; vPhi2[0] += DeltavPhi;

	cfg.get(HvR          , "HvR"          , split<dvec>("0 0 0 0 0"));	haloMeans[0] = &HvR;
	cfg.get(HvPhi        , "HvPhi"        , split<dvec>("0 0 0 0 0"));	haloMeans[1] = &HvPhi;
	cfg.get(HvZ          , "HvZ"          , split<dvec>("0 0 0 0 0"));	haloMeans[2] = &HvZ;
	cfg.get(HsigmaRR     , "HsigmaRR"     , split<dvec>("135 0 0 0 0"));	haloEllip[0] = &HsigmaRR;
	cfg.get(HsigmaRPhi   , "HsigmaRPhi"   , split<dvec>("0 0 0 0 0"));	haloEllip[1] = &HsigmaRPhi;
	cfg.get(HsigmaRZ     , "HsigmaRZ"     , split<dvec>("0 0 0 0 0"));	haloEllip[2] = &HsigmaRZ;
	cfg.get(HsigmaPhiPhi , "HsigmaPhiPhi" , split<dvec>("85 0 0 0 0"));	haloEllip[3] = &HsigmaPhiPhi;
	cfg.get(HsigmaZPhi   , "HsigmaZPhi"   , split<dvec>("0 0 0 0 0"));	haloEllip[4] = &HsigmaZPhi;
	cfg.get(HsigmaZZ     , "HsigmaZZ"     , split<dvec>("85 0 0 0 0"));	haloEllip[5] = &HsigmaZZ;

	// some info
	MLOG(verb1) << "Disk gaussian normalizations: " << fk << " : " << (1-fk);
	MLOG(verb1) << "Second disk gaussian offset:  " << DeltavPhi;

	MLOG(verb1) << "vR coefficients:              " << vR;
	MLOG(verb1) << "vZ coefficients:              " << vZ;
	MLOG(verb1) << "sigmaRR coefficients:         " << sigmaRR;
	MLOG(verb1) << "sigmaRPhi coefficients:       " << sigmaRPhi;
	MLOG(verb1) << "sigmaRZ coefficients:         " << sigmaRZ;
	MLOG(verb1) << "sigmaPhiPhi1 coefficients:    " << sigmaPhiPhi1;
	MLOG(verb1) << "sigmaPhiPhi2 coefficients:    " << sigmaPhiPhi2;
	MLOG(verb1) << "sigmaZPhi coefficients:       " << sigmaZPhi;
	MLOG(verb1) << "sigmaZZ coefficients:         " << sigmaZZ;

	MLOG(verb1) << "HvR coefficients:             " << HvR;
	MLOG(verb1) << "HvZ coefficients:             " << HvZ;
	MLOG(verb1) << "HsigmaRR coefficients:        " << HsigmaRR;
	MLOG(verb1) << "HsigmaRPhi coefficients:      " << HsigmaRPhi;
	MLOG(verb1) << "HsigmaRZ coefficients:        " << HsigmaRZ;
	MLOG(verb1) << "HsigmaPhiPhi coefficients:    " << HsigmaPhiPhi;
	MLOG(verb1) << "HsigmaZPhi coefficients:      " << HsigmaZPhi;
	MLOG(verb1) << "HsigmaZZ coefficients:        " << HsigmaZZ;

	return true;
}



// convert input absolute/apparent magnitudes to ugriz colors
class os_photometry : public osink
{
	protected:
		std::string bandset;			// name of this filter set
		std::string bband;			// band off which to bootstrap other bands, using color relations. Must be supplied by other modules.
		std::string absbband;			// Absolute magnitude band for which the datafile gives col(absmag,FeH) values. By default, it's equal to "abs$bband". Must be supplied by other modules.
		std::string photoFlags;			// Name of the photometric flags field
		int bidx;				// index of bootstrap band in bnames
		size_t offset_absmag, offset_mag;	// sstruct offsets to bootstrap apparent and absolute magnitudes [input]
		size_t offset_mags;			// sstruct offset to magnitudes in computed bands [output]
		size_t offset_photoflags;		// sstruct offset to photometric flags [outout]
		std::vector<std::string> bnames;	// band names (e.g., LSSTr, LSSTg, SDSSr, V, B, R, ...)
		std::vector<std::vector<float> > clt;
		typedef char cbool;			// to avoid the special vector<bool> semantics, while maintaining a smaller memory footprint than vector<int>
		std::vector<std::vector<cbool> > eclt;	// extrapolation flags
		int nMr, nFeH;
		double Mr0, Mr1, dMr;
		double FeH0, FeH1, dFeH;
		int ncolors;

		float color(int ic, double FeH, double Mr, bool *e = NULL)
		{
			ASSERT(ic >= 0 && ic < clt.size());
			ASSERT(Mr0 <= Mr && Mr <= Mr1) { std::cerr << Mr0 << " <= " << Mr << " <= " << Mr1 << "\n"; }
			ASSERT(FeH0 <= FeH && FeH <= FeH1) { std::cerr << FeH0 << " <= " << FeH << " <= " << FeH1 << "\n"; }

			int f = (int)((FeH - FeH0) / dFeH);
			int m = (int)((Mr  -  Mr0) / dMr);
			int idx = m*nFeH + f;
			if(e) { *e = eclt[ic][idx]; }
			return clt[ic][idx];
		}
	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual const std::string &name() const { static std::string s("photometry"); return s; }

		os_photometry() : osink(), offset_absmag(-1), offset_mag(-1), offset_mags(-1)
		{
			req.insert("FeH");
		}
};

struct Mr2col
{
	static const size_t maxcolors = 20;	// ug,gr,ri,iz,zy

	double Mr;
	double c[maxcolors];

	bool operator <(const Mr2col &a) const { return Mr < a.Mr; }
};

struct colsplines
{
	double Mrmin, Mrmax;
	spline s[Mr2col::maxcolors];
	spline &operator [](const size_t i) { return s[i]; }
};

bool os_photometry::init(const Config &cfg)
{
	// band definitions
	std::string tmp, bname;
	cfg.get(bandset,   "bandset",   "LSSTugrizy");

	cfg.get(tmp,   "bands",   "LSSTu LSSTg LSSTr LSSTi LSSTz LSSTy");
	std::istringstream ss(tmp);
	while(ss >> bname)
	{
		bnames.push_back(bname);
	}
	ncolors = bnames.size()-1;
	if(ncolors >= Mr2col::maxcolors)
	{
		THROW(EAny, "This module cannot handle more than " + str(Mr2col::maxcolors+1) + " bands.");
	}

	// check if the bandset has the number of bands appended
	if(bandset.find('[') == std::string::npos || bandset.find(']') == std::string::npos)
	{
		bandset += '[' + str(bnames.size()) + ']';
	}
	prov.insert(bandset);

	// deduce photoFlags name
	std::string bandsetname = bandset.substr(0, bandset.find('['));
	photoFlags = bandsetname + "PhotoFlags{i|class=flags,fmt=%5d}";
	//std::cerr << bandset << " " << bandsetname << " " << photoFlags << "\n"; exit(0);
	prov.insert(photoFlags);

	// bootstap band setup
	cfg.get(bband,   "bootstrap_band",   "LSSTr");
	cfg.get(absbband,   "absband",   "abs"+bband);
	req.insert(bband);
	req.insert(absbband);
	typeof(bnames.begin()) b = std::find(bnames.begin(), bnames.end(), bband);
	if(b == bnames.end())
	{
		THROW(EAny, "Bootstrap band must be listed in the 'bands' keyword.");
	}
	bidx = b - bnames.begin();

	// sampling parameters
	cfg.get(tmp,   "absmag_grid",   "3 15 0.01");
	ss.clear(); ss.str(tmp);
	if(!(ss >> Mr0 >> Mr1 >> dMr))
	{
		THROW(EAny, "Error reading Mr field from config file. Expect Mr = <Mr0> <Mr1> <dMr>, got Mr = " + tmp);
	}
	cfg.get(tmp,   "FeH_grid",   "-3 0.5 0.01");
	ss.clear(); ss.str(tmp);
	if(!(ss >> FeH0 >> FeH1 >> dFeH))
	{
		THROW(EAny, "Error reading FeH field from config file. Expect FeH = <FeH0> <FeH1> <dFeH>, got FeH = " + tmp);
	}

	cfg.get(tmp,   "file",   "");
	if(tmp == "") { // Try to find internal definitions for the photometric system
		tmp = datadir() + "/" + bandset + ".photosys.txt";
		if(!file_exists(tmp))
		{
			THROW(EAny, "Default photosys. definition file " + tmp + " for " + bandset + " not found. Specify it explicitly using the file=xxx keyword.");
		}
	}
	text_input_or_die(in, tmp);

	MLOG(verb1) << "Generating " << bandset << " photometry.";
	MLOG(verb1) << bandset << ": Generating " << bnames.size() << " bands: " << bnames;
	MLOG(verb1) << bandset << ": Using color(" << absbband << ", FeH) table from " << tmp << ".";
	MLOG(verb1) << bandset << ": Using " << bband << " to bootstrap appmags from colors.";
	MLOG(verb1) << bandset << ": Resampling color table to fast lookup grid:";
	MLOG(verb1) << bandset << ":    " << absbband << "0, " << absbband << "1, d(" << absbband << ") = " << Mr0 << ", " << Mr1 << ", " << dMr << ".";
	MLOG(verb1) << bandset << ":    FeH0, FeH1, dFeH = " << FeH0 << ", " << FeH1 << ", " << dFeH << ".";

	// load the data
	std::map<double, std::vector<Mr2col> > v;
	std::set<double> uFeH;
	Mr2col mc; double FeH;
	bind(in, mc.Mr,0, FeH,1);
	FOR(0, ncolors) { bind(in, mc.c[i], i+2); }
	while(in.next())
	{
		v[FeH].push_back(mc);
		uFeH.insert(FeH);
	}
	std::vector<double> vFeH(uFeH.begin(), uFeH.end());

	// construct col(Mr) splines for each FeH line present in the input
	std::map<double, colsplines> s;
	FOREACH(v)
	{
		double FeH = i->first;
		std::vector<Mr2col> &m2c = i->second;
		colsplines &ss = s[FeH];		// Splines that will return col(Mr) for this FeH

		std::sort(m2c.begin(), m2c.end());
		std::vector<double> Mr(m2c.size()), c(m2c.size());
		FOR(0, m2c.size()) { Mr[i] = m2c[i].Mr; }	// Construct (sorted) array of Mr
		FOR(0, ncolors) // for each color i, construct its col_i(Mr) spline
		{
			FORj(j, 0, m2c.size()) { c[j] = m2c[j].c[i]; }
			ss[i].construct(Mr, c);
//			std::cerr << FeH << " " << ss[i].f << "\n";
		}
		// remember beginning/end for test of extrapolation
		ss.Mrmin = Mr.front();
		ss.Mrmax = Mr.back();
	}

/*	FOREACH(s) {
		std::cerr << i->first << ": ";
		FORj(j,0,ncolors) {
			std::cerr << i->second[j].f << "  ";
		}
		std::cerr << "\n";
	}*/
	
	// allocate memory for output tables
	nMr  = (int)((Mr1 -Mr0) /dMr  + 1);
	nFeH = (int)((FeH1-FeH0)/dFeH + 1);
	clt.resize(ncolors);  FOREACH(clt)  { i->resize(nMr*nFeH); }
	eclt.resize(ncolors); FOREACH(eclt) { i->resize(nMr*nFeH); }

	// thread in Fe/H direction, constructing col(FeH) spline for each given Mr,
	// using knots derived from previously calculated col(Mr) splines for FeHs given in the input file
	std::vector<double> vcol(vFeH.size());	// col knots
	std::vector<double> vecol(vFeH.size());	// 1 if corresponding col knot was extrapolated, 0 otherwise
	FORj(ic, 0, ncolors)
	{
		FORj(m, 0, nMr)
		{
			double Mr = Mr0 + m*dMr;

			// construct col(FeH) splines for given Mr
			int k=0;
//			std::cerr << "Mr = " << Mr << "\n";
			FOREACH(s)
			{
				vcol[k]  = i->second[ic](Mr);
				vecol[k] = i->second.Mrmin > Mr || Mr > i->second.Mrmax;
//				std::cerr << i->first << " " << vcol[k] << " " << vecol[k] << "\n";
				k++;
			}
			spline s, es;
			s.construct(vFeH, vcol);
			es.construct(vFeH, vecol);	// zero where neither adjacent knot was extrapolated, nonzero otherwise

			FORj(f, 0, nFeH) // compute col(FeH, Mr)
			{
				double FeH = FeH0 + f*dFeH;
				int idx = m*nFeH + f;

				clt[ic][idx] = s(FeH);
				eclt[ic][idx] = vFeH.front() > FeH || FeH > vFeH.back() || es(FeH) != 0.;
//				std::cerr << Mr << " " << FeH << " : " << clt[ic][idx] << " " << eclt[ic][idx] << "\n";
			}
		}
	}

	std::vector<double> nextrap(ncolors);
	FOR(0, ncolors)
	{
		nextrap[i] = (double)count_if(eclt[i].begin(), eclt[i].end(), _1 != 0) / eclt[i].size();
	}
	MLOG(verb1) << bandset << ":    grid size = " << nMr << " x " << nFeH << " (" << clt[0].size() << ").";
	MLOG(verb1) << bandset << ":    extrapolation fractions = " << nextrap;

	return true;
}

size_t os_photometry::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	const size_t nbands = bnames.size();
	if(offset_mag == -1)
	{
		offset_mag    = sstruct::factory.getOffset(bband);
		offset_absmag = sstruct::factory.getOffset(absbband);
		offset_mags   = sstruct::factory.getOffset(bandset);
		offset_photoflags = sstruct::factory.getOffset(photoFlags);
	}

	// ASSUMPTIONS:
	//	- Fe/H exists in input
	//	- Apparent and absolute magnitude in the requested band exist in input
	for(size_t i=0; i != count; i++)
	{
		sstruct &s = in[i];

		// calculate magnitudes, absolute magnitudes and distance
		const float Mr   = s.get<float>(offset_absmag);
		const float bmag = s.get<float>(offset_mag);
		const float FeH  = s.FeH();

		// construct colors given the absolute magnitude and metallicity
		int &flags = *s.getptr<int>(offset_photoflags);
		flags = 0;
		bool ex;
		float c[ncolors];
		FOR(0, ncolors)
		{
			c[i] = color(i, FeH, Mr, &ex);
			flags |= ex << i;
		}
		// ug gr ri iz zy
		// convert colors to apparent magnitudes
		float *mags = s.getptr<float>(offset_mags);
		FORj(b, 0, nbands)
		{
			if(b < bidx) { FOR(b, bidx) { mags[b] += c[i]; } }
			if(b > bidx) { FOR(bidx, b) { mags[b] -= c[i]; } }
			mags[b] += bmag;
		}
#if 0
		FORj(b,0,ncolors)
		{
			std::cerr << mags[b]-mags[b+1] << " =?= " << c[b] << "\n";
		}
#endif
/*		mags[0] = r + ug + gr;
		mags[1] = r + gr;
		mags[2] = r;
		mags[3] = r - ri;
		mags[4] = r - ri - iz;*/
	}

	return nextlink->push(in, count, rng);
}

// convert input absolute/apparent magnitudes to ugriz colors
class os_ugriz : public osink
{
	protected:
		spline Mr2gi, gi2ri;
		std::vector<double> FeHDeltaC;

		void Mr2gi_build(const std::string &gi2Mr_str);
		void gi2ri_build();

		float MrFeH2gi(float Mr, float FeH);
		float grFeH2ug(float gr, float FeH);
		float ri2iz(float ri);
	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual const std::string &name() const { static std::string s("ugriz"); return s; }

		os_ugriz() : osink()
		{
			prov.insert("SDSSugriz[5]");
			req.insert("SDSSr");
			req.insert("absSDSSr");
			req.insert("FeH");
		}
};

float os_ugriz::grFeH2ug(float gr, float FeH)
{
	double y = gr;

	// TODO: we don't know how to set u band for g-r > 0.6
	// so if g-r > 0.6, just set it to 0.6 for now
	if(y > 0.6) {
		y = 0.6;
	}

	//// DR7 best-fit values
	// fitParabola23mixed:
	// ABCDEFGHI = -13.12526674 14.09473031 28.03659028
	//              -5.505543181 -5.902217618 -58.6842534
	//               9.135926568 -20.60655771 58.19755486
	//   differences are in the range -0.7295560609 to 0.4241155876
	//   mean = 7.66452796e-10, rms = 0.08717628066
	// this rounded up version returns metallicity
	// within 0.001 and with rms=0.062 min/max=-0.32/0.32
	// from the training sample
	const static double A = -13.13;
	const static double B =  14.09;
	const static double C =  28.04;
	const static double D =  -5.51;
	const static double E =  -5.90;
	const static double F = -58.68;
	const static double G =   9.14;
	const static double H = -20.61;
	const static double I =  58.20;

	const double y2 = y*y, y3 = y2*y;
	const double a = E + G * y;
	const double b = B + D * y + H * y2;
	const double c = A + C * y + F * y2 + I * y3 - FeH;
	double x0, x1;
	if(gsl_poly_solve_quadratic(a, b, c, &x0, &x1) == 0)
	{
	
		std::stringstream ss; ss << "Could not find roots for ug given FeH=" << FeH << " and gr=" << gr;
		MLOG(verb2) << ss.str();
		x0 = x1 = 1.;
		//THROW(EAny, ss.str());
	}

	bool x0q = 0.5 < x0 && x0 < 1.8;
	bool x1q = 0.5 < x1 && x1 < 1.8;
	if(x0q == x1q) {
		if(x0q) { DLOG(verb3) << "WARNING: Real u-g root ambiguous: x0,x1 = " << x0 << ", " << x1; }
		else    { DLOG(verb3) << "WARNING: Both roots probably wrong: x0,x1 = " << x0 << ", " << x1; }
		x0q = fabs(x0 - 1.05) < fabs(x1 - 1.05) ? true : false;	// return the one closer to 1.05
		x1q = !x0q;
	}
	return x0q ? x0 : x1;
#if 0
	// add Gaussian noise
	randomizeG Zraw 0 0.1 $3
	unset Zraw
#endif
}

void os_ugriz::Mr2gi_build(const std::string &gi2Mr_str)
{
	if(0) {
		// coefficients of Mr(g-i) relation
		// double coeff[] = {-0.56, 14.32, -12.97, 6.127, -1.267, 0.0967};
		std::string fn = datadir() + "/gi2Mr.FeH=0.txt";
		FILE_EXISTS_OR_THROW(fn);
		MLOG(verb2) << "Mr2gi locus: " << fn;

		text_input_or_die(f, fn);
		std::vector<double> gi, Mr;
		load(f, gi, 0, Mr, 1);
		Mr2gi.construct(Mr, gi);
	} else {
		std::vector<double> giarr, Mrarr, gi2Mr;
		split(gi2Mr, gi2Mr_str);
		for(double gi=-0.5; gi <= 5.; gi+=0.01)
		{
			double Mr = gsl_poly_eval(&gi2Mr[0], gi2Mr.size(), gi);
			Mrarr.push_back(Mr);
			giarr.push_back(gi);
		}
		Mr2gi.construct(Mrarr, giarr);

		MLOG(verb2) << "Constructed gi(Mr) using Mr(gi) coeffs = " << gi2Mr;
	}
}

float os_ugriz::MrFeH2gi(float Mr, float FeH)
{
	// offset due to metallicity
//	double FeHDelta = FeH*(-1.11 - 0.18*FeH);
	double FeHDelta = gsl_poly_eval(&FeHDeltaC[0], FeHDeltaC.size(), FeH);

	return Mr2gi(Mr - FeHDelta);
}
void os_ugriz::gi2ri_build()
{
	std::string fn = datadir() + "/gi2ri.spline.txt";
	FILE_EXISTS_OR_THROW(fn);
	MLOG(verb2) << "gi2ri locus: " << fn;

	text_input_or_die(f, fn);
	std::vector<double> gi, ri;
	load(f, gi, 0, ri, 1);
	gi2ri.construct(gi, ri);
}
float os_ugriz::ri2iz(float ri)
{
	return 0.53*ri;
}

size_t os_ugriz::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	// ASSUMPTIONS:
	//	- Fe/H exists in input
	//	- paralax returns the absolute magnitude in SDSS r band
	//	- all stars are main sequence
	for(size_t i=0; i != count; i++)
	{
		sstruct &s = in[i];

		// calculate magnitudes, absolute magnitudes and distance
		const double Mr = s.absSDSSr();
		const float FeH = s.FeH();

		// construct SDSS colors given the absolute magnitude and metallicity
		float gi = MrFeH2gi(Mr, FeH);
		float ri = gi2ri(gi);
		float gr = gi - ri;
		float ug = grFeH2ug(gr, FeH);
		float iz = ri2iz(ri);

		// convert colors to apparent magnitudes
		float r  = s.SDSSr();
		float *mags = s.SDSSugriz();
		mags[0] = ug + gr + r;
		mags[1] = gr + r;
		mags[2] = r;
		mags[3] = r - ri;
		mags[4] = r - ri - iz;
	}

	return nextlink->push(in, count, rng);
}

bool os_ugriz::init(const Config &cfg)
{
	std::string gi2Mr_str, FeH2dMr_str;
	cfg.get(gi2Mr_str,   "gi2Mr",   "-0.56 14.32 -12.97 6.127 -1.267 0.0967");
	cfg.get(FeH2dMr_str, "FeH2dMr", "0.0 -1.11 -0.18");

	split(FeHDeltaC, FeH2dMr_str);
	MLOG(verb2) << "DeltaMr(Fe/H) coeffs = " << FeHDeltaC;

	Mr2gi_build(gi2Mr_str);
	gi2ri_build();

	return true;
}


// convert velocities to proper motions
class os_vel2pm : public osink
{
public:
	int coordsys;
	static const int GAL = 0;
	static const int EQU = 1;

	float vLSR,		// Local standard of rest velocity
 	      u0, v0, w0;	// Solar peculiar motion

public:
	virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
	virtual bool init(const Config &cfg);
	virtual const std::string &name() const { static std::string s("vel2pm"); return s; }

	os_vel2pm() : osink(), coordsys(GAL)
	{
		req.insert("lb[2]");
		req.insert("XYZ[3]");
		req.insert("vcyl[3]");
	}
};

// convert cartesian galactocentric velocities to vl,vr,vb wrt. the observer
void vel_xyz2lbr(float &vl, float &vr, float &vb, const float vx, const float vy, const float vz, const double l, const double b)
{
	double cl, sl, cb, sb;
	cl = cos(l);
	sl = sin(l);
	cb = cos(b);
	sb = sin(b);

	float tmp;
	vl  =  vx*sl  - vy*cl;
	tmp =  vx*cl  + vy*sl;
	vr  = -cb*tmp + vz*sb;
	vb  =  sb*tmp + vz*cb;
}

void vel_cyl2xyz(float &vx, float &vy, float &vz, const float vr, const float vphi, const float vz0, const float X, const float Y)
{
	// convert galactocentric cylindrical to cartesian velocities
	float rho = sqrt(X*X + Y*Y);
	float cphi = X / rho;
	float sphi = Y / rho;

	vx = -sphi * vr + cphi * vphi;
	vy =  cphi * vr + sphi * vphi;
	vz = vz0;
}

template<typename T>
void array_copy(T *dest, const T *src, const size_t n)
{
	FOR(0, n) { dest[i] = src[i]; }
}

size_t os_vel2pm::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	// ASSUMPTIONS:
	//	vcyl() velocities are in km/s, XYZ() distances in parsecs
	//
	// OUTPUT:
	//	Proper motions in mas/yr for l,b directions in pm[0], pm[1]
	//	Radial velocity in km/s in pm[2]
	for(size_t i=0; i != count; i++)
	{
		sstruct &s = in[i];

		// fetch prerequisites
		const double *lb0 = s.lb(); double lb[2];
		lb[0] = rad(lb0[0]);
		lb[1] = rad(lb0[1]);
		const float *vcyl = s.vcyl();
		const float *XYZ = s.XYZ();

		// convert the velocities from cylindrical to galactocentric cartesian system
		float pm[3];
		vel_cyl2xyz(pm[0], pm[1], pm[2],   vcyl[0], vcyl[1], vcyl[2],   XYZ[0], XYZ[1]);

		// switch to Solar coordinate frame
		pm[0] -= u0;
		pm[1] -= v0 + vLSR;
		pm[2] -= w0;

		// convert to velocities wrt. the observer
		vel_xyz2lbr(pm[0], pm[1], pm[2],   pm[0], pm[1], pm[2],   lb[0], lb[1]);

		// convert to proper motions
		float D = sqrt(sqr(XYZ[0]) + sqr(XYZ[1]) + sqr(XYZ[2]));
		pm[0] /= 4.74 * D*1e-3;	// proper motion in mas/yr (4.74 km/s @ 1kpc is 1mas/yr)
		pm[1] /= 4.74 * D*1e-3;

		// rotate to output coordinate system
		switch(coordsys)
		{
		case GAL:
			array_copy(s.pmlb(), pm, 3);
			break;
		case EQU:
			THROW(EAny, "Output in equatorial system not implemented yet.");
			array_copy(s.pmradec(), pm, 3);
			break;
		default:
			THROW(EAny, "Unknown coordinate system [id=" + str(coordsys) + "] requested");
			break;
		}
	}

	return nextlink->push(in, count, rng);
}

bool os_vel2pm::init(const Config &cfg)
{
	//if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	std::string cs;
	cfg.get(cs, "coordsys", "gal");
	     if(cs == "gal") { coordsys = GAL; prov.insert("pmlb[3]"); }
	else if(cs == "equ") { coordsys = EQU; prov.insert("pmradec[3]"); }
	else { THROW(EAny, "Unknown coordinate system (" + cs + ") requested."); }

	// LSR and peculiar Solar motion
	cfg.get(vLSR, "vLSR", -220.0f);
	cfg.get(u0,   "u0",    -10.0f);
	cfg.get(v0,   "v0",     -5.3f);
	cfg.get(w0,   "w0",      7.2f);

	return true;
}


/////////////////////////////////////////////////////////////

// convert between coordinate systems
class os_gal2other : public osink
{
public:
	int coordsys;
	static const int GAL = 0;
	static const int EQU = 1;

public:
	virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
	virtual bool init(const Config &cfg);
	virtual const std::string &name() const { static std::string s("gal2other"); return s; }

	os_gal2other() : osink(), coordsys(GAL)
	{
		req.insert("lb[2]");
	}
};

size_t os_gal2other::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	// ASSUMPTIONS:
	//	vcyl() velocities are in km/s, XYZ() distances in parsecs
	//
	// OUTPUT:
	//	Proper motions in mas/yr for l,b directions in pm[0], pm[1]
	//	Radial velocity in km/s in pm[2]
	if(coordsys != GAL)
	{
		for(size_t i=0; i != count; i++)
		{
			sstruct &s = in[i];

			// fetch prerequisites
			const double *lb0 = s.lb(); double lb[2];
			lb[0] = rad(lb0[0]);
			lb[1] = rad(lb0[1]);

			// rotate to output coordinate system
			double *out;
			switch(coordsys)
			{
			case EQU:
				out = s.radec();
				galequ(lb[0], lb[1], out[0], out[1]);
				break;
			default:
				THROW(EAny, "Unknown coordinate system [id=" + str(coordsys) + "] requested");
				break;
			}

			// convert to degrees
			out[0] /= ctn::d2r;
			out[1] /= ctn::d2r;
		}
	}

	return nextlink->push(in, count, rng);
}

bool os_gal2other::init(const Config &cfg)
{
	//if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	std::string cs;
	cfg.get(cs, "coordsys", "gal");
	     if(cs == "gal") { coordsys = GAL; /* -- noop -- */ }
	else if(cs == "equ") { coordsys = EQU; prov.insert("radec[2]"); }
	else { THROW(EAny, "Unknown coordinate system (" + cs + ") requested."); }

	return true;
}


/////////////////////////////////////////////////////////////
#if 0
// mix in photometric errors
class os_photoErrors : public osink
{
public:
	struct photoerr_t
	{
		spline sigma;

		photoerr_t() {}

		float draw(const float mag, gsl_rng *rng)
		{
			double s = sigma(mag);
			float err = gsl_ran_gaussian(rng, s);
			return err;
		}
	};
	std::map<std::string, photoerr_t> photoerrs;	// band -> error definitions

public:
	std::map<size_t, photoerr_t *> photoerrsI;	// sstruct idx -> error definitions (optimization, this is populated on init)

public:
	virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
	virtual bool init(const Config &cfg);
	virtual const std::string &name() const { static std::string s("photoErrors"); return s; }

	os_photoErrors() : osink()
	{
		req.insert("lb[2]");
	}
};

size_t os_photoErrors::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	// ASSUMPTIONS:
	//	vcyl() velocities are in km/s, XYZ() distances in parsecs
	//
	// OUTPUT:
	//	Proper motions in mas/yr for l,b directions in pm[0], pm[1]
	//	Radial velocity in km/s in pm[2]
	for(size_t i=0; i != count; i++)
	{
		sstruct &s = in[i];

		// fetch prerequisites
		const double *lb0 = s.lb(); double lb[2];
		lb[0] = rad(lb0[0]);
		lb[1] = rad(lb0[1]);

		// rotate to output coordinate system
		double *out;
		switch(coordsys)
		{
		case EQU:
			out = s.radec();
			galequ(lb[0], lb[1], out[0], out[1]);
			break;
		default:
			THROW(EAny, "Unknown coordinate system [id=" + str(coordsys) + "] requested");
			break;
		}

		// convert to degrees
		out[0] /= ctn::d2r;
		out[1] /= ctn::d2r;
	}

	return nextlink->push(in, count, rng);
}

bool os_photoErrors::init(const Config &cfg)
{
	//if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	cfg.get(cs, "coordsys", "gal");

	return true;
}
#endif

// in/out ends of the chain
class os_textout : public osink
{
	protected:
//		std::ofstream out;
//		peyton::io::gzstream::ofstream out;
		flex_output out;

		bool headerWritten;
		ticker tick;

	protected:
		// map of field name -> formatter string, for output
		std::map<std::string, std::string> outputs;
		// map of field index -> formatter string (optimization)
		std::map<size_t, std::string> outputsI;

	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual int priority() { return 1000000; }	// ensure this stage has the least priority
		virtual const std::string &name() const { static std::string s("textout"); return s; }
		virtual const std::string &type() const { static std::string s("output"); return s; }

		os_textout() : osink(), headerWritten(false), tick(-1)
		{
		}
};

size_t os_textout::push(sstruct *&data, const size_t count, gsl_rng *rng)
{
	if(tick.step <= 0) { tick.open("Processing", 10000); }

	if(!headerWritten)
	{
		out.out() << "# " << sstruct::factory << "\n";
		headerWritten = true;
	}

	size_t cnt = 0;
	while(cnt < count && (out.out() << data[cnt] << "\n")) { cnt++; tick.tick(); }

	if(!out.out()) { THROW(EIOException, "Error outputing data"); }

	delete [] data;
	data = NULL;

	return count;
}

bool os_textout::init(const Config &cfg)
{
	if(!cfg.count("filename")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
//	out.open(cfg["filename"].c_str());
	out.open(cfg["filename"].c_str());

	// slurp up any output/formatting information
	// the formats are written in the configuration file as:
	//	format.<fieldname> = <formatter_fmt_string>
	std::set<std::string> keys;
	cfg.get_matching_keys(keys, "format\\.[a-zA-Z0-9_]+$");
	FOREACH(keys)
	{
		std::string name = i->substr(7);
		outputs[name] = cfg[*i];
	}

	return out.out();
}

class osource : public opipeline_stage
{
	public:
		virtual int priority() { return -1000000; } // ensure highest priority for this stage

	public:
		osource() : opipeline_stage()
		{
			prov.insert("_source");
		}
};

class os_textin : public osource
{
	protected:
//		std::ifstream in;
		flex_input in;

	public:
		virtual bool init(const Config &cfg);
		virtual int run(gsl_rng *rng);
		virtual const std::string &name() const { static std::string s("textin"); return s; }
		virtual const std::string &type() const { static std::string s("input"); return s; }

		os_textin() {};
};

bool os_textin::init(const Config &cfg)
{
	const char *fn = cfg["filename"].c_str();
	in.open(fn);
	if(!in.in()) { THROW(EFile, "Failed to open '" + (std::string)fn + "' for input."); }

	in.in() >> sstruct::factory;
	return in.in();
}

int os_textin::run(gsl_rng *rng)
{
	size_t block = 10000;

	size_t total = 0;
	do {
		sstruct *data = sstruct::create(block);
		size_t cnt = 0;
		while(cnt < block && in.in() >> data[cnt]) { cnt++; }
		total += nextlink->push(data, cnt, rng);
	} while(in.in());

	return total;
}

boost::shared_ptr<opipeline_stage> opipeline_stage::create(const std::string &name)
{
	boost::shared_ptr<opipeline_stage> s;

	if(name == "textin") { s.reset(new os_textin); }
	else if(name == "textout") { s.reset(new os_textout); }
	else if(name == "ugriz") { s.reset(new os_ugriz); }
	else if(name == "photometry") { s.reset(new os_photometry); }
	else if(name == "FeH") { s.reset(new os_FeH); }
	else if(name == "fixedFeH") { s.reset(new os_fixedFeH); }
	else if(name == "vel2pm") { s.reset(new os_vel2pm); }
	else if(name == "gal2other") { s.reset(new os_gal2other); }
//	else if(name == "photoErrors") { s.reset(new os_photoErrors); }
	else if(name == "kinTMIII") { s.reset(new os_kinTMIII); }
	else { THROW(EAny, "Module " + name + " unknown."); }

	ASSERT(name == s->name());

	return s;
}

class opipeline
{
	public:
		std::list<boost::shared_ptr<opipeline_stage> > stages;

	public:
		void add(boost::shared_ptr<opipeline_stage> pipe) { stages.push_back(pipe); }
		int run(gsl_rng *rng);
};

int opipeline::run(gsl_rng *rng)
{
	// construct the pipeline based on requirements and provisions
	std::set<std::string> haves;
	sstruct::factory.gettags(haves);

	// form a priority queue of requested stages. The priorities ensure that _output stage
	// will end up last (as well as allow some control of which stage executes first if
	// multiple stages have all prerequisites satisfied)
	std::set<std::pair<int, opipeline_stage *> > stages;
	FOREACH(this->stages)
	{
		stages.insert(std::make_pair((*i)->priority(), (*i).get()));
	}

	std::list<opipeline_stage *> pipeline;
	std::string which;
	while(!stages.empty())
	{
		std::ostringstream ss;
		FOREACH(haves) { ss << *i << " "; };
		DLOG(verb2) << "haves: " << ss.str();

		bool foundOne = false;
		FOREACH(stages)
		{
			opipeline_stage &s = *i->second;
			if(!s.satisfied_with(haves)) { continue; }

			// check for collisions
			if(s.provides_any_of(haves, which))
			{
				THROW(EAny, "Another module already provides " + which + ", that " + s.name() + " is trying to provide");
			}

			// append to pipeline
			pipeline.push_back(&s);
			haves.insert(s.provides().begin(), s.provides().end());

			// use tags which the stage will provide
			FOREACHj(j, s.provides())
			{
				if(j->at(0) == '_') { continue; }
				sstruct::factory.useTag(*j, true);
			}

			// erase from the list of pending stages
			stages.erase(i);
			foundOne = true;
			break;
		}
		if(!foundOne)
		{
			std::stringstream ss;
			std::string sep;
			FOREACH(stages)
			{
				opipeline_stage &s = *i->second;
				ss << sep << s.name();
				sep = ", ";
			}
			THROW(EAny, "Module(s) '" + ss.str() + "' require one or more fields that no other module (or input) provides.");
		}
	}

	// chain the constructed pipeline
	opipeline_stage *last, *source = NULL;
	std::stringstream ss;
	FOREACH(pipeline)
	{
		if(source == NULL) { last = source = *i; ss << last->name(); continue; }
		osink *next = dynamic_cast<osink*>(*i);
		ASSERT(next);
		last->chain(next);
		last = next;
		ss << " -> " << last->name();
	}
	MLOG(verb1) << "Postprocessing module chain: " << ss.str();

	return source->run(rng);
}

void observe_catalog(const std::string &conffn, const std::string &input, const std::string &output)
{
	Config cfg; cfg.load(conffn);

	int seed;
	std::string inmod, outmod;
	cfg.get(seed,	  "seed", 	  42);
	cfg.get(inmod,	  "input_module",  "textin");
	cfg.get(outmod,	  "output_module", "textout");

	opipeline pipe;

	std::set<std::string> keys;
	cfg.get_matching_keys(keys, "module\\[.+\\]");
	FOREACH(keys)
	{
		std::pair<std::string, std::string> moddef = cfg[*i];

		Config modcfg;
		if(!moddef.second.empty())
		{
			modcfg.load(moddef.second);
		}
		if(moddef.first ==  inmod) { modcfg.insert(make_pair("filename", input)); }	// some special handling for input/output modules
		if(moddef.first == outmod) { modcfg.insert(make_pair("filename", output)); }

		boost::shared_ptr<opipeline_stage> stage( opipeline_stage::create(moddef.first) );
		if(!stage.get()) { THROW(EAny, "Module " + moddef.first + " unknown or failed to load."); }

		if(!stage->init(modcfg)) { THROW(EAny, "Failed to initialize output pipeline stage '" + moddef.first + "'"); }

		pipe.add(stage);
	}

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, seed);

	int nstars = pipe.run(rng);
	MLOG(verb1) << "Observing pipeline generated " << nstars << " point sources.";

	gsl_rng_free(rng);
#if 0
	cfg.get(Ar,	"Ar", 		0.);	// extinction

	// load magnitude error map
					if(cfg.count("photo_error_map"))
			{
					std::string fn = cfg["photo_error_map"];
					std::ifstream iin(fn.c_str()); ASSERT(iin) { std::cerr << "Could not open magnitude error map " << fn << "\n"; }
					io::ibstream in(iin);

					if(!(in >> magerrs_cvm)) { ASSERT(in >> magerrs); }
					magerrs.initialize(magerrs_cvm);

					cfg.get(magerrs.use_median_beam_for_all, "magerr_only_use_median", false);
					if(magerrs.use_median_beam_for_all) { std::cerr << "Using median error only.\n"; }
}

					cfg.get(constant_photo_error, "constant_photo_error", 0.);
					cfg.get(paralax_dispersion, "paralax_dispersion", 0.);

	// load various flags
					int flag;
					cfg.get(flag, "apply_photo_errors", 0); 	if(flag) flags |= APPLY_PHOTO_ERRORS;
	
	// dump short status
					std::cerr << "constant_photo_error = " << constant_photo_error << "\n";
					std::cerr << "paralax_dispersion = " << paralax_dispersion << "\n";
					std::cerr << "magerrs.use_median_beam_for_all = " << magerrs.use_median_beam_for_all << "\n";
					std::cerr << "magerrs.empty() = " << magerrs.empty() << "\n";
					std::cerr << "Flags: "
					<< (flags & APPLY_PHOTO_ERRORS ? "APPLY_PHOTO_ERRORS" : "")
					<< "\n";
#endif
}


void postprocess_catalog(const std::string &conffn, const std::string &input, const std::string &output, std::vector<std::string> modules)
{
	Config cfg; cfg.load(conffn);

	int seed;
	std::string inmod, outmod;
	cfg.get(seed,	  "seed", 	  42);

	// merge-in any modules included from the config file via the 'modules' keyword
	// modules keyword may either contain the module name (in which case the config
	// will be read from the current config file's module.<module_name>.XXXX keys)
	// or a filename with module configuration. (*********DEPRECATED********)
	std::string name;
	std::ostringstream msg;
	if(cfg.count("modules"))
	{
		std::istringstream ss(cfg["modules"]);
		while(ss >> name)
		{
			msg << " " << name;
			modules.push_back(name);
		}
	}

	// merge in any modules with module.<module_name>.XXXX present and
	// module.<module_name>.enabled != 0. Configuration will be read from
	// module.<module_name>.XXXX keys.
	std::set<std::string> keys, smodules;
	cfg.get_matching_keys(keys, "module\\.[a-zA-Z0-9_}{]+\\..*$");
	FOREACH(keys)
	{
		const std::string &key = *i;
		name = key.substr(key.find('.')+1);
		name = name.substr(0, name.find('.'));
		//std::cerr << "key = " << key << " name=" << name << "\n";

		std::string enkey = "module." + name + ".enabled";
		if(cfg.count(enkey) && !cfg[enkey].vbool()) { continue; } // skip the disabled modules

		if(!smodules.count(name)) { msg << " " << name; }
		smodules.insert(name);
	}
	modules.insert(modules.end(), smodules.begin(), smodules.end());
	MLOG(verb2) << "Adding modules from config file:" << msg.str();

	// merge-in modules with options given in the config file
	opipeline pipe;

	FOREACH(modules)
	{
		const std::string &cffn = *i;

		Config modcfg;
		if(file_exists(cffn))
		{
			// load from filename
			modcfg.load(cffn);
			if(!modcfg.count("module")) { THROW(EAny, "Configuration file " + cffn + " does not specify the module name"); }
			name = modcfg["module"];
		}
		else
		{
			// load from subkeys
			cfg.get_subset(modcfg, "module." + cffn + ".", true);
			modcfg.insert(make_pair("module", cffn));

			// get module name, in case the user is using module.name{uniqident}.xxx syntax
			name = cffn.substr(0, cffn.find('{'));
		}

		boost::shared_ptr<opipeline_stage> stage( opipeline_stage::create(name) );
		if(!stage) { THROW(EAny, "Module " + name + " unknown or failed to load."); }

		if(stage->type() == "input")  { modcfg.insert(make_pair("filename", input)); }
		if(stage->type() == "output") { modcfg.insert(make_pair("filename", output)); }

		if(!stage->init(modcfg)) { THROW(EAny, "Failed to initialize output pipeline stage '" + name + "'"); }
		MLOG(verb2) << "Loaded " << name << " (module type: " << stage->type() << ")";

		pipe.add(stage);
	}

	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, seed);

	int nstars = pipe.run(rng);
	MLOG(verb1) << "Observing pipeline generated " << nstars << " point sources.";

	gsl_rng_free(rng);
}
#endif
