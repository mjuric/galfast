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

#include "gpc_cpp.h"

#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <sstream>
#include <iostream>
#include <fstream>

#include "simulate.h"
#include "projections.h"
#include "model.h"
#include "paralax.h"
#include "analysis.h"
//#include "dm.h"
#include "io.h"
#include "gpu.h"

#include "simulate_base.h"

#include <vector>
#include <map>
#include <string>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <astro/io/format.h>
#include <astro/system/fs.h>
#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

using namespace boost::lambda;
namespace ct = column_types;

#include "observe.h"

bool opipeline_stage::runtime_init(/*const std::list<opipeline_stage *> &pipeline, */otable &t)
{
	// test if otable has all the necessary prerequisites
	FOREACH(req)
	{
		if(!t.using_column(*i))
		{
			DLOG(verb2) << "Failed on: " << *i;
			return false;
		}
	}

	// use tags which the stage will provide
	FOREACH(prov)
	{
		if(i->at(0) == '_') { continue; }
		t.use_column(*i);	// just touch the column to initialize it
	}

	return true;
}

#if 0
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
#endif



#if 1
// add photometric errors information
class os_photometricErrors : public osink
{
protected:
	struct errdef
	{
		std::string trueBandset;
		std::string obsBandset;
		int bandIdx;
		const spline *sgma;	// spline giving gaussian sigma of errors given true magnitude

		errdef(const std::string &obsBandset_, const std::string &trueBandset_, int bandIdx_, const spline &bandErrors)
			: obsBandset(obsBandset_), trueBandset(trueBandset_), bandIdx(bandIdx_), sgma(&bandErrors) {}
		float sigma(float mag) { return (*sgma)(mag); }
	};

protected:
	std::map<std::string, std::map<std::string, spline> > availableErrors;
	std::vector<errdef> columnsToTransform;

	void addErrorCurve(const std::string &bandset, const std::string &band, const std::string &file);

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
	virtual bool runtime_init(otable &t);

	virtual const std::string &name() const { static std::string s("photometricErrors"); return s; }
	virtual int priority() { return PRIORITY_INSTRUMENT; }	// ensure this stage has the least priority

	os_photometricErrors() : osink()
	{
	}
};

size_t os_photometricErrors::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// mix-in gaussian error, with sigma drawn from preloaded splines
	FOREACH(columnsToTransform)
	{
		ct::cfloat::host_t magObs  = in.col<float>(i->obsBandset);
		ct::cfloat::host_t magTrue = in.col<float>(i->trueBandset);
		int bandIdx = i->bandIdx;

		for(size_t row=begin; row <= end; row++)
		{
			float mag = magTrue(row, bandIdx);
			magObs(row, bandIdx) = mag + rng.gaussian(i->sigma(mag));
		}
	}

	return nextlink->process(in, begin, end, rng);
}

bool os_photometricErrors::runtime_init(otable &t)
{
	// Search the configuration for all photometric tags that are defined.
	// Note that the priority of this module ensures it's run after any module
	// that may generate photometric information has already run.
	std::set<std::string> bandset;
	if(!t.get_used_columns_by_class(bandset, "magnitude"))
	{
		MLOG(verb1) << "WARNING: Not mixing in photometric errors, as no photometric information is being generated.";
		return true;
	}

	FOREACH(bandset)
	{
		if(availableErrors.count(*i) == 0 ) { continue; }
		std::map<std::string, spline> &errors = availableErrors[*i];

		otable::columndef &cdef = t.getColumn(*i);
		if(column_type_traits::get<float>() != cdef.type())
		{
			THROW(EAny, "Photometry errors module expects all photometric information to be stored as single-precision floats, and " + *i + " is not.");
		}

		const std::string &trueBandset = *i;	// e.g. obsSDSSugriz
		std::string obsBandset  = "obs" + *i;	// e.g. obsSDSSugriz

		std::set<std::string> bands;
		cdef.getFieldNames(bands);
		FOREACH(bands)
		{
			if(!errors.count(*i)) { continue; }			// don't have errors for this band
			spline &bandErrors = errors[*i];

			if(!t.using_column(obsBandset))
			{
				std::map<int, std::string> fieldNames;
				t.getColumn(trueBandset).getFieldNames(fieldNames);
				FOREACH(fieldNames)
				{
					fieldNames[i->first] = "obs" + i->second;
					//std::cerr << fieldNames[i->first];
				}
				//abort();
				t.use_column_by_cloning(obsBandset, trueBandset, &fieldNames);
			}

			int bandIdx = cdef.getFieldIndex(*i);
			columnsToTransform.push_back(errdef(obsBandset, trueBandset, bandIdx, bandErrors));

			MLOG(verb2) << "Adding photometric errors to " << trueBandset << "." << *i << " (output in " << obsBandset << "." << *i << ")";
		}
	}
	return true;
}

void os_photometricErrors::addErrorCurve(const std::string &bandset, const std::string &band, const std::string &file)
{
	text_input_or_die(in, file);
	std::vector<double> mag, sigma;
	load(in, mag, 0, sigma, 1);

	availableErrors[bandset][band].construct(mag, sigma);
	MLOG(verb2) << "Loaded photometric errors for " << bandset << ", " << band << " band";
}

bool os_photometricErrors::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// Expected configuration format:
	// 	<photosys>.<bandidx>.file = banderrs.txt
	// If instead of a filename the keyword 'internal' is specified, built-in files will be used.
	// Built-in filenames are of the form $datadir/<photosys>.<bandname>.photoerr.txt
	//
	// Example:
	//   SDSSugriz.SDSSu.file = SDSSugriz.SDSSu.photoerr.txt
	//
	// For convenience, this is also allowed:
	//
	//   SDSSugriz.file = internal
	//
	// where all files of the form SDSSugriz.*.photoerr.txt will be looked up and loaded.
	//
	// Expected banderrs.txt format:
	//   <mag>   <sigma(mag)>
	
	// convenience -- 'SDSSugriz.file = internal' slurps up anything with
	// SDSSugriz.*.photoerr.txt from data directory
	std::set<std::string> skeys;
	cfg.get_matching_keys(skeys, "[a-zA-Z0-9_]+\\.file");
	FOREACH(skeys)
	{
		size_t p1 = i->find('.');
		std::string bandset = i->substr(0, p1);

		std::string path = datadir() + "/" + bandset + ".";
		size_t pos = path.size();
		path += "*.photoerr.txt";
		peyton::io::dir dir(path);
		
		MLOG(verb2) << "Looking for error definition files for " << bandset << " (" << path << ")";

		FOREACH(dir)
		{
			std::string band = i->substr(pos, i->find('.', pos)-pos);
			addErrorCurve(bandset, band, *i);
		}
	}

	// parse per-band keys
	std::set<std::string> keys;
	cfg.get_matching_keys(keys, "[a-zA-Z0-9_]+\\.[a-zA-Z0-9_]+\\.file");
	FOREACH(keys)
	{
		size_t p1 = i->find('.');
		std::string bandset = i->substr(0, p1);
		p1++;
		std::string band = i->substr(p1, i->find('.', p1)-p1);

		std::string file = cfg[*i];
		if(file == "internal")
		{
			file = datadir() + "/" + bandset + "." + band + ".photoerr.txt";
		}

		addErrorCurve(bandset, band, file);
	}

	return true;
}
#endif


#if 1
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

	void draw(gsl_vector *y, rng_t &rng, bool zero = false)
	{
		gsl_vector_set(Z, 0, rng.gaussian(1.));
		gsl_vector_set(Z, 1, rng.gaussian(1.));
		gsl_vector_set(Z, 2, rng.gaussian(1.));

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
class os_kinTMIII_OLD : public osink
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
		void add_dispersion(double v[3], double Rsquared, double Z, dvec *ellip[6], rng_t &rng);
		void compute_means(double v[3], double Rsquared, double Z, dvec *means[3]);

		void get_disk_kinematics(double v[3], double Rsquared, double Z, rng_t &rng, bool &firstGaussian);
		void get_halo_kinematics(double v[3], double Rsquared, double Z, rng_t &rng);

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool init(const Config &cfg, otable &t);
		virtual const std::string &name() const { static std::string s("kinTMIII"); return s; }

		os_kinTMIII_OLD() : osink()
		{
			prov.insert("vcyl");
			req.insert("comp");
			req.insert("XYZ");
		}
};

#if 0
void test_kin()
{
	return;
// 	Radians l, b;
// 	equgal(0, ctn::pi/2., l, b);
// 	printf("%.10f\n", deg(l));
// 	exit(-1);

	os_kinTMIII_OLD o;

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
#endif

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

size_t os_kinTMIII_OLD::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input
	double tmp[3]; bool firstGaussian;
	ct::cint::host_t   comp = in.col<int>("comp");
	ct::cfloat::host_t XYZ  = in.col<float>("XYZ");
	ct::cfloat::host_t vcyl = in.col<float>("vcyl");

	// ASSUMPTIONS:
	//	- Fe/H exists in input
	//	- Apparent and absolute magnitude in the requested band exist in input
	for(size_t row=begin; row != end; row++)
	{
		// fetch prerequisites
		const int component = comp[row];
		float X = XYZ(row, 0);
		float Y = XYZ(row, 1);
		float Zpc = XYZ(row, 2);
		const double Rsquared = 1e-6 * (sqr(X) + sqr(Y));
		const double Z = 1e-3 * Zpc;

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
		vcyl(row, 0) = tmp[0];
		vcyl(row, 1) = tmp[1];
		vcyl(row, 2) = tmp[2];
	}

	return nextlink->process(in, begin, end, rng);
}



void os_kinTMIII_OLD::add_dispersion(double v[3], double Rsquared, double Z, dvec *ellip[6], rng_t &rng)
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

void os_kinTMIII_OLD::compute_means(double v[3], double Rsquared, double Z, dvec *means[3])
{
	// returns means in v[3]
	FOR(0, 3)
	{
		dvec &p = *means[i];
		v[i] = modfun(Rsquared, Z, p[0], p[1], p[2], p[3], p[4]);
	}
}

void os_kinTMIII_OLD::get_disk_kinematics(double v[3], double Rsquared, double Z, rng_t &rng, bool &firstGaussian)
{
	// set up which gaussian are we drawing from
	double p = rng.uniform();
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

void os_kinTMIII_OLD::get_halo_kinematics(double v[3], double Rsquared, double Z, rng_t &rng)
{
	compute_means(v, Rsquared, Z, haloMeans);
	add_dispersion(v, Rsquared, Z, haloEllip, rng);
}

template<typename T> inline OSTREAM(const std::vector<T> &v) { FOREACH(v) { out << *i << " "; }; return out; }

bool os_kinTMIII_OLD::init(const Config &cfg, otable &t)
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
	MLOG(verb2) << "Disk gaussian normalizations: " << fk << " : " << (1-fk);
	MLOG(verb2) << "Second disk gaussian offset:  " << DeltavPhi;

	MLOG(verb2) << "vR coefficients:              " << vR;
	MLOG(verb2) << "vZ coefficients:              " << vZ;
	MLOG(verb2) << "sigmaRR coefficients:         " << sigmaRR;
	MLOG(verb2) << "sigmaRPhi coefficients:       " << sigmaRPhi;
	MLOG(verb2) << "sigmaRZ coefficients:         " << sigmaRZ;
	MLOG(verb2) << "sigmaPhiPhi1 coefficients:    " << sigmaPhiPhi1;
	MLOG(verb2) << "sigmaPhiPhi2 coefficients:    " << sigmaPhiPhi2;
	MLOG(verb2) << "sigmaZPhi coefficients:       " << sigmaZPhi;
	MLOG(verb2) << "sigmaZZ coefficients:         " << sigmaZZ;

	MLOG(verb2) << "HvR coefficients:             " << HvR;
	MLOG(verb2) << "HvZ coefficients:             " << HvZ;
	MLOG(verb2) << "HsigmaRR coefficients:        " << HsigmaRR;
	MLOG(verb2) << "HsigmaRPhi coefficients:      " << HsigmaRPhi;
	MLOG(verb2) << "HsigmaRZ coefficients:        " << HsigmaRZ;
	MLOG(verb2) << "HsigmaPhiPhi coefficients:    " << HsigmaPhiPhi;
	MLOG(verb2) << "HsigmaZPhi coefficients:      " << HsigmaZPhi;
	MLOG(verb2) << "HsigmaZZ coefficients:        " << HsigmaZZ;

	return true;
}

#endif


// convert input absolute/apparent magnitudes to ugriz colors
class os_photometry : public osink
{
	protected:
		std::string bandset2;			// name of this filter set
		//std::string bband;			// band off which to bootstrap other bands, using color relations. Must be supplied by other modules.
		std::string absbband;			// Absolute magnitude band for which the datafile gives col(absmag,FeH) values. By default, it's equal to "abs$bband". Must be supplied by other modules.
		std::string photoFlagsName;		// Name of the photometric flags field
		int bidx;				// index of bootstrap band in bnames
//		size_t offset_absmag, offset_mag;	// sstruct offsets to bootstrap apparent and absolute magnitudes [input]
//		std::vector<ct::cfloat*> mags;		// sstruct offset to magnitudes in computed bands [output]
		size_t offset_photoflags;		// sstruct offset to photometric flags [outout]
		std::vector<std::string> bnames;	// band names (e.g., LSSTr, LSSTg, SDSSr, V, B, R, ...)
///		std::vector<std::vector<float> > clt;
		std::vector<xptrng::tptr<float> > isochrones;	// A rectangular, fine-grained, (Mr,FeH) -> colors map
		std::vector<xptrng::tptr<uint> > eflags;	// Flags noting if a pixel in an isochrone was extrapolated
/*		typedef char cbool;			// to avoid the special vector<bool> semantics, while maintaining a smaller memory footprint than vector<int>
		std::vector<std::vector<cbool> > eclt;	// extrapolation flags*/
		int nMr, nFeH;
		double Mr0, Mr1, dMr;
		double FeH0, FeH1, dFeH;
		int ncolors;

		float color(int ic, double FeH, double Mr, int *e = NULL)
		{
///			ASSERT(ic >= 0 && ic < clt.size()) { std::cerr << "ic = " << ic << "\nclt.size() = " << clt.size() << "\n"; }
			ASSERT(ic >= 0 && ic < isochrones.size()) { std::cerr << "ic = " << ic << "\nclt.size() = " << isochrones.size() << "\n"; }
			ASSERT(Mr0 <= Mr && Mr <= Mr1) { std::cerr << Mr0 << " <= " << Mr << " <= " << Mr1 << "\n"; }
			ASSERT(FeH0 <= FeH && FeH <= FeH1) { std::cerr << FeH0 << " <= " << FeH << " <= " << FeH1 << "\n"; }

			int f = (int)((FeH - FeH0) / dFeH);
			int m = (int)((Mr  -  Mr0) / dMr);
//			std::cerr << "fm = " << f << " " << m << "   " << ((FeH - FeH0) / dFeH) << " " << ((Mr  -  Mr0) / dMr) << "\n";
//			int idx = m*nFeH + f;
//			if(e) { *e = eclt[ic][idx]; }
			if(e) { *e = eflags[ic].elem(f, m); }
//			return clt[ic][idx];
			return isochrones[ic].elem(f, m);
		}
	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
		virtual bool runtime_init(otable &t);
		virtual const std::string &name() const { static std::string s("photometry"); return s; }

		os_photometry() : osink() /*, offset_absmag(-1), offset_mag(-1)*/
		{
			req.insert("FeH");
		}
};

bool os_photometry::runtime_init(otable &t)
{
	return osink::runtime_init(t);
}

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

bool os_photometry::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// load bandset name
	std::string tmp, bname;
	cfg.get(bandset2,   "bandset",   "LSSTugrizy");

	// load band names and construct field definition
	cfg.get(tmp,   "bands",   "LSSTu LSSTg LSSTr LSSTi LSSTz LSSTy");
	std::istringstream ss(tmp);
	std::ostringstream sbnames;
	while(ss >> bname)
	{
		if(bnames.size()) { sbnames << ","; }
		sbnames << bnames.size() << ":" << bname;

		bnames.push_back(bname);
	}

	const size_t nbands = bnames.size();
	std::string bfield = bandset2 + "[" + str(nbands) + "]{class=magnitude;fieldNames=" + sbnames.str() + ";}";
	prov.insert(bfield);

	// determine the number of colors
	ncolors = bnames.size()-1;
	if(ncolors >= Mr2col::maxcolors)
	{
		THROW(EAny, "This module cannot handle more than " + str(Mr2col::maxcolors+1) + " bands.");
	}

	// deduce photoFlags name
	std::string bandsetname = bandset2;
	photoFlagsName = bandsetname + "PhotoFlags";
	prov.insert(photoFlagsName + "{class=flags;}");

	// bootstap band setup
	std::string bband;
	cfg.get(bband,   "bootstrap_band",   "LSSTr");
	cfg.get(absbband,   "absband",   "abs"+bband);
	//req.insert(bband);
	req.insert("DM");
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
		tmp = datadir() + "/" + bandset2 + ".photosys.txt";
		if(!file_exists(tmp))
		{
			THROW(EAny, "Default photosys. definition file " + tmp + " for " + bandset2 + " not found. Specify it explicitly using the file=xxx keyword.");
		}
	}
	text_input_or_die(in, tmp);

	MLOG(verb1) << "Photometry: Generating " << bandset2 << " ( " << bnames << ")";
	MLOG(verb2) << bandset2 << ": Generating " << bnames.size() << " bands: " << bnames;
	MLOG(verb2) << bandset2 << ": Input absolute magnitude assumed to be in " << bband << " band.";
	MLOG(verb2) << bandset2 << ": Using color(" << absbband << ", FeH) table from " << tmp << ".";
	MLOG(verb2) << bandset2 << ": Resampling color table to fast lookup grid:";
	MLOG(verb2) << bandset2 << ":    " << absbband << "0, " << absbband << "1, d(" << absbband << ") = " << Mr0 << ", " << Mr1 << ", " << dMr << ".";
	MLOG(verb2) << bandset2 << ":    FeH0, FeH1, dFeH = " << FeH0 << ", " << FeH1 << ", " << dFeH << ".";

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

	// allocate memory for output tables
	nMr  = (int)((Mr1 -Mr0) /dMr  + 1);
	nFeH = (int)((FeH1-FeH0)/dFeH + 1);
///	clt.resize(ncolors);  FOREACH(clt)  { i->resize(nMr*nFeH); }
//	eclt.resize(ncolors); FOREACH(eclt) { i->resize(nMr*nFeH); }
	isochrones.resize(ncolors); FOREACH(isochrones)  { i->realloc(nFeH, nMr); }
	eflags.resize(ncolors);     FOREACH(eflags)      { i->realloc(nFeH, nMr); }

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

///				clt[ic][idx] = s(FeH);
				isochrones[ic].elem(f, m) = s(FeH);
//				eclt[ic][idx] = (vFeH.front() > FeH || FeH > vFeH.back() || es(FeH) != 0.) << ic;
				eflags[ic].elem(f, m) = (vFeH.front() > FeH || FeH > vFeH.back() || es(FeH) != 0.) << ic;
//				if(eflags[ic](f, m) && ic > 1) { std::cerr << ic << " " << Mr << " " << FeH << " : " << isochrones[ic](f, m) << " " << eflags[ic](f,m) << "\n"; }
			}
		}
	}

	std::vector<double> nextrap(ncolors);
	FOR(0, ncolors)
	{
//		nextrap[i] = (double)count_if(eclt[i].begin(), eclt[i].end(), _1 != 0) / eclt[i].size();
		nextrap[i] = 0;
		xptrng::hptr<uint> ptr = eflags[i];
/*		FOREACHj(cf, ptr)
		{
			if(*cf != 0) { nextrap[i] += 1; }
		}*/
		FORj(x, 0, nFeH)
		{
			FORj(y, 0, nMr)
			{
				if(ptr(x, y) != 0) { nextrap[i] += 1; }
			}
		}
		nextrap[i] /= nFeH*nMr;
	}
///	MLOG(verb1) << bandset2 << ":    grid size = " << nMr << " x " << nFeH << " (" << clt[0].size() << ").";
	MLOG(verb2) << bandset2 << ":    grid size = " << nFeH << " x " << nMr << " (" << isochrones[0].size() << ").";
	MLOG(verb2) << bandset2 << ":    extrapolation fractions = " << nextrap;

#if HAVE_CUDA
	// upload lookup table to CUDA textures
	
#endif

#if 0
	std::ofstream ff("dump.0.txt");
	for(double Mr = Mr0; Mr < Mr1; Mr += dMr*2.)
	{
		for(double FeH = FeH0; FeH < FeH1; FeH += dFeH*2.)
		{
			ff << Mr << " " << FeH << " " << color(1, FeH, Mr) << "\n";
		}
	}
#endif
	return true;
}


#if 0
size_t os_photometry::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	ct::cint  &flags  = in.col<int>(photoFlagsName);
	ct::cfloat &bmag  = in.col<float>(bband);
	ct::cfloat &Mr    = in.col<float>(absbband);
	ct::cfloat &mags  = in.col<float>(bandset2);
	ct::cfloat &FeH   = in.col<float>("FeH");
	const size_t nbands = bnames.size();

	// ASSUMPTIONS:
	//	- Fe/H exists in input
	//	- Apparent and absolute magnitude in the requested band exist in input
	for(size_t row=begin; row != end; row++)
	{
		// construct colors given the absolute magnitude and metallicity
		int ex;
		float c[ncolors];

		int f = 0;
		FOR(0, ncolors)
		{
			c[i] = color(i, FeH[row], Mr[row], &ex);
//			if(ex != 0) { std::cerr << i << " " << ex << "\n"; }
//			f |= ex << i;
			f |= ex;
		}
		flags[row] = f;

		// convert colors to apparent magnitudes
		FORj(b, 0, nbands)
		{
			float mag = bmag[row];
			if(b < bidx) { FOR(b, bidx) { mag += c[i]; } }
			if(b > bidx) { FOR(bidx, b) { mag -= c[i]; } }
			mags(row, b) = mag;
		}

/*		if(row == 30)
		{
			for(int i=0; i != ncolors; i++) { std::cerr << "c[" << i << "]=" << c[i] << "\n"; }
			for(int i=0; i != ncolors+1; i++) { std::cerr << "mags[" << i << "]=" << mags(row, i) << "\n"; }
			exit(0);
		}*/
	}

	return nextlink->process(in, begin, end, rng);
}
#else
typedef ct::cfloat::gpu_t gcfloat;
typedef ct::cint::gpu_t gcint;
DECLARE_KERNEL(os_photometry_kernel(otable_ks ks, os_photometry_data lt, gcint flags, gcfloat DM, gcfloat Mr, int nabsmag, gcfloat mags, gcfloat FeH));
void os_photometry_set_isochrones(const char *id, std::vector<xptrng::tptr<float> > *loc, std::vector<xptrng::tptr<uint> > *flgs);

size_t os_photometry::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	ct::cint  &flags  = in.col<int>(photoFlagsName);
	ct::cfloat &DM    = in.col<float>("DM");
	ct::cfloat &mags  = in.col<float>(bandset2);
	ct::cfloat &FeH   = in.col<float>("FeH");

	std::string absmagSys = absbband + "Sys";
	ct::cfloat &Mr    = in.using_column(absmagSys) ?
				in.col<float>(absmagSys) :
				in.col<float>(absbband);

	os_photometry_data lt = { ncolors, bidx, FeH0, dFeH, Mr0, dMr };
	os_photometry_set_isochrones(getUniqueId().c_str(), &isochrones, &eflags);
	CALL_KERNEL(os_photometry_kernel, otable_ks(begin, end, -1, sizeof(float)*ncolors), lt, flags, DM, Mr, Mr.width(), mags, FeH);

	return nextlink->process(in, begin, end, rng);
}
#endif

#if 0
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
		virtual bool init(const Config &cfg, otable &t);
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

bool os_ugriz::init(const Config &cfg, otable &t)
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
#endif


/////////////////////////////////////////////////////////////

// convert between coordinate systems
class os_gal2other : public osink
{
public:
	int coordsys;

public:
	size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("gal2other"); return s; }

	os_gal2other() : osink(), coordsys(GAL)
	{
		req.insert("lb");
	}
};

DECLARE_KERNEL(os_gal2other_kernel(otable_ks ks, int coordsys, ct::cdouble::gpu_t lb0, ct::cdouble::gpu_t out));

size_t os_gal2other::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	using namespace column_types;
	cdouble &lb   = in.col<double>("lb");

	if(coordsys == EQU)
	{
		cdouble &out = in.col<double>("radec");
		CALL_KERNEL(os_gal2other_kernel, otable_ks(begin, end), coordsys, lb, out);
	}

	return nextlink->process(in, begin, end, rng);
}

bool os_gal2other::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	//if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	std::string cs;
	cfg.get(cs, "coordsys", "gal");
	     if(cs == "gal") { coordsys = GAL; /* -- noop -- */ }
	else if(cs == "equ") { coordsys = EQU; prov.insert("radec[2]"); }
	else { THROW(EAny, "Unknown coordinate system (" + cs + ") requested."); }

	return true;
}

#if 0
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
	virtual bool init(const Config &cfg, otable &t);
	virtual const std::string &name() const { static std::string s("photoErrors"); return s; }

	os_photoErrors() : osink()
	{
		req.insert("lb");
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

bool os_photoErrors::init(const Config &cfg, otable &t)
{
	//if(!cfg.count("FeH")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	cfg.get(cs, "coordsys", "gal");

	return true;
}
#endif

#endif
// in/out ends of the chain
class os_textout : public osink
{
	protected:
		flex_output out;

		bool headerWritten;
		ticker tick;

#if 0 // not implemented yet
	protected:
		// map of field name -> formatter string, for output
		std::map<std::string, std::string> outputs;
		// map of field index -> formatter string (optimization)
		std::map<size_t, std::string> outputsI;
#endif
	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
		virtual int priority() { return PRIORITY_OUTPUT; }	// ensure this stage has the least priority
		virtual const std::string &name() const { static std::string s("textout"); return s; }
		virtual const std::string &type() const { static std::string s("output"); return s; }

		os_textout() : osink(), headerWritten(false), tick(-1)
		{
		}
};


struct mask_output : otable::mask_functor
{
	column_types::cint::host_t hidden;
	ticker &tick;
	mask_output(column_types::cint::host_t &h, ticker &tck) : hidden(h), tick(tck) {}

	virtual bool shouldOutput(int row) const
	{
		tick.tick();
		return !hidden[row];
	}
};

size_t os_textout::process(otable &t, size_t from, size_t to, rng_t &rng)
{
//	if(tick.step <= 0) { tick.open("Writing output", 10000); }
	ticker tick("Writing output", (int)ceil((to-from)/50.));

	if(!headerWritten)
	{
		out.out() << "# ";
		t.serialize_header(out.out());
		out.out() << "\n";
		headerWritten = true;
	}

	swatch.start();

	size_t nserialized = 0;
#if 1
	if(t.using_column("hidden"))
	{
		column_types::cint::host_t   hidden = t.col<int>("hidden");
		nserialized = t.serialize_body(out.out(), from, to, mask_output(hidden, tick));
	}
	else
	{
		nserialized = t.serialize_body(out.out(), from, to);
	}
#endif
	swatch.stop();
	static bool firstTime = true; if(firstTime) { swatch.reset(); kernelRunSwatch.reset(); firstTime = false; }

	if(!out.out()) { THROW(EIOException, "Error outputing data"); }

// 	// quick hack to limit the number of stars that can be generated
// 	static int ntotal = 0;
// 	ntotal += nserialized;
// 	const static int nmax = 100*1000*1000;
// 	if(ntotal > nmax)
// 	{
// 		THROW(EAny, "Output currently limited to not (much) more than " + str(ntotal) + " stars.");
// 	}

	return nserialized;
}

bool os_textout::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	if(!cfg.count("filename")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	out.open(cfg["filename"].c_str());

#if 0 // not implemented
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
#endif
	return out.out();
}

class os_textin : public osource
{
	protected:
		flex_input in;

	public:
		virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
		virtual bool runtime_init(otable &t);
		virtual size_t run(otable &t, rng_t &rng);
		virtual const std::string &name() const { static std::string s("textin"); return s; }
		virtual const std::string &type() const { static std::string s("input"); return s; }

		os_textin() {};
};

bool os_textin::runtime_init(otable &t)
{
	// this unserializes the header and fills the prov vector with columns this module will provide
	t.unserialize_header(in.in(), &prov);

	osource::runtime_init(t);
}

bool os_textin::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	const char *fn = cfg["filename"].c_str();
	in.open(fn);
	if(!in.in()) { THROW(EFile, "Failed to open '" + (std::string)fn + "' for input."); }

	return in.in();
}

size_t os_textin::run(otable &t, rng_t &rng)
{
	size_t total = 0;
	do {
		swatch.start();
		t.clear();
		t.unserialize_body(in.in());
		swatch.stop();
		if(t.size() > 0)
		{
			static bool firstTime = true; if(firstTime) { swatch.reset(); kernelRunSwatch.reset(); firstTime = false; }
			total += nextlink->process(t, 0, t.size(), rng);
		}
	} while(in.in());

	return total;
}

boost::shared_ptr<opipeline_stage> opipeline_stage::create(const std::string &name)
{
	boost::shared_ptr<opipeline_stage> s;

	if(name == "textin") { s.reset(new os_textin); }
	else if(name == "skygen") { s.reset(new os_skygen); }
	else if(name == "textout") { s.reset(new os_textout); }
	else if(name == "unresolvedMultiples") { s.reset(new os_unresolvedMultiples); }
	else if(name == "FeH") { s.reset(new os_FeH); }
	else if(name == "fixedFeH") { s.reset(new os_fixedFeH); }
	else if(name == "photometry") { s.reset(new os_photometry); }
	else if(name == "photometricErrors") { s.reset(new os_photometricErrors); }
	else if(name == "clipper") { s.reset(new os_clipper); }
#if 0
	else if(name == "ugriz") { s.reset(new os_ugriz); }
#endif
	else if(name == "vel2pm") { s.reset(new os_vel2pm); }
	else if(name == "gal2other") { s.reset(new os_gal2other); }
//	else if(name == "photoErrors") { s.reset(new os_photoErrors); }
	else if(name == "kinTMIII") { s.reset(new os_kinTMIII); }
	else { THROW(EAny, "Module " + name + " unknown."); }

	ASSERT(name == s->name());

	return s;
}

// construct the pipeline based on requirements and provisions
size_t opipeline::run(otable &t, rng_t &rng)
{
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
		// find next pipeline stage that is satisfied with the available tags
		bool foundOne = false;
		FOREACH(stages)
		{
			opipeline_stage &s = *i->second;

			// initialize this pipeline stage (this typically adds and uses the columns
			// this stage will add)
			if(!s.runtime_init(t)) { continue; }

			// append to pipeline
			pipeline.push_back(&s);

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
		ss << " | " << last->name();
	}
	MLOG(verb1) << "Pipeline: " << ss.str();

#if 0
	//
	bool dryRun = true;
	if(dryRun)
	{
		std::cout << "###GENERATED_COLUMNS: ";
		t.serialize_header(std::cout);
		std::cout << "\n";
		return 0;
	}
#endif

	int ret = source->run(t, rng);

	DLOG(verb2) << "Postprocessing times:";
	FOREACH(pipeline)
	{
		DLOG(verb2) << io::format("  %10s: %f") << (*i)->name() << (*i)->getProcessingTime();
	}
	DLOG(verb2) << "Pure kernel run time: " << kernelRunSwatch.getTime();

	return ret;
}

bool opipeline::has_module_of_type(const std::string &type) const
{
	FOREACH(stages)
	{
		if((*i)->type() == type) { return true; }
	}
	return false;
}

void postprocess_catalog(const std::string &conffn, const std::string &input, const std::string &output, std::set<std::string> modules)
{
	Config cfg; cfg.load(conffn);

	int seed;
	std::string inmod, outmod;
	cfg.get(seed,	  "seed", 	  42);
	rng_gsl_t rng(seed);

	// output table
//	static const size_t Kbatch = 99999;
//	static const size_t Kbatch = 500000;
//	static const size_t Kbatch = 2500000 / 2;
//	static const size_t Kbatch = 735000;
	static const size_t Kbatch = 5000000;
	DLOG(verb1) << "Postprocessing in batches of " << Kbatch << " objects";
	otable t(Kbatch);

	std::string name;
	std::ostringstream msg;

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
	modules.insert(smodules.begin(), smodules.end());
	//DLOG(verb2) << "Adding modules from config file:" << msg.str();

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

		stage->setUniqueId(modcfg["module"]);

		if(stage->type() == "input")  { modcfg.insert(make_pair("filename", input)); }
		if(stage->type() == "output") { modcfg.insert(make_pair("filename", output)); }

		if(!stage->construct(modcfg, t, pipe)) { THROW(EAny, "Failed to initialize output pipeline stage '" + name + "'"); }
		DLOG(verb2) << "postprocessing module loaded: " << name << " (type: " << stage->type() << ")";

		pipe.add(stage);
	}

	// set default I/O, if not overriden by othe modules
	if(!pipe.has_module_of_type("input"))
	{
		name = "textin";
		Config modcfg;
		boost::shared_ptr<opipeline_stage> stage( opipeline_stage::create(name) );
		modcfg.insert(make_pair("filename", input));
		if(!stage->construct(modcfg, t, pipe)) { THROW(EAny, "Failed to initialize output pipeline stage '" + name + "'"); }
		DLOG(verb2) << "postprocessing module loaded: " << name << " (type: " << stage->type() << ")";
		pipe.add(stage);
	}
	if(!pipe.has_module_of_type("output"))
	{
		name = "textout";
		Config modcfg;
		boost::shared_ptr<opipeline_stage> stage( opipeline_stage::create(name) );
		modcfg.insert(make_pair("filename", output));
		if(!stage->construct(modcfg, t, pipe)) { THROW(EAny, "Failed to initialize output pipeline stage '" + name + "'"); }
		DLOG(verb2) << "postprocessing module loaded: " << name << " (type: " << stage->type() << ")";
		pipe.add(stage);
	}

	int nstars = pipe.run(t, rng);
	MLOG(verb2) << "Observing pipeline generated " << nstars << " point sources.";
}
