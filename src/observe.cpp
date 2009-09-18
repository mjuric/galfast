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

#include "observe.h"
#include "spline.h"
#include "analysis.h"
#include "io.h"
#include "gpu.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <astro/io/format.h>
#include <astro/system/log.h>
#include <astro/useall.h>

#include <fstream>

bool opipeline_stage::runtime_init(otable &t)
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
	void addErrorCurve(const std::string &bandset, const std::string &band, const std::vector<double> &mag, const std::vector<double> &sigma);

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
		cfloat_t::host_t magObs  = in.col<float>(i->obsBandset);
		cfloat_t::host_t magTrue = in.col<float>(i->trueBandset);
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

void os_photometricErrors::addErrorCurve(const std::string &bandset, const std::string &band, const std::vector<double> &mag, const std::vector<double> &sigma)
{
	availableErrors[bandset][band].construct(mag, sigma);
	MLOG(verb2) << "Acquired photometric errors for " << bandset << ", " << band << " band";
}

void os_photometricErrors::addErrorCurve(const std::string &bandset, const std::string &band, const std::string &file)
{
	text_input_or_die(in, file);
	std::vector<double> mag, sigma;
	load(in, mag, 0, sigma, 1);

	addErrorCurve(bandset, band, mag, sigma);
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

class os_modelPhotoErrors : public os_photometricErrors
{
protected:
	struct model_t
	{
		virtual double eval(double m) = 0;
		virtual std::istream &load(std::istream &in) = 0;
	};
	struct lsst_t : public model_t
	{
		double PHsys, m5, fGamma, Nobs;

		double eval(double m)
		{
			double x = pow10(0.4*(m-m5));
			double PHrand = sqrt((0.04-fGamma)*x + fGamma*sqr(x)) / sqrt(Nobs);
			double e = sqrt(sqr(PHsys) + sqr(PHrand));
			return e;
		}
		std::istream &load(std::istream &in)
		{
			return in >> PHsys >> m5 >> fGamma >> Nobs;
		}
	};
public:
	virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);

	virtual const std::string &name() const { static std::string s("modelPhotoErrors"); return s; }

	os_modelPhotoErrors() : os_photometricErrors()
	{
	}
};

bool os_modelPhotoErrors::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// Computes photometric errors according to:
	//
	//         ## LSST errors (eqs.4-6 from astroph/0805.2366)
	//              x = 10^(0.4*(m-m5))
	//         PHrand = sqrt((0.04-fGamma)*x + fGamma*x^2) / sqrt(Nobs)
	//            err = sqrt(PHsys^2 + PHrand^2)
	//
	//	Input parameters:
	//		<PHsys>	-- systematic error
	//		<m5>	-- 5 sigma depth for a point source
	//		<fGamma>-- random error behavior term
	//		<Nobs>	-- number of observations
	//
	// Expected configuration format:
	// 	<photosys>.<bandidx> = lsst <PHsys> <m5> <fGamma> <Nobs>
	//

	std::vector<double> mag, err;
	for(double m = 0; m < 40; m += 0.01)
	{
		mag.push_back(m);
	}
	err.resize(mag.size());

	std::set<std::string> skeys;
	cfg.get_matching_keys(skeys, "[a-zA-Z0-9_]+\\.[a-zA-Z0-9_]+");
	FOREACH(skeys)
	{
		// parse the photometric system and band name
		size_t pos = i->find('.');
		std::string bandset = i->substr(0, pos);
		std::string band = i->substr(pos+1);

		// parse the error model parameters
		std::istringstream ss(cfg[*i]);

		std::string model;
		ss >> model;
		std::auto_ptr<model_t> m;
		if(model == "lsst")	{ m.reset(new lsst_t); }
		else 			{ THROW(EAny, "Unknown photometric error model '" + model + "'"); }

		if(!m->load(ss))
		{
			MLOG(verb1) <<  "Problem loading photometric error model parameters. Problematic input: " << cfg[*i];
			THROW(EAny, "Failed to load photometric error model parameters. Aborting.");
		}

		// sample the error curve
		FOR(0, mag.size())
		{
			err[i] = m->eval(mag[i]);
		}

		addErrorCurve(bandset, band, mag, err);
	}
	return true;
}





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

inline double modfun(double Rsquared, double Z, double a, double b, double c, double d, double e)
{
	return a + b*pow(fabs(Z), c) + d*pow(Rsquared, 0.5*e);
}

template<typename T> inline OSTREAM(const std::vector<T> &v) { FOREACH(v) { out << *i << " "; }; return out; }


// convert input absolute/apparent magnitudes to ugriz colors
class os_photometry : public osink, public os_photometry_data
{
	protected:
		std::string bandset2;			// name of this filter set
		std::string absbband;			// Absolute magnitude band for which the datafile gives col(absmag,FeH) values. By default, it's equal to "abs$bband". Must be supplied by other modules.
		std::string photoFlagsName;		// Name of the photometric flags field
		std::vector<std::string> bnames;	// band names (e.g., LSSTr, LSSTg, SDSSr, V, B, R, ...)
		std::vector<cuxSmartPtr<float> > isochrones;	// A rectangular, fine-grained, (Mr,FeH) -> colors map
		std::vector<cuxSmartPtr<uint> > eflags;	// Flags noting if a pixel in an isochrone was extrapolated
		std::pair<cuxTexture<float, 3>, cuxTexture<float, 3> > extinction;	// north/south extinction maps (for bootstrap band)

 		int nMr, nFeH;
		float Mr1, FeH1;

		float color(int ic, double FeH, double Mr, int *e = NULL)
		{
			ASSERT(ic >= 0 && ic < isochrones.size()) { std::cerr << "ic = " << ic << "\nclt.size() = " << isochrones.size() << "\n"; }
			ASSERT(Mr0 <= Mr && Mr <= Mr1) { std::cerr << Mr0 << " <= " << Mr << " <= " << Mr1 << "\n"; }
			ASSERT(FeH0 <= FeH && FeH <= FeH1) { std::cerr << FeH0 << " <= " << FeH << " <= " << FeH1 << "\n"; }

			int f = (int)((FeH - FeH0) / dFeH);
			int m = (int)((Mr  -  Mr0) / dMr);
//			std::cerr << "fm = " << f << " " << m << "   " << ((FeH - FeH0) / dFeH) << " " << ((Mr  -  Mr0) / dMr) << "\n";
//			int idx = m*nFeH + f;
//			if(e) { *e = eclt[ic][idx]; }
			if(e) { *e = eflags[ic](f, m); }
//			return clt[ic][idx];
			return isochrones[ic](f, m);
		}
	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
		virtual bool runtime_init(otable &t);
		virtual const std::string &name() const { static std::string s("photometry"); return s; }

		os_photometry() : osink()
		{
			comp0 = 0;
			comp1 = 0xffffffff;
			FOR(0, N_REDDENING) { reddening[i] = 1.f; }

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
	// load component IDs to which this module will apply
	cfg.get(comp0,   "comp0",   0U);
	cfg.get(comp1,   "comp1",   0xffffffff);

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

	// load reddening coefficients, if given
	cfg.get(tmp, "reddening", "");
	ss.clear(); ss.str(tmp);
	FOR(0, nbands)
	{
		if(!(ss >> reddening[i])) { break; }
	}

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
	isochrones.resize(ncolors); FOREACH(isochrones)  { *i = cuxSmartPtr<float>(nFeH, nMr); }
	eflags.resize(ncolors);     FOREACH(eflags)      { *i = cuxSmartPtr<uint>(nFeH, nMr); }

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

//				clt[ic][idx] = s(FeH);
				isochrones[ic](f, m) = s(FeH);
//				eclt[ic][idx] = (vFeH.front() > FeH || FeH > vFeH.back() || es(FeH) != 0.) << ic;
				eflags[ic](f, m) = (vFeH.front() > FeH || FeH > vFeH.back() || es(FeH) != 0.) << ic;
//				if(eflags[ic](f, m) && ic > 1) { std::cerr << ic << " " << Mr << " " << FeH << " : " << isochrones[ic](f, m) << " " << eflags[ic](f,m) << "\n"; }
			}
		}
	}

	std::vector<double> nextrap(ncolors);
	FOR(0, ncolors)
	{
		nextrap[i] = 0;
		hptr<uint> ptr = eflags[i];
		FORj(x, 0, nFeH)
		{
			FORj(y, 0, nMr)
			{
				if(ptr(x, y) != 0) { nextrap[i] += 1; }
			}
		}
		nextrap[i] /= nFeH*nMr;
	}
	MLOG(verb2) << bandset2 << ":    grid size = " << nFeH << " x " << nMr << " (" << isochrones[0].size() << ").";
	MLOG(verb2) << bandset2 << ":    extrapolation fractions = " << nextrap;

	return true;
}


typedef cfloat_t::gpu_t gcfloat;
typedef cint_t::gpu_t gcint;
DECLARE_KERNEL(os_photometry_kernel(otable_ks ks, gcfloat Am, gcint flags, gcfloat DM, gcfloat Mr, int nabsmag, gcfloat mags, gcfloat FeH, gcint comp));
void os_photometry_set_isochrones(const char *id, std::vector<cuxSmartPtr<float> > *loc, std::vector<cuxSmartPtr<uint> > *flgs);
void os_photometry_cleanup_isochrones(const char *id, std::vector<cuxSmartPtr<float> > *loc, std::vector<cuxSmartPtr<uint> > *flgs);

size_t os_photometry::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	cint_t &comp     = in.col<int>("comp");
	cint_t &flags    = in.col<int>(photoFlagsName);
	cfloat_t &DM     = in.col<float>("DM");
	cfloat_t &mags   = in.col<float>(bandset2);
	cfloat_t &FeH    = in.col<float>("FeH");
	cfloat_t &Am     = in.col<float>("Am");

	std::string absmagSys = absbband + "Sys";
	cfloat_t &Mr    = in.using_column(absmagSys) ?
				in.col<float>(absmagSys) :
				in.col<float>(absbband);

/*	os_photometry_data lt = { ncolors, bidx, FeH0, dFeH, Mr0, dMr, comp0, comp1 };*/
	cuxUploadConst("os_photometry_params", static_cast<os_photometry_data&>(*this));
	os_photometry_set_isochrones(getUniqueId().c_str(), &isochrones, &eflags);

	CALL_KERNEL(os_photometry_kernel, otable_ks(begin, end, -1, sizeof(float)*ncolors), Am, flags, DM, Mr, Mr.width(), mags, FeH, comp);

	os_photometry_cleanup_isochrones(getUniqueId().c_str(), &isochrones, &eflags);

	return nextlink->process(in, begin, end, rng);
}


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

DECLARE_KERNEL(os_gal2other_kernel(otable_ks ks, int coordsys, cdouble_t::gpu_t lb0, cdouble_t::gpu_t out));

size_t os_gal2other::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	cdouble_t &lb   = in.col<double>("lb");

	if(coordsys == EQU)
	{
		cdouble_t &out = in.col<double>("radec");
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
	cint_t::host_t hidden;
	ticker &tick;
	mask_output(cint_t::host_t &h, ticker &tck) : hidden(h), tick(tck) {}

	virtual bool shouldOutput(int row) const
	{
		tick.tick();
		return !hidden(row);
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
		cint_t::host_t   hidden = t.col<int>("hidden");
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
	const char *fn = cfg.count("filename") ? cfg["filename"].c_str() : "sky.obs.txt";
	out.open(fn);
	MLOG(verb1) << "Output file: " << fn << " (text)\n";

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


/////////////////////////////

#if HAVE_LIBCFITSIO
#include "fitsio2.h"

// in/out ends of the chain
class os_fitsout : public osink
{
	public:
		struct coldef
		{
			char *data;
			int width;
			int elementSize;
			size_t pitch;
		};

	protected:
		fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
		std::vector<iteratorCol> data;
		std::vector<coldef> columns;

		bool headerWritten;
		std::string header_def;
		ticker tick;

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
		virtual int priority() { return PRIORITY_OUTPUT; }	// ensure this stage has the least priority
		virtual const std::string &name() const { static std::string s("fitsout"); return s; }
		virtual const std::string &type() const { static std::string s("output"); return s; }

		os_fitsout() : osink(), headerWritten(false), tick(-1), fptr(NULL)
		{
		}
		~os_fitsout();
};

struct write_fits_rows_state
{
	cint_t::host_t hidden;
	os_fitsout::coldef *columns;
	int row, to, rowswritten;

	write_fits_rows_state(os_fitsout::coldef *columns_, int from_, int to_)
	{
		hidden.reset(); // make_hptr2D<int>(NULL, 0);
		rowswritten = 0;

		columns = columns_;
		row = from_;
		to = to_;
	}
};

int write_fits_rows(long totaln, long offset, long firstn, long nvalues, int narrays, iteratorCol *data,  void *userPointer)
{
	write_fits_rows_state &st = *(write_fits_rows_state *)userPointer;
	os_fitsout::coldef *columns = st.columns;
	firstn--;		// adjust to 0-based convention
	firstn -= offset;	// refer to the first row of the otable we're dumping

	// for each column...
	int orow = st.row;
	FOR(0, narrays)
	{
		os_fitsout::coldef &c = columns[i];
		char *f = (char *)fits_iter_get_array(data+i);
		ASSERT(c.width == data[i].repeat);

		// tell cfitsio we have no NULLs
		memset(f, 0, c.elementSize);
		f += c.elementSize;

		// for each row...
		orow = st.row;
		FORj(row, 0, nvalues)
		{
			if(st.hidden)
			{
				while(orow != st.to && st.hidden(orow)) { orow++; }	// find next non-hidden row
			}
			if(orow == st.to) { break; }				// have we reached the end of the table?

			char *from = c.data + c.elementSize * orow;	// 0th element in this row
			char *to = f + c.width*c.elementSize*row;
			// for each vector element in row...
			FORj(elem, 0, c.width)
			{
				char *elemto = to + c.elementSize*elem;
				char *elemfrom = from + c.pitch*elem;
				memcpy(
					elemto,
					elemfrom,
					c.elementSize
				);
#if 0
				memset(
					elemto,
					'A'+i,
					c.elementSize
				);
				elemto[0] = '0'+(row%10);
#endif
#if 0
				switch(data[i].datatype)
				{
					case TFLOAT:
						std::cerr << *(float*)elemfrom << " ";
						break;
					case TDOUBLE:
						std::cerr << *(double*)elemfrom << " ";
						break;
					case TINT:
						std::cerr << *(int*)elemfrom << " ";
						break;
				}
#endif
			}

			if(i == 0) { st.rowswritten++; }
			orow++;
		}
	}
	st.row = orow;
	return 0;
}

size_t os_fitsout::process(otable &t, size_t from, size_t to, rng_t &rng)
{
	ticker tick("Writing output", (int)ceil((to-from)/50.));

	// fetch columns we're going to write
	std::vector<const otable::columndef *> cols;
	t.getSortedColumnsForOutput(cols);
	const int tfields = cols.size();

	if(!headerWritten)
	{
		// collect header metadata
		std::ostringstream hdr;
		t.serialize_header(hdr);
		header_def = hdr.str();

		// create the output table
		char *ttype[tfields], *tform[tfields];
		data.resize(tfields);
		FOR(0, tfields)
		{
			coldef c;
			(const_cast<otable::columndef *>(cols[i]))->rawdataptr(c.elementSize, c.width, c.pitch);

			ttype[i] = strdup(cols[i]->getPrimaryName().c_str());
	 		asprintf(&tform[i], "%d%c", c.width, cols[i]->type()->fits_tform());

			//std::cerr << ttype[i] << " " << tform[i] << "\n";
		}

		int status = 0;
		fits_create_tbl(fptr, BINARY_TBL, 0, tfields, ttype, tform, NULL, "CATALOG", &status);
		ASSERT(status == 0) { fits_report_error(stderr, status); }

		// construct array for cfitsio/Iterator routines
		columns.resize(tfields);
		FOR(0, tfields)
		{
			int dtype;
			switch(cols[i]->type()->fits_tform())
			{
				case 'A': dtype = TSTRING; break;
				case 'J': dtype = TINT; break;
				case 'E': dtype = TFLOAT; break;
				case 'D': dtype = TDOUBLE; break;
				default: ASSERT(0);
			}

			fits_iter_set_by_num(&data[i], fptr, i+1, dtype,  OutputCol);

			coldef &c = columns[i];
			c.data = (char*)(const_cast<otable::columndef *>(cols[i]))->rawdataptr(c.elementSize, c.width, c.pitch);

			free(ttype[i]);
			free(tform[i]);
		}

		headerWritten = true;
	}

	swatch.start();

//	if(t.using_column("hidden"))
	// call cfitsio Iterator
	int status = 0;
	long nrows;
	fits_get_num_rows(fptr, &nrows, &status);		ASSERT(status == 0) { fits_report_error(stderr, status); }
	fits_insert_rows(fptr, nrows, to-from, &status);	ASSERT(status == 0) { fits_report_error(stderr, status); }

	write_fits_rows_state st(&columns[0], from, to);
	if(t.using_column("hidden"))
	{
		st.hidden = t.col<int>("hidden");
	}
	fits_iterate_data(tfields, &data[0], nrows, 0, write_fits_rows, &st, &status);		ASSERT(status == 0) { fits_report_error(stderr, status); }
	fits_delete_rows(fptr, nrows + st.rowswritten + 1, to-from-st.rowswritten, &status);	ASSERT(status == 0) { fits_report_error(stderr, status); }

	swatch.stop();
	static bool firstTime = true; if(firstTime) { swatch.reset(); kernelRunSwatch.reset(); firstTime = false; }

	if(status != 0) { fits_report_error(stderr, status); }

	return st.rowswritten;
}

bool os_fitsout::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	const char *fn = cfg.count("filename") ? cfg["filename"].c_str() : "sky.fits";

	int status = 0;         /* initialize status before calling fitsio routines */
	unlink(fn);
	fits_create_file(&fptr, fn, &status);   /* create new file */
	MLOG(verb1) << "Output file: " << fn << " (FITS)\n";
	if(status) { abort(); }

	return true;
}

os_fitsout::~os_fitsout()
{
	if(fptr)
	{
		int status = 0;
		if(!header_def.empty())
		{
			int len = header_def.size();

			// create an additional extension with a single column exactly wide enough to store
			// our header
			char *ttype = "HEADER";
			char *tform;
			asprintf(&tform, "%dA", len);
			fits_create_tbl(fptr, BINARY_TBL, 0, 1, &ttype, &tform, NULL, "METADATA", &status);
			ASSERT(status == 0) { fits_report_error(stderr, status); }
			free(tform);

			// write header
			fits_insert_rows(fptr, 0, 1, &status);
			ASSERT(status == 0) { fits_report_error(stderr, status); }
			const char *hstr = header_def.c_str();
			fits_write_col(fptr, TSTRING, 1, 1, 1, 1, &hstr, &status);
			ASSERT(status == 0) { fits_report_error(stderr, status); }
		}
		fits_close_file(fptr, &status);
	}
}
#endif

/////////////////////////////


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
	const char *fn = cfg.count("filename") ? cfg["filename"].c_str() : "sky.cat.txt";
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
#if HAVE_LIBCFITSIO
	else if(name == "fitsout") { s.reset(new os_fitsout); }
#endif
	else if(name == "modelPhotoErrors") { s.reset(new os_modelPhotoErrors); }
	else if(name == "unresolvedMultiples") { s.reset(new os_unresolvedMultiples); }
	else if(name == "FeH") { s.reset(new os_FeH); }
	else if(name == "fixedFeH") { s.reset(new os_fixedFeH); }
	else if(name == "photometry") { s.reset(new os_photometry); }
	else if(name == "photometricErrors") { s.reset(new os_photometricErrors); }
	else if(name == "clipper") { s.reset(new os_clipper); }
	else if(name == "vel2pm") { s.reset(new os_vel2pm); }
	else if(name == "gal2other") { s.reset(new os_gal2other); }
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
//	static const size_t Kbatch = 5000000;
//	static const size_t Kbatch = 23;

	// HACK: Kbatch should be read from skygen.conf, or auto-computed to maximize memory use otherwise
	size_t Kbatch = 5000000;
	EnvVar kb("KBATCH");
	if(kb) { Kbatch = (int)atof(kb.c_str()); } // atof instead of atoi to allow shorthands such as 1e5

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

		if(stage->type() == "input" && !input.empty())  { modcfg.insert(make_pair("filename", input)); }
		if(stage->type() == "output" && !output.empty()) { modcfg.insert(make_pair("filename", output)); }

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
		if(!input.empty()) modcfg.insert(make_pair("filename", input));
		if(!stage->construct(modcfg, t, pipe)) { THROW(EAny, "Failed to initialize output pipeline stage '" + name + "'"); }
		DLOG(verb2) << "postprocessing module loaded: " << name << " (type: " << stage->type() << ")";
		pipe.add(stage);
	}
	if(!pipe.has_module_of_type("output"))
	{
		name = "textout";
		Config modcfg;
		boost::shared_ptr<opipeline_stage> stage( opipeline_stage::create(name) );
		if(!output.empty()) modcfg.insert(make_pair("filename", output));
		if(!stage->construct(modcfg, t, pipe)) { THROW(EAny, "Failed to initialize output pipeline stage '" + name + "'"); }
		DLOG(verb2) << "postprocessing module loaded: " << name << " (type: " << stage->type() << ")";
		pipe.add(stage);
	}

	int nstars = pipe.run(t, rng);
	MLOG(verb2) << "Observing pipeline generated " << nstars << " point sources.";
}
