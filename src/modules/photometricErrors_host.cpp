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

#include "../pipeline.h"
#include "FeH_gpu.cu.h"

#include "spline.h"
#include "analysis.h"
#include <fstream>
#include "io.h"

#include <astro/system/config.h>
#include <astro/useall.h>

// add photometric errors
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
	//virtual int priority() { return PRIORITY_INSTRUMENT; }	// ensure this stage has the least priority
	virtual double ordering() const { return ord_detector; }

	os_photometricErrors() : osink()
	{
	}
};
extern "C" opipeline_stage *create_module_photometricerrors() { return new os_photometricErrors(); }	// Factory; called by opipeline_stage::create()

// add photometric errors, givan an analytic model
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
extern "C" opipeline_stage *create_module_modelphotoerrors() { return new os_modelPhotoErrors(); }	// Factory; called by opipeline_stage::create()

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

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
	// Built-in filenames are of the form $/<photosys>.<bandname>.photoerr.txt
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

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

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

