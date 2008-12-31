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
#include <sstream>

#include "simulate.h"
#include "projections.h"
#include "model.h"
#include "paralax.h"
#include "analysis.h"
#include "dm.h"

#include <vector>
#include <map>
#include <string>

#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

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
		static boost::shared_ptr<opipeline_stage> create(const std::string &name, Config &cfg);
		virtual const std::string &name() = 0;

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
	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual const std::string &name() { static std::string s("FeH"); return s; }

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

		// assign metallicity
		static double A[]     = { 0.63, 0.37 };
		static double sigma[] = { 0.20, 0.20,  0.30};
		static double offs[]  = { 0.00, 0.14, -1.46};

		static double Hmu = 0.5;
		static double muInf = -0.82;
		static double DeltaMu = 0.55;

		switch(component)
		{
			case BahcallSoneira_model::THIN:
			case BahcallSoneira_model::THICK: {
				// choose the gaussian to draw from
				double p = gsl_rng_uniform(rng)*(A[0]+A[1]);
				int i = p < A[0] ? 0 : 1;

				// calculate mean
				double muD = muInf + DeltaMu*exp(-abs(Z/1000)/Hmu);		// Bond et al. A2
				double aZ = muD - 0.067;

				// draw
				FeH = gsl_ran_gaussian(rng, sigma[i]);
				FeH += aZ + offs[i];
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
	return true;
}



// add kinematic information
class os_kin_TMIII : public osink
{
	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual const std::string &name() { static std::string s("kin_TMIII"); return s; }

		os_kin_TMIII() : osink()
		{
			prov.insert("vPhivRvZ[3]");
			req.insert("comp");
			req.insert("XYZ[3]");
		}
};

void getDiskKinematics(float &vPhi, float &vR, float &vZ, const double Z, gsl_rng *rng)
{
	// given XYZ=($1-$3), in pc, generate vPhi, vR, vZ = ($4-$6) for DISK
	// stars. The velocities are physical, i.e. in Galactic coordinate frame
	// (km/s), and based on a vPhi model from Tomography II (astro-ph/0804.3850)
	// with vR and vZ models from Tomography III

	//// Phi component ////
	// we assume that vPhi is a function of Z, as given in TomographyII
	// the distribution shape of vPhi is a sum of 2 (fixed) gaussians
	// with the 3:1 normalization ratio; the first Gaussian is offset
	// by 23 km/s from vPhiMedian(Z) to more negative values, and the
	// second Gaussian is offset by 11 km/s to more positive values;
	double vPhiMedian = -225. + 20. + 19.2*pow(Z/1000, 1.25);   // =0 at Z=6.5 kpc
	if(vPhiMedian > 0) { vPhiMedian = 0; }

	// choose gaussian to draw from
	static const double norm1 = 0.25;  // fraction of stars in the first gaussian
	int i = (gsl_rng_uniform(rng) < norm1) ? 0 : 1;

	// this is eq.13 from Tomography II, with a -220 offset coming from
	// accounting for the vLSR (220 km/s), and 5 km/s from Solar peculiar
	// motion (which was not included in TomographyII)
	static const double offs[]  = { -23., +11. };
	static const double sigPhi[] = { 12.0, 34.0 }; // the intrinsic widths of the Gaussians are sigPhi1 and sigPhi2,
	vPhi = (vPhiMedian + offs[i]) + gsl_ran_gaussian(rng, sigPhi[i]);

#if 0
	// for dispersion gradient (to be added in quadrature)
	define aPhi  0.0
	define bPhi 10.0
	define cPhi  1.0
	set sigPhiZ = $aPhi + $bPhi * (Z/1000)**$cPhi
	set sigPhi1z = sqrt($sigPhi1**2 + sigPhiZ**2)
	set sigPhi2z = sqrt($sigPhi2**2 + sigPhiZ**2)
#endif

	//// R and Z components ////
	// for R and Z, the distribution shape is Gaussian, with spatially
	// invariant mean value, while the dispersions can increase with Z
	static double vRmean = 0.0;
	static double vZmean = 0.0;
	// for dispersions:
	static double aR = 35.0;
	static double bR = 10.0;
	static double cR =  1.0;
	static double aZ = 20.0;
	static double bZ = 10.0;
	static double cZ =  1.0;
	double sigR = aR + bR * pow(Z/1000, cR);
	double sigZ = aZ + bZ * pow(Z/1000, cZ);

	// and now generate vR and vZ
	vR = vRmean + gsl_ran_gaussian(rng, sigR);
	vZ = vZmean + gsl_ran_gaussian(rng, sigZ);
}

void getHaloKinematics(float &vPhi, float &vR, float &vZ, const double Z, gsl_rng *rng)
{
	// given XYZ=($1-$3), in pc, generate vPhi, vR, vZ = ($4-$6) for HALO
	// stars. The velocities are physical, i.e. in Galactic coordinate frame
	// (km/s), and based on a vPhi model from Tomography II (astro-ph/0804.3850)
	// with vR and vZ models from Tomography III

	// it is assumed that the halo kinematics are spatially invariant
	static double sigPhi = 100.0;
	static double sigR   = 110.0;
	static double sigZ   =  90.0;
	static double vPhimean = 0.0;
	static double vRmean   = 0.0;
	static double vZmean   = 0.0;

	vPhi = vPhimean + gsl_ran_gaussian(rng, sigPhi);
	vR = vRmean + gsl_ran_gaussian(rng, sigR);
	vZ = vZmean + gsl_ran_gaussian(rng, sigZ);
}

size_t os_kin_TMIII::push(sstruct *&in, const size_t count, gsl_rng *rng)
{
	// ASSUMPTIONS:
	//	- Bahcall-Soneira component tags exist in input
	//	- galactocentric XYZ coordinates exist in input
	for(size_t i=0; i != count; i++)
	{
		sstruct &s = in[i];

		// fetch prerequisites
		const int component = s.component();
		float Z = s.XYZ()[2];
		float *v = s.vPhivRvZ();

		switch(component)
		{
			case BahcallSoneira_model::THIN:
			case BahcallSoneira_model::THICK:
				getDiskKinematics(v[0], v[1], v[2], Z, rng);
				break;
			case BahcallSoneira_model::HALO:
				getHaloKinematics(v[0], v[1], v[2], Z, rng);
				break;
			default:
				THROW(ENotImplemented, "We should have never gotten here");
		}
	}

	return nextlink->push(in, count, rng);
}

bool os_kin_TMIII::init(const Config &cfg)
{
	return true;
}




// convert input absolute/apparent magnitudes to ugriz colors
class os_ugriz : public osink
{
	protected:
		spline Mr2gi, gi2ri;

		void Mr2gi_build();
		void gi2ri_build();

		float MrFeH2gi(float Mr, float FeH);
		float grFeH2ug(float gr, float FeH);
		float ri2iz(float ri);
	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual const std::string &name() { static std::string s("ugriz"); return s; }

		os_ugriz() : osink()
		{
			prov.insert("ugriz[5]");
			req.insert("sdss_r");
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
	
		std::stringstream ss; ss << "Could not find roots for ug given FeH=" << FeH << " and gr=" << gr << "\n";
		THROW(EAny, ss.str());
	}

	bool x0q = 0.5 < x0 && x0 < 1.8;
	bool x1q = 0.5 < x1 && x1 < 1.8;
	if(x0q == x1q) {
		if(x0q) { DLOG(verb3) << "WARNING: Real u-g root ambiguous: x0,x1 = " << x0 << ", " << x1; }
		else    { DLOG(verb3) << "WARNING: Both roots probably wrong: x0,x1 = " << x0 << ", " << x1; }
		x0q = abs(x0 - 1.05) < abs(x1 - 1.05) ? true : false;	// return the one closer to 1.05
		x1q = !x0q;
	}
	return x0q ? x0 : x1;
#if 0
	// add Gaussian noise
	randomizeG Zraw 0 0.1 $3
	unset Zraw
#endif
}

void os_ugriz::Mr2gi_build()
{
	// coefficients of Mr(g-i) relation
	// double coeff[] = {-0.56, 14.32, -12.97, 6.127, -1.267, 0.0967};
	std::string fn = datadir() + "/gi2Mr.FeH=0.txt";
	FILE_EXISTS_OR_THROW(fn);
	MLOG(verb2) << "Mr2gi locus: " << fn;

	text_input_or_die(f, fn);
	std::vector<double> gi, Mr;
	load(f, gi, 0, Mr, 1);
	Mr2gi.construct(Mr, gi);
}

float os_ugriz::MrFeH2gi(float Mr, float FeH)
{
	// offset due to metallicity
	double FeHDelta = FeH*(-1.11 - 0.18*FeH);

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
		const double Mr = s.absmag();
		const float FeH = s.FeH();

		// construct SDSS colors given the absolute magnitude and metallicity
		float gi = MrFeH2gi(Mr, FeH);
		float ri = gi2ri(gi);
		float gr = gi - ri;
		float ug = grFeH2ug(gr, FeH);
		float iz = ri2iz(ri);

		// convert colors to apparent magnitudes
		float r  = s.sdss_r();
		float *mags = s.sdss_mag();
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
	Mr2gi_build();
	gi2ri_build();

	return true;
}



// in/out ends of the chain
class os_textout : public osink
{
	protected:
		std::ofstream out;
		bool headerWritten;
		ticker tick;

	public:
		virtual size_t push(sstruct *&data, const size_t count, gsl_rng *rng);
		virtual bool init(const Config &cfg);
		virtual int priority() { return 1000000; }	// ensure this stage has the least priority
		virtual const std::string &name() { static std::string s("textout"); return s; }

		os_textout() : osink(), headerWritten(false), tick(-1)
		{
		}
};

size_t os_textout::push(sstruct *&data, const size_t count, gsl_rng *rng)
{
	if(tick.step <= 0) { tick.open("Processing", 10000); }

	if(!headerWritten)
	{
		out << "# " << sstruct::factory << "\n";
		headerWritten = true;
	}

	size_t cnt = 0;
	while(cnt < count && (out << data[cnt] << "\n")) { cnt++; tick.tick(); }

	if(!out) { THROW(EIOException, "Error outputing data"); }

	delete [] data;
	data = NULL;

	return count;
}

bool os_textout::init(const Config &cfg)
{
	if(!cfg.count("filename")) { THROW(EAny, "Keyword 'filename' must exist in config file"); }
	out.open(cfg["filename"].c_str());

	return out;
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
		std::ifstream in;

	public:
		virtual bool init(const Config &cfg);
		virtual int run(gsl_rng *rng);
		virtual const std::string &name() { static std::string s("textin"); return s; }

		os_textin() {};
};

bool os_textin::init(const Config &cfg)
{
	const char *fn = cfg["filename"].c_str();
	in.open(fn);
	if(!in) { THROW(EFile, "Failed to open '" + (std::string)fn + "' for input."); }

	in >> sstruct::factory;
	return in;
}

int os_textin::run(gsl_rng *rng)
{
	size_t block = 10000;

	size_t total = 0;
	do {
		sstruct *data = sstruct::create(block);
		size_t cnt = 0;
		while(cnt < block && in >> data[cnt]) { cnt++; }
		total += nextlink->push(data, cnt, rng);
	} while(in);

	return total;
}

boost::shared_ptr<opipeline_stage>
	opipeline_stage::create(const std::string &name, Config &cfg)
{
	boost::shared_ptr<opipeline_stage> s;

	if(name == "textin") { s.reset(new os_textin); }
	else if(name == "textout") { s.reset(new os_textout); }
	else if(name == "ugriz") { s.reset(new os_ugriz); }
	else if(name == "FeH") { s.reset(new os_FeH); }
	else if(name == "kin_TMIII") { s.reset(new os_kin_TMIII); }
	else { THROW(EAny, "Module " + name + " unknown."); }

	ASSERT(name == s->name());

	bool succ = s->init(cfg);
	if(!succ) { THROW(EAny, "Failed to initialize output pipeline stage '" + name + "'"); }
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
				sstruct::factory.useTag(*j);
			}

			// erase from the list of pending stages
			stages.erase(i);
			foundOne = true;
			break;
		}
		if(!foundOne)
		{
			std::stringstream ss;
			FOREACH(stages)
			{
				opipeline_stage &s = *i->second;
				ss << s.name() << ", ";
			}
			THROW(EAny, "Modules " + ss.str() + " depend on fields that no other modules (or input) provide.");
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
	MLOG(verb1) << "Output module chain: " << ss.str();

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

		boost::shared_ptr<opipeline_stage> stage( opipeline_stage::create(moddef.first, modcfg) );

		ASSERT(stage) { std::cerr << "Module " << moddef.first << " unknown or failed to load.\n"; }
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
#endif
