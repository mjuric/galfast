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
#include "photometry.h"

#include "spline.h"
#include "analysis.h"
#include "io.h"
#include <fstream>

#include <astro/system/config.h>
#include <astro/useall.h>

// compute magnitudes in the desired photometric system
class os_photometry : public osink, public os_photometry_data
{
protected:
	/* Aux structure used while loading the isochrones */
	struct Mr2col
	{
		static const size_t maxcolors = N_REDDENING-1;

		double Mr;
		double c[maxcolors];

		bool operator <(const Mr2col &a) const { return Mr < a.Mr; }
	};

	/* Aux structure used while resampling isochrones to textures */
	struct colsplines
	{
		double Mrmin, Mrmax;
		spline s[Mr2col::maxcolors];
		spline &operator [](const size_t i) { return s[i]; }
	};

protected:
	std::string bandset2;			// name of this filter set
	std::string absbband;			// Absolute magnitude band for which the datafile gives col(absmag,FeH) values. By default, it's equal to "abs$bband". Must be supplied by other modules.
	std::string photoFlagsName;		// Name of the photometric flags field
	std::vector<std::string> bnames;	// band names (e.g., LSSTr, LSSTg, SDSSr, V, B, R, ...)
	std::vector<cuxTexture<float4, 2> > isochrones;	// Isochrone (stellar loci) texture
	std::vector<cuxTexture<float4, 2> > eflags;	// Flags
	std::pair<cuxTexture<float, 3>, cuxTexture<float, 3> > extinction;	// north/south extinction maps (for bootstrap band)

protected:
	typedef boost::shared_ptr<cuxTextureBinder> tbptr;
	void bind_isochrone(std::list<tbptr> &binders, cuxTextureReferenceInterface &texc, cuxTextureReferenceInterface &texf, int idx);

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("photometry"); return s; }

	os_photometry() : osink()
	{
		comp0 = 0;
		comp1 = 0xffffffff;
		FOR(0, N_REDDENING) { reddening[i] = 1.f; }

		req.insert("FeH");
	}
};
extern "C" opipeline_stage *create_module_photometry() { return new os_photometry(); }	// Factory; called by opipeline_stage::create()

template<typename T> inline OSTREAM(const std::vector<T> &v) { FOREACH(v) { out << *i << " "; }; return out; }

bool os_photometry::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// load component IDs to which this module will apply
	cfg.get(comp0,   "comp0",   0U);
	cfg.get(comp1,   "comp1",   0xffffffff);

	// load bandset name
	std::string tmp, bname;
	cfg.get(bandset2,   "bandset",   "LSSTugrizy");

	// load band names and construct otable field definition
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
	ncolors = bnames.size()-1;

	std::string bfield = bandset2 + "[" + str(nbands) + "]{class=magnitude;fieldNames=" + sbnames.str() + ";}";
	prov.insert(bfield);

	// check the number of colors doesn't exceed the maximum
	if(ncolors >= Mr2col::maxcolors)
	{
		THROW(EAny, "This module cannot handle more than " + str(Mr2col::maxcolors+1) + " bands.");
	}

	// construct photoFlags otable field definition
	std::string bandsetname = bandset2;
	photoFlagsName = bandsetname + "PhotoFlags";
	prov.insert(photoFlagsName + "{class=flags;}");

	// convert the bootstrap band name into band index
	std::string bband;
	cfg.get(   bband,   "bootstrap_band",   "LSSTr");
	cfg.get(absbband,   "absband",      "abs"+bband);
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

	// load isochrone texture sampling parameters
	float FeH0, FeH1, dFeH, Mr0, Mr1, dMr;
	cfg.get(tmp,   "absmag_grid",   "3 15 0.01"); ss.clear(); ss.str(tmp);
	if(!(ss >> Mr0 >> Mr1 >> dMr))
	{
		THROW(EAny, "Error reading Mr field from config file. Expect Mr = <Mr0> <Mr1> <dMr>, got Mr = " + tmp);
	}
	cfg.get(tmp,   "FeH_grid",   "-3 0.5 0.01"); ss.clear(); ss.str(tmp);
	if(!(ss >> FeH0 >> FeH1 >> dFeH))
	{
		THROW(EAny, "Error reading FeH field from config file. Expect FeH = <FeH0> <FeH1> <dFeH>, got FeH = " + tmp);
	}

	// find and open the definition file for the isochrones
	cfg.get(tmp,   "file",   "");
	if(tmp == "")
	{
		// Try to find internal definitions for the photometric system
		tmp = datadir() + "/" + bandset2 + ".photosys.txt";
		if(!file_exists(tmp))
		{
			THROW(EAny, "Default photosys. definition file " + tmp + " for " + bandset2 + " not found. Specify it explicitly using the file=xxx keyword.");
		}
	}

	// load the isochrone table into a map indexed by FeH
	text_input_or_die(in, tmp);
	std::map<double, std::vector<Mr2col> > v; // map of the form v[FeH] -> vec(Mr,colors)
	Mr2col mc; double FeH;
	bind(in, mc.Mr,0, FeH,1);
	FOR(0, ncolors) { bind(in, mc.c[i], i+2); }
	while(in.next())
	{
		v[FeH].push_back(mc);
	}

	// compute the needed texture size, and the texture coordinates
	int nFeH = (int)((FeH1-FeH0)/dFeH + 1);	// +1 pads it a little, to ensure FeH1 gets covered
	int nMr  = (int)((Mr1 -Mr0) /dMr  + 1);
	float2 tcFeH = texcoord_from_range(0, nFeH, FeH0, FeH0 + nFeH*dFeH);
	float2 tcMr  = texcoord_from_range(0, nMr,   Mr0,  Mr0 +  nMr*dMr);

	MLOG(verb1) << "Photometry: Generating " << bandset2 << " ( " << bnames << ")";
	MLOG(verb2) << bandset2 << ": Generating " << bnames.size() << " bands: " << bnames;
	MLOG(verb2) << bandset2 << ": Input absolute magnitude assumed to be in " << bband << " band.";
	MLOG(verb2) << bandset2 << ": Using color(" << absbband << ", FeH) table from " << tmp << ".";
	MLOG(verb2) << bandset2 << ": Resampling color table to fast lookup grid:";
	MLOG(verb2) << bandset2 << ":    " << absbband << "0, " << absbband << "1, d(" << absbband << ") = " << Mr0 << ", " << Mr1 << ", " << dMr << ".";
	MLOG(verb2) << bandset2 << ":    FeH0, FeH1, dFeH = " << FeH0 << ", " << FeH1 << ", " << dFeH << ".";

	// construct col(Mr) splines for each FeH line present in the input.
	// The map s will hold nband splines giving color(Mr) at fixed FeH (FeH is the key).
	std::map<double, colsplines> s;	// s[FeH] -> (colors)=spline(Mr)
	std::vector<double> vFeH;	// keys of s
	vFeH.reserve(v.size());
	FOREACH(v)
	{
		double FeH = i->first;
		vFeH.push_back(FeH);

		std::vector<Mr2col> &m2c = i->second;
		colsplines &ss = s[FeH];			// Splines that will return col(Mr) for this FeH

		std::sort(m2c.begin(), m2c.end());		// Sort m2c array by Mr
		std::vector<double> 	Mr(m2c.size()),
					 c(m2c.size());
		FOR(0, m2c.size()) { Mr[i] = m2c[i].Mr; }	// A sorted array of Mr
		FOR(0, ncolors) 				// for each color i, construct and store its color_i(Mr) spline
		{
			FORj(j, 0, m2c.size()) { c[j] = m2c[j].c[i]; }
			ss[i].construct(Mr, c);
		}

		// store beginning/end for test of extrapolation
		ss.Mrmin = Mr.front();
		ss.Mrmax = Mr.back();
	}

	// thread in Fe/H direction through the texture, constructing a col(FeH)
	// spline at given Mr using knots computed from previously precomputed
	// col(Mr) splines. Resample the resulting spline to texture grid.
	std::vector<double>  vcol(s.size());			// this will hold the knots of col(FeH) spline
	std::vector<double> vecol(s.size());			// will be set to 1 if the corresponding color knot was extrapolated, 0 otherwise
	int compidx = -1;
	std::vector<double> fraction_extrapolated(ncolors, 0.);	// collecting statistics: how many points in (FeH,Mr) grid had to be extrapolated
	FORj(ic, 0, ncolors)
	{
		//
		// Since the cost of a texture fetch of float4 is equal to that of just a float (I _think_! TODO: check)
		// we pack up to four colors into a single float4 texture, instead of having each
		// (Mr, FeH)->color mapping from a separate texture.
		//
		compidx++; compidx %= 4;	// ic % 4; the index of the (float/uint)4 component within the texture where this color will be stored
		if(compidx == 0)
		{
			// add new texture
			isochrones.push_back( cuxTexture<float4, 2>(nFeH, nMr, tcFeH, tcMr) );
			    eflags.push_back( cuxTexture<float4, 2>(nFeH, nMr, tcFeH, tcMr) );
		}

		cuxTexture<float4, 2> texc = isochrones.back();
		cuxTexture<float4, 2> texf  = eflags.back();

		// go throught all texture pixels, Mr first
		FORj(m, 0, nMr)
		{
			double Mr = Mr0 + m*dMr;

			// construct col(FeH) and eflag(FeH) spline at given Mr
			int k=0;
			FOREACH(s)
			{
				 vcol[k] = i->second[ic](Mr);					// compute color(FeH|Mr). FeH is i->first (and in vFeH)
				vecol[k] = i->second.Mrmin > Mr || Mr > i->second.Mrmax;	// 1 if extrapolation was needed to compute vcol[k]
				k++;
			}

			spline FeH_to_color(vFeH, vcol);
			spline FeH_to_eflag(vFeH, vecol);	// zero where neither adjacent knot was extrapolated, nonzero otherwise
			FORj(f, 0, nFeH)
			{
				// compute color(FeH, Mr)
				double FeH = FeH0 + f*dFeH;

				float color   = FeH_to_color(FeH);
				bool  eflag   = vFeH.front() > FeH || FeH > vFeH.back()		// extrapolated in FeH direction ?
						|| FeH_to_eflag(FeH) != 0.;			// extrapolated in Mr direction ?

				// pack into float4 texture
				float *texcval = (float *)&(texc(f, m));
				float *texfval = (float *)&(texf(f, m));

				texcval[compidx] = color;
				texfval[compidx] = eflag;

				// compute some summary statistics
				if(eflag) { fraction_extrapolated[ic]++; }
			}
		}

		fraction_extrapolated[ic] /= nFeH*nMr;
	}

	MLOG(verb2) << bandset2 << ":    grid size = " << nFeH << " x " << nMr << " (" << isochrones[0].memsize() << " bytes).";
	MLOG(verb2) << bandset2 << ":    extrapolation fractions = " << fraction_extrapolated;

	return true;
}

// helper for os_photometry::process()
inline void os_photometry::bind_isochrone(std::list<tbptr> &binders, cuxTextureReferenceInterface &texc, cuxTextureReferenceInterface &texf, int idx)
{
	if(idx >= isochrones.size()) { return; }

	binders.push_back(tbptr(new cuxTextureBinder(texc, isochrones[idx])));
	binders.push_back(tbptr(new cuxTextureBinder(texf, eflags[idx])));
}

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

	{
		// Bind all used textures. The list will be autodeallocated on exit
		// from the block, triggering cuxTextureBinder destructors and
		// unbinding the textures.
		std::list<tbptr> binders;
		bind_isochrone(binders, color0, cflags0, 0);
		bind_isochrone(binders, color1, cflags1, 1);
		bind_isochrone(binders, color2, cflags2, 2);
		bind_isochrone(binders, color3, cflags3, 3);

		cuxUploadConst("os_photometry_params", static_cast<os_photometry_data&>(*this));

		CALL_KERNEL(os_photometry_kernel, otable_ks(begin, end, -1, sizeof(float)*ncolors), Am, flags, DM, Mr, Mr.width(), mags, FeH, comp);
	}

	return nextlink->process(in, begin, end, rng);
}
