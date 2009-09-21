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

#ifndef __observe_h
#define __observe_h

#include "simulate_base.h"
#include "gpu.h"
#include "gpc_cpp.h"
#include <astro/system/config.h>
#include "projections.h"
#include <vector>

#include "pipeline.h"

// Clips out stars not within the requested observation area
// Currently only used as support to os_skygen
class partitioned_skymap;
class os_clipper : public osink
{
protected:
	// hemisphere -- pixelized north/south sky (aux class).
	struct hemisphere
	{
		peyton::math::lambert 	proj;	// projection
		partitioned_skymap 	*sky;	// skymap of the hemisphere polygon

		hemisphere() : sky(NULL) {}
		~hemisphere() { delete sky; }
	private:
		// disallow copy, copy constructable
		hemisphere &operator=(const hemisphere&);
		hemisphere(const hemisphere&);
	};

	hemisphere hemispheres[2];	// pixelized northern and southern sky

public:
	struct pixel
	{
		double l, b; // in Radians
		double X, Y; // lambert coordinates
		int projIdx; // index of the projection for (X,Y)<->(l,b) transform (hemispheres[projIdx].proj)
		float pixelArea, coveredArea;	// area of the nominal pixel, area covered by the footprint within the pixel
		
		pixel(double l_, double b_, double X_, double Y_, int projIdx_, float pixelArea_, float coveredArea_)
			: X(X_), Y(Y_), l(l_), b(b_), projIdx(projIdx_), pixelArea(pixelArea_), coveredArea(coveredArea_) {}
	};

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe); // NOTE: overriden as it's abstract, but should NEVER be called directly. Use construct_from_hemispheres() instead.

	virtual const std::string &name() const { static std::string s("clipper"); return s; }
	virtual int priority() { return PRIORITY_INSTRUMENT; } // ensure this module is placed near the end of the pipeline

	// constructs the clipper object from north/south hemispheres in projection proj
	void construct_from_hemispheres(float dx, const peyton::math::lambert &nproj, const std::pair<gpc_polygon, gpc_polygon> &sky);

	int getPixelCenters(std::vector<os_clipper::pixel> &pix) const;			// returns the centers of all pixels
	int getProjections(std::vector<std::pair<double, double> > &ppoles) const;	// returns the poles of all used projections

	os_clipper() : osink()
	{
		req.insert("projIdx");
		req.insert("projXY");
		prov.insert("hidden");
	}
};

// GPU generator input
class os_skygen : public osource
{
protected:
	std::vector<boost::shared_ptr<skyConfigInterface> > kernels;
	size_t nstarLimit;	// maximum number of stars to generate
	float nstars;		// the mean number of stars to generate (if nstars=0, the number will be determined by the model)
	cuxTexture<float, 3>	ext_north, ext_south;	// north/south extinction maps

public:
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
	virtual size_t run(otable &t, rng_t &rng);
	virtual const std::string &name() const { static std::string s("skygen"); return s; }
	virtual const std::string &type() const { static std::string s("input"); return s; }

protected:
	const os_clipper &load_footprints(const std::string &footprints, float dx, opipeline &pipe);
	skyConfigInterface *create_kernel_for_model(const std::string &model);
	int load_models(otable &t, skygenConfig &sc, const std::string &model_cfg_list, const os_clipper &clipper);
	void load_pdf(float &dx, skygenConfig &sc, otable &t, const std::string &cfgfn);
	void load_extinction_maps(const std::string &econf);
};

#endif // __observe_h
