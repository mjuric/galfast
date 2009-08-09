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
#include "model.h"
#include "gpu.h"
#include "gpc_cpp.h"
#include <astro/system/config.h>
#include <projections.h>
#include <vector>

class opipeline
{
	public:
		std::list<boost::shared_ptr<opipeline_stage> > stages;

	public:
		void add(boost::shared_ptr<opipeline_stage> pipe) { stages.push_back(pipe); }
		virtual size_t run(otable &t, rng_t &rng);

		bool has_module_of_type(const std::string &type) const;
};

class osource : public opipeline_stage
{
	public:
		virtual int priority() { return PRIORITY_INPUT; } // ensure highest priority for this stage

	public:
		osource() : opipeline_stage() {}
};

// Clips out stars not within the requested observation area
// Currently only used as support to os_skygen
class partitioned_skymap;
class os_clipper : public osink
{
protected:
//	float dx/*, dA*/;			// linear scale and angular area of each pixel (rad, rad^2)
//	float bmin;				// zone of avoidance around the Galactic plane

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

	hemisphere hemispheres[2];

public:
	struct pixel
	{
		double l, b; // in Radians
		int projIdx;
		float pixelArea, coveredArea;	// area of the nominal pixel, area covered by the footprint within the pixel
		
		pixel(double l_, double b_, int projIdx_, float pixelArea_, float coveredArea_)
			: l(l_), b(b_), projIdx(projIdx_), pixelArea(pixelArea_), coveredArea(coveredArea_) {}
	};

public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);

	virtual const std::string &name() const { static std::string s("clipper"); return s; }
	virtual int priority() { return PRIORITY_INSTRUMENT; } // ensure this is placed near the end of the pipeline

	// constructs the clipper object from north/south hemispheres in projection proj
	void construct_from_hemispheres(float dx, const peyton::math::lambert &nproj, const std::pair<gpc_polygon, gpc_polygon> &sky);

	int getPixelCenters(std::vector<os_clipper::pixel> &pix) const;			// returns the centers of all pixels
	int getProjections(std::vector<std::pair<double, double> > &ppoles) const;	// returns the poles of all used projections

	os_clipper() : osink()
	{
		req.insert("lb");
		prov.insert("hidden");
	}
};

// GPU generator input
class os_skygen : public osource
{
protected:
	std::vector<boost::shared_ptr<skyConfigInterface> > kernels;
	size_t nstarLimit;	// maximum number of stars to generate

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

// 	bool init(
// 		const peyton::system::Config &cfg,
// 		const peyton::system::Config &pdf_cfg,
// 		const peyton::system::Config &foot_cfg,
// 		const peyton::system::Config &model_cfg,
// 		otable &t,
// 		opipeline &pipe);

/*	os_skygen() : skygen(NULL) {};
	virtual ~os_skygen() { delete skygen; }*/
};

// add Fe/H information
class os_FeH : public osink, os_FeH_data
{
public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("FeH"); return s; }

	os_FeH() : osink()
	{
		prov.insert("FeH");
		req.insert("comp");
		req.insert("XYZ");
	}
};

// add Fe/H information
class os_fixedFeH : public osink
{
	protected:
		float fixedFeH;

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
		virtual const std::string &name() const { static std::string s("fixedFeH"); return s; }

		os_fixedFeH() : osink(), fixedFeH(0)
		{
			prov.insert("FeH");
		}
};

// unresolved multiple system creator
class os_unresolvedMultiples : public osink
{
	protected:
		std::string absmagSys;
		multiplesAlgorithms::algo algo;			// algorithm for magnitude assignment to secondaries

	public:
		virtual bool runtime_init(otable &t);
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
		virtual const std::string &name() const { static std::string s("unresolvedMultiples"); return s; }
		virtual int priority() { return PRIORITY_STAR; } // ensure this is placed near the beginning of the pipeline

		os_unresolvedMultiples() : osink()
		{
			req.insert("absmag");
		}
};

// convert velocities to proper motions
class os_vel2pm : public osink , public os_vel2pm_data
{
protected:
	std::string output_col_name;
public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
	virtual const std::string &name() const { static std::string s("vel2pm"); return s; }

	os_vel2pm() : osink()
	{
		coordsys=GAL;
		req.insert("lb");
		req.insert("XYZ");
		req.insert("vcyl");
	}
};
	

// add kinematic information
class os_kinTMIII : public osink, os_kinTMIII_data
{	
	float DeltavPhi;
	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe);
		virtual const std::string &name() const { static std::string s("kinTMIII"); return s; }

		os_kinTMIII() : osink()
		{
			prov.insert("vcyl");
			req.insert("comp");
			req.insert("XYZ");
		}
};

#endif // __observe_h
