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
#include <astro/system/config.h>

class osource : public opipeline_stage
{
	public:
		virtual int priority() { return PRIORITY_INPUT; } // ensure highest priority for this stage

	public:
		osource() : opipeline_stage()
		{
			prov.insert("_source");
		}
};

// GPU generator input
class os_skygen : public osource
{
	public:
		virtual bool init(const peyton::system::Config &cfg, otable &t);
		virtual size_t run(otable &t, rng_t &rng);
		virtual const std::string &name() const { static std::string s("skygen"); return s; }
		virtual const std::string &type() const { static std::string s("input"); return s; }

		os_skygen() {};
};

// add Fe/H information
class os_FeH : public osink, os_FeH_data
{
public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool init(const peyton::system::Config &cfg, otable &t);
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
		virtual bool init(const peyton::system::Config &cfg, otable &t);
		virtual const std::string &name() const { static std::string s("fixedFeH"); return s; }

		os_fixedFeH() : osink(), fixedFeH(0)
		{
			prov.insert("FeH");
		}
};

// convert velocities to proper motions
class os_vel2pm : public osink , public os_vel2pm_data
{	
public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool init(const peyton::system::Config &cfg, otable &t);
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
		virtual bool init(const peyton::system::Config &cfg, otable &t);
		virtual const std::string &name() const { static std::string s("kinTMIII"); return s; }

		os_kinTMIII() : osink()
		{
			prov.insert("vcyl");
			req.insert("comp");
			req.insert("XYZ");
		}
};

#endif // __observe_h
