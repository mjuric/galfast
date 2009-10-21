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

#ifndef pipeline_h__
#define pipeline_h__

#include "otable.h"

#include <list>

namespace peyton { namespace system { class Config; }};

class osink;
class opipeline;
class opipeline_stage
{
	protected:
		std::set<std::string> prov, req;	// add here the fields required/provided by this modules, if using stock runtime_init() implementation
//		std::string uniqueId;			// a string uniquely identifying this module instance
		stopwatch swatch;			// times how long it takes to process() this stage

		osink *nextlink;
	public:
		void chain(osink *nl) { nextlink = nl; }
		float getProcessingTime() { return swatch.getTime(); }

//		void setUniqueId(const std::string &uid) { uniqueId = uid; }
//		const std::string &getUniqueId() const { return uniqueId; }

	public:
		static boost::shared_ptr<opipeline_stage> create(const std::string &name);
		virtual const std::string &name() const = 0;
		virtual const std::string &type() const { static std::string s("stage"); return s; }

	public:
		virtual bool construct(const peyton::system::Config &cfg, otable &t, opipeline &pipe) = 0;
		virtual bool runtime_init(otable &t);
		virtual size_t run(otable &t, rng_t &rng) = 0;
		virtual ~opipeline_stage() {};

		static const int PRIORITY_INPUT      = -10000;
		static const int PRIORITY_STAR       =      0;
		static const int PRIORITY_SPACE      =    100;
		static const int PRIORITY_INSTRUMENT =   1000;
		static const int PRIORITY_OUTPUT     =  10000;
		virtual int priority() { return PRIORITY_SPACE; }

	public:
		opipeline_stage() : nextlink(NULL)
		{
		}
};

class osink : public opipeline_stage
{
	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng) = 0;

	public:
		virtual size_t run(otable &t, rng_t &rng) { THROW(peyton::exceptions::ENotImplemented, "We should have never gotten here"); } // we should never reach this place

	public:
		osink() : opipeline_stage() {}
};

// Class encapsulating the pipeline of modules that add/modify
// properties of generated objects, as they move from generation
// to output.
//
// Constructed and run by ::postprocess_catalog()
class opipeline
{
	public:
		bool dryrun;	// whether to gerate (draw) the catalog, or stop after computing the expected number of stars

	public:
		std::list<boost::shared_ptr<opipeline_stage> > stages;	// the pipeline (an ordered list of stages)

	public:
		void add(const boost::shared_ptr<opipeline_stage> &pipe) { stages.push_back(pipe); }
		bool create_and_add(
			peyton::system::Config &modcfg, otable &t,
			size_t maxstars, size_t nstars,
			const std::string &models, const std::string &foots, const std::string &extmaps,
			const  std::string &input, const std::string &output);

		virtual size_t run(otable &t, rng_t &rng);

		bool has_module_of_type(const std::string &type) const;
	public:
		opipeline(bool dryrun_) : dryrun(dryrun_) {}
};

//
// Specialization for input (source) modules. os_skygen overrides this
// class.
//
class osource : public opipeline_stage
{
	public:
		virtual int priority() { return PRIORITY_INPUT; } // ensure highest priority for this stage

	public:
		osource() : opipeline_stage() {}
};

#endif // pipeline_h__
