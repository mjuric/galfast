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

#include "analysis.h"
#include "io.h"
#include "pipeline.h"

#include <fstream>

#include <astro/io/format.h>
#include <astro/system/log.h>
#include <astro/system/config.h>
#include <astro/useall.h>

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


// in/out ends of the chain
class os_textout : public osink
{
	protected:
		flex_output out;

		bool headerWritten;
		ticker tick;

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
	swatch.start();

	ticker tick("Writing output", (int)ceil((to-from)/50.));

	if(!headerWritten)
	{
		out.out() << "# ";
		t.serialize_header(out.out());
		out.out() << "\n";
		headerWritten = true;
	}


	size_t nserialized = 0;

	if(t.using_column("hidden"))
	{
		cint_t::host_t   hidden = t.col<int>("hidden");
		nserialized = t.serialize_body(out.out(), from, to, mask_output(hidden, tick));
	}
	else
	{
		nserialized = t.serialize_body(out.out(), from, to);
	}

	if(!out.out()) { THROW(EIOException, "Error outputing data"); }

	swatch.stop();
	//static bool firstTime = true; if(firstTime) { swatch.reset(); kernelRunSwatch.reset(); firstTime = false; }

	return nserialized;
}

bool os_textout::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	const char *fn = cfg.count("filename") ? cfg["filename"].c_str() : "sky.obs.txt";
	out.open(fn);
	MLOG(verb1) << "Output file: " << fn << " (text)\n";

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

		void createOutputTable(otable &t);

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
		hidden.reset();
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

void os_fitsout::createOutputTable(otable &t)
{
	if(headerWritten) { return; }

	// fetch columns we're going to write
	std::vector<const otable::columndef *> cols;
	t.getSortedColumnsForOutput(cols);
	const int tfields = cols.size();

	// collect header metadata -- we'll write this out into a separate table
	// in the destructor
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

	FOR(0, tfields)
	{
		free(ttype[i]);
		free(tform[i]);
	}

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
	}

	headerWritten = true;
}

size_t os_fitsout::process(otable &t, size_t from, size_t to, rng_t &rng)
{
	ticker tick("Writing output", (int)ceil((to-from)/50.));

	createOutputTable(t);

	// Get data pointers from output the table
	std::vector<const otable::columndef *> cols;
	t.getSortedColumnsForOutput(cols);
	FOR(0, columns.size())
	{
		coldef &c = columns[i];
		c.data = (char*)(const_cast<otable::columndef *>(cols[i]))->rawdataptr(c.elementSize, c.width, c.pitch);
	}

	swatch.start();

	// append the (maximum) number of rows we're going to write
	int status = 0;
	long nrows;
	fits_get_num_rows(fptr, &nrows, &status);		ASSERT(status == 0) { fits_report_error(stderr, status); }
	fits_insert_rows(fptr, nrows, to-from, &status);	ASSERT(status == 0) { fits_report_error(stderr, status); }

	// call cfitsio Iterator
	write_fits_rows_state st(&columns[0], from, to);
	if(t.using_column("hidden"))
	{
		st.hidden = t.col<int>("hidden");
	}
	fits_iterate_data(columns.size(), &data[0], nrows, 0, write_fits_rows, &st, &status);	ASSERT(status == 0) { fits_report_error(stderr, status); }

	// truncate any extra rows
	fits_delete_rows(fptr, nrows + st.rowswritten + 1, to-from-st.rowswritten, &status);	ASSERT(status == 0) { fits_report_error(stderr, status); }

	swatch.stop();
	//static bool firstTime = true; if(firstTime) { swatch.reset(); kernelRunSwatch.reset(); firstTime = false; }

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
			//static bool firstTime = true; if(firstTime) { swatch.reset(); kernelRunSwatch.reset(); firstTime = false; }
			total += nextlink->process(t, 0, t.size(), rng);
		}
	} while(in.in());

	return total;
}

#include <dlfcn.h>

typedef opipeline_stage *(*moduleFactory_t)();

// TODO: Replace all explicit instantiations with calls to factory functions
boost::shared_ptr<opipeline_stage> opipeline_stage::create(const std::string &name)
{
	boost::shared_ptr<opipeline_stage> s;

	if(name == "textin") { s.reset(new os_textin); }
//	else if(name == "skygen") { s.reset(new os_skygen); }
	else if(name == "textout") { s.reset(new os_textout); }
#if HAVE_LIBCFITSIO
	else if(name == "fitsout") { s.reset(new os_fitsout); }
#endif
//	else if(name == "modelPhotoErrors") { s.reset(new os_modelPhotoErrors); }
//	else if(name == "unresolvedMultiples") { s.reset(new os_unresolvedMultiples); }
//	else if(name == "FeH") { s.reset(new os_FeH); }
//	else if(name == "fixedFeH") { s.reset(new os_fixedFeH); }
//	else if(name == "photometry") { s.reset(new os_photometry); }
//	else if(name == "photometricErrors") { s.reset(new os_photometricErrors); }
//	else if(name == "clipper") { s.reset(new os_clipper); }
//	else if(name == "vel2pm") { s.reset(new os_vel2pm); }
//	else if(name == "gal2other") { s.reset(new os_gal2other); }
//	else if(name == "kinTMIII") { s.reset(new os_); }
	else
	{
		// try loading using a factory function
		void *me = dlopen(NULL, RTLD_LAZY);
		if(me == NULL)
		{
			const char *err = dlerror();
			THROW(EAny, err);
		}

		std::string factory_name = "create_module_" + normalizeKeyword(name);

		DLOG(verb2) << "Looking for " << factory_name << " factory function (for module '" << name << "')";
		moduleFactory_t factory = (moduleFactory_t)dlsym(me, factory_name.c_str());
		if(factory)
		{
			s.reset(factory());
		}
		else
		{
			THROW(EAny, "Module " + name + " unknown.");
		}
	}

//	ASSERT(name == s->name());

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

	int ret = source->run(t, rng);

	MLOG(verb2) << "Module runtimes:";
	FOREACH(pipeline)
	{
		MLOG(verb2) << io::format("  %17s: %f") << (*i)->name() << (*i)->getProcessingTime();
	}
	MLOG(verb2) << "GPU kernels runtime: " << kernelRunSwatch.getTime();

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

bool opipeline::create_and_add(
	Config &modcfg, otable &t,
	size_t maxstars, size_t nstars,
	const std::string &models, const std::string &foots, const std::string &extmaps,
	const std::string &input, const std::string &output
)
{
	// create the module
	std::string module = modcfg["module"];
	boost::shared_ptr<opipeline_stage> stage( opipeline_stage::create( module ) );
	if(!stage) { THROW(EAny, "Module " + module + " unknown or failed to load."); }

	/**
		Convenience: allow the user to override the default input/output
		specified in input/output module configuration files from the
		command line.
	*/
	if(stage->type() == "input")
	{
		if(!input.empty()) { modcfg.insert(make_pair("filename", input)); } // override only if explicitly given
		modcfg.insert(make_pair("model", models));
		modcfg.insert(make_pair("foot", foots));
		modcfg.insert(make_pair("maxstars", str(maxstars)));
		modcfg.insert(make_pair("nstars", str(nstars)));
		modcfg.insert(make_pair("extmaps", extmaps));
		modcfg.insert(make_pair("dryrun", str(this->dryrun)));
	}
	if(stage->type() == "output" && !output.empty())
	{
		modcfg.insert(make_pair("filename", output));
	}

	// allow the module to construct itself from the command line
	if(!stage->construct(modcfg, t, *this)) { THROW(EAny, "Failed to initialize module '" + module + "'"); }
	DLOG(verb2) << "module loaded: " << module << " (type: " << stage->type() << ")";

	add(stage);
}

void apply_definitions(const std::string &def_modules)
{
	if(def_modules.empty()) { return; }

	// concatenate all files
	std::string deffn, merged_config_text;
	std::istringstream in(def_modules.c_str());
	while(in >> deffn)
	{
		std::ifstream cfin(deffn.c_str());
		std::string line;
		while(std::getline(cfin, line))
		{
			merged_config_text += line + "\n";
		}
	}

	// load them as config (and expand any variables)
	std::istringstream cfgstrm(merged_config_text);
	Config cfg;
	cfg.load(cfgstrm);

	// set environment variables corresponding to every entry
	// in the config file
	FOREACH(cfg)
	{
		const std::string var = i->first, value = i->second;
		//std::cerr << "POTENTIAL ENVVAR [" << var << "] = [" << value << "]\n";
		if(var == "module") { continue; }

		EnvVar(var).set(value);
	}
}

void generate_catalog(int seed, size_t maxstars, size_t nstars, const std::set<std::string> &modules, const std::string &input, const std::string &output, bool dryrun)
{
	rng_gsl_t rng(seed);

	// find and load (== set corresponding EnvVars) definitions from
	// definition module(s) before doing anything else
	std::string defs;
	FOREACH(modules)
	{
		const std::string &cffn = *i;

		if(!file_exists(cffn))
		{
			THROW(EAny, "Module configuration file " + cffn + " is inaccessible or doesn't exist.");
		}

		// load from file
		Config cfg;
		cfg.load(cffn, false);
		std::string module = normalizeKeyword(cfg["module"]);

		if(module == "definitions")
		{
			if(!defs.empty()) { defs += " "; }
			defs += cffn;
		}
	}
	apply_definitions(defs);

	// output table setup
	// HACK: Kbatch should be read from skygen.conf, or auto-computed to maximize memory use otherwise
	size_t Kbatch = 5000000;
	EnvVar kb("KBATCH");
	if(kb) { Kbatch = (int)atof(kb.c_str()); } // atof instead of atoi to allow shorthands such as 1e5
	DLOG(verb1) << "Processing in batches of " << Kbatch << " objects";
	otable t(Kbatch);

	// load all module config files and detect and set aside
	// models and footprints
	std::list<Config> module_configs;
	std::string models, foots, extmaps;
	FOREACH(modules)
	{
		const std::string &cffn = *i;

		if(!file_exists(cffn))
		{
			THROW(EAny, "Module configuration file " + cffn + " is inaccessible or doesn't exist.");
		}

		// load from file
		Config modcfg(cffn);

		if(!modcfg.count("module")) { THROW(EAny, "Configuration file " + cffn + " does not specify the module name"); }
		std::string module = normalizeKeyword(modcfg["module"]);

		// set aside "special" modules
		if(module == "model")
		{
			if(!models.empty()) { models += " "; }
			models += cffn;
		}
		else if(module == "footprint")
		{
			if(!foots.empty()) { foots += " "; }
			foots += cffn;
		}
		else if(module == "extinction")
		{
			if(!extmaps.empty()) { extmaps += " "; }
			extmaps += cffn;
		}
		else if(module == "definitions")
		{
			// do nothing
		}
		else
		{
			module_configs.push_back(modcfg);
		}
	}
	MLOG(verb1) << "Definitions: " << (defs.empty() ? "<none>" : defs);
	MLOG(verb1) << "Models: " << (models.empty() ? "<none>" : models);
	MLOG(verb1) << "Footprints: " << (foots.empty() ? "<none>" : foots);
	MLOG(verb1) << "Extinction maps: " << (extmaps.empty() ? "<none>" : extmaps);

	// Create the modules and construct the pipeline
	opipeline pipe(dryrun);
	FOREACH(module_configs)
	{
		pipe.create_and_add(*i, t, maxstars, nstars, models, foots, extmaps, input, output);
	}

	// Add default I/O modules, if no I/O modules were found above
	if(!pipe.has_module_of_type("input"))
	{
		Config modcfg;
		modcfg.insert(std::make_pair("module", "textin"));
		pipe.create_and_add(modcfg, t, maxstars, nstars, models, foots, extmaps, input, output);
	}
	if(!pipe.has_module_of_type("output"))
	{
		Config modcfg;
		modcfg.insert(std::make_pair("module", "textout"));
		pipe.create_and_add(modcfg, t, maxstars, nstars, models, foots, extmaps, input, output);
	}

	// execute the pipeline
	int nstarsGenerated = pipe.run(t, rng);
}
