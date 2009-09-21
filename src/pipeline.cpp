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

		std::string factory_name = "create_module_";
		FOREACH(name)
		{
			if(!isalnum(*i)) { continue; }
			factory_name += tolower(*i);
		}

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
