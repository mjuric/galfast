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

#include "galfast_config.h"
#include "galfast_version.h"

#include "analysis.h"
#include "io.h"
#include "pipeline.h"
#include "projections.h"

#include <fstream>

#include <astro/io/format.h>
#include <astro/system/log.h>
#include <astro/system/config.h>
#include <astro/useall.h>

#include <boost/thread.hpp>

///////////////////////////////////////////////////////////////////////
// Component map implementation

componentMap_t componentMap;

// return an index corresponding to compID, creating one if it doesn't exist
uint32_t componentMap_t::seqIdx(uint32_t compID)
{
	// check if we already know of this component
	if(!comp2seq.count(compID))
	{
		comp2seq[compID] = seq2comp.size();
		seq2comp.push_back(compID);

		DLOG(verb1) << "Mapped compID=" << compID << " to seqIdx=" << comp2seq[compID] << "\n";
	}

	return comp2seq[compID];
}

uint32_t componentMap_t::compID(uint32_t seqIdx)
{
	if(seqIdx >= seq2comp.size()) { return 0xffffffff; }
	return seq2comp[seqIdx];
}

// convert a list of (closed!) ranges to a bitmap, and upload it to the GPU
bit_map::bit_map(const interval_list &cl)
{
	// set all bits in the bit_map that are in any of the intervals
	// in the interval list. We assume the intervals are closed
	// (that is, [from,to]).
	set_all(0);
	FOREACH(cl)
	{
		std::pair<uint32_t, uint32_t> iv = *i;
		FOREACHj(comp, componentMap.comp2seq)
		{
			if(iv.first > comp->first || comp->first > iv.second) { continue; }
			set(comp->second);
		}
	}

	//FOR(0, bit_map::maxbits()) { std::cerr << (isset(i) ? "1" : "0"); if((i+1)%10 == 0) std::cerr << " "; } std::cerr << "\n";
}

///////////////////////////////////////////////////////////////////////

int opipeline_stage::instanceId()
{
	// FIXME: This function should really be made const, and assignment of component Id should somehow be moved to the constructor
	// FIXME: This is not thread safe (not that it has to be, but just to keep it in mind...)

	if(m_instanceId == -1)
	{
		// construct-on-first-use idiom
		static std::map<std::string, int> instanceIdMap;
		m_instanceId = instanceIdMap[name()]++;
	}

	return m_instanceId;
}

std::string opipeline_stage::instanceName()
{
	std::ostringstream ss;
	ss << name() << "[" << instanceId() << "]";
	return ss.str();
}

bool opipeline_stage::runtime_init(otable &t)
{
//	// components we care about
//	bit_map comps = componentMap.get(applyToComponents);

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

void opipeline_stage::read_component_map(interval_list &applyToComponents, const peyton::system::Config &cfg, const std::string &compCfgKey, uint32_t compFirst, uint32_t compLast)
{
	applyToComponents.clear();

	// new, flexible, syntax
	if(!compCfgKey.empty() && cfg.count(compCfgKey))
	{
		// parse the line of the form
		// comp = 1,2, 4-8, 9-12
		std::string line = cfg.get(compCfgKey);
		int at = 0;
		do {
			int range_sep = 0;
			int at0 = at;
			while(at != line.size() && line[at] != ',')
			{
				if(line[at] != '-') { at++; continue; }
				if(range_sep || at0 == at) { THROW(EAny, "Syntax error in the value of 'applyToComponents' key: " + line); }
				range_sep = at;
				at++;
			}

			if(range_sep)
			{
				compFirst = (uint32_t)atof(line.substr(at0, at - range_sep).c_str());
				range_sep++;
				compLast =   (uint32_t)atof(line.substr(range_sep, at - range_sep).c_str());
			}
			else
			{
				compFirst = compLast = (uint32_t)atof(line.substr(at0, at - at0).c_str());
			}

			applyToComponents.push_back(std::make_pair(compFirst, compLast));

			at++;
		} while(at < line.size());
	}
	else
	{
#if 0
		// backwards compatibility
		if(compCfgKey.empty())
		{
			// backwards-compatible syntax, where keys are allowed to be missing
			if(cfg.count("compFirst"))
			{
				cfg.get(compFirst, "compFirst", 0U);
				cfg.get(compLast, "compLast", 0xffffffff);
			}
			applyToComponents.push_back(std::make_pair(compFirst, compLast));
		}
		else
		{
#endif
			// explode otherwise -- the key _must_ exist
			if(!cfg.count(compCfgKey)) { THROW(EAny, "Value for " + compCfgKey + " missing in configuration file."); }
#if 0
		}
#endif
	}
}

/////////////////////////////

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
		//virtual int priority() { return PRIORITY_OUTPUT; }	// ensure this stage has the least priority
		virtual double ordering() const { return ord_output; }
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

void osink::transformComponentIds(otable &t, size_t begin, size_t end)
{
	// transform 'comp' column from seqIdx to compID
	cint_t::host_t comp = t.col<int>("comp");
	FORj(row, begin, end)
	{
		int cmp = comp(row);
		comp(row) = componentMap.compID(cmp);
	}
}

size_t os_textout::process(otable &t, size_t from, size_t to, rng_t &rng)
{
	swatch.start();

	ticker tick("Writing output", (int)ceil((to-from)/50.));

	transformComponentIds(t, from, to);

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

class os_countsMap : public osink
{
	protected:
		flex_output out;
		std::string xy_column, z_column;		// column names for xgitude/yitude/distance
		float x0, dx;
		float y0, dy;
		float z0, dz;
		int n_x, n_y, n_z,	// the number of x, y, z bins
		    z_mag_width,	// the width of the z vector (if z is a scalar, 1)
		    z_width;		// == z_mag_width + coadd_coeffs.size() (the number of elements in the output record)
		int dense_output;	// whether to output all bins, or only those containing the data (boolean)

		std::vector<std::vector<float> > coadd_offsets;	// coefficients used to compute coadd depths
		std::vector<std::string>         coadd_names;   // the names of the coadds (used for headings in output)

		struct beam
		{
			int X, Y, map;
			beam(int X_ = 0, int Y_ = 0, int map_ = 0) : X(X_), Y(Y_), map(map_) {}
			bool operator <(const beam &a) const
			{
				return X < a.X ||
					X == a.X && Y < a.Y ||
					X == a.X && Y == a.Y && map < a.map;
			}
		};
		typedef std::map<beam, std::vector<int> > counts_map_t;
		counts_map_t countsX;
		std::vector<int> &get_z_array(counts_map_t &counts, const beam &b) const;

		long long m_total;
		bool m_equalarea;

		typedef peyton::math::lambert lambert;
		lambert proj[2];

		void get_columns(otable &t, cdouble_t::host_t &xy, cfloat_t::host_t &z, cint_t::host_t &hidden) const;
		size_t threaded_process(counts_map_t &counts, long long &total, 
			cdouble_t::host_t xy, cfloat_t::host_t z, cint_t::host_t hidden,
			size_t from, size_t to, rng_t &rng, int startoffs, int stride) const;
		friend class mt_binner;

		void create_dense_output_array();
	public:
		virtual bool construct(const Config &cfg, otable &t, opipeline &pipe);
		virtual bool runtime_init(otable &t);
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		//virtual int priority() { return PRIORITY_OUTPUT; }	// ensure this stage has the least priority
		virtual double ordering() const { return ord_output; }
		virtual const std::string &name() const { static std::string s("countsMap"); return s; }
		virtual const std::string &type() const { static std::string s("output"); return s; }

		~os_countsMap();
		os_countsMap() : osink()
		{
			m_total = 0;
			dense_output = 0;
			proj[0] = lambert(rad(90), rad(90));
			proj[1] = lambert(rad(-90), rad(-90));
		}
};

extern "C" opipeline_stage *create_module_countsmap() { return new os_countsMap; }

// Bin n bins of width dx, with the first bin centered on x0+dx/2
// Points that fall below/above the range are binned into bin 0 and n-1, respectively
//   (similar to how SM does it)
inline int find_bin(double x, double x0, double dx, int n)
{
	int b = (int)((x - x0) / dx + 0.5);	// find the bin
	b = std::max(0, b);
	b = std::min(n - 1, b);
	return b;
}

void os_countsMap::create_dense_output_array()
{
	// pre-initialize the output map with zeros, when we want
	// the output to contain them even if there's no data
	// (requested by ZI so that he can easily add results from
	// different simulations in SM)

	countsX.clear();

	beam b;
	for(b.X = 0; b.X != n_x; b.X++)
	{
		for(b.Y = 0; b.Y != n_y; b.Y++)
		{
			// this will implicitly initialize the array to zero
			get_z_array(countsX, b);
		}
	}
}	

std::vector<int> &os_countsMap::get_z_array(counts_map_t &counts, const beam &b) const
{
	std::vector<int> &c = counts[b];
	bool first = c.empty();
	if(c.empty())
	{
		c.resize(n_z*z_width);
		memset(&c[0], 0, sizeof(*&c[0])*c.size());
	}
	return c;
}

void os_countsMap::get_columns(otable &t, cdouble_t::host_t &xy, cfloat_t::host_t &z, cint_t::host_t &hidden) const
{
	// get the needed columns
	xy = t.col<double>(xy_column);
	z   = t.col<float>(z_column);

	// find out if we're using the 'hidden' flag column (TODO: this column should be made mandatory to avoid each output module having to check for its existence)
	if(t.using_column("hidden"))
	{
		hidden = t.col<int>("hidden");
	}
	else
	{
		hidden.reset();
	}
};

size_t os_countsMap::threaded_process(counts_map_t &counts, long long &total, 
	cdouble_t::host_t xy,
	cfloat_t::host_t z,
	cint_t::host_t hidden,
	size_t from, size_t to, rng_t &rng, int startoffs, int stride) const
{
	// bin
	size_t nserialized = 0;
	for(size_t row = from + startoffs; row < to; row += stride)
	{
		if(hidden && hidden(row)) { continue; }

		beam b;
		double x = xy(row, 0);
		double y = xy(row, 1);
		if(m_equalarea)
		{
			// lambert-project the coordinates to north/south hemispheres
			// and set the output map accordingly
			b.map = y > 0 ? 0 : 1;
			const lambert &proj = this->proj[b.map];
			proj.project(x, y, rad(x), rad(y));
			
		}

		b.X = find_bin(x, x0, dx, n_x);
		b.Y = find_bin(y, y0, dy, n_y);

		// fetch the correct bin (note there are n_d*d_width elements of the pencil beam,
		// where d_width is usually the number of bands)
		std::vector<int> &c = get_z_array(counts, b);

		// bin
		for(int i = 0; i != z_mag_width; i++)
		{
			int Z = find_bin(z(row, i), z0, dz, n_z);
			int idx = Z*z_width + i;
			assert(idx >= 0 && idx < c.size());
			c[idx]++;
			total++;
		}
		
		// compute coadd bins, if any
		FORj(coadd, 0, coadd_offsets.size())
		{
			const std::vector<float> &coeffs = coadd_offsets[coadd];
			assert(coeffs.size() >= z_mag_width);

			// Find the maximum r-band depth to which this object is detected in _any_ band
			float magCut = std::numeric_limits<float>::infinity();
			for(int i = 0; i != z_mag_width; i++)
			{
				float mag = z(row, i) - coeffs[i];
//				std::cerr << z(row, i) << " -> " << mag << "\n";
				magCut = std::min(magCut, mag);
			}

			// This object is to be measured in all stacks deeper than magCut
			int bin = find_bin(magCut, z0, dz, n_z);
//			std::cerr << bin << " " << z0 << " " << dz << "\n";
			for(int Z = bin; Z < n_z; Z++)
			{
				int idx = Z*z_width + z_mag_width + coadd;
				assert(idx >= 0 && idx < c.size());
				c[idx]++;
			}

			// ---
//			std::cerr << "magCut = " << magCut << "\n";
//			for(int i = 0; i != z_mag_width; i++)
//			{
//				std::cerr << z(row, i) << " ";
//			}
//			std::cerr << "\n";
		}

		nserialized++;
	}

	return nserialized;
}

struct mt_binner
{
	const os_countsMap &ctmap;
	
	// will fill out and returns
	os_countsMap::counts_map_t counts;
	long long m_total;
	size_t nserialized;

	// input parameters
	size_t from, to;
	rng_t &rng;
	int start, stride;
	cdouble_t::host_t xy;
	cfloat_t::host_t z;
	cint_t::host_t hidden;

	mt_binner(const os_countsMap &ctmap_, otable &t_, size_t from_, size_t to_, rng_t &rng_, int start_, int stride_)
		: ctmap(ctmap_), from(from_), to(to_), rng(rng_), start(start_), stride(stride_)
	{
		m_total = 0;
		nserialized = 0;
		ctmap.get_columns(t_, xy, z, hidden);
	}

	void operator()()
	{
		nserialized = ctmap.threaded_process(counts, m_total, xy, z, hidden, from, to, rng, start, stride);
	}
};

size_t os_countsMap::process(otable &t, size_t from, size_t to, rng_t &rng)
{
	swatch.start();

#if 1
	// launch NCORE threads to do the binning
	typedef boost::shared_ptr<mt_binner> mtbp_t;
	std::vector<mtbp_t> kernels(boost::thread::hardware_concurrency());
	boost::thread_group threads;

	// bin in threads
	FOR(0, kernels.size())
	{
		kernels[i] = mtbp_t(new mt_binner(*this, t, from, to, rng, i, kernels.size()));
		threads.create_thread(boost::ref(*kernels[i]));
	}
	threads.join_all();
	
	// accumulate the results
	int nserialized = 0;
	FOREACHj(kp, kernels)
	{
		mt_binner &k = **kp;
		FOREACH(k.counts)
		{
			std::vector<int> &c = get_z_array(countsX, i->first);
			FORj(j, 0, c.size())
			{
				c[j] += i->second[j];
			}
		}
		m_total += k.m_total;
		nserialized += k.nserialized;
	}
#else
	cdouble_t::host_t xy;
	cfloat_t::host_t z;
	cint_t::host_t hidden;
	get_columns(t, xy, z, hidden);
	int nserialized = threaded_process(countsX, m_total, xy, z, hidden, from, to, rng, 0, 1);
#endif
	swatch.stop();

	return nserialized;
}

bool os_countsMap::runtime_init(otable &t)
{
	if(!osink::runtime_init(t)) { return false; }

	z_mag_width = t.col<float>(z_column).width();
	z_width = z_mag_width + coadd_offsets.size();

	if(dense_output)
		create_dense_output_array();

	return true;
}

bool os_countsMap::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// output filename
	std::string fn = cfg.get("filename");
	out.open(fn);
	MLOG(verb1) << "Output file: " << fn << " (text)\n";

	float tmp;

	// see if equal-area output has been requested
	cfg.get(m_equalarea, "equalarea", false);
	double x1, y1, z1;
	if(m_equalarea)
	{
		// reasonable defaults for lambert-projected x/y
		x0 = -2; x1 = 2; dx = rad(1.);
		y0 = -2; y1 = 2; dy = rad(1.);
	}
	else
	{
		// reasonable defaults for full-sky x/y
		x0 =   0; x1 = 360; dx = 1;
		y0 = -90; y1 =  90; dy = 1;
	}
	// reasonable defaults for z-coordinate (usually the magnitude)
	z0 = 0; z1 = 40; dz = 1;

	xy_column = cfg.get("xy");
	req.insert(xy_column);

	// x coordinate (usually longitude) range setup
	cfg.get(x0, "x0", x0);
	cfg.get(x1, "x1", x1);
	cfg.get(dx, "dx", dx);
	tmp = (x1 - x0) / dx;
	n_x = (int)(fabs(tmp - round(tmp)) < 1e-4 ? round(tmp) : ceil(tmp)) + 1; // adding 1 to ensure x1 is fully covered (because bins are centered on x0)

	// y coordinate (usually latitude) range setup
	cfg.get(y0, "y0", y0);
	cfg.get(y1, "y1", y1);
	cfg.get(dy, "dy", dy);
	tmp = (y1 - y0) / dy;
	n_y = (int)(fabs(tmp - round(tmp)) < 1e-4 ? round(tmp) : ceil(tmp)) + 1;

	// z coordinate (usually magnitude) range setup
	z_column = cfg.get("z");
	req.insert(z_column);
	cfg.get(z0, "z0", z0);
	cfg.get(z1, "z1", z1);
	cfg.get(dz, "dz", dz);

	tmp = (z1 - z0) / dz;
	n_z = (int)(fabs(tmp - round(tmp)) < 1e-4 ? round(tmp) : ceil(tmp)) + 1;

	// coadds setup
	std::set<std::string> skeys;
	cfg.get_matching_keys(skeys, "coadd\\..*");
	FOREACHj(key, skeys)
	{
		std::vector<float> offsets;
		std::istringstream ss(cfg[*key]);

		float coeff;
		while(ss >> coeff)
		{
			offsets.push_back(coeff);
		}

		coadd_names.push_back(*key);
		coadd_offsets.push_back(offsets);
	}
	MLOG(verb1) << "Num. coadd columns: " << coadd_names.size() << "\n";

	cfg.get(dense_output, "output all bins", dense_output);
	if(dense_output)
		MLOG(verb1) << "Output type: Writing out all bins (dense output)\n";
	else
		MLOG(verb1) << "Output type: Writing out bins with data (sparse output)\n";
//	if(coadd_names.size())
//	{
//		FOR(0, coadd_names.size())
//		{
//			std::cerr << coadd_names[i] << ":";
//			FOREACHj(j, coadd_offsets[i])
//			{
//				std::cerr << *j << " ";
//			}
//			std::cerr << "\n";
//		}
//	}

	return out.out();
}

os_countsMap::~os_countsMap()
{
	// header
	if(m_equalarea)
	{
		out.out() << "# Binned in lambert equal area projection, using the following poles:\n";
		for(int i = 0; i != 2; i++)
		{
			out.out() << "#\tmap = " << i << "   :  l0 = " << deg(proj[i].l0) << ",  phi1 = " << deg(proj[i].phi1) << "\n";
		}
		out.out() << "#\n";
		out.out() << "# map\tlam_X\tlam_Y\t" << z_column;
	}
	else
	{
		out.out() << "# map\t" << xy_column << "[0]\t" << xy_column << "[1]\t" << z_column;
	}
	for(int k = 0; k != z_mag_width; k++) { out.out() << "\tN(" << z_column << "[" << k << "])"; }
	FOREACH(coadd_names) { out.out() << "\t" << *i; }
	if(m_equalarea)
	{
		out.out() << "\t" << xy_column << "[0]\t" << xy_column << "[1]";
	}
	else
	{
		out.out() << "\t" << "dA";
	}
	out.out() << "\n#\n";

	// write out the (sparse) binned array into the text output file
	long long total = 0;
	FOREACHj(XY, countsX)
	{
		double x = x0 + dx * XY->first.X;
		double y = y0 + dy * XY->first.Y;
		int map = XY->first.map;
		std::vector<int> &zbins = XY->second;
		double lon, lat;
		if(m_equalarea)
		{
			proj[map].deproject(lon, lat, x, y);
		}
		for(int i = 0; i != zbins.size(); i += z_width)
		{
			bool hasDatum = false;
			for(int k = 0; k != z_width; k++)
			{
				hasDatum |= zbins[i+k] != 0;
			}
//			if(!hasDatum) { continue; }

			float z = z0 + dz * ((float)i / z_width);
			out.out() << map << "\t" << x << "\t" << y << "\t" << z;
			for(int k = 0; k != z_width; k++)
			{
				out.out() << "\t" << zbins[i+k];
				if(k < z_mag_width) { total += zbins[i+k]; }
			}
			if(m_equalarea)
			{
				// print out deprojected coordinates
				out.out() << "\t" << deg(lon) << "\t" << deg(lat);
			}
			else
			{
				// print out the area of the pixel, assuming xy are spherical longitude and latitutde, in degrees
				// Exact formula: area = (lon2-lon1) * (sin(lat2) - sin(lat1))
				double area = dx * deg( sin(rad(std::min(y + 0.5*dy, 90.))) - sin(rad(std::max(y - 0.5*dy, -90.))) );
				out.out() << "\t" << area;
			}
			out.out() << "\n";
		}
	}
	// simple sanity check
	if(total != m_total)
	{
		THROW(EAny, "Bug: total != m_total (" + str(total) + " != " + str(m_total) + "). Notify the authors.");
	}
}

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
		//virtual int priority() { return PRIORITY_OUTPUT; }	// ensure this stage has the least priority
		virtual double ordering() const { return ord_output; }
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

	transformComponentIds(t, from, to);
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
			const char *ttype = "HEADER";
			char *tform;
			asprintf(&tform, "%dA", len);
			fits_create_tbl(fptr, BINARY_TBL, 0, 1, (char**)&ttype, &tform, NULL, "METADATA", &status);
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
	else if(name == "fitsout") { s.reset(new os_fitsout); }
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
	std::set<std::pair<double, opipeline_stage *> > stages;
	FOREACH(this->stages)
	{
		stages.insert(std::make_pair((*i)->ordering(), (*i).get()));
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
			//std::cerr << s.name() << "\n";

			// initialize this pipeline stage (this typically adds and uses the columns
			// this stage will add)
			if(!s.runtime_init(t)) { /*continue;*/ std::cerr << s.name() << "\n"; assert(0); } // TODO: I switched from dependency tracking, to explicit ordering

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
	std::vector<boost::shared_ptr<std::stringstream> > ss(componentMap.size());
	FOREACH(ss) { i->reset(new std::stringstream); }
	FOREACH(pipeline)
	{
		if(source == NULL) { last = source = *i; FOREACH(ss) { **i << last->name(); }; continue; }
		osink *next = dynamic_cast<osink*>(*i);
		ASSERT(next);
		last->chain(next);
		last = next;

		bit_map comps = last->getAffectedComponents();
		std::string name = last->instanceName();
		std::string empty(name.size(), ' ');
		FOR(0, ss.size())
		{
			std::string res = (comps.isset(i) ? name : empty);
			*ss[i] << " | " << res;
		}
	}
	FOREACH(componentMap.comp2seq)
	{
		MLOG(verb1) << "Pipeline [comp=" << i->first << "] : " << ss[i->second]->str() << "\n";
	}

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

#if 1
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

	// load them as global config (and expand any variables)
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
#endif
}

void generate_catalog(int seed, size_t maxstars, size_t nstars, const std::set<Config::filespec> &modules, const std::string &input, const std::string &output, bool dryrun)
{
	rng_gsl_t rng(seed);

	// find and load into Config::globals all definitions from
	// definition module(s) before doing anything else
	std::string defs;
	FOREACH(modules)
	{
		const std::string &cffn = *i;
		const std::string fn = cffn.substr(0, cffn.find('{'));

		if(!file_exists(fn))
		{
			THROW(EAny, "Module configuration file " + fn + " is inaccessible or doesn't exist.");
		}

		// load from file
		Config cfg;
		cfg.load(cffn, false);
		std::string module = normalizeKeyword(cfg["module"]);

		if(module == "definitions")
		{
			// load into globals, without expanding the variables
			Config::globals.load(cffn, false);
			Config::globals.erase("module");
			
			// this is just informational (printed out for the user, below...)
			if(!defs.empty()) { defs += " "; }
			defs += cffn;
		}
	}
	Config::globals.expandVariables();
// 	apply_definitions(defs);

	// output table setup
	// HACK: Kbatch should be read from skygen.conf, or auto-computed to maximize memory use otherwise
	size_t Kbatch = 5000000;
	EnvVar kb("KBATCH");
	if(kb) { Kbatch = (int)atof(kb.c_str()); } // atof instead of atoi to allow shorthands such as 1e5
	if(kb) { MLOG(verb1) << "Batch size: " << Kbatch << " objects"; }
	DLOG(verb1) << "Processing in batches of " << Kbatch << " objects";

	std::string ver;
	EnvVar nover("NOVERSION");
	if(!nover || atoi(nover.c_str()) == 0) { ver = VERSION_STRING; }
	else { MLOG(verb1) << "WARNING: Not recording software version in output files."; }
	otable t(Kbatch, ver);

	// load all module config files and detect and set aside
	// models and footprints
	std::list<Config> module_configs;
	std::string models, foots, extmaps;
	FOREACH(modules)
	{
		const std::string &cffn = *i;
		const std::string fn = cffn.substr(0, cffn.find('{'));

		if(!file_exists(fn))
		{
			THROW(EAny, "Module configuration file " + fn + " is inaccessible or doesn't exist.");
		}

		// load from file
		Config modcfg(cffn);

		#if 0
		std::cout << "# ======== " << cffn << "\n";
		FOREACH(modcfg)
		{
			if(Config::globals.count(i->first)) { continue; }
			std::cout << i->first << " = " << i->second << "\n";
		}
		#endif

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
	MLOG(verb2) << "Definitions: " << (defs.empty() ? "<none>" : defs);
	MLOG(verb2) << "Models: " << (models.empty() ? "<none>" : models);
	MLOG(verb2) << "Footprints: " << (foots.empty() ? "<none>" : foots);
	MLOG(verb2) << "Extinction maps: " << (extmaps.empty() ? "<none>" : extmaps);

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
