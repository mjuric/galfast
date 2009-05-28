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

#include "observe.h"
#include "io.h"
#include "skygen.h"
#include "analysis.h"
#include "gpc_cpp.h"
#include <iomanip>
#include <fstream>
#include <vector>

#include <astro/useall.h>

/***********************************************************************/

lfParams lfTextureManager::loadConstant(float val)
{
	const int nlf = 2;
	float lfp[nlf] = { val, val };
	float lfM0 = -100, lfM1 = 100, lfdM = (lfM1 - lfM0) / (nlf-1);
	return set(&lfp[0], nlf, lfM0, lfM1, lfdM);
}

lfParams lfTextureManager::load(const char *fn)
{
	// load the luminosity function and normalize to m.rho0.
	// rho0 is assumed to contain the number of stars per cubic parsec
	// per 1mag of r-i
	text_input_or_die(lfin, fn);
	std::vector<double> M, phi;
	::load(lfin, M, 0, phi, 1);
	spline lf;
	lf.construct(M, phi);

	// resample the LF to texture
	const int nlf = 1024;
	std::vector<float> lfp(nlf);
	float lfM0 = M.front(), lfM1 = M.back(), lfdM = (lfM1 - lfM0) / (nlf-1);
	for(int i=0; i != nlf; i++)
	{
		lfp[i] = lf(lfM0 + i*lfdM);
	}
	return set(&lfp[0], nlf, lfM0, lfM1, lfdM);
}

void expModel::load(const peyton::system::Config &cfg)
{
	// density distribution parameters
	rho0  = cfg.get("rho0");
	l     = cfg.get("l");
	h     = cfg.get("h");
	z0    = cfg.get("z0");
	f     = cfg.get("f");
	lt    = cfg.get("lt");
	ht    = cfg.get("ht");
	fh    = cfg.get("fh");
	q     = cfg.get("q");
	n     = cfg.get("n");

	// luminosity function
	if(cfg.count("lumfunc"))
	{
		lf = texLFMgr.load(cfg["lumfunc"].c_str());
	}
	else
	{
		lf = texLFMgr.loadConstant(1.f);
	}

	// cutoff radius (default: 1Mpc)
	cfg.get(r_cut2,  "rcut",   1e6f);
	r_cut2 *= r_cut2;
}

/***********************************************************************/

template<typename T>
skyConfig<T>::skyConfig()
{
	cpu_pixels = NULL;
	cpu_hist = NULL;
	cpu_maxCount = NULL;
	cpu_state = NULL;
	rng = NULL;

	this->pixels = 0;
	this->lock = 0;
	this->counts = 0;
	this->nstars = 0;

	this->ks.constructor();
}

template<typename T>
skyConfig<T>::~skyConfig()
{
	delete [] cpu_pixels;
	delete [] cpu_hist;
	delete [] cpu_maxCount;
	delete [] cpu_state;
	delete rng;

	this->lock.free();
	this->pixels.free();
	this->counts.free();
	this->nstars.free();

	this->ks.destructor();
}

template<typename T>
int skyConfig<T>::bufferSafetyMargin()
{
	// Compute the extra padding necessary to make buffer overfilling
	// exceedingly unlikely

	if(cpu_maxCount == NULL) { abort(); }

	float maxCount = std::accumulate(cpu_maxCount, cpu_maxCount + this->nthreads, 0.f);
	int bufsize = (int)ceil(maxCount + 7.*sqrt(maxCount))	// for maxCount > 100, P(k > k+7sqrt(k)) ~ nthreads*7e-11
			 + 3*this->nthreads;			// and this part is to cover the selection effect that we only break if k>1

#if 0
	// some diagnostic info
	sort(cpu_maxCount, cpu_maxCount + this->nthreads, std::greater<float>());
	for(int i = 0; i != 3 && i != this->nthreads; i++)
	{
		std::cerr << "rho" << i << " = " << cpu_maxCount[i] << "\n";
	}
	std::cerr << "Total of max density bins: " << maxCount << "\n";
	std::cerr << "Buffer safety margin: " << bufsize << "\n";
#endif

	return bufsize;
}

template<typename T>
void skyConfig<T>::upload_generic(bool draw)
{
	// Upload pixels to be processed
	this->pixels.upload(cpu_pixels, this->npixels);

	// prepare locks and outputs
	int zero = 0;

	if(!draw)
	{
		this->lock.upload(&zero, 1);

		this->rhoHistograms.alloc(this->nthreads*this->nhistbins);
		cudaMemset(this->rhoHistograms.ptr, 0, this->nthreads*this->nhistbins*4);

		this->maxCount.alloc(this->nthreads);
		cudaMemset(this->maxCount.ptr, 0, this->nthreads*4);

		this->counts.alloc(this->nthreads);
		cudaMemset(this->counts.ptr, 0, this->nthreads*4);
	}
	else
	{
		this->nstars.upload(&zero, 1);

		this->stopstars = Kbatch - bufferSafetyMargin();
		assert(this->stopstars > 0);
	}

	cuxUploadConst("Rg_gpu", this->Rg);
	cuxUploadConst("rng", *this->rng);
	cuxUploadConst("proj", this->proj);
}

template<typename T>
void skyConfig<T>::upload(bool draw)
{
	assert(0); // Generic version of upload should _never_ be called
}

template<>
void skyConfig<expModel>::upload(bool draw)
{
	this->upload_generic(draw);
	cuxUploadConst("expModelSky", (skyConfigGPU<expModel>&)*this); // upload thyself (note: this MUST come last)
}

template<typename T>
void skyConfig<T>::loadConfig()
{
	// Galactic center distance (in pc)
// 	Rg = 8000.;

	// Galactic model initialization
//	this->model.setmodel(2150., 245., 25,  0.13, 3261, 743,   0.0051, 1./0.64, 2.77);
//	this->model.setmodel(2150., 300., 25.,  0.13, 2150, 900,   0.0051, 1./0.64, 2.77);

	// luminosity function initialization
// 	               //  0      4      8    12    16    20   24
// 	const int nlf = 7;
// 	float cpu_lf[] = { 0.00, 0.03, 0.05, 0.07, 0.04, 0.02, 0.00 };
// 	float lfM0 = 0, lfM1 = 24, lfdM = (lfM1 - lfM0) / (nlf-1);
// 	this->model.lf = texLFMgr.set(cpu_lf, nlf, lfM0, lfM1, lfdM);

	// density histogram setup
// 	this->lrho0 = -3.5f;
// 	this->dlrho = 1.0f;

	// integration limits
/*	float deltam = 0.01;*/
//	this->M0 = 10.; this->M1 = 10.05;  this->dM = deltam;
// 	this->M0 =  4.; this->M1 = 15.;  this->dM = deltam;
// 	this->m0 = 15.; this->m1 = 24.;  this->dm = deltam;
// 	this->nM = (int)round((this->M1 - this->M0) / this->dM);
// 	this->nm = (int)round((this->m1 - this->m0) / this->dm);
// 	this->m0 += 0.5*this->dm;
// 	this->M1 -= 0.5*this->dM;

	// generate sky pixels
// 	this->proj.init(rad(90.), rad(90.));
// 	this->npixels = 1; //2000;
// 	this->dx = rad(1.);
// 	this->dA = sqr(this->dx);
// 	cpu_pixels = new direction[this->npixels];
// 	Radians l0 = rad(44), b0 = rad(85.), dx = rad(1), l, b;
// //	Radians l0 = rad(0.), b0 = rad(5.), dx = rad(1), l, b;
// 	for(int i=0; i != this->npixels; i++)
// 	{
// 		l = l0 + i*dx;
// 		b = b0;
// 		cpu_pixels[i] = direction(l, b);
// 	}

	// compute/set kernel execution parameters
/*	blockDim.x = 64; //256;
	gridDim.x = 120; // 30;
	this->nthreads = blockDim.x * blockDim.y * blockDim.z * gridDim.x * gridDim.y * gridDim.z;
	std::cout << "nthreads=" << this->nthreads << "\n";
	shb = 12 * blockDim.x; // for RNG*/
}

gpc_polygon makeBeamMap(const peyton::system::Config &cfg, peyton::math::lambert &proj);
gpc_polygon makeRectMap(const peyton::system::Config &cfg, peyton::math::lambert &proj);

gpc_polygon make_footprint(const Config &cfg, peyton::math::lambert &proj)
{
	if(!cfg.count("type")) { THROW(EAny, "The footprint configuration file must contain the 'type' keyword."); }

	if(cfg["type"] == "beam") { return makeBeamMap(cfg, proj); }
	if(cfg["type"] == "rect") { return makeRectMap(cfg, proj); }

	THROW(EAny, "The requested footprint type " + cfg["type"] + " is unknown.");
}

// Calculate the 'skymap', for fast point-in-polygon lookups of the observed
// sky footprint.
partitioned_skymap *make_skymap(Radians dx, gpc_polygon sky)
{
	using boost::shared_ptr;
	using std::cout;

	partitioned_skymap *skymap = new partitioned_skymap;
	skymap->dx = dx;

 	poly_bounding_box(skymap->x0, skymap->x1, skymap->y0, skymap->y1, sky);

	int X = 0;
	for(double x = skymap->x0; x < skymap->x1; x += skymap->dx, X++) // loop over all x values in the bounding rectangle
	{
		int Y = 0;
		double xa = x, xb = x+skymap->dx;
		for(double y = skymap->y0; y < skymap->y1; y += skymap->dx, Y++) // loop over all y values for a given x
		{
			double ya = y, yb = y+skymap->dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon poly;
			gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
			if(poly.num_contours == 0) continue; // if there are no observations in this direction

			// store the polygon into a fast lookup map
			partitioned_skymap::pixel_t &pix = skymap->skymap[std::make_pair(X, Y)];
			pix.poly = poly;
			pix.area = polygon_area(poly);
		}
	}
	return skymap;
}


template<typename T>
bool skyConfig<T>::init(
			const peyton::system::Config &cfg,
			const peyton::system::Config &foot_cfg,
			const peyton::system::Config &model_cfg,
			otable &t,
			opipeline &pipe)
{
	// Galactic center distance (in pc)
	Rg = 8000.;

#if 0
	// Footprint polygon computation
	cfg.get(this->dx,	"dx", 	1.f);
	this->dx = rad(this->dx);
	this->dA = sqr(this->dx);

	peyton::math::lambert proj;
	gpc_polygon sky = make_footprint(foot_cfg, proj);
	partitioned_skymap *skymap = make_skymap(this->dx, sky);

	this->npixels = skymap->skymap.size();
	cpu_pixels = new direction[this->npixels];
	this->proj.init(proj.l0, proj.phi1);
	int cnt = 0;
	FOREACH(skymap->skymap)
	{
		Radians x, y, l, b;
		x = skymap->x0 + skymap->dx*(i->first.first  + 0.5);
		y = skymap->y0 + skymap->dx*(i->first.second + 0.5);
		proj.inverse(x, y, l, b);
		cpu_pixels[cnt++] = direction(l, b);
	}
#else
	cfg.get(this->dx,	"dx", 	1.f);
	this->dx = rad(this->dx);
	this->dA = sqr(this->dx);

	boost::shared_ptr<opipeline_stage> foot( opipeline_stage::create("clipper") );
	Config fcfg(foot_cfg);
	fcfg.insert(std::make_pair("dx", cfg.get("dx")));
	if(!foot->construct(fcfg, t, pipe))
	{
		return false;
	}
	pipe.add(foot);

	std::vector<std::pair<double, double> > lb;	// pixels
	double l0, b0;					// projection pole
	this->npixels = ((os_clipper*)foot.get())->getPixelCenters(lb, l0, b0);
	this->proj.init(l0, b0);
	cpu_pixels = new direction[this->npixels];
	int cnt = 0;
	FOREACH(lb)
	{
		cpu_pixels[cnt++] = direction(i->first, i->second);
	}
#endif

	// Density binning parameters
	this->M0 = cfg.get_any_of("M0", "ri0");
	this->M1 = cfg.get_any_of("M1", "ri1");
	this->dM = cfg.get_any_of("dM", "dri");
	this->m0 = cfg.get_any_of("m0", "r0");
	this->m1 = cfg.get_any_of("m1", "r1");
	this->dm = cfg.get("dm");
	assert(this->dM == this->dm);

	this->nM = (int)round((this->M1 - this->M0) / this->dM);
	this->nm = (int)round((this->m1 - this->m0) / this->dm);
	this->m0 += 0.5*this->dm;
	this->M1 -= 0.5*this->dM;

	// catalog generation results
	t.use_column("lb");
	t.use_column("XYZ");
	t.use_column("comp");

	std::string 	colorName, absmagName, magName;
	colorName = model_cfg.get("color");
	magName = model_cfg.get("band");
	absmagName = "abs" + magName;
	assert(colorName == absmagName);

	t.use_column(colorName);  t.alias_column(colorName, "color");
	t.use_column(absmagName); t.alias_column(absmagName, "absmag");
	t.use_column(magName);	  t.alias_column(magName, "mag");

	this->model.load(model_cfg);

	// For debugging/stats
	this->lrho0 = -3.5f;
	this->dlrho = 1.0f;

	// GPU kernel execution setup
	blockDim.x = 64; //256;
	gridDim.x = 120; // 30;
	this->nthreads = blockDim.x * blockDim.y * blockDim.z * gridDim.x * gridDim.y * gridDim.z;
	std::cout << "nthreads=" << this->nthreads << "\n";
	shb = gpu_rng_t::state_bytes() * blockDim.x; // for RNG

	std::cout << "nm=" << this->nm << " nM=" << this->nM << " npixels=" << this->npixels << "\n";

	return true;
}

template<typename T>
size_t skyConfig<T>::run(otable &in, osink *nextlink, rng_t &cpurng)
{
	// initialize rng
	rng = new gpu_rng_t(cpurng);

	//
	// First pass: compute total expected starcounts
	//
	upload(false);
	compute(false);
	download(false);
	
	std::ostringstream ss;
	ss << "Histogram:     log(rho) |";
	for(int i=0; i != this->nhistbins; i++)
	{
		ss << std::setw(8) << this->lrho0+i*this->dlrho << "|";
	}
	MLOG(verb1) << ss.str();
	ss.str("");
	ss << "Histogram: rho=" << pow10f(this->lrho0 - this->dlrho*0.5) << "-->|";
	for(int i=0; i != this->nhistbins; i++)
	{
		ss << std::setw(8) << this->cpu_hist[i] << "|";
	}
	ss << "<--" << pow10f(this->lrho0+(this->nhistbins-0.5)*this->dlrho);
	MLOG(verb1) << ss.str();

	MLOG(verb1) << "Total expected star count: " << std::setprecision(9) << this->nstarsExpected;
	MLOG(verb1) << "Kernel runtime: " << this->swatch.getAverageTime();

	//
	// Second pass: draw stars
	//
	this->Kbatch = in.capacity();

	this->ks.alloc(this->nthreads);
	uint64_t total = 0;
	do
	{
		// setup output destination
		this->stars.lb    = in.col<double>("lb");
		this->stars.XYZ   = in.col<float>("XYZ");
		this->stars.comp  = in.col<int>("comp");
		this->stars.M     = in.col<float>("absmag");
		this->stars.m     = in.col<float>("mag");
		//this->stars.color = in.col<float>("color"); -- not implemented yet

		this->upload(true);
		this->compute(true);
		this->download(true);
		in.set_size(this->stars_generated);

		MLOG(verb1) << "Skygen generated " << this->stars_generated << " stars (" << this->nstarsExpected - total << " expected).";
		MLOG(verb1) << "Kernel runtime: " << this->swatch.getAverageTime();

		if(in.size())
		{
			total += nextlink->process(in, 0, in.size(), cpurng);
		}
#if 0
		// write out where each thread stopped
		for(int i=0; i != this->nthreads; i++)
		{
			printf("Thread %d: %d %d %d\n", i, this->cpu_state[i].x, this->cpu_state[i].y, this->cpu_state[i].z);
		}
#endif

#if 0
		// write out the results
		static bool first = true;
		FILE *fp = fopen("result.txt", first ? "w" : "a");
		if(first)
		{
			first = false;
			fprintf(fp, "# lb[2] absSDSSr{alias=absmag;alias=color;} SDSSr{alias=mag;} XYZ[3] comp\n");
		}
		column_types::cdouble::host_t lb = in.col<double>("lb");
		column_types::cfloat::host_t XYZ = in.col<float>("XYZ");
		column_types::cint::host_t comp  = in.col<int>("comp");
		column_types::cfloat::host_t M   = in.col<float>("absmag");
		column_types::cfloat::host_t m   = in.col<float>("mag");
		for(int row=0; row != in.size(); row++)
		{
			fprintf(fp, "%13.8f %13.8f %6.3f %6.3f %10.2f %10.2f %10.2f %3d\n",
				lb(row, 0), lb(row, 1), M[row], m[row], XYZ(row, 0), XYZ(row, 1), XYZ(row, 2), comp[row]);
		}
		fclose(fp);
#endif
	} while(this->stars_generated >= this->stopstars);
	MLOG(verb1) << "Generated " << total << " stars (" << this->nstarsExpected << " expected in _unclipped_ volume).";

	#if _EMU_DEBUG
	long long voxelsVisitedExpected = this->npixels;
	voxelsVisitedExpected *= this->nm*this->nM;
	std::cout << "voxels visited = " << voxelsVisited << " (expected: " << voxelsVisitedExpected << ").";
	#endif

	return total;
}

////////////////////////////////////////////////////////////////////////////

bool os_skygen::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// Load density model configuration and instantiate
	// the correct skyConfig kernel
	std::string modelfn = cfg.get("model");
	Config cfgModel(modelfn);
	std::string model = cfgModel.get("model");
	if(model == "BahcallSoneira")
	{
		skygen = new skyConfig<expModel>();
	}
	else
	{
		THROW(EAny, "Unknow density model '" + model + "'");
	}

	// Load footprint configuration
	std::string footfn = cfg.get("foot");
	Config cfgFoot(footfn);

	if(cfg.count("pdf"))
	{
		Config cfgPDF(cfg["pdf"]);
		return skygen->init(cfgPDF, cfgFoot, cfgModel, t, pipe);
	}
	else
	{
		return skygen->init(cfg, cfgFoot, cfgModel, t, pipe);
	}
}

size_t os_skygen::run(otable &in, rng_t &rng)
{
	return skygen->run(in, nextlink, rng);
}

////////////////////////////////////////////////////////////////////////////

bool os_clipper::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	cfg.get(this->dx,	"dx", 	1.f);
	this->dx = rad(this->dx);
	this->dA = sqr(this->dx);

	gpc_polygon sky = make_footprint(cfg, proj);
	skymap = make_skymap(this->dx, sky);

	return true;
}

int os_clipper::getPixelCenters(std::vector<std::pair<double, double> > &lb, double &l0, double &b0)
{
	l0 = proj.l0;
	b0 = proj.phi1;

	lb.clear();
	lb.reserve(skymap->skymap.size());
	FOREACH(skymap->skymap)
	{
		Radians x, y, l, b;
		x = skymap->x0 + skymap->dx*(i->first.first  + 0.5);
		y = skymap->y0 + skymap->dx*(i->first.second + 0.5);
		proj.inverse(x, y, l, b);
		lb.push_back(std::make_pair(l, b));
	}
	return lb.size();
}

size_t os_clipper::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// fetch prerequisites
	using namespace column_types;
	cdouble::host_t lb     = in.col<double>("lb");
	cint::host_t	hidden = in.col<int>("hidden");

	for(size_t row=begin; row < end; row++)
	{
		// clip everything outside the footprint polygon
		Radians l = rad(lb(row, 0));
		Radians b = rad(lb(row, 1));
		Radians x, y;
		proj.convert(l, b, x, y);
		std::pair<int, int> XY;
		XY.first = (int)((x - skymap->x0) / skymap->dx);
		XY.second = (int)((y - skymap->y0) / skymap->dx);

		typeof(skymap->skymap.begin()) it = skymap->skymap.find(XY);
		if(it == skymap->skymap.end()) { continue; }
// 		ASSERT(skymap->skymap.count(XY))
// 		{
// 			std::cerr << "Skymap (x0,x1,y0,y1): " << skymap->x0 << " " << skymap->x1 << " " << skymap->y0 << " " << skymap->y1;
// 			std::cerr << "(x,y) = " << x << " " << y << "\n";
// 			std::cerr << "(X,Y) = " << XY.first << " " << XY.second << "\n";
// 		}

		// check that the star is inside survey footprint, reject if it's not
		gpc_polygon &poly = it->second.poly;
		gpc_vertex vtmp = { x, y };

		bool infoot = in_polygon(vtmp, poly);
		hidden[row] = !infoot;
	}

	return nextlink->process(in, begin, end, rng);
}
