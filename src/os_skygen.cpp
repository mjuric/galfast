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
#include "simulate.h"
#include "io.h"
#include "skygen.h"
#include "analysis.h"
#include <iomanip>
#include <fstream>
#include <vector>
#include <ext/numeric>	// for iota

#include <astro/useall.h>

/***********************************************************************/

#if 0
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
#endif

xptrng::tptr<float> load_constant_texture(float2 &texCoords, float val, float X0 = -100, float X1 = 100)
{
	xptrng::tptr<float> tex(2);
	MLOG(verb2) << tex.desc->memsize() << "\n";
	tex(0U) = val;
	tex(1U) = val;
	texCoords.x = X0;
	texCoords.y = 1./(X1 - X0);
	return tex;
}

xptrng::tptr<float> load_and_resample_1D_texture(float2 &texCoords, const char *fn, int nsamp = 1024)
{
	// load the points from the file, and construct
	// a spline to resample from
	text_input_or_die(txin, fn);
	std::vector<double> X, phi;
	::load(txin, X, 0, phi, 1);
	spline tx;
	tx.construct(X, phi);

	// resample to texture
	xptrng::tptr<float> tex(nsamp, 0, 0);
	float X0 = X.front(), X1 = X.back(), dX = (X1 - X0) / (nsamp-1);
	for(int i=0; i != nsamp; i++)
	{
		tex(i) = tx(X0 + i*dX);
	}

	// construct 
	texCoords.x = X0;
	texCoords.y = 1./dX;
	return tex;
}

void expModel::load(host_state_t &hstate, const peyton::system::Config &cfg)
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
		hstate.lf = load_and_resample_1D_texture(lf, cfg["lumfunc"].c_str());
	}
	else
	{
		hstate.lf = load_constant_texture(lf, 1.f);
	}

	// cutoff radius (default: 1Mpc)
	cfg.get(r_cut2,  "rcut",   1e6f);
	r_cut2 *= r_cut2;
}

void expDisk::load(host_state_t &hstate, const peyton::system::Config &cfg)
{
	// density distribution parameters
	l     = cfg.get("l");
	h     = cfg.get("h");
	z0    = cfg.get("z0");
	f     = cfg.get("f");
	comp  = cfg.get("comp");

	// luminosity function
	if(cfg.count("lumfunc"))
	{
		hstate.lf = load_and_resample_1D_texture(lf, cfg["lumfunc"].c_str());
	}
	else
	{
		hstate.lf = load_constant_texture(lf, 1.f);
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
	cpurng = NULL;

	this->pixels = 0;
	this->lock = 0;
	this->counts = 0;
	this->countsCovered = 0;
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
	this->countsCovered.free();
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
		this->countsCovered.alloc(this->nthreads);
		cudaMemset(this->countsCovered.ptr, 0, this->nthreads*4);
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

	this->model.prerun(model_host_state, draw);
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

template<>
void skyConfig<expDisk>::upload(bool draw)
{
	this->upload_generic(draw);

	cuxUploadConst("expDiskSky", (skyConfigGPU<expDisk>&)*this); // upload thyself (note: this MUST come last)
}

int add_lonlat_rect(sph_polygon &poly, Radians l0, Radians b0, Radians l1, Radians b1, Radians dl, gpc_op op = GPC_UNION);
Radians clip_zone_of_avoidance(std::list<sph_polygon> &poly, int npoly, const Config &cfg)
{
	Radians bmin;
	cfg.get(bmin,	"bmin", 0.);
	bmin = rad(bmin);

	if(bmin)
	{
		MLOG(verb1) << "Zone of avoidance: clipping |b| < " << deg(bmin) << " deg";

		// add a "clip" contour to last npoly polygons
		REVEACH(poly)
		{
			if(!npoly--) { break; }

			const double dl = rad(.1);	// TODO: Make this user settable?
			add_lonlat_rect(*i, 0, -bmin, 0, bmin, dl, GPC_DIFF);
		}
	}

	return bmin;
}

int load_footprint_rect(std::list<sph_polygon> &foot, const peyton::system::Config &cfg);
int load_footprint_beam(std::list<sph_polygon> &foot, const peyton::system::Config &cfg);

int load_footprint(std::list<sph_polygon> &foot, const Config &cfg)
{
	if(!cfg.count("type")) { THROW(EAny, "The footprint configuration file must contain the 'type' keyword."); }

	int k = -1;
	     if(cfg["type"] == "beam") { k = load_footprint_beam(foot, cfg); }
	else if(cfg["type"] == "rect") { k = load_footprint_rect(foot, cfg); }

	if(k >= 0)
	{
		clip_zone_of_avoidance(foot, k, cfg);
		return k;
	}

	THROW(EAny, "The requested footprint type " + cfg["type"] + " is unknown.");
}

// Calculate the 'skymap', for fast point-in-polygon lookups of the observed
// sky footprint.
partitioned_skymap *make_skymap_slow(Radians dx, gpc_polygon sky)
{
	using boost::shared_ptr;
	using std::cout;

	partitioned_skymap *skymap = new partitioned_skymap;
	skymap->dx = dx;

 	poly_bounding_box(skymap->x0, skymap->x1, skymap->y0, skymap->y1, sky);
	std::cerr << "Creating fast skymap lookup (this may take a while) ... ";
	//std::cerr << skymap->x0 << " " <<  skymap->x1 << " " <<  skymap->y0 << " " <<  skymap->y1 << "\n";

	int X = 0;
	for(double x = skymap->x0; x < skymap->x1; x += skymap->dx, X++) // loop over all x values in the bounding rectangle
	{
		int Y = 0;
		double xa = x, xb = x+skymap->dx;
		for(double y = skymap->y0; y < skymap->y1; y += skymap->dx, Y++) // loop over all y values for a given x
		{
//			std::cerr << "Doing " << X << " " << Y << "\n";
			double ya = y, yb = y+skymap->dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon poly;
			gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
			if(poly.num_contours == 0) continue; // if there are no observations in this direction

			// store the polygon into a fast lookup map
			partitioned_skymap::pixel_t &pix = skymap->skymap[std::make_pair(X, Y)];
			pix.poly = poly;
			pix.coveredArea = polygon_area(poly);
			pix.pixelArea = sqr(skymap->dx);
		}
	}

	std::cerr << " done.\n";
	return skymap;
}

void make_skymap_piece(gpc_polygon sky, partitioned_skymap *skymap, int X0, int X1, int Y0, int Y1)
{
	// test for intersection
	double xa = skymap->x0 + X0*skymap->dx;
	double xb = skymap->x0 + X1*skymap->dx;
	double ya = skymap->y0 + Y0*skymap->dx;
	double yb = skymap->y0 + Y1*skymap->dx;
	gpc_polygon r = poly_rect(xa, xb, ya, yb);

//	std::cerr << "Testing (" << X0 << "," << X1 << "," << Y0 << "," << Y1 << ") == ("
//		  << xa << "," << xb << "," << ya << "," << yb << " answer=";
	gpc_polygon poly;
	gpc_polygon_clip(GPC_INT, &sky, &r, &poly);
	if(poly.num_contours == 0)
	{
//		std::cerr << "no.\n";
		return; // if there are no observations in this region
	}

//	std::cerr << "yes.";
	int DX = X1-X0, DY = Y1-Y0;
	if(DX == 1 && DY == 1) // leaf
	{
//		std::cerr << " Leaf node.\n";
		partitioned_skymap::pixel_t &pix = skymap->skymap[std::make_pair(X0, Y0)];
		pix.poly = poly;
		pix.coveredArea = polygon_area(poly);
		pix.pixelArea = sqr(skymap->dx);
//		std::cerr << pix.pixelArea/sqr(ctn::pi/180) << " " << pix.coveredArea/sqr(ctn::pi/180) << "\n";
		return;
	}

	// subdivide further
//	std::cerr << " Subdividing furthed.\n";
	make_skymap_piece(poly, skymap, X0 +    0, X0 + DX/2, Y0 +    0, Y0 + DY/2);
	make_skymap_piece(poly, skymap, X0 + DX/2, X0 + DX,   Y0 +    0, Y0 + DY/2);
	make_skymap_piece(poly, skymap, X0 + DX/2, X0 + DX,   Y0 + DY/2, Y0 +   DY);
	make_skymap_piece(poly, skymap, X0 +    0, X0 + DX/2, Y0 + DY/2, Y0 +   DY);
	gpc_free_polygon(&poly);
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
//	std::cerr << "Creating fast skymap lookup. ";
	//std::cerr << skymap->x0 << " " <<  skymap->x1 << " " <<  skymap->y0 << " " <<  skymap->y1 << "\n";

	int NX = (int)ceil((skymap->x1 - skymap->x0) / skymap->dx);
	int NY = (int)ceil((skymap->y1 - skymap->y0) / skymap->dx);
	// round up to nearest power-of-two
	NX = 1 << (int)ceil(log2(NX));
	NY = 1 << (int)ceil(log2(NY));
	if(NX == 1) NX++;
	if(NY == 1) NY++;
	NX = NY = std::max(NX, NY);	// make things really simple...
//	std::cerr << "NX,NY = " << NX << " " << NY << "\n";
	assert(skymap->x0 + NX*skymap->dx >= skymap->x1);
	assert(skymap->y0 + NY*skymap->dx >= skymap->y1);

	// hierarchically subdivide
	make_skymap_piece(sky, skymap,    0, NX/2,    0, NY/2);
	make_skymap_piece(sky, skymap, NX/2, NX,      0, NY/2);
	make_skymap_piece(sky, skymap, NX/2, NX,   NY/2,   NY);
	make_skymap_piece(sky, skymap,    0, NX/2, NY/2,   NY);

	return skymap;
}



template<typename T>
bool skyConfig<T>::init(
		otable &t,
		const peyton::system::Config &cfg,	// model cfg file
		const skygenConfig &sc,
		const skypixel *pixels
)
{
	// Galactic center distance (in pc)
	Rg = 8000.;

	// load the model
	this->model.load(model_host_state, cfg);

	// setup config & pixels
	(skygenConfig &)*this = sc;
	cpu_pixels = new skypixel[this->npixels];
	FOR(0, this->npixels) { cpu_pixels[i] = pixels[i]; }

	// For debugging/stats
	this->lrho0 = -3.5f;
	this->dlrho = 1.0f;

	// GPU kernel execution setup (TODO: should I load this through skygenConfig?)
	blockDim.x = 64; //256;
	gridDim.x = 120; // 30;
	this->nthreads = blockDim.x * blockDim.y * blockDim.z * gridDim.x * gridDim.y * gridDim.z;
	shb = gpu_rng_t::state_bytes() * blockDim.x; // for RNG

	DLOG(verb1) << "nthreads=" << this->nthreads;
	DLOG(verb1) << "nm=" << this->nm << " nM=" << this->nM << " npixels=" << this->npixels;

	return true;
}

struct star_comp
{
	cdouble_t::host_t	lb;
	cint_t::host_t		projIdx;
	cfloat_t::host_t	XYZ;
	cint_t::host_t		comp;
	cfloat_t::host_t	M;
	cfloat_t::host_t	DM;

	star_comp(otable &in)
	{
		lb      = in.col<double>("lb");
		projIdx = in.col<int>("projIdx");
		XYZ     = in.col<float>("XYZ");
		comp    = in.col<int>("comp");
		M       = in.col<float>("absmag");
		DM      = in.col<float>("DM");
	}

	bool operator()(const size_t a, const size_t b) const	// less semantics
	{
		return	lb(a, 0) <  lb(b, 0)	||	lb(a, 0) == lb(b, 0) && (
			lb(a, 1) <  lb(b, 1)	||	lb(a, 1) == lb(b, 1) && (
			DM(a)    <  DM(b)	||	DM(a)    == DM(b) && (
			M(a)     <  M(b)
			)));
	}
};

template<typename T>
void skyConfig<T>::initRNG(rng_t &cpurng)	// initialize the random number generator from CPU RNG
{
	// initialize rng
	this->cpurng = &cpurng;
	rng = new gpu_rng_t(cpurng);
}

template<typename T>
double skyConfig<T>::integrateCounts()
{
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
	MLOG(verb2) << ss.str();
	ss.str("");
	ss << "Histogram: rho=" << pow10f(this->lrho0 - this->dlrho*0.5) << "-->|";
	for(int i=0; i != this->nhistbins; i++)
	{
		ss << std::setw(8) << this->cpu_hist[i] << "|";
	}
	ss << "<--" << pow10f(this->lrho0+(this->nhistbins-0.5)*this->dlrho);
	MLOG(verb2) << ss.str();

	MLOG(verb1) << "Expected starcount: " << std::setprecision(9) << this->nstarsExpected;
	MLOG(verb2) << "Total expected star count: " <<
		std::setprecision(9) << this->nstarsExpected <<
		" (" << this->nstarsExpectedToGenerate << " in pixelized area)";
	DLOG(verb1) << "Skygen kernel runtime: " << this->swatch.getAverageTime();

	return this->nstarsExpected;
}

template<typename T>
size_t skyConfig<T>::run(otable &in, osink *nextlink)
{
	//
	// Second pass: draw stars
	//
	this->Kbatch = in.capacity();

	this->ks.alloc(this->nthreads);
	uint64_t total = 0;
	do
	{
		// setup output destination
		this->stars.lb      = in.col<double>("lb");
		this->stars.projIdx = in.col<int>("projIdx");
		this->stars.XYZ     = in.col<float>("XYZ");
		this->stars.comp    = in.col<int>("comp");
		this->stars.M       = in.col<float>("absmag");
		this->stars.DM      = in.col<float>("DM");

		this->upload(true);
		this->compute(true);
		this->download(true);
		in.set_size(this->stars_generated);

		{
			// sort the generated stars by l,b,DM,M
			std::vector<size_t> s(this->stars_generated);
			__gnu_cxx::iota(s.begin(), s.end(), size_t(0));
			std::vector<size_t> h2a(s), a2h(s);
			star_comp sc(in);
			std::sort(s.begin(), s.end(), sc);
			
			FOR(0, s.size())
			{
				#define  SWAP(arr)    std::swap(sc.arr(i),    sc.arr(j))
				#define SWAP1(arr, k) std::swap(sc.arr(i, k), sc.arr(j, k))

				int hi = a2h[i];
				int hj = s[i];
				int j = h2a[hj];

				SWAP1(lb, 0);	SWAP1(lb, 1);
				SWAP(projIdx);
				SWAP1(XYZ, 0);	SWAP1(XYZ, 1);	SWAP1(XYZ, 2);
				SWAP(comp);
				SWAP(M);
				SWAP(DM);

				std::swap(h2a[hi], h2a[hj]);
				std::swap(a2h[i],  a2h[j]);

				#undef SWAP
				#undef SWAP1
			}
		}

		DLOG(verb1) << "Skygen generated " << this->stars_generated << " stars ";
		DLOG(verb1) << "Kernel runtime: " << this->swatch.getAverageTime();

		if(in.size())
		{
			total += nextlink->process(in, 0, in.size(), *cpurng);
		}
		
		double pctdone = 100. * total / this->nstarsExpected;
		char pcts[50]; sprintf(pcts, "%.0f", pctdone);
		MLOG(verb1) << "Skygen progress: " << pcts << "% done.";

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

	double sigma = (total - this->nstarsExpected) / sqrt(this->nstarsExpected);
	char sigmas[50]; sprintf(sigmas, "% 4.1f", sigma);
	MLOG(verb1) << "Skygen completed: " << total << " stars, " << sigmas << " sigma from input model mean (" << this->nstarsExpected << ").";

	#if _EMU_DEBUG
	long long voxelsVisitedExpected = this->npixels;
	voxelsVisitedExpected *= this->nm*this->nM;
	std::cout << "voxels visited = " << voxelsVisited << " (expected: " << voxelsVisitedExpected << ").";
	#endif

	return total;
}

////////////////////////////////////////////////////////////////////////////

//void makeHemisphereMaps(gpc_polygon &nsky, gpc_polygon &ssky, peyton::math::lambert &sproj, const peyton::math::lambert &nproj, gpc_polygon allsky, Radians dx, Radians margin);
//gpc_polygon clip_zone_of_avoidance(gpc_polygon &sky, double bmin, const peyton::math::lambert &proj);
std::pair<gpc_polygon, gpc_polygon> project_to_hemispheres(const std::list<sph_polygon> &foot, const peyton::math::lambert &proj, Radians dx);

const os_clipper &os_skygen::load_footprints(const std::string &footprints, float dx, opipeline &pipe)
{
	// Load footprints. The result is the list of spherical polygons on the sky;
	std::list<sph_polygon> foot;

	std::vector<std::string> footstr;
	split(footstr, footprints);
	FOREACH(footstr)
	{
		Config cfg(*i);							// load footprint config
		load_footprint(foot, cfg);
	}

	// Project the footprint onto north/south hemispheres
	peyton::math::lambert proj(rad(90), rad(90));
//	peyton::math::lambert proj(rad(33), rad(22));
	std::pair<gpc_polygon, gpc_polygon> sky = project_to_hemispheres(foot, proj, dx);

	// setup clipper for the footprint
	boost::shared_ptr<opipeline_stage> clipper_s(opipeline_stage::create("clipper"));		// clipper for this footprint
	os_clipper &clipper = *static_cast<os_clipper*>(clipper_s.get());
	clipper.construct_from_hemispheres(dx, proj, sky);
	pipe.add(clipper_s);

	gpc_free_polygon(&sky.first);
	gpc_free_polygon(&sky.second);

	return clipper;
}

skyConfigInterface *os_skygen::create_kernel_for_model(const std::string &model)
{
	if(model == "BahcallSoneira")
		return new skyConfig<expModel>();
	if(model == "ExponentialDisk")
		return new skyConfig<expDisk>();

	THROW(EAny, "Unknow density model '" + model + "'");
	return NULL;
}

int os_skygen::load_models(otable &t, skygenConfig &sc, const std::string &model_cfg_list, const os_clipper &clipper)
{
	// prepare the skypixels to be processed by subsequently loaded models
	std::vector<os_clipper::pixel> pix;
	sc.npixels = clipper.getPixelCenters(pix);
	std::vector<skypixel> skypixels(sc.npixels);
	FOR(0, sc.npixels)
	{
		os_clipper::pixel &p = pix[i];
		float coveredFraction = p.coveredArea / p.pixelArea;
		skypixels[i] = skypixel(p.l, p.b, p.projIdx, sqrt(p.pixelArea), coveredFraction);
	}

	// projections the pixels are bound to
	std::vector<std::pair<double, double> > lb0;	// projection poles
	int nproj = clipper.getProjections(lb0);
	assert(nproj == 2);
	for(int i=0; i != lb0.size(); i++)
	{
		sc.proj[i].init(lb0[i].first, lb0[i].second);
	}

	// Load density model configurations, instantiate
	// and configure the correct skyConfig kernels
	std::vector<std::string> models;
	split(models, model_cfg_list);
	FOREACH(models)
	{
		Config cfg(*i);		// model configuration file
		boost::shared_ptr<skyConfigInterface> kernel(create_kernel_for_model(cfg.get("model")));

		// setup the sky pixels to be processed by this model
		kernel->init(t, cfg, sc, &skypixels[0]);

		kernels.push_back(kernel);
	}
}

void os_skygen::load_pdf(float &dx, skygenConfig &sc, otable &t, const std::string &cfgfn)
{
	// load PDF config
	Config cfg(cfgfn);

	// sky binning scale
	cfg.get(dx,	"dx", 1.f);
	dx = rad(dx);

	// Density binning parameters
	sc.M0 = cfg.get_any_of("M0", "ri0");
	sc.M1 = cfg.get_any_of("M1", "ri1");
	sc.dM = cfg.get_any_of("dM", "dri");
	sc.m0 = cfg.get_any_of("m0", "r0");
	sc.m1 = cfg.get_any_of("m1", "r1");
	sc.dm = cfg.get("dm");
	assert(sc.dM == sc.dm);

	sc.nM = (int)round((sc.M1 - sc.M0) / sc.dM);
	sc.nm = (int)round((sc.m1 - sc.m0) / sc.dm);
	sc.m0 += 0.5*sc.dm;
	sc.M1 -= 0.5*sc.dM;

	// prepare the table for output
	t.use_column("lb");
	t.use_column("projIdx");
	t.use_column("XYZ");
	t.use_column("comp");

	// prepare the PDF-related table columns for output
	std::string 	absmagName, bandName;
	bandName = cfg.get("band");
	absmagName = "abs" + bandName;

	t.use_column(absmagName); t.alias_column(absmagName, "absmag");
				  t.alias_column(absmagName, "M1");
				  t.set_column_property("absmag", "band", bandName);
	t.use_column("DM");
}

bool os_skygen::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	// maximum number of stars skygen is allowed to generate (0 for unlimited)
	cfg.get(nstarLimit, "nstarlimit", (size_t)100*1000*1000);

	// load PDF/sky pixelization config and prepare the table for output
	skygenConfig sc;
	float dx;
	load_pdf(dx, sc, t, cfg.get("pdf"));

	// load footprints and construct the clipper
	const os_clipper &clipper = load_footprints(cfg.get("foot"), dx, pipe);

	// load models
	load_models(t, sc, cfg.get("model"), clipper);

	return true;
}

size_t os_skygen::run(otable &in, rng_t &rng)
{
	size_t nstarsExpected = 0;
	FOREACH(kernels)
	{
		(*i)->initRNG(rng);
		nstarsExpected += (*i)->integrateCounts();
	}

	//
	// Stop if the expected number of stars is more than the limit
	//
	if(nstarsExpected > nstarLimit && nstarLimit != 0)
	{
		THROW(EAny, "The expected number of generated stars (" + str(nstarsExpected) + ") "
			"exceeds the configuration limit (" + str(nstarLimit) + "). Increase "
			"nstarlimit or set to zero for unlimited."
			);
	}

	size_t starsGenerated = 0;
	FOREACH(kernels)
	{
		starsGenerated += (*i)->run(in, nextlink);
	}
	return starsGenerated;
}

////////////////////////////////////////////////////////////////////////////

bool os_clipper::construct(const Config &cfg, otable &t, opipeline &pipe)
{
	THROW(EAny, "Module 'clipper' must never be instantiated directly.");

	return true;
}

void os_clipper::construct_from_hemispheres(float dx, const peyton::math::lambert &proj, const std::pair<gpc_polygon, gpc_polygon> &sky)
{
	// set the north hemisphere projection to input map projection
	hemispheres[0].proj = proj;
	hemispheres[1].proj = peyton::math::lambert(modulo(proj.l0 + ctn::pi, ctn::twopi), -proj.phi1);

	// construct north/south skymaps
	hemispheres[0].sky = make_skymap(dx, sky.first);
	hemispheres[1].sky = make_skymap(dx, sky.second);

	DLOG(verb1) << "Sky pixels in the north: " << hemispheres[0].sky->skymap.size();
	DLOG(verb1) << "Sky pixels in the south: " << hemispheres[1].sky->skymap.size();
}

int os_clipper::getProjections(std::vector<std::pair<double, double> > &ppoles) const	// returns the poles of all used projections
{
	ppoles.clear();
	for(int i = 0; i != 2; i++)
	{
		const peyton::math::lambert &proj = hemispheres[i].proj;
		ppoles.push_back(std::make_pair(proj.l0, proj.phi1));
	}
	return ppoles.size();
}

int os_clipper::getPixelCenters(std::vector<os_clipper::pixel> &pix) const		// returns the centers of all pixels
{
	pix.clear();

	FORj(projIdx, 0, 2)
	{
		partitioned_skymap *skymap = hemispheres[projIdx].sky;
		FOREACH(skymap->skymap)
		{
			Radians x, y, l, b;
			x = skymap->x0 + skymap->dx*(i->first.first  + 0.5);
			y = skymap->y0 + skymap->dx*(i->first.second + 0.5);
			hemispheres[projIdx].proj.deproject(l, b, x, y);
			pix.push_back(pixel(l, b, projIdx, i->second.pixelArea, i->second.coveredArea));
		}
	}

	return pix.size();
}

size_t os_clipper::process(otable &in, size_t begin, size_t end, rng_t &rng)
{
	// fetch prerequisites
	cdouble_t::host_t lb      = in.col<double>("lb");
	cint_t::host_t pIdx       = in.col<int>("projIdx");
	cint_t::host_t	hidden  = in.col<int>("hidden");

	// debugging statistics
	int nstars[2] = { 0, 0 };
	for(size_t row=begin; row < end; row++)
	{
		// clip everything outside the footprint polygon
		Radians l = rad(lb(row, 0));
		Radians b = rad(lb(row, 1));
		int projIdx = pIdx(row);
		nstars[projIdx]++;

// 		// Galactic plane zone of avoidance
// 		if(bmin && fabs(b) <= bmin)
// 		{
// 			hidden[row] = 1;
// 			continue;
// 		}

		Radians x, y;
		hemispheres[projIdx].proj.project(x, y, l, b);
		
		// immediately reject if in the southern hemisphere (for this projection)
		if(sqr(x) + sqr(y) > 2.)
		{
			hidden(row) = 1;
			continue;
		}

		partitioned_skymap *skymap = hemispheres[projIdx].sky;
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
		hidden(row) = !infoot;
	}

	DLOG(verb1) << "nstars north: " << nstars[0];
	DLOG(verb1) << "nstars south: " << nstars[1];

	return nextlink->process(in, begin, end, rng);
}
