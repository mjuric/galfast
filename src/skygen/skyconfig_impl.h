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

#include "skygen.h"
#include "otable.h"
#include "../pipeline.h"
#include "analysis.h"

#include <iomanip>
#include <ext/numeric>	// for iota

#include <astro/useall.h>

/***********************************************************************/

// aux. class for sorting the rows. Used in skygenHost<T>::drawSources().
struct star_comp
{
	cdouble_t::host_t	lb;
	cint_t::host_t		projIdx;
	cfloat_t::host_t	XYZ, projXY;
	cint_t::host_t		comp;
	cfloat_t::host_t	M, Am;
	cfloat_t::host_t	DM;

	star_comp(otable &in)
	{
		lb      = in.col<double>("lb");
		projIdx = in.col<int>("projIdx");
		projXY  = in.col<float>("projXY");
		XYZ     = in.col<float>("XYZ");
		comp    = in.col<int>("comp");
		M       = in.col<float>("absmag");
		Am      = in.col<float>("Am");
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
skygenHost<T>::skygenHost()
{
	cpu_pixels = NULL;
	cpu_hist = NULL;
	cpu_maxCount = NULL;
	cpu_state = NULL;
	rng = NULL;
	cpurng = NULL;

	this->pixels = 0;
	this->nstars = 0;
	this->counts = 0;
	this->countsCovered = 0;
	this->rhoHistograms = 0;
	this->maxCount = 0;

	this->norm = 1.f;
	this->ks.constructor();
}

template<typename T>
skygenHost<T>::~skygenHost()
{
	delete [] cpu_pixels;
	delete [] cpu_hist;
	delete [] cpu_maxCount;
	delete [] cpu_state;
	delete rng;

	this->pixels.free();
	this->counts.free();
	this->countsCovered.free();
	this->nstars.free();

	this->ks.destructor();
}

//
// Download skygen data to GPU. Flag 'draw' denotes if this call was preceeded
// by launch of a kernel to draw stars, or to compute the overal normalization.
//
template<typename T>
void skygenHost<T>::download(bool draw, int pixfrom, int pixto)
{
	if(draw)
	{
		this->nstars.download(&stars_generated, 1);
		stars_generated = std::min(stars_generated, this->stopstars);

		// this is for debugging purposes mostly
		int *ilb = new int[this->nthreads];
		int *im = new int[this->nthreads];
		int *iM = new int[this->nthreads];
		delete [] cpu_state;
		cpu_state = new int3[this->nthreads];
		this->ks.ilb.download(ilb, this->nthreads);
		this->ks.im.download(im,   this->nthreads);
		this->ks.iM.download(iM,   this->nthreads);
		for(int i=0; i != this->nthreads; i++)
		{
			cpu_state[i] = make_int3(ilb[i], im[i], iM[i]);
		}
		delete [] ilb; delete [] im, delete [] iM;
	}
	else
	{
		double *cpu_counts = new double[this->nthreads];
		double *cpu_countsCovered = new double[this->nthreads];
		this->counts.download(cpu_counts, this->nthreads);
		this->countsCovered.download(cpu_countsCovered, this->nthreads);
		// sum up the total expected number of stars
		if(pixfrom == 0)
		{
			nstarsExpectedToGenerate = 0;
			nstarsExpected = 0;
		}
		double counts = 0, countsCovered = 0;
		for(int i=0; i != this->nthreads; i++)
		{
			nstarsExpectedToGenerate += cpu_counts[i];
			nstarsExpected += cpu_countsCovered[i];
		}
		nstarsExpectedToGenerate += counts;
		nstarsExpected += countsCovered;
		delete [] cpu_counts;
		delete [] cpu_countsCovered;

		int *cpu_rhoHistograms = new int[this->nthreads*this->nhistbins];
		this->rhoHistograms.download(cpu_rhoHistograms, this->nthreads*this->nhistbins);

		// sum up the total
		if(pixfrom == 0)
		{
			delete [] cpu_hist;
			cpu_hist = new int[this->nhistbins];
			memset(cpu_hist, 0, sizeof(float)*this->nhistbins);
		}
		for(int i=0; i != this->nthreads; i++)
		{
			for(int j=0; j != this->nhistbins; j++)
			{
				cpu_hist[j] += cpu_rhoHistograms[this->nthreads*j + i];
			}
		}
		delete [] cpu_rhoHistograms;

		// Download list of maximum densities encountered by each thread
		delete [] cpu_maxCount;
		cpu_maxCount = new float[this->nthreads];
		this->maxCount.download(cpu_maxCount, this->nthreads);
	}

	this->model.postrun(model_host_state, draw);
}

#if !SKYGEN_ON_GPU
extern gpuRng::constant rng;	// GPU RNG
extern lambert proj[2];		// Projection definitions for the two hemispheres
#endif

//
// Upload skygen data to GPU. Flag 'draw' denotes if this call will be followed
// by a launch of a kernel to draw stars, or to compute the overal normalization.
//
template<typename T>
void skygenHost<T>::upload(bool draw, int pixfrom, int pixto)
{
	// Upload pixels to be processed
	this->npixels = pixto - pixfrom;
	this->pixels.free();
	this->pixels.upload(cpu_pixels + pixfrom, this->npixels);

	if(!draw)
	{
		cpu_countsCoveredPerBeam = cuxSmartPtr<float>(this->nthreads, this->npixels);
		//std::cerr << cpu_countsCoveredPerBeam.size() << " " << this->nthreads << " " << this->npixels << "\n";
		FOREACH(cpu_countsCoveredPerBeam) { *i = 0.; }
		this->countsCoveredPerBeam = cpu_countsCoveredPerBeam;

		this->rhoHistograms.realloc(this->nthreads*this->nhistbins);
		cudaMemset(this->rhoHistograms.ptr, 0, this->nthreads*this->nhistbins*4);

		this->maxCount.realloc(this->nthreads);
		cudaMemset(this->maxCount.ptr, 0, this->nthreads*4);

		this->counts.realloc(this->nthreads);
		cudaMemset(this->counts.ptr, 0, this->nthreads*4);
		this->countsCovered.realloc(this->nthreads);
		cudaMemset(this->countsCovered.ptr, 0, this->nthreads*4);
	}
	else
	{
		int zero = 0;
		this->nstars.upload(&zero, 1);

		this->stopstars = output_table_capacity;
		assert(this->stopstars > 0);
	}

#if SKYGEN_ON_GPU
	cuxUploadConst("rng", *this->rng);
	cuxUploadConst("proj", this->proj);
#else
	::rng = *this->rng;
	::proj[0] = this->proj[0];
	::proj[1] = this->proj[1];
#endif

	this->model.prerun(model_host_state, draw);
	this->upload_self(draw);
}

template<typename T>
bool skygenHost<T>::init(
		otable &t,
		const peyton::system::Config &cfg,	// model cfg file
		const skygenParams &sc,
		const pencilBeam *pixels
)
{
	// Galactic center distance (in pc)
//	Rg = 8000.;

	// load the model
	this->model.load(model_host_state, cfg);

	// setup config & pixels
	(skygenParams &)*this = sc;
	cpu_pixels = new pencilBeam[this->npixels];
	FOR(0, this->npixels) { cpu_pixels[i] = pixels[i]; }

	// For debugging/stats
	this->lrho0 = -3.5f;
	this->dlrho = 1.0f;

#if 1	// should be moved elsewhere

	// GPU kernel execution setup (TODO: should I load this through skygenParams? Or autodetect based on the GPU?)
	blockDim.x = 64; //256;
	gridDim.x = 120; // 30;

	// HACK: This should be settable through config files
	char *blk = getenv("SKYGEN_KCONF"); // expect the form of "gx gy gz bx by bz"
	if(blk != NULL)
	{
		assert(sscanf(blk, "%d %d %d %d %d %d", &gridDim.x, &gridDim.y, &gridDim.z, &blockDim.x, &blockDim.y, &blockDim.z) == 6);
		DLOG(verb1) << "Reading config from SKYGEN_KCONF; "
			"grid=(" << gridDim.x  << ", " << gridDim.y  << ", " << gridDim.z << ") "
			"block=(" << blockDim.x << ", " << blockDim.y << ", " << blockDim.z << ")\n";
	}

	this->nthreads = blockDim.x * blockDim.y * blockDim.z * gridDim.x * gridDim.y * gridDim.z;
	shb = gpu_rng_t::state_bytes() * blockDim.x; // for RNG

	DLOG(verb1) << "nthreads=" << this->nthreads;
	DLOG(verb1) << "nm=" << this->nm << " nM=" << this->nM << " npixels=" << this->npixels;
#endif

	return true;
}

template<typename T>
void skygenHost<T>::initRNG(rng_t &cpurng)	// initialize the random number generator from CPU RNG
{
	// initialize rng
	this->cpurng = &cpurng;
	rng = new gpu_rng_t(cpurng);
}

template<typename T>
void skygenHost<T>::setDensityNorm(float norm_)	// set the density normalization of the model
{
	this->norm = norm_;
	if(cpu_maxCount != NULL)
	{
		FOR(0, this->nthreads)
		{
			cpu_maxCount[i] *= norm_;
		}
	}
	this->nstarsExpected *= norm_;
	this->nstarsExpectedToGenerate *= norm_;
}

template<typename T>
double skygenHost<T>::integrateCounts(float &runtime, const char *denmapPrefix)
{
	//
	// First pass: compute total expected starcounts
	//
	swatch.reset();
	swatch.start();

	FILE *fp = NULL;
	if(denmapPrefix && denmapPrefix[0] != 0)
	{
		std::string skyDensityMapFile(denmapPrefix);
		skyDensityMapFile += "." + str(this->model.component()) + ".txt";

		fp = fopen(skyDensityMapFile.c_str(), "w");
		ASSERT(fp != NULL);
	}

	int lastpix = this->npixels;
	std::vector<double> den(lastpix);

	// Process in blocks of PIXBLOCK pixels, where PIXBLOCK is computed
	// so that countsCoveredPerBeam array has about ~10M entries max
	const int PIXBLOCK = 1024*1024*2 / this->nthreads;
	//const int PIXBLOCK = this->nthreads;
	ASSERT(PIXBLOCK >= 1);
	//std::cerr << "npix, PIXELBLOCK = " << lastpix << " " << PIXBLOCK << "\n";

	int step = (int)ceil(((float)lastpix/PIXBLOCK)/50);
	ticker tick("Integrating", step);

	int startpix = 0;
	while(startpix < lastpix)
	{
		tick.tick();
		int endpix = std::min(startpix + PIXBLOCK, lastpix);

		upload(false, startpix, endpix);
		compute(false);
		download(false, startpix, endpix);

		// sum up and print out the on-sky densities, if requested
		if(fp != NULL)
		{
			int npix = endpix - startpix;
			for(int px = startpix; px != endpix; px++)
			{
				double rho = 0;
				for(int th = 0; th != this->nthreads; th++)
				{
					rho += cpu_countsCoveredPerBeam(th, px-startpix);
				}
				den[px] = rho;
				//std::cout << cpu_countsCoveredPerBeam(100, px-startpix) << " " << rho << "\n";
				
				Radians l, b;
				pencilBeam pb = cpu_pixels[px];
				this->proj[pb.projIdx].deproject(l, b, pb.X, pb.Y);
				l *= 180./dbl_pi; b *= 180./dbl_pi;
				fprintf(fp, "% 9.4f % 8.4f %12.2f   % 8.4f % 8.4f %d\n", l, b, rho, pb.X, pb.Y, pb.projIdx);
			}
		}

		startpix = endpix;
	}
	this->npixels = lastpix;
	tick.close();

	if(fp)
	{
		fclose(fp);

		double total = accumulate(den.begin(), den.end(), 0.);
		ASSERT(fabs(total / this->nstarsExpected - 1.) < 1e-5)
		{
			std::cerr << "The expected number of stars, computed in two different ways is inconsistent.\n";
			std::cerr << "   nstarsExpected=" << this->nstarsExpected << "\n";
			std::cerr << "   total=" << total << "\n";
		}
	}

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

	MLOG(verb1) << "Comp. " << this->model.component() << " starcount : " << std::setprecision(9) << this->nstarsExpected
		<< " (" << this->nstarsExpectedToGenerate << " in pixelized area)";

	swatch.stop();
	runtime = swatch.getTime();
	DLOG(verb1) << "integrateCounts() runtime: " << runtime << "s";
	
	return this->nstarsExpected;
}

//
// Draw the catalog
//
template<typename T>
size_t skygenHost<T>::drawSources(otable &in, osink *nextlink, float &runtime)
{
	swatch.reset();
	swatch.start();

	this->output_table_capacity = in.capacity();

	this->ks.alloc(this->nthreads);
	uint64_t total = 0;
	do
	{
		// setup output destination
		this->stars.lb      = in.col<double>("lb");
		this->stars.projIdx = in.col<int>("projIdx");
		this->stars.projXY  = in.col<float>("projXY");
		this->stars.Am      = in.col<float>("Am");
		this->stars.XYZ     = in.col<float>("XYZ");
		this->stars.comp    = in.col<int>("comp");
		this->stars.M       = in.col<float>("absmag");
		this->stars.DM      = in.col<float>("DM");

		// call the generation kernel
		upload(true, 0, this->npixels);
		compute(true);
		download(true, 0, this->npixels);

		in.set_size(this->stars_generated);	// truncate the table to the number of rows that were actually generated

		{
			//
			// sort the generated stars by l,b,DM,M
			//
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
				SWAP1(projXY, 0);	SWAP1(projXY, 1);
				SWAP(Am);
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
		//DLOG(verb1) << "Kernel runtime: " << swatch.getAverageTime();

		swatch.stop();
		if(in.size())
		{
			total += nextlink->process(in, 0, in.size(), *cpurng);
		}
		swatch.start();
		
		double pctdone = 100. * total / this->nstarsExpected;
		char pcts[50]; sprintf(pcts, "%.0f", pctdone);
		MLOG(verb1) << "Comp. "<< this->model.component() << " progress: " << pcts << "% done.";

	} while(this->stars_generated >= this->stopstars);

	double sigma = (total - this->nstarsExpected) / sqrt(this->nstarsExpected);
	char sigmas[50]; sprintf(sigmas, "% 4.1f", sigma);
	MLOG(verb1) << "Comp. "<< this->model.component() << " completed: " << total << " stars, " << sigmas << " sigma from input model mean (" << this->nstarsExpected << ").";

	swatch.stop();
	runtime = swatch.getTime();
	DLOG(verb1) << "run() runtime: " << runtime << "s";

	return total;
}
