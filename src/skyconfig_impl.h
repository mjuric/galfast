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
#include "observe.h"
#include "simulate.h"
#include <iomanip>
#include <ext/numeric>	// for iota

#include <astro/useall.h>

/***********************************************************************/

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
	this->norm = 1.f;

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

#if 0
template<typename T>
int skyConfig<T>::bufferSafetyMargin()
{
	// Compute the extra padding necessary to make buffer overfilling
	// exceedingly unlikely

	if(cpu_maxCount == NULL) { abort(); }

	float maxCount = std::accumulate(cpu_maxCount, cpu_maxCount + this->nthreads, 0.f);
	int bufsize = (int)ceil(maxCount + 7.*sqrt(maxCount))	// for maxCount > 100, P(k > k+7sqrt(k)) ~ nthreads*7e-11
			 + 3*this->nthreads;			// and this part is to cover the selection effect that we only break if k>1

#if 1
	// some diagnostic info
	DLOG(verb1) << "Maxcount diagnostics:";
	sort(cpu_maxCount, cpu_maxCount + this->nthreads, std::greater<float>());
	for(int i = 0; i != 3 && i != this->nthreads; i++)
	{
		DLOG(verb1) << "rho" << i << " = " << cpu_maxCount[i] << "\n";
	}
	DLOG(verb1) << "Total of max density bins: " << maxCount << "\n";
	DLOG(verb1) << "Buffer safety margin: " << bufsize << "\n";
#endif

	return bufsize;
}
#endif

template<typename T>
void skyConfig<T>::download(bool draw)
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
		float *cpu_counts = new float[this->nthreads];
		float *cpu_countsCovered = new float[this->nthreads];
		this->counts.download(cpu_counts, this->nthreads);
		this->countsCovered.download(cpu_countsCovered, this->nthreads);

		// sum up the total expected number of stars
		nstarsExpectedToGenerate = 0;
		nstarsExpected = 0;
		for(int i=0; i != this->nthreads; i++)
		{
			nstarsExpectedToGenerate += cpu_counts[i];
			nstarsExpected += cpu_countsCovered[i];
		}
		delete [] cpu_counts;
		delete [] cpu_countsCovered;

		int *cpu_rhoHistograms = new int[this->nthreads*this->nhistbins];
		this->rhoHistograms.download(cpu_rhoHistograms, this->nthreads*this->nhistbins);

		// sum up the total
		cpu_hist = new int[this->nhistbins];
		memset(cpu_hist, 0, sizeof(float)*this->nhistbins);
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

template<typename T>
void skyConfig<T>::upload(bool draw)
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

//		this->stopstars = Kbatch - bufferSafetyMargin();
		this->stopstars = output_table_capacity;
		assert(this->stopstars > 0);
	}

	cuxUploadConst("Rg_gpu", this->Rg);
	cuxUploadConst("rng", *this->rng);
	cuxUploadConst("proj", this->proj);

	this->model.prerun(model_host_state, draw);
	this->upload_self(draw);
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

#if 1	// should be moved elsewhere

	// GPU kernel execution setup (TODO: should I load this through skygenConfig? Or autodetect based on the GPU?)
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
void skyConfig<T>::initRNG(rng_t &cpurng)	// initialize the random number generator from CPU RNG
{
	// initialize rng
	this->cpurng = &cpurng;
	rng = new gpu_rng_t(cpurng);
}

template<typename T>
void skyConfig<T>::setDensityNorm(float norm_)	// set the density normalization of the model
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
	this->output_table_capacity = in.capacity();

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

