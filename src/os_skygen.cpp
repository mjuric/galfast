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
#include <iomanip>

#include <astro/useall.h>


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

template<>
void skyConfig<expModel>::upload(bool draw)
{
	// Upload pixels to be processed
	pixels.upload(cpu_pixels, npixels);

	// prepare locks and outputs
	int zero = 0;

	if(!draw)
	{
		lock.upload(&zero, 1);

		rhoHistograms.alloc(nthreads*nhistbins);
		cudaMemset(rhoHistograms.ptr, 0, nthreads*nhistbins*4);

		maxCount.alloc(nthreads);
		cudaMemset(maxCount.ptr, 0, nthreads*4);

		counts.alloc(nthreads);
		cudaMemset(counts.ptr, 0, nthreads*4);
	}
	else
	{
		nstars.upload(&zero, 1);

		stopstars = Kbatch - bufferSafetyMargin();
		assert(stopstars > 0);

		ks.alloc(nthreads);
	}

	cuxUploadConst("Rg_gpu", this->Rg);
	cuxUploadConst("rng", *this->rng);
	cuxUploadConst("proj", this->proj);
	cuxUploadConst("expModelSky", (skyConfigGPU<expModel>&)*this); // upload thyself (note: this MUST come last)
}

template<typename T>
void skyConfig<T>::loadConfig(rng_t &seeder)
{
	// Galactic center distance (in pc)
	Rg = 8000.;

	// Galactic model initialization
//	this->model.setmodel(2150., 245., 25,  0.13, 3261, 743,   0.0051, 1./0.64, 2.77);
	this->model.setmodel(2150., 300., 25.,  0.13, 2150, 900,   0.0051, 1./0.64, 2.77);

	// density histogram setup
	this->lrho0 = -3.5f;
	this->dlrho = 1.0f;

	// luminosity function initialization
	               //  0      4      8    12    16    20   24
	const int nlf = 7;
	float cpu_lf[] = { 0.00, 0.03, 0.05, 0.07, 0.04, 0.02, 0.00 };
	float lfM0 = 0, lfM1 = 24, lfdM = (lfM1 - lfM0) / (nlf-1);
	this->model.lf = texLFMgr.set(cpu_lf, nlf, lfM0, lfM1, lfdM);

	// integration limits
	float deltam = 0.01;
//	this->M0 = 10.; this->M1 = 10.05;  this->dM = deltam;
	this->M0 =  4.; this->M1 = 15.;  this->dM = deltam;
	this->m0 = 15.; this->m1 = 24.;  this->dm = deltam;
	this->nM = (int)round((this->M1 - this->M0) / this->dM);
	this->nm = (int)round((this->m1 - this->m0) / this->dm);
	this->m0 += 0.5*this->dm;
	this->M1 -= 0.5*this->dM;

	// generate sky pixels
	this->proj.init(rad(90.), rad(90.));
	this->npixels = 1; //2000;
	this->dx = rad(1.);
	this->dA = sqr(this->dx);
	cpu_pixels = new direction[this->npixels];
	Radians l0 = rad(44), b0 = rad(85.), dx = rad(1), l, b;
//	Radians l0 = rad(0.), b0 = rad(5.), dx = rad(1), l, b;
	for(int i=0; i != this->npixels; i++)
	{
		l = l0 + i*dx;
		b = b0;
		cpu_pixels[i] = direction(l, b);
	}

	// compute/set kernel execution parameters
//	this->Kbatch = 3000000;
	blockDim.x = 64; //256;
	gridDim.x = 120; // 30;
	this->nthreads = blockDim.x * blockDim.y * blockDim.z * gridDim.x * gridDim.y * gridDim.z;
	std::cout << "nthreads=" << this->nthreads << "\n";
	shb = 12 * blockDim.x; // for RNG

	// initialize rng
#if 1
	seed = 43;
#else
	FILE *fp;
	fp = fopen("/dev/urandom", "r");
	fread(&seed, sizeof(seed), 1, fp);
	fclose(fp);
	std::cerr << "Using seed " << seed << "\n";
#endif
//	rng = new gpuRng(seed, this->nthreads);
	rng = new gpu_rng_t(seeder);
}

////////////////////////////////////////////////////////////////////////////

bool os_skygen::init(const Config &cfg1, otable &t)
{
//	const char *fn = cfg1["filename"].c_str();
//	Config cfg(fn);

	// this fills the prov vector with columns this module will provide
/*	prov.insert("lb");
	prov.insert("XYZ");
	prov.insert("comp");*/
	t.use_column("lb");
	t.use_column("XYZ");
	t.use_column("comp");

	std::string 	colorName = "absSDSSr",
			absmagName = "absSDSSr",
			magName = "SDSSr";

	t.use_column(colorName);  t.alias_column(colorName, "color");
	t.use_column(absmagName); t.alias_column(absmagName, "absmag");
	t.use_column(magName);	  t.alias_column(magName, "mag");

	return true;
}

size_t os_skygen::run(otable &in, rng_t &rng)
{
	skyConfig<expModel> sky;
	sky.loadConfig(rng);
	std::cout << "nm=" << sky.nm << " nM=" << sky.nM << " npixels=" << sky.npixels << "\n";

	//
	// First pass: compute total expected starcounts
	//
	sky.upload(false);
	sky.compute(false);
	sky.download(false);
	std::cout << "Histogram:     log(rho) |";
	for(int i=0; i != sky.nhistbins; i++)
	{
		std::cout << std::setw(8) << sky.lrho0+i*sky.dlrho << "|";
	}
	std::cout << "\n";
	std::cout << "Histogram: rho=" << pow10f(sky.lrho0 - sky.dlrho*0.5) << "-->|";
	for(int i=0; i != sky.nhistbins; i++)
	{
		std::cout << std::setw(8) << sky.cpu_hist[i] << "|";
	}
	std::cout << "<--" << pow10f(sky.lrho0+(sky.nhistbins-0.5)*sky.dlrho) << "\n";

	std::cout << "Total expected star count: " << std::setprecision(9) << sky.nstarsExpected << "\n";
	std::cout << "Kernel runtime: " << sky.swatch.getAverageTime() << "\n";

	//
	// Second pass: draw stars
	//
	sky.Kbatch = in.capacity();

	// setup output destination
	sky.stars.lb    = in.col<double>("lb");
	sky.stars.XYZ   = in.col<float>("XYZ");
	sky.stars.comp  = in.col<int>("comp");
	sky.stars.M     = in.col<float>("absmag");
	sky.stars.m     = in.col<float>("mag");
	//sky.stars.color = in.col<float>("color"); -- not implemented yet
	sky.upload(true);

	uint64_t total = 0;
	do
	{
		int zero = 0;
		sky.nstars.upload(&zero, 1);
		sky.compute(true);
		sky.download(true);
		in.set_size(sky.stars_generated);

		std::cout << "Generated " << sky.stars_generated << " stars (" << sky.nstarsExpected - total << " expected).\n";
		std::cout << "Kernel runtime: " << sky.swatch.getAverageTime() << "\n";

		total += nextlink->process(in, 0, in.size(), rng);
#if 0
		// write out where each thread stopped
		for(int i=0; i != sky.nthreads; i++)
		{
			printf("Thread %d: %d %d %d\n", i, sky.cpu_state[i].x, sky.cpu_state[i].y, sky.cpu_state[i].z);
		}
#endif

#if 0
		// write out the results
		FILE *fp = fopen("result.txt", total ? "a" : "w");
		if(total == 0)
		{
			fprintf(fp, "# lb[2] absSDSSr{alias=absmag;alias=color;} SDSSr{alias=mag;} XYZ[3] comp\n");
		}
		for(int i=0; i != sky.stars_generated; i++)
		{
			star &s = sky.cpu_stars[i];
			fprintf(fp, "%13.8f %13.8f %6.3f %6.3f %10.2f %10.2f %10.2f %3d\n",
				deg(s.l), deg(s.b), s.M, s.m, s.pos.x, s.pos.y, s.pos.z, s.comp);
		}
		fclose(fp);
#endif
	} while(sky.stars_generated >= sky.stopstars);
	std::cout << "Generated " << total << " stars (" << sky.nstarsExpected << " expected).\n";

	#if _EMU_DEBUG
	long long voxelsVisitedExpected = sky.npixels;
	voxelsVisitedExpected *= sky.nm*sky.nM;
	std::cout << "voxels visited = " << voxelsVisited << " (expected: " << voxelsVisitedExpected << ")\n";
	#endif

	return total;
}
