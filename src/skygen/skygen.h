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

#ifndef skygen__h_
#define skygen__h_

#include <cuda_runtime_api.h>
#include "gpu.h"
#include "column.h"
#include "module_lib.h"
#include <astro/types.h>
#include <string>

typedef prngs::gpu::mwc gpuRng;
using peyton::Radians;

// some forward declarations
namespace peyton { namespace system { class Config; }};
class otable;
class osink;
struct skygenParams;
struct pencilBeam;

//
// Abstract interface to mock catalog generator for a model. For each density model,
// an instance of skygenHost<> (that derives from this class) is constructed
// and returned to os_skygen. os_skygen calls its methods (notably, integrateCounts()
// and run()) to draw the catalog.
//
struct ALIGN(16) skygenInterface
{
	virtual double integrateCounts(float &runtime, const char *denmappfix) = 0;	// computes the expected source count
	virtual size_t drawSources(otable &in, osink *nextlink, float &runtime) = 0;	// draws the sources, stores the output in the table, and invokes the rest of the pipeline

	virtual bool init(								// initialize the catalog generator for this model
		otable &t,
		const peyton::system::Config &model_cfg,
		const skygenParams &sc,
		const pencilBeam *pixels) = 0;
	virtual void initRNG(rng_t &rng) = 0;		// initialize the random number generator from CPU RNG
	virtual void setDensityNorm(float norm) = 0;	// explicitly set the overall density normalization of the model.
	virtual ~skygenInterface() {};
};

//
// These must be overloaded by individual the density model classes.
//
struct modelConcept
{
	struct host_state_t {};		// the state of the model that needs to be kept on the host
	struct state {};		// the state of the model that needs to be kept on the GPU

	void load(host_state_t &hstate, const peyton::system::Config &cfg);	// loads the density model configuration
	void prerun(host_state_t &hstate, bool draw);				// called before executing the model generation kernel
	void postrun(host_state_t &hstate, bool draw);				// called after executing the model generation kernel

	__device__ void setpos(state &s, float x, float y, float z) const;	// set the 3D position which will be implied in subsequent calls to rho()
	__device__ float rho(state &s, float M) const;				// return the number density at the position set by setpos, and absolute magnitude M

	__device__ int component() const;					// return the component ID of this model
};

//
// A rectangular (in projected coordinates) pencil beam
// on the sky, pointing in direction <direction>
//
struct ALIGN(16) pencilBeam : public direction
{
	float coveredFraction;	// the fraction of the pixel that is covered by the footprint being observed
	float dx, dA;		// linear size (in radians), and area (radians^2) of the pencil beam

	int projIdx;		// index of the lambert projection in which this pixel was constructed
	float X, Y;		// lambert coordinates of this pencil beam (in projIdx projection)

	pencilBeam() {}
	pencilBeam(Radians l_, Radians b_, float X_, float Y_, int projIdx_, float dx_, float coveredFraction_)
	: direction(l_, b_), X(X_), Y(Y_), projIdx(projIdx_), dx(dx_), dA(dx_*dx_), coveredFraction(coveredFraction_)
	{ }
};

//
// Lambert projection/deprojection class
//
class ALIGN(16) lambert
{
public:
	Radians l0, b0;
	struct { float cl, cb, sl, sb; } p;	// projection pole

public:
	void init(Radians l0_, Radians b0_)
	{
		l0 = l0_; b0 = b0_;
		p.cl = cos(l0); p.cb = cos(b0);
		p.sl = sin(l0); p.sb = sin(b0);
	}

	__device__ void project(float &x, float &y, const direction &d) const
	{
		double cll0 = d.cl*p.cl + d.sl*p.sl; // cos(l - l0);
		double sll0 = p.cl*d.sl - d.cl*p.sl; // sin(l - l0);
		double kp = 1.414213562373095 * sqrt(1./(1.+p.sb*d.sb+p.cb*d.cb*cll0));
		x = kp*d.cb*sll0;
		y = kp*(p.cb*d.sb-p.sb*d.cb*cll0);
	}
	
	__device__ void deproject(double &l, double &b, const float x, const float y) const
	{
		double r2 = x*x + y*y;
		if(r2 == 0.)
		{
			b = b0;
			l = 0.;
			return;
		}
		double ir = 1./sqrt(r2);			// ir == 1/r
		double sc = 1./ir * sqrt(1. - 0.25*r2);		// sin(c) == r*Sqrt(1 - Power(r,2)/4.)		where c = 2*asin(0.5*r)
		double cc = 1. - 0.5*r2;			// cos(c) == 1 - Power(r,2)/2			where c = 2*asin(0.5*r)

		b = asin(cc*p.sb + y*sc*p.cb*ir);
		l = l0 + atan2(x*sc, p.cb*cc/ir - y*p.sb*sc);
	}
};

//
// Configuration parameters of catalog generator.
//
struct ALIGN(16) skygenParams
{
	float m0, m1, dm;		// apparent magnitude range
	float M0, M1, dM;		// absolute magnitude range
	int npixels, nm, nM;		// number of sky pixels (pencil beams), number of magnitude bins, number of absmag bins

	int nthreads;			// total number of threads processing the sky
	int stopstars;			// stop after this many stars have been generated

	lambert proj[2];		// north/south sky lambert projections
};

//
// Columns used by the catalog generator
//
struct ALIGN(16) ocolumns
{
	cdouble_t::gpu_t	lb;
	cint_t::gpu_t		projIdx;
	cfloat_t::gpu_t		projXY, DM, M, XYZ, Am;
	cint_t::gpu_t		comp;
};

//
// Runtime state of a skygen kernel. Think of this as a "stack dump"
// for the kernel, if the kernel needs to exit early (e.g., because
// the output table is filled), but expects to resume later.
//
template<typename Model>
struct ALIGN(16) runtime_state
{
	typedef typename Model::state ms_t;

	cuxDevicePtr<int> cont, ilb, im, iM, k, bc, ndraw;
	cuxDevicePtr<float3> pos;
	cuxDevicePtr<float> D, Am;
	cuxDevicePtr<pencilBeam> pix;
	cuxDevicePtr<ms_t> ms;

	void alloc(int nthreads)
	{
		cont.alloc(nthreads); ilb.alloc(nthreads); im.alloc(nthreads); iM.alloc(nthreads); k.alloc(nthreads); bc.alloc(nthreads); ndraw.alloc(nthreads);
		pos.alloc(nthreads); D.alloc(nthreads); Am.alloc(nthreads); pix.alloc(nthreads);
		ms.alloc(nthreads);

		cudaMemset(cont.ptr, 0, nthreads*4); // not continuing a previous run
	}
	void free()
	{
		cont.free(); ilb.free(); im.free(); iM.free(); k.free(); bc.free(); ndraw.free();
		pos.free(); D.free(); Am.free(); pix.free();
		ms.free();
	}
	void constructor()	// as CUDA doesn't allow real constructors
	{
		cont = ilb = im = iM = k = bc = ndraw = 0;
		pos = 0; D = 0; Am = 0; pix = 0;
		ms = 0;
	}
	void destructor()
	{
		free();
	}
	void reset()
	{
		free();
	}

	__device__ void load(int &tid, int &ilb, int &im, int &iM, int &k, int &bc, float3 &pos, float &D, pencilBeam &pix, float &Am, typename Model::state &ms, int &ndraw) const
	{
		ilb = this->ilb(tid);
		im  = this->im(tid);
		iM  = this->iM(tid);
		k   = this->k(tid);
		bc  = this->bc(tid);
		pos = this->pos(tid);
		pix = this->pix(tid);
		D   = this->D(tid);
		Am  = this->Am(tid);
		ms  = this->ms(tid);
		ndraw  = this->ndraw(tid);
	}

	__device__ void store(int tid, int ilb, int im, int iM, int k, int bc, float3 pos, float D, pencilBeam pix, float Am, typename Model::state ms, int ndraw) const
	{
		this->ilb(tid)   = ilb;
		this->im(tid)    = im;
		this->iM(tid)    = iM;
		this->k(tid)     = k;
		this->bc(tid)    = bc;
		this->pos(tid)   = pos;
		this->pix(tid)   = pix;
		this->D(tid)     = D;
		this->Am(tid)    = Am;
		this->ms(tid)    = ms;
		this->ndraw(tid) = ndraw;

		this->cont(tid) = 1;
	}
	
	bool continuing(int tid) const { return this->cont(tid); }
};

//
// Device functions and data for skygen (instantiated/defined in skygen.cu.h)
// This piece gets uploaded as a __constant__ to the GPU. It differs from
// skygenParams in that it contains pointers to output as well as model
// configuration, but no configuration parameters beyond those provided in
// skygenParams.
//
template<typename Model>
struct ALIGN(16) skygenGPU : public skygenParams
{
	typedef Model Model_t;

	Model_t model;

	cuxDevicePtr<pencilBeam> pixels;	// pixels on the sky to process

	cuxDevicePtr<int> nstars;
	cuxDevicePtr<double> counts, countsCovered;
	gptr<float, 2> countsCoveredPerBeam;
	ocolumns stars;

	static const int nhistbins = 10;
	cuxDevicePtr<int> rhoHistograms;	// nthreads*nbins sized array
	float lrho0, dlrho;		// histogram array start (bin midpoint), bin size

	cuxDevicePtr<float> maxCount;	// [nthreads] sized array, returning the maximum density found by each thread

	runtime_state<Model> ks;
	float norm;			// normalization of overall density (usually 1.f)

	template<int draw> __device__ void kernel() const;
	__device__ float3 compute_pos(float &D, float &Am, float M, const int im, const pencilBeam &dir) const;
	__device__ bool advance(int &ilb, int &i, int &j, pencilBeam &pix, const int x, const int y) const;
	__device__ void draw_stars(int &ndraw, const float &M, const int &im, const pencilBeam &pix) const;
};

//
// Host interface to the mock catalog generator. Derived from skygenGPU part (which
// is uploaded to the device), and the abstract skygenInterface that allows instantiations
// of this template to be manipulated by os_skygen in a uniform manner.
//
class opipeline;
template<typename Model>
class ALIGN(16) skygenHost : public skygenGPU<Model>, public skygenInterface
{
protected:
	// Any state the model needs to maintain on the host (e.g., loaded textures)
	typename Model::host_state_t	model_host_state;

//	float Rg;		// distance to the galactic center

	gpuRng *rng;
	rng_t *cpurng;
	unsigned seed;
	pencilBeam *cpu_pixels;
	cuxSmartPtr<float> cpu_countsCoveredPerBeam;

	std::string skyDensityMapFile;

	// return
	double nstarsExpectedToGenerate;	// expected number of stars in the pixelized footprint
	double nstarsExpected;			// expected number of stars in the actual footprint
	int stars_generated;
	int *cpu_hist;
	float *cpu_maxCount;
	int3 *cpu_state;

	// kernel execution configuration
	dim3 gridDim, blockDim;		// CUDA grid dimension, block dimension
	int shb;			// shared memory per block needed by skygen kernel
	int output_table_capacity;	// number of rows in the output table
	stopwatch swatch;		// runtime of this model (measured in integrateCounts() and run()).

	void upload_self(bool draw = false);

	void compute(bool draw = false);
	void upload(bool draw, int pixfrom, int pixto);
	void download(bool draw, int pixfrom, int pixto);

public:
	skygenHost();
	~skygenHost();

	// external interface (skygenInterface)
	virtual double integrateCounts(float &runtime, const char *denmapPrefix);	// return the expected starcounts contributed by this model
	virtual size_t drawSources(otable &in, osink *nextlink, float &runtime);

	virtual void initRNG(rng_t &rng);		// initialize the random number generator from CPU RNG
	virtual void setDensityNorm(float norm);	// explicitly set the overall density normalization of the model.
	virtual bool init(
		otable &t,
		const peyton::system::Config &cfg,	// model cfg file
		const skygenParams &sc,
		const pencilBeam *pixels);
};

//
// Some generic stub code to instantiate the GPU kernels, the class,
// and a specialized upload method for the model in question
//
// Active only when compiled with nvcc
//
#define SKYGEN_ON_GPU	0
#ifdef __CUDACC__
	#if SKYGEN_ON_GPU
		#define MODEL_IMPLEMENTATION(name) \
			__device__ __constant__ skygenGPU<name> name##Sky; \
			\
			template<> __global__ void integrate_counts_kernel<name>() { name##Sky.kernel<0>(); } \
			template<> __global__ void     draw_sources_kernel<name>() { name##Sky.kernel<1>(); } \
			\
			template<> \
			void skygenHost<name>::upload_self(bool draw) \
			{ \
				cuxUploadConst(name##Sky, (skygenGPU<name>&)*this); \
			} \
			\
			template class skygenHost<name>
	#else
		#define MODEL_IMPLEMENTATION(name)
	#endif
#else
	#if BUILD_FOR_CPU && !SKYGEN_ON_GPU
		#define MODEL_IMPLEMENTATION(name) \
			__device__ __constant__ skygenGPU<name> name##Sky; \
			\
			template<> __global__ void integrate_counts_kernel<name>() { name##Sky.kernel<0>(); } \
			template<> __global__ void     draw_sources_kernel<name>() { name##Sky.kernel<1>(); } \
			\
			template<> \
			void skygenHost<name>::upload_self(bool draw) \
			{ \
				cuxUploadConst(name##Sky, (skygenGPU<name>&)*this); \
			} \
			\
			template class skygenHost<name>
	#else
		#define MODEL_IMPLEMENTATION(name)
	#endif
#endif

#endif
