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

#define _EMU_DEBUG (_DEBUG && __DEVICE_EMULATION__)

#include <cuda_runtime_api.h>
#include "gpu.h"
#include "column.h"
#include "simulate_base.h"
#include <astro/types.h>

typedef prngs::gpu::mwc gpuRng;
using peyton::Radians;

#ifdef __CUDACC__
	// distance to the Galactic center
	__device__ __constant__ float Rg_gpu;
	__device__ inline float Rg() { return Rg_gpu; }
#else
	#include "paralax.h"
#endif

// these must be overloaded by models
struct modelConcept
{
	struct host_state_t {};
	struct state {};

	void load(host_state_t &hstate, const peyton::system::Config &cfg);
	void prerun(host_state_t &hstate, bool draw);
	void postrun(host_state_t &hstate, bool draw);

	__device__ void setpos(state &s, float x, float y, float z) const;
	__device__ float rho(state &s, float M) const;
	__device__ int component(float x, float y, float z, float M, gpuRng::constant &rng) const;
};

static const double dbl_pi  = 3.14159265358979323846264338;
static const double dbl_d2r = 0.01745329251994329576923691; // (pi/180.0)
static const double dbl_r2d = 57.2957795130823208767981548; // (180.0/pi)
static const float flt_pi  = 3.14159265358979323846264338f;
static const float flt_d2r = 0.01745329251994329576923691f; // (pi/180.0)
static const float flt_r2d = 57.2957795130823208767981548f; // (180.0/pi)

namespace cudacc
{
	inline __device__ float deg(float rd) { return flt_r2d * rd; }
	inline __device__ float rad(float dg) { return flt_d2r * dg; }
	template<typename T>
		__host__ __device__ inline __device__ T sqr(const T x) { return x*x; }
}

// inline double rad(double dgr) { return dgr * dbl_d2r; }
// inline double deg(double rd)  { return rd  / dbl_d2r; }

struct direction
{
	float cl, cb, sl, sb;

	direction() {}
	direction(Radians l_, Radians b_)
	: cl(cos(l_)), cb(cos(b_)), sl(sin(l_)), sb(sin(b_))
	{ }
};

struct ALIGN(16) skypixel : public direction
{
	int projIdx;
	float coveredFraction;
	float dx, dA;
	float X, Y;		// lambert coordinates of this direction (in projIdx projection)

	skypixel() {}
	skypixel(Radians l_, Radians b_, float X_, float Y_, int projIdx_, float dx_, float coveredFraction_)
	: direction(l_, b_), X(X_), Y(Y_), projIdx(projIdx_), dx(dx_), dA(dx_*dx_), coveredFraction(coveredFraction_)
	{ }
};

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

	__device__ void convert(const direction &d, float &x, float &y) const
	{
		double cll0 = d.cl*p.cl + d.sl*p.sl; // cos(l - l0);
		double sll0 = p.cl*d.sl - d.cl*p.sl; // sin(l - l0);
		double kp = 1.414213562373095 * sqrt(1./(1.+p.sb*d.sb+p.cb*d.cb*cll0));
		x = kp*d.cb*sll0;
		y = kp*(p.cb*d.sb-p.sb*d.cb*cll0);
	}
	
	__device__ void inverse(const float x, const float y, double &l, double &b) const
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

__device__ inline float3 position(const direction &p, const float d)
{
	float3 ret;

	ret.x = Rg() - d*p.cl*p.cb;
	ret.y =      - d*p.sl*p.cb;
	ret.z =        d*p.sb;

	//rcyl = sqrt(x*x + y*y);
	return ret;
};

struct ALIGN(16) ocolumns
{
	cdouble_t::gpu_t	lb;
	cint_t::gpu_t		projIdx;
	cfloat_t::gpu_t		projXY, DM, M, XYZ, Am;
	cint_t::gpu_t		comp;
};

template<typename Model>
struct ALIGN(16) runtime_state
{
	typedef typename Model::state ms_t;

	cuxDevicePtr<int> cont, ilb, im, iM, k, bc, ndraw;
	cuxDevicePtr<float3> pos;
	cuxDevicePtr<float> D, Am;
	cuxDevicePtr<skypixel> pix;
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
		pos = 0; D = 0; pix = 0;
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

	__device__ void load(int &tid, int &ilb, int &im, int &iM, int &k, int &bc, float3 &pos, float &D, skypixel &pix, float &Am, typename Model::state &ms, int &ndraw) const
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

	__device__ void store(int tid, int ilb, int im, int iM, int k, int bc, float3 pos, float D, skypixel pix, float Am, typename Model::state ms, int ndraw) const
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

struct ALIGN(16) skygenConfig
{
	float m0, m1, dm;		// apparent magnitude range
	float M0, M1, dM;		// absolute magnitude range
	int npixels, nm, nM;
	int nthreads;			// total number of threads processing the sky
	int stopstars;			// stop after this many stars have been generated

	lambert proj[2];		// north/south sky lambert projections
};

//
// Device functions and data for skygen (instantiated/defined in skygen.cu.h)
// This piece gets uploaded as a __constant__ to the GPU
//
template<typename Model>
struct ALIGN(16) skyConfigGPU : public skygenConfig
{
	typedef Model Model_t;

	Model_t model;

	cuxDevicePtr<skypixel> pixels;	// pixels on the sky to process

	cuxDevicePtr<int> lock;
	cuxDevicePtr<int> nstars;
	cuxDevicePtr<float> counts, countsCovered;
	ocolumns stars;

	static const int nhistbins = 10;
	cuxDevicePtr<int> rhoHistograms;	// nthreads*nbins sized array
	float lrho0, dlrho;		// histogram array start (bin midpoint), bin size

	cuxDevicePtr<float> maxCount;	// [nthreads] sized array, returning the maximum density found by each thread

	runtime_state<Model> ks;
	float norm;			// normalization of overall density (usually 1.f)

	template<int draw> __device__ void kernel() const;
	__device__ float3 compute_pos(float &D, float &Am, float M, const int im, const skypixel &dir) const;
	__device__ bool advance(int &ilb, int &i, int &j, skypixel &pix, const int x, const int y) const;
	__device__ void draw_stars(int &ndraw, const float &M, const int &im, const skypixel &pix) const;
};

class opipeline;
template<typename Model>
struct ALIGN(16) skyConfig : public skyConfigGPU<Model>, public skyConfigInterface
{
	// Any state the model needs to maintain on the host (e.g., loaded textures)
	typename Model::host_state_t	model_host_state;

	float Rg;		// distance to the galactic center

	gpuRng *rng;
	rng_t *cpurng;
	unsigned seed;
	skypixel *cpu_pixels;

	// return
	float nstarsExpectedToGenerate;	// expected number of stars in the pixelized footprint
	float nstarsExpected;		// expected number of stars in the actual footprint
	int stars_generated;
	int *cpu_hist;
	float *cpu_maxCount;
	int3 *cpu_state;

	// kernel execution configuration
	dim3 gridDim, blockDim;		// CUDA grid dimension, block dimension
	int shb;			// shared memory per block needed by skygen kernel
	int output_table_capacity;	// number of rows in the output table
	stopwatch swatch;

	int bufferSafetyMargin();
	void upload_self(bool draw = false);

	void compute(bool draw = false);
	void upload(bool draw = false);
	void download(bool draw = false);

	skyConfig();
	~skyConfig();

	// external interface
	virtual void initRNG(rng_t &rng);	// initialize the random number generator from CPU RNG
	virtual double integrateCounts();	// return the expected starcounts contributed by this model
	virtual void setDensityNorm(float norm);// set the density normalization of the model
	virtual bool init(
		otable &t,
		const peyton::system::Config &cfg,	// model cfg file
		const skygenConfig &sc,
		const skypixel *pixels);
	virtual size_t run(otable &in, osink *nextlink);
};

#if _EMU_DEBUG
extern long long voxelsVisited;
#endif

//
// Some generic stub code to instantiate the GPU kernels, the class,
// and a specialized upload method for the model in question
//
// Active only when compiled with nvcc
//
#ifdef __CUDACC__
	#define MODEL_IMPLEMENTATION(name) \
		__device__ __constant__ skyConfigGPU<name> name##Sky; \
		\
		template<> __global__ void compute_sky<name>() { name##Sky.kernel<0>(); } \
		template<> __global__ void    draw_sky<name>() { name##Sky.kernel<1>(); } \
		\
		template<> \
		void skyConfig<name>::upload_self(bool draw) \
		{ \
			cuxUploadConst(name##Sky, (skyConfigGPU<name>&)*this); \
		} \
		\
		template class skyConfig<name>
#else
	#define MODEL_IMPLEMENTATION(name)
#endif

#endif
