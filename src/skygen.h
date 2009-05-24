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
#include <astro/types.h>

typedef prngs::gpu::mwc gpuRng;
using peyton::Radians;

__device__ __constant__ float Rg_gpu;
#ifdef __CUDACC__
// distance to the Galactic center
__device__ inline float Rg() { return Rg_gpu; }
#endif

struct lfParams
{
	float M0, M1, inv_dM;
	float padding;		// to get this properly aligned both for CUDA and GCC (there has to be a smarter way to do this!)

	template<typename T>
	__device__ float sample(T &texref, float M) const
	{
		float val, iMr;
/*		for(M=M0; M <= M1; M += 1./(inv_dM*5.))
		{*/
			iMr  = (M  -  M0) * inv_dM + 0.5;
			val = tex1D(texref, iMr);
/*			printf("%f %f %f\n", M, iMr, val);
		}*/
		return val;
	}
};
inline lfParams make_lfParams(float M0, float M1, float dM)
{
	lfParams lfp = { M0, M1, 1.f/dM };
	return lfp;
}

struct lfTextureManager
{
	lfParams par;
	textureReference &texref;
	cudaArray *lfArray;

	lfTextureManager(textureReference &tr) : texref(tr), lfArray(0) { par = make_lfParams(0, 0, 1.); }
	~lfTextureManager() { free(); }

	lfParams set(float *cpu_lf, int lfLen, float M0, float M1, float dM);
	void free();
};

// luminosity function texture and manager
#ifdef __CUDACC__
texture<float, 1, cudaReadModeElementType> texLF(false, cudaFilterModeLinear, cudaAddressModeClamp);
#else
extern lfTextureManager texLFMgr;
#endif

// exponential model
struct ALIGN(16) expModel
{
	struct state {
		float rho;
	};
	float l, h, z0, f, lt, ht, fh, q, n;
	lfParams lf;

	// Management functions
	void setmodel(float l_, float h_, float z0_, float f_, float lt_, float ht_, float fh_, float q_, float n_)
	{
		l = l_; h = h_; z0 = z0_; f = f_; lt = lt_; ht = ht_; fh = fh_; q = q_; n = n_;
	}

 	// Model functions
 	__device__ float sqr(float x) const { return x*x; }
	__device__ float halo_denom(float r, float z) const { return sqr(r) + sqr(q*(z + z0)); }

	__device__ float rho_thin(float r, float z)  const
	{
		float rho = expf((Rg()-r)/l  + (fabsf(z0) - fabsf(z + z0))/h);
		//fprintf(stderr, "rho=%f\n", rho);
		return rho;
	}
	__device__ float rho_thick(float r, float z) const
	{
		float rho = f * expf((Rg()-r)/lt + (fabsf(z0) - fabsf(z + z0))/ht);
		//fprintf(stderr, "rho=%f\n", rho);
		return rho;
	}
	__device__ float rho_halo(float r, float z)  const
	{
		float rho = fh * powf(Rg()/sqrtf(halo_denom(r,z)),n);
		//fprintf(stderr, "rho=%f\n", rho);
		return rho;
	}
	__device__ float rho(float r, float z)       const {
		float rho = rho_thin(r, z) + rho_thick(r, z) + rho_halo(r, z);
		//fprintf(stderr, "rho=%f\n", rho);
		return rho;
	}

	__device__ float rho(float x, float y, float z, float M) const
	{
		float rh = rho(sqrtf(x*x + y*y), z);
		//fprintf(stderr, "rho=%f\n", rh);
		return rh;
	}

#ifdef __CUDACC__
	__device__ void setpos(expModel::state &s, float x, float y, float z) const
	{
//		((float*)shmem)[threadIdx.x] = rho(x, y, z, 0.f);
		s.rho = rho(x, y, z, 0.f);
	}

	__device__ float rho(expModel::state &s, float M) const
	{
//		return 1.f;
//		return 1.f * ((float*)shmem)[threadIdx.x];
//		return 0.05f * s.rho;
		return lf.sample(texLF, M) * s.rho;
	}

	static const int THIN  = 0;
	static const int THICK = 1;
	static const int HALO  = 2;
	__device__ int component(float x, float y, float z, float M, gpuRng &rng) const
	{
		float r = sqrt(x*x + y*y);

		float thin = rho_thin(r, z);
		float thick = rho_thick(r, z);
		float halo = rho_halo(r, z);
		float rho = thin+thick+halo;

		float pthin  = thin / rho;
		float pthick = (thin + thick) / rho;

		float u = rng.uniform();
		if(u < pthin) { return THIN; }
		else if(u < pthick) { return THICK; }
		else { return HALO; }
	}
#endif
};

static const double dbl_pi  = 3.14159265358979323846264338;
static const double dbl_d2r = 0.01745329251994329576923691; // (pi/180.0)
static const double flt_pi  = 3.14159265358979323846264338;
static const double flt_d2r = 0.01745329251994329576923691; // (pi/180.0)

inline double rad(double dgr) { return dgr * dbl_d2r; }
inline double deg(double rd)  { return rd  / dbl_d2r; }

struct direction
{
	float cl, cb, sl, sb;

	direction() {}
	direction(Radians l_, Radians b_)
	: cl(cos(l_)), cb(cos(b_)), sl(sin(l_)), sb(sin(b_))
	{ }
};

class lambert
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
		float cll0 = d.cl*p.cl + d.sl*p.sl; // cos(l - l0);
		float sll0 = p.cl*d.sl - d.cl*p.sl; // sin(l - l0);
		float kp = 1.414213562373095f * sqrtf(1./(1+p.sb*d.sb+p.cb*d.cb*cll0));
		x = kp*d.cb*sll0;
		y = kp*(p.cb*d.sb-p.sb*d.cb*cll0);
	}

	__device__ void inverse(const float x, const float y, double &l, double &b) const
	{
		float r2 = x*x + y*y;
		if(r2 == 0.f)
		{
			b = b0;
			l = 0.;
			return;
		}
		float ir = 1.f/sqrtf(r2);			// ir == 1/r
		float sc = 1.f/ir * sqrtf(1.f - 0.25f*r2);	// sin(c) == r*Sqrt(1 - Power(r,2)/4.)		where c = 2*asin(0.5*r)
		float cc = 1.f - 0.5*r2;			// cos(c) == 1 - Power(r,2)/2			where c = 2*asin(0.5*r)

		b = asin(cc*p.sb + y*sc*p.cb*ir);
		l = l0 + atan2(x*sc, p.cb*cc/ir - y*p.sb*sc);
	}
};

__device__ float3 position(const direction &p, const float d)
{
	float3 ret;

	ret.x = Rg() - d*p.cl*p.cb;
	ret.y =      - d*p.sl*p.cb;
	ret.z =        d*p.sb;

	//rcyl = sqrt(x*x + y*y);
	return ret;
};

/*struct star
{
	Radians l, b;
	float m, M;
	float D;
	float3 pos;
	int comp;
	int alignment_dummy;
};
*/
struct ALIGN(16) ocolumns
{
	column_types::cdouble::gpu_t	lb;
	column_types::cfloat::gpu_t	m, M, XYZ;
	column_types::cint::gpu_t	comp;
};

template<typename Model>
struct ALIGN(16) runtime_state
{
	typedef typename Model::state ms_t;

	cux_ptr<int> cont, ilb, im, iM, k, bc;
	cux_ptr<float3> pos;
	cux_ptr<float> D;
	cux_ptr<direction> dir;
	cux_ptr<ms_t> ms;

	void alloc(int nthreads)
	{
		cont.alloc(nthreads); ilb.alloc(nthreads); im.alloc(nthreads); iM.alloc(nthreads); k.alloc(nthreads); bc.alloc(nthreads);
		pos.alloc(nthreads); D.alloc(nthreads); dir.alloc(nthreads);
		ms.alloc(nthreads);

		cudaMemset(cont.ptr, 0, nthreads*4); // not continuing a previous run
	}
	void free()
	{
		cont.free(); ilb.free(); im.free(); iM.free(); k.free(); bc.free();
		pos.free(); D.free(); dir.free();
		ms.free();
	}
	void constructor()	// as CUDA doesn't allow real constructors
	{
		cont = ilb = im = iM = k = bc = 0;
		pos = 0; D = 0; dir = 0;
		ms = 0;
	}
	void reset()
	{
		free();
	}

	__device__ void load(int &tid, int &ilb, int &im, int &iM, int &k, int &bc, float3 &pos, float &D, direction &dir, typename Model::state &ms) const
	{
		ilb = this->ilb[tid];
		im  = this->im[tid];
		iM  = this->iM[tid];
		k   = this->k[tid];
		bc  = this->bc[tid];
		pos = this->pos[tid];
		dir = this->dir[tid];
		D   = this->D[tid];
		ms  = this->ms[tid];
	}

	__device__ void store(int tid, int ilb, int im, int iM, int k, int bc, float3 pos, float D, direction dir, typename Model::state ms) const
	{
		this->ilb[tid] = ilb;
		this->im[tid]  = im;
		this->iM[tid]  = iM;
		this->k[tid]   = k;
		this->bc[tid]  = bc;
		this->pos[tid] = pos;
		this->dir[tid] = dir;
		this->D[tid]   = D;
		this->ms[tid]  = ms;

		this->cont[tid] = 1;
	}
	
	bool continuing(int tid) const { return this->cont[tid]; }
};

template<typename Model>
struct ALIGN(16) skyConfigGPU
{
	typedef Model Model_t;

	Model_t model;

	cux_ptr<direction> pixels;	// pixels on the sky to process
	float dx, dA;			// linear scale and angular area of each pixel (rad, rad^2)
	float m0, m1, dm;		// apparent magnitude range
	float M0, M1, dM;		// absolute magnitude range
	int npixels, nm, nM;
	int nthreads;			// total number of threads processing the sky
	int stopstars;			// stop after this many stars have been generated

	cux_ptr<int> lock;
	//cux_ptr<star> stars;
	cux_ptr<int> nstars;
	cux_ptr<float> counts;
	ocolumns stars;

	static const int nhistbins = 10;
	cux_ptr<int> rhoHistograms;	// nthreads*nbins sized array
	float lrho0, dlrho;		// histogram array start (bin midpoint), bin size

	cux_ptr<float> maxCount;		// [nthreads] sized array, returning the maximum density found by each thread

	runtime_state<Model> ks;

	template<int draw> __device__ void kernel() const;
	__device__ float3 compute_pos(float &D, float M, const int im, const direction &dir) const;
	__device__ bool advance(int &ilb, int &i, int &j, direction &dir, const int x, const int y) const;
};

int computeSafetyBuffer(float *cpu_maxCount, int nthreads);

template<typename Model>
struct skyConfig : public skyConfigGPU<Model>
{
	float Rg;	// distance to the galactic center

	gpuRng *rng;
	unsigned seed;
	direction *cpu_pixels;
	lambert proj;

	// return
	float nstarsExpected;
/*	star *cpu_stars;*/
	int stars_generated;
	int *cpu_hist;
	float *cpu_maxCount;
	int3 *cpu_state;

	// kernel execution configuration
	dim3 gridDim, blockDim;
	int shb;
	stopwatch swatch;
	int Kbatch;

	int bufferSafetyMargin();

	void compute(bool draw = false);
	void upload(bool draw = false);
	void download(bool draw = false);
	void loadConfig(rng_t &seeder);

	skyConfig()  : rng(NULL), cpu_state(NULL), cpu_maxCount(NULL), /*cpu_stars(NULL), */cpu_pixels(NULL), cpu_hist(NULL) { this->pixels = 0; this->lock = 0; this->counts = 0; this->nstars = 0; this->ks.constructor(); }
	~skyConfig()
	{
		delete [] cpu_pixels;
/*		delete [] cpu_stars;*/
		delete [] cpu_hist;
		delete [] cpu_maxCount;
		delete rng;

		this->lock.free();
		//this->stars.free();
		this->pixels.free();
		this->counts.free();
	}
};

#if _DEBUG
extern long long voxelsVisited;
#endif

#endif
