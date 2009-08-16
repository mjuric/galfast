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
	//inline float Rg() { return Rg_cpu; }
	#include "paralax.h"
#endif

#if 0
struct ALIGN(16) lfParams
{
	float M0, M1, inv_dM;
	int isLoaded;

	template<typename T>
	__device__ float sample(T &texref, float M) const
	{
		float val, iMr;
/*		if(!isLoaded)
		{
			val = 1.f;
		}
		else
		{
*/			iMr  = (M  -  M0) * inv_dM + 0.5;
			val = tex1D(texref, iMr);
//		}
		return val;
	}
};
inline lfParams make_lfParams(float M0, float M1, float dM)
{
	int isLoaded = 0; //M0 < M1;
	lfParams lfp = { M0, M1, 1.f/dM, isLoaded };
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
	lfParams load(const char *fn);
	lfParams loadConstant(float val);
	void bind();
	void free();
};

// luminosity function texture and manager
#ifdef __CUDACC__
texture<float, 1, cudaReadModeElementType> texLF(false, cudaFilterModeLinear, cudaAddressModeClamp);
#else
extern lfTextureManager texLFMgr;
#endif
#endif


	////////////////////////////////////////////////
	//  Texturing support, next generation
	////////////////////////////////////////////////

	// host interface
	template<typename T, int dim = 1>
	class cuxTexture
	{
	public:
		xptrng::xptr<T> m_data;
		float x0[dim], dx[dim];

		void init(T *data, int n0, float x0, float dx0, int n1 = 0, float x1 = 0, float dx1 = 0, int n2 = 0, float x2 = 0, float dx2 = 0);
	};

	template<typename T, int dim>
	void cuxTexture<T, dim>::init(T *data, int n0, float x0a, float dx0, int n1, float x1, float dx1, int n2, float x2, float dx2)
	{
		if(dim >= 1) { x0[0] = x0a; dx[0] = dx0; }
		if(dim >= 2) { x0[1] = x1;  dx[1] = dx1; }
		if(dim >= 3) { x0[2] = x2;  dx[2] = dx2; }

		dim3 d = {n0, n1, n2};
		m_data = xptrng::xptr<T>(d);
		m_data.copyFrom(data);
	};

#if !__CUDACC__
	template<typename T, int dim = 1, enum cudaTextureReadMode mode = cudaReadModeElementType>
	class texture
	{
	};
#endif

/*
	// GPU interface
	template<typename T, int dim = 1, enum cudaTextureReadMode mode = cudaReadModeElementType>
	class cuxTextureReference
	{
	public:
		typedef texture<T, dim, mode> gpu_t;

	protected:
		float x0[dim], dx[dim];

	public:
		// binds the texture, uploads texture coordinate offset and dx parameters
		bool bind(cuxTexture<T, dim> &texref);

		__device__ T sample(reference &r, float x)
		{
			float xx = x0[0] + dx[0]*x;
			return tex1D(r, xx);
		}
	};
*/

/*
	Work around CUDA defficiency with some built-in struct alignments.

	CUDA header files declare some structs (float2 being an example) with
	__builtin_align() attribute that resolves to __align__ only when using
	CUDACC. This makes those structure's alignments different in nvcc compiled
	object code, compared to GCC's. Example: with nvcc, float2 is 8-byte
	aligned; on gcc, it's 4-byte aligned (given its members are all 4-byte
	aligned floats). Therefore, a struct that has a float2 member may be
	packed differently on gcc and nvcc. Example: struct { float a; float2 b; };
	On nvcc, &b = &a + 8 (in bytes). On gcc, &b = &a + 4 (bytes).

	This cludge works around the problem by deriving an aligned type from
	the problematic CUDA type. It should be used instead of the CUDA type
	in structures where this problem may occur.
*/
struct ALIGN(8) afloat2 : public float2 {};

template<typename T, int dim, enum cudaTextureReadMode mode>
	static inline __device__ T sample(texture<T, dim, mode> &r, float x, float2 xx)
	{
		float xi = (x - xx.x) * xx.y + 0.5;
		T y = tex1D(r, xi);
		//T y = 0.01f;
#if __DEVICE_EMULATION__
//		printf("phi=%f\n", y);
#endif
		return y;
	}



#ifdef __CUDACC__
texture<float, 1, cudaReadModeElementType> expModelLF(false, cudaFilterModeLinear, cudaAddressModeClamp);
#endif

// double-exponential+powerlaw Halo model
struct ALIGN(16) expModel
{
	struct ALIGN(16) host_state_t
	{
		xptrng::xptr<float> lf;
	};
	struct state
	{
		float rho;
	};
	float rho0, l, h, z0, f, lt, ht, fh, q, n;
	float r_cut2;
	afloat2 lf;

	// Management functions
/*	void setmodel(float l_, float h_, float z0_, float f_, float lt_, float ht_, float fh_, float q_, float n_)
	{
		l = l_; h = h_; z0 = z0_; f = f_; lt = lt_; ht = ht_; fh = fh_; q = q_; n = n_;
	}*/
	void load(host_state_t &hstate, const peyton::system::Config &cfg);
	void prerun(host_state_t &hstate, bool draw);
	void postrun(host_state_t &hstate, bool draw);

protected:
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
	__device__ float rho(float r, float z)       const
	{
		if(sqr(r) + sqr(z) > r_cut2) { return 0.f; }

		float rho = rho0 * (rho_thin(r, z) + rho_thick(r, z) + rho_halo(r, z));
		//fprintf(stderr, "rho=%f\n", rho);
		return rho;
	}

public:
	__device__ float rho(float x, float y, float z, float M) const
	{
		float rh = rho(sqrtf(x*x + y*y), z);
		//fprintf(stderr, "rho=%f\n", rh);
		return rh;
	}

#ifdef __CUDACC__
	__device__ void setpos(state &s, float x, float y, float z) const
	{
//		x = 7720.f; y= -770.1f; z = 2252.f;
//		((float*)shmem)[threadIdx.x] = rho(x, y, z, 0.f);
		s.rho = rho(x, y, z, 0.f);
#if __DEVICE_EMULATION__
//		printf("rho_A=%f\n", s.rho);
#endif
	}

	__device__ float rho(state &s, float M) const
	{
//		return 1.f;
//		return 1.f * ((float*)shmem)[threadIdx.x];
//		return 0.05f * s.rho;
//		M = 5.80;
		float phi = sample(expModelLF, M, lf);
#if __DEVICE_EMULATION__
//		printf("phi=%f rho=%f\n", phi, phi*s.rho);
#endif
		return phi * s.rho;
	}

	static const int THIN  = 0;
	static const int THICK = 1;
	static const int HALO  = 2;
	__device__ int component(float x, float y, float z, float M, gpuRng::constant &rng) const
	{
		float r = sqrtf(x*x + y*y);

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


// exponential disk model
#ifdef __CUDACC__
texture<float, 1, cudaReadModeElementType> expDiskLF(false, cudaFilterModeLinear, cudaAddressModeClamp);
#endif

struct ALIGN(16) expDisk
{
public:
	struct ALIGN(16) host_state_t
	{
		xptrng::xptr<float> lf;
	};

protected:
	float f, l, h, z0;
	float r_cut2;

	// luminosity function
	afloat2 lf;

	int comp;

public:
	struct state
	{
		float rho;
	};
	void load(host_state_t &hstate, const peyton::system::Config &cfg);
	void prerun(host_state_t &hstate, bool draw);
	void postrun(host_state_t &hstate, bool draw);

protected:
	__device__ float rho(float x, float y, float z, float M) const
	{
		float r2 = x*x + y*y;
		if(r2 + z*z > r_cut2) { return 0.f; }

		float r = sqrtf(r2);
		float rho = f * expf((Rg()-r)/l  + (fabsf(z0) - fabsf(z + z0))/h);

		return rho;
	}

public:
#ifdef __CUDACC__
	__device__ void setpos(state &s, float x, float y, float z) const
	{
		s.rho = rho(x, y, z, 0.f);
	}

	__device__ float rho(state &s, float M) const
	{
		float phi = sample(expDiskLF, M, lf);
		return phi * s.rho;
	}

	__device__ int component(float x, float y, float z, float M, gpuRng::constant &rng) const
	{
		return comp;
	}
#endif
};

static const double dbl_pi  = 3.14159265358979323846264338;
static const double dbl_d2r = 0.01745329251994329576923691; // (pi/180.0)
static const double dbl_r2d = 57.2957795130823208767981548; // (180.0/pi)
static const double flt_pi  = 3.14159265358979323846264338;
static const double flt_d2r = 0.01745329251994329576923691; // (pi/180.0)

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

	skypixel() {}
	skypixel(Radians l_, Radians b_, int projIdx_, float dx_, float coveredFraction_)
	: direction(l_, b_), projIdx(projIdx_), dx(dx_), dA(dx_*dx_), coveredFraction(coveredFraction_)
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
	cfloat_t::gpu_t		DM, M, XYZ;
	cint_t::gpu_t		comp;
};

template<typename Model>
struct ALIGN(16) runtime_state
{
	typedef typename Model::state ms_t;

	cux_ptr<int> cont, ilb, im, iM, k, bc;
	cux_ptr<float3> pos;
	cux_ptr<float> D;
	cux_ptr<skypixel> pix;
	cux_ptr<ms_t> ms;

	void alloc(int nthreads)
	{
		cont.alloc(nthreads); ilb.alloc(nthreads); im.alloc(nthreads); iM.alloc(nthreads); k.alloc(nthreads); bc.alloc(nthreads);
		pos.alloc(nthreads); D.alloc(nthreads); pix.alloc(nthreads);
		ms.alloc(nthreads);

		cudaMemset(cont.ptr, 0, nthreads*4); // not continuing a previous run
	}
	void free()
	{
		cont.free(); ilb.free(); im.free(); iM.free(); k.free(); bc.free();
		pos.free(); D.free(); pix.free();
		ms.free();
	}
	void constructor()	// as CUDA doesn't allow real constructors
	{
		cont = ilb = im = iM = k = bc = 0;
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

	__device__ void load(int &tid, int &ilb, int &im, int &iM, int &k, int &bc, float3 &pos, float &D, skypixel &pix, typename Model::state &ms) const
	{
		ilb = this->ilb(tid);
		im  = this->im(tid);
		iM  = this->iM(tid);
		k   = this->k(tid);
		bc  = this->bc(tid);
		pos = this->pos(tid);
		pix = this->pix(tid);
		D   = this->D(tid);
		ms  = this->ms(tid);
	}

	__device__ void store(int tid, int ilb, int im, int iM, int k, int bc, float3 pos, float D, skypixel pix, typename Model::state ms) const
	{
		this->ilb(tid) = ilb;
		this->im(tid)  = im;
		this->iM(tid)  = iM;
		this->k(tid)   = k;
		this->bc(tid)  = bc;
		this->pos(tid) = pos;
		this->pix(tid) = pix;
		this->D(tid)   = D;
		this->ms(tid)  = ms;

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

template<typename Model>
struct ALIGN(16) skyConfigGPU : public skygenConfig
{
	typedef Model Model_t;

	Model_t model;

	cux_ptr<skypixel> pixels;	// pixels on the sky to process
//	float dx, dA;			// linear scale and angular area of each pixel (rad, rad^2)

	cux_ptr<int> lock;
	cux_ptr<int> nstars;
	cux_ptr<float> counts, countsCovered;
	ocolumns stars;

	static const int nhistbins = 10;
	cux_ptr<int> rhoHistograms;	// nthreads*nbins sized array
	float lrho0, dlrho;		// histogram array start (bin midpoint), bin size

	cux_ptr<float> maxCount;		// [nthreads] sized array, returning the maximum density found by each thread

	runtime_state<Model> ks;

	template<int draw> __device__ void kernel() const;
	__device__ float3 compute_pos(float &D, float M, const int im, const direction &dir) const;
	__device__ bool advance(int &ilb, int &i, int &j, skypixel &pix, const int x, const int y) const;
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
	dim3 gridDim, blockDim;
	int shb;
	stopwatch swatch;
	int Kbatch;

	int bufferSafetyMargin();
	void upload_generic(bool draw = false);

	void compute(bool draw = false);
	void upload(bool draw = false);
	void download(bool draw = false);

	skyConfig();
	~skyConfig();

	// external interface
	virtual void initRNG(rng_t &rng);	// initialize the random number generator from CPU RNG
	virtual double integrateCounts();	// return the expected starcounts contributed by this model
	virtual bool init(
		otable &t,
		const peyton::system::Config &cfg,	// model cfg file
		const skygenConfig &sc,
		const skypixel *pixels);
/*	virtual bool init(
			const peyton::system::Config &cfg,
			const peyton::system::Config &pdf_cfg,
			const peyton::system::Config &foot_cfg,
			const peyton::system::Config &model_cfg,
			otable &t,
			opipeline &pipe);*/
	virtual size_t run(otable &in, osink *nextlink);
};

#if _EMU_DEBUG
extern long long voxelsVisited;
#endif

#endif
