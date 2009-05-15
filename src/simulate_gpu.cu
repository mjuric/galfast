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

#ifndef __CONFIG_H
#define __CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include <math.h>

#include <astro/constants.h>

#include "column.h"
#include "gpu2.h"

#include "simulate_base.h"

namespace ct = column_types;

//======================================================================
//======================================================================
//    FeH
//======================================================================
//======================================================================

KERNEL(
	ks, 3*4,
	os_FeH_kernel(
		otable_ks ks, os_FeH_data par, gpu_rng_t rng, 
		ct::cint::gpu_t comp, 
		ct::cfloat::gpu_t XYZ, 
		ct::cfloat::gpu_t FeH),
	os_FeH_kernel,
	(ks, par, rng, comp, XYZ, FeH)
)
{
	rng.load(ks);
	uint32_t tid = threadID();
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
//		FeH[row] = -((threadID() % 300) / 100.f);
// 		float temp = -rng.uniform();
// 		FeH[row] = temp;
// 		continue; 

//		feh = -comp[row]; FeH[row] = feh; continue;
		float feh;
		int component = comp[row];
#if 1
/*		if (component==0) //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			feh=0;
		else*/
		if(component < 2)
		{
			// choose the gaussian to draw from
			float p = rng.uniform()*(par.A[0]+par.A[1]);
			int i = p < par.A[0] ? 0 : 1;

			// calculate mean
			float muD = par.muInf + par.DeltaMu*exp(-fabs(XYZ(row, 2))/par.Hmu);		// Bond et al. A2
			float aZ = muD - 0.067f;

			// draw
			feh = rng.gaussian(par.sigma[i]) + aZ + par.offs[i];			
		}
		else if(component == 2)
		{
			feh = par.offs[2] + rng.gaussian(par.sigma[2]);
		}
		else
		{
			feh = -9999.f;
		}
#else
		/*
			This chunk of code appears to compile incorrectly with CUDA 2.1 nvcc (V0.2.1221).
			Expected result of running with below: feh=1 when comp=1, 0 otherwise.
			Actual result: feh=comp when comp=1,2
		*/
		feh = float(component);
		//FeH[row] = feh;
		//continue;
		switch(component)
		{
/*			case 2: //BahcallSoneira_model::HALO:
				feh = -2.f;
				break;*/
			case 3: // BahcallSoneira_model::THIN:
				feh = -3.0f;
				break;
			case 1: // BahcallSoneira_model::THICK:
				feh = component-1;
				break;
			case 4: // BahcallSoneira_model::THIN:
				feh = -4.0f;
				break;
/*			default:
				//THROW(ENotImplemented, "We should have never gotten here");
				feh = -9.f;
				break;*/
		}
#endif
		FeH[row] = feh;
	}
	rng.store(ks);
}

//======================================================================
//======================================================================
//    vel2pm
//======================================================================
//======================================================================

inline __device__ double sqrDouble(double x) { return x*x;}

static const double degToRadCoef=0.0174532925;//PI/180
inline __device__ double radDouble(double degrees) { return degrees*degToRadCoef; }

inline __device__ void vel_cyl2xyz(
		float &vx, float &vy, float &vz, 
		const float vr, const float vphi, const float vz0, 
		const float X, const float Y)
{
	// convert galactocentric cylindrical to cartesian velocities
	float rho = sqrt(X*X + Y*Y);
	float cphi = X / rho;
	float sphi = Y / rho;

	vx = -sphi * vr + cphi * vphi;
	vy =  cphi * vr + sphi * vphi;
	vz = vz0;
}

// convert cartesian galactocentric velocities to vl,vr,vb wrt. the observer
inline __device__ void vel_xyz2lbr(
		float &vl, float &vr, float &vb, 
		const float vx, const float vy, const float vz, 
		const double l, const double b)
{
	double cl, sl, cb, sb;
	cl = cos(l);
	sl = sin(l);
	cb = cos(b);
	sb = sin(b);

	float tmp;
	vl  =  vx*sl  - vy*cl;
	tmp =  vx*cl  + vy*sl;
	vr  = -cb*tmp + vz*sb;
	vb  =  sb*tmp + vz*cb;
}

KERNEL(
	ks, 0, 
	os_vel2pm_kernel(
		otable_ks ks, os_vel2pm_data par, gpu_rng_t rng, 
		ct::cdouble::gpu_t lb0, 
		ct::cfloat::gpu_t XYZ,
		ct::cfloat::gpu_t vcyl,  
		ct::cfloat::gpu_t pmlb),
	os_vel2pm_kernel,
	(ks, par, rng, lb0, XYZ, vcyl,pmlb)
)
{

	for(uint32_t row=ks.row_begin(); row < ks.row_end(); row++)
	{
		// fetch prerequisites
		double l = radDouble(lb0(row, 0));
		double b = radDouble(lb0(row, 1));
		float X = XYZ(row, 0);
		float Y = XYZ(row, 1);
		float Z = XYZ(row, 2);
		float vx = vcyl(row, 0);
		float vy = vcyl(row, 1);
		float vz = vcyl(row, 2);

		// convert the velocities from cylindrical to galactocentric cartesian system
		float pm[3];
		vel_cyl2xyz(pm[0], pm[1], pm[2],   vx, vy, vz,   X, Y);

		// switch to Solar coordinate frame
		pm[0] -= par.u0;
		pm[1] -= par.v0 + par.vLSR;
		pm[2] -= par.w0;

		// convert to velocities wrt. the observer
		vel_xyz2lbr(pm[0], pm[1], pm[2],   pm[0], pm[1], pm[2],  l, b);

		// convert to proper motions
		float D = sqrt(sqrDouble(X) + sqrDouble(Y) + sqrDouble(Z));
		pm[0] /= 4.74 * D*1e-3;	// proper motion in mas/yr (4.74 km/s @ 1kpc is 1mas/yr)
		pm[1] /= 4.74 * D*1e-3;

		// rotate to output coordinate system
		switch(par.coordsys)
		{
		case GAL:
//			array_copy(s.pmlb(), pm, 3);
			pmlb(row, 0) = pm[0];
			pmlb(row, 1) = pm[1];
			pmlb(row, 2) = pm[2];
			break;
		case EQU:
			//THROW(EAny, "Output in equatorial system not implemented yet.");
			//array_copy(s.pmradec(), pm, 3);
			break;
		default:
			//THROW(EAny, "Unknown coordinate system [id=" + str(coordsys) + "] requested");
			break;
		}
	}
}


//======================================================================
//======================================================================
//    kinTMIII
//======================================================================
//======================================================================
float __device__ sqrf(float x) { return x*x;}	

#if 1

#define K_IO
#define K_OUT

void __device__ trivar_gaussK_draw(float y[3], float s11, float s12, float s13, float s22, float s23, float s33, gpu_rng_t &rng)
{		
	// NOTE: ADDS THE RESULT TO v, DOES NOT ZERO v BEFORE !!	
	float b10=s12;      float b11=sqrf(s22); 
	float b20=s13;      float b21=s23;      float b22=sqrf(s33);
	
	float a00=s11; float oneDivA00=1.0f/a00;
		float a10=( b10 ) * oneDivA00;	
	
	float a11=sqrtf( b11 - a10*a10 );
		float a20=( b20 ) * oneDivA00;
		float a21=(b21 - a20*a10)/a11;		

	float a22=sqrtf( b22 - a20*a20 - a21*a21);				
		
	float z0=rng.gaussian(1.0f);
	float z1=rng.gaussian(1.0f);
	float z2=rng.gaussian(1.0f);

	y[0]+=a00*z0;
	y[1]+=a10*z0+a11*z1;
	y[2]+=a20*z0+a21*z1+a22*z2;
}


__device__ inline float modfun(float Rsquared, float Z, farray5 val)
{
	return val[0] + val[1]*powf(fabs(Z), val[2]) + val[3]*powf(Rsquared, 0.5*val[4]);
}

__device__ void add_dispersion(K_IO float v[3], float Rsquared, float Z, farray5 ellip[6],K_IO gpu_rng_t &rng)
{
	// compute velocity dispersions at this position, and draw from trivariate gaussian
	// NOTE: ADDS THE RESULT TO v, DOES NOT ZERO v BEFORE !!
	float sigma[6];
	for (int i=0;i<6;i++)		
		sigma[i] = modfun(Rsquared, Z, ellip[i]);	
			
	trivar_gaussK_draw(v, sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5], rng);
}

__device__ void compute_means(K_OUT float v[3], float Rsquared, float Z, farray5 means[3])
{
	// returns means in v[3]
	for (int i=0;i<3;i++) 	
		v[i] = modfun(Rsquared, Z, means[i]);	
}

__device__ void get_disk_kinematics(K_OUT float v[3], float Rsquared, float Z, gpu_rng_t &rng,
			os_kinTMIII_data& par, farray5 diskMeans[3], farray5 diskEllip[6] )
{
	// set up which gaussian are we drawing from
	float p = rng.uniform();
	if(p < par.fk)
	{
		// first gaussian
		diskMeans[1] = par.vPhi1;
		diskEllip[3] = par.sigmaPhiPhi1;
	}
	else
	{
		// second gaussian
		diskMeans[1] = par.vPhi2;
		diskEllip[3] = par.sigmaPhiPhi2;
	}

	compute_means(v, Rsquared, Z, diskMeans);
	// truncate v_phi > 0
	if(v[1] > 0.) { v[1] = 0.; }

	add_dispersion(v, Rsquared, Z, diskEllip, rng);
}

__device__ void get_halo_kinematics(K_OUT float v[3], float Rsquared, float Z, K_IO gpu_rng_t &rng, farray5 haloMeans[3], farray5 haloEllip[6])
{
	compute_means(v, Rsquared, Z, haloMeans);
	add_dispersion(v, Rsquared, Z, haloEllip, rng);
}

__device__ void iarray_to_farray(farray5& fa, iarray5& ia)
{
	for (int i=0;i<5;i++)
		fa[i]=ia[i]/100.0f;
}

__device__ void i8array_to_farray(farray5& fa, i8array5& ia)
{
	for (int i=0;i<5;i++)
		fa[i]=ia[i]*10.0f;
}

__device__ __constant__ os_kinTMIII_data os_kinTMIII_par;

KERNEL(
	ks, 3*4,
	os_kinTMIII_kernel(
		otable_ks ks, gpu_rng_t rng, 
		ct::cint::gpu_t comp, 
		ct::cfloat::gpu_t XYZ, 
		ct::cfloat::gpu_t vcyl),
	os_kinTMIII_kernel,
	(ks, rng, comp, XYZ, vcyl)
)
{
	rng.load(ks); 
#define par os_kinTMIII_par
//	os_kinTMIII_data par;

// 	iarray_to_farray(par.vR, par_int.vR);
//     iarray_to_farray(par.vPhi1, par_int.vPhi1);	
//     iarray_to_farray(par.vZ, par_int.vZ);
//     iarray_to_farray(par.sigmaRR, par_int.sigmaRR);
//     iarray_to_farray(par.sigmaRPhi, par_int.sigmaRPhi);
//     iarray_to_farray(par.sigmaRZ, par_int.sigmaRZ);
//     iarray_to_farray(par.sigmaPhiPhi1, par_int.sigmaPhiPhi1);
//     iarray_to_farray(par.sigmaPhiPhi2, par_int.sigmaPhiPhi2);
//     iarray_to_farray(par.sigmaZPhi, par_int.sigmaZPhi);
//     iarray_to_farray(par.sigmaZZ, par_int.sigmaZZ);

	// v2 is v1 + DeltavPhi, which is what this does.
// 	par.vPhi2 = par.vPhi1;
// 	par.vPhi2[0] += par.DeltavPhi;

//   	i8array_to_farray(par.HvR, par_int.HvR);
// 	i8array_to_farray(par.HvPhi, par_int.HvPhi);
// 	i8array_to_farray(par.HvZ, par_int.HvZ);
// 	i8array_to_farray(par.HsigmaRR, par_int.HsigmaRR);
// 	i8array_to_farray(par.HsigmaRPhi, par_int.HsigmaRPhi);
// 	i8array_to_farray(par.HsigmaRZ, par_int.HsigmaRZ);
// 	i8array_to_farray(par.HsigmaPhiPhi, par_int.HsigmaPhiPhi);
// 	i8array_to_farray(par.HsigmaZPhi, par_int.HsigmaZPhi);
// 	i8array_to_farray(par.HsigmaZZ, par_int.HsigmaZZ);
/*
	par.fk=par_int.fk;*/

	farray5 diskEllip[6], haloEllip[6], diskMeans[3], haloMeans[3];

	diskMeans[0] = par.vR;
	diskMeans[1] = par.vPhi1;
	diskMeans[2] = par.vZ;
	diskEllip[0] = par.sigmaRR;
	diskEllip[1] = par.sigmaRPhi;
	diskEllip[2] = par.sigmaRZ;
	diskEllip[3] = par.sigmaPhiPhi1;
	diskEllip[4] = par.sigmaZPhi;
	diskEllip[5] = par.sigmaZZ;

	haloMeans[0] = par.HvR;
	haloMeans[1] = par.HvPhi;
	haloMeans[2] = par.HvZ;
	haloEllip[0] = par.HsigmaRR;
	haloEllip[1] = par.HsigmaRPhi;
	haloEllip[2] = par.HsigmaRZ;
	haloEllip[3] = par.HsigmaPhiPhi;
	haloEllip[4] = par.HsigmaZPhi;
	haloEllip[5] = par.HsigmaZZ;

	float tmp[3]; 
	uint32_t tid = threadID();
		// ASSUMPTIONS:
	//	- Fe/H exists in input
	//	- Apparent and absolute magnitude in the requested band exist in input
	for(size_t row=ks.row_begin(); row < ks.row_end(); row++)
	{
		// fetch prerequisites
		const int component = comp[row];
		float X = XYZ(row, 0);
		float Y = XYZ(row, 1);
		float Zpc = XYZ(row, 2);
		const float Rsquared = 1e-6 * (X*X + Y*Y);
		const float Z = 1e-3 * Zpc;

		switch(component) 
		{
			case BahcallSoneira_model_THIN:
			case BahcallSoneira_model_THICK:
				get_disk_kinematics(tmp, Rsquared, Z, rng, par, diskMeans, diskEllip);
				break;
			case BahcallSoneira_model_HALO:		
				get_halo_kinematics(tmp, Rsquared, Z, rng, haloMeans, haloEllip);
				break;
			//default:
				//THROW(ENotImplemented, "We should have never gotten here");
		}
		vcyl(row, 0) = tmp[0];
		vcyl(row, 1) = tmp[1];
		vcyl(row, 2) = tmp[2];
	}
	rng.store(ks);
#undef par
}
#endif

//======================================================================
//======================================================================
//    equgal
//======================================================================
//======================================================================

// equgal - Equatorial to Galactic coordinates
using namespace peyton;
typedef double Radians;
static const double angp = ctn::d2r * 192.859508333; //  12h 51m 26.282s (J2000)
//static const double dngp = ctn::d2r * 27.128336111;  // +27d 07' 42.01" (J2000)
static const double l0 = ctn::d2r * 32.932;
static const double ce = 0.88998740217659689; // cos(dngp)
static const double se = 0.45598511375586859; // sin(dngp)

inline __device__ double2 galequ(const double2 lb)
{
	const double cb = cos(lb.y);
	const double sb = sin(lb.y);
	const double cl = cos(lb.x-l0);
	const double sl = sin(lb.x-l0);

//	// TODO: These should be precomputed constants
//	const double ce = cos(dngp);
//	const double se = sin(dngp);

	double2 r;
	r.x = atan2(
			cb*cl,
			sb*ce-cb*se*sl
		) + angp;
	r.y = asin(cb*ce*sl + sb*se);

	while(r.x < 0.) { r.x += ctn::pi2; }
	return r;
}

KERNEL(
	ks, 0,
	os_gal2other_kernel(otable_ks ks, int coordsys, ct::cdouble::gpu_t lb0, ct::cdouble::gpu_t out),
	os_gal2other_kernel,
	(ks, coordsys, lb0, out)
)
{
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		double2 lb, ret;

		// convert to radians
		lb.x = lb0(row, 0) * ctn::d2r;
		lb.y = lb0(row, 1) * ctn::d2r;

		// rotate to output coordinate system
		switch(coordsys)
		{
		case EQU:
			ret = galequ(lb);
			break;
		default:
			ret.x = ret.y = -9999.;
			break;
		}

		// convert to degrees
		out(row, 0) = ret.x / ctn::d2r;
		out(row, 1) = ret.y / ctn::d2r;
	}
}

KERNEL(
	ks, 0,
	os_fixedFeH_kernel(otable_ks ks, float fixedFeH, ct::cfloat::gpu_t FeH),
	os_fixedFeH_kernel,
	(ks, fixedFeH, FeH)
)
{
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		FeH[row] = fixedFeH;
	}
}

#include <vector>

#if BUILD_FOR_CPU
extern __TLS std::vector<tptr<float> > *locuses;
extern __TLS std::vector<tptr<uint> >   *flags;
#else
__TLS std::vector<tptr<float> > *locuses;
__TLS std::vector<tptr<uint> >   *flags;
#endif

#if HAVE_CUDA && !BUILD_FOR_CPU
texture<float4, 2, cudaReadModeElementType> color0;
texture<float4, 2, cudaReadModeElementType> color1;
texture<float4, 2, cudaReadModeElementType> color2;
texture<float4, 2, cudaReadModeElementType> color3;
texture<uint4, 2, cudaReadModeElementType> cflags0;
texture<uint4, 2, cudaReadModeElementType> cflags1;
texture<uint4, 2, cudaReadModeElementType> cflags2;
texture<uint4, 2, cudaReadModeElementType> cflags3;

texture<float4, 2, cudaReadModeElementType> *colorTextures[] = { &color0, &color1, &color2, &color3 };
texture<uint4, 2, cudaReadModeElementType> *cflagsTextures[] = { &cflags0, &cflags1, &cflags2, &cflags3 };

__device__ uint fill(float *&colors, uint &flags, const float4 clr, const uint4 f, int &ncolors)
{
	*colors = clr.x; colors++; flags |= f.x; if(--ncolors == 0) return flags;
	*colors = clr.y; colors++; flags |= f.y; if(--ncolors == 0) return flags;
	*colors = clr.z; colors++; flags |= f.z; if(--ncolors == 0) return flags;
	*colors = clr.w; colors++; flags |= f.w; --ncolors; return flags;
}

__device__ uint sampleColors(float *colors, float FeH, float Mr, int ncolors)
{
	float4 clr;
	uint4 f;
	uint flags;

	clr = tex2D(color0, FeH, Mr); f = tex2D(cflags0, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) return flags;
	clr = tex2D(color1, FeH, Mr); f = tex2D(cflags1, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) return flags;
	clr = tex2D(color2, FeH, Mr); f = tex2D(cflags2, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) return flags;
	clr = tex2D(color3, FeH, Mr); f = tex2D(cflags3, FeH, Mr); fill(colors, flags, clr, f, ncolors); if(ncolors == 0) return flags;
	return 0xFFFFFFFF;
}

#else
uint sampleColors(float *colors, float FeH, float Mr, int ncolors)
{
	int f = (int)FeH;
	int m = (int)Mr;
//	std::cerr << "fm = " << f << " " << m << "   " << FeH << " " << Mr << "\n";

	uint fl = 0;
	for(int ic=0; ic != ncolors; ic++)
	{
		colors[ic] = (*locuses)[ic](f, m);
		fl |= (*flags)[ic](f, m);
	}
	return fl;
}
#endif

#if HAVE_CUDA && BUILD_FOR_CPU
#include <map>
std::map<std::string, xptrng::tptr<float4> > os_photometry_tex_c;
std::map<std::string, xptrng::tptr<uint4> >  os_photometry_tex_f;
void os_photometry_tex_get(const char *id, xptrng::tptr<float4> &c, xptrng::tptr<uint4> &f)
{
	c = os_photometry_tex_c[id];
	f = os_photometry_tex_f[id];
}
void os_photometry_tex_set(const char *id, xptrng::tptr<float4> &c, xptrng::tptr<uint4> &f)
{
	os_photometry_tex_c[id] = c;
	os_photometry_tex_f[id] = f;
}
#endif

#if !BUILD_FOR_CPU || !HAVE_CUDA

void os_photometry_tex_get(const char *id, xptrng::tptr<float4> &c, xptrng::tptr<uint4> &f);
void os_photometry_tex_set(const char *id, xptrng::tptr<float4> &c, xptrng::tptr<uint4> &f);

void os_photometry_set_isochrones(const char *id, std::vector<tptr<float> > *loc, std::vector<tptr<uint> > *flgs)
{
	locuses = loc;
	flags = flgs;

#if HAVE_CUDA
	if(gpuGetActiveDevice() < 0) { return; }

	size_t width  = (*loc)[0].width();
	size_t height = (*loc)[0].height();

	xptrng::tptr<float4> texc;
	xptrng::tptr<uint4>  texf;
	cudaError err;
	int texid = 0;
	for(int i=0; i != loc->size(); i += 4)
	{
		// Get the pre-built arrays cached across kernel calls
		char idx[50];
		sprintf(idx, "%s%d", id, i);
		os_photometry_tex_get(idx, texc, texf);
		if(!texc || !texf)
		{
			texc = xptrng::tptr<float4>(width, height);
			texf =  xptrng::tptr<uint4>(width, height);

			// Pack the lookups to float4
			for(int y=0; y != height; y++)
			{
				for(int x=0; x != width; x++)
				{
					texc(x, y).x =                     (*loc)[i+0](x, y)    ;
					texc(x, y).y = i+1 < loc->size() ? (*loc)[i+1](x, y) : 0;
					texc(x, y).z = i+2 < loc->size() ? (*loc)[i+2](x, y) : 0;
					texc(x, y).w = i+3 < loc->size() ? (*loc)[i+3](x, y) : 0;

					texf(x, y).x =                      (*flgs)[i+0](x, y)    ;
					texf(x, y).y = i+1 < flgs->size() ? (*flgs)[i+1](x, y) : 0;
					texf(x, y).z = i+2 < flgs->size() ? (*flgs)[i+2](x, y) : 0;
					texf(x, y).w = i+3 < flgs->size() ? (*flgs)[i+3](x, y) : 0;
				}
			}
			os_photometry_tex_set(idx, texc, texf);
		}
		printf("%f %f %f %f\n%f %f %f %f\n",
			texc(317, 28).x,
			texc(317, 28).y,
			texc(317, 28).z,
			texc(317, 28).w,
			(*loc)[i+0](317, 28),
			(*loc)[i+1](317, 28),
			(*loc)[i+2](317, 28),
			(*loc)[i+3](317, 28));

		// Bind isochrone array to texture reference
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);
		//cudaArray* cu_array = gpuMMU.mapToCUDAArray(texc, channelDesc);
		cudaArray* cu_array = texc.getCUDAArray(channelDesc);
		// set texture parameters & bind the array to the texture
		colorTextures[texid]->addressMode[0] = cudaAddressModeClamp;
		colorTextures[texid]->addressMode[1] = cudaAddressModeClamp;
		colorTextures[texid]->filterMode = cudaFilterModeLinear;
		colorTextures[texid]->normalized = false;    // access with normalized texture coordinates
		err = cudaBindTextureToArray( *colorTextures[texid], cu_array, channelDesc);
		CUDA_ASSERT(err);

		// Bind flags array to texture reference
		channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindUnsigned);
		//cu_array = gpuMMU.mapToCUDAArray(texf, channelDesc);
		cu_array = texf.getCUDAArray(channelDesc);
		// set texture parameters & bind the array to the texture
		cflagsTextures[texid]->addressMode[0] = cudaAddressModeClamp;
		cflagsTextures[texid]->addressMode[1] = cudaAddressModeClamp;
		cflagsTextures[texid]->filterMode = cudaFilterModePoint;
		cflagsTextures[texid]->normalized = false;    // access with normalized texture coordinates
		err = cudaBindTextureToArray( *cflagsTextures[texid], cu_array, channelDesc);
		CUDA_ASSERT(err);
	}
#endif
}
#endif

typedef ct::cfloat::gpu_t gcfloat;
typedef ct::cint::gpu_t gcint;

#if !__CUDACC__
#include <iostream>
#endif

KERNEL(
	ks, 0,
	os_photometry_kernel(otable_ks ks, os_photometry_data lt, gcint flags, gcfloat bmag, gcfloat Mr, gcfloat mags, gcfloat FeH),
	os_photometry_kernel,
	(ks, lt, flags, bmag, Mr, mags, FeH)
)
{
	float *c = ks.sharedMemory<float>();
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		// construct colors given the absolute magnitude and metallicity
		float fFeH = FeH[row];
		float fMr = Mr[row];
		float iFeH = (fFeH - lt.FeH0) / lt.dFeH;
		float iMr  = (fMr  -  lt.Mr0) / lt.dMr;
		flags[row] = sampleColors(c, iFeH, iMr, lt.ncolors);

		// convert colors to apparent magnitudes
		float mag0 = bmag[row];
		for(int b = 0; b <= lt.ncolors; b++)
		{
			float mag = mag0;
			if(b < lt.bidx) { for(int i=b;       i != lt.bidx; i++) { mag += c[i]; } }
			if(b > lt.bidx) { for(int i=lt.bidx; i != b;    i++) { mag -= c[i]; } }
			mags(row, b) = mag;
		}

/*		if(row == 30)
		{
			for(int i=0; i != lt.ncolors; i++) { std::cerr << "c[" << i << "]=" << c[i] << "\n"; }
			for(int i=0; i != lt.ncolors+1; i++) { std::cerr << "mags[" << i << "]=" << mags(row, i) << "\n"; }
			exit(0);
		}*/
	}
}
