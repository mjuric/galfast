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
#include "gpu.h"
#include "skygen.h"

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
	uint32_t tid = threadID();
	rng.load(tid);
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
	rng.store(threadID());
}

//======================================================================
//======================================================================
//    vel2pm
//======================================================================
//======================================================================

// template<typename T>
// __device__ inline __device__ T sqr(const T x) { return x*x; }

using namespace peyton;

__device__ inline double rad(double deg) { return 0.017453292519943295769  * deg; } // convert deg to rad
__device__ inline float  radf(float deg) { return 0.017453292519943295769f * deg; } // convert deg to rad

/* convert cartesian to celestial coordinates */
__device__ void xyz2lbr(Radians &l, Radians &b, float &r, float x, float y, float z)
{
	r = sqrt(x*x + y*y + z*z);
	b = asin(z / r);
	l = atan2(y, x);
	if(l < 0) l += ctn::twopi;
}

/* convert celestial to cartesian coordinates */
__device__ void lbr2xyz(float &x, float &y, float &z, Radians l, Radians b, float r)
{
	x = r*cos(l)*cos(b);
	y = r*sin(l)*cos(b);
	z = r*sin(b);
}

// convert cartesian _Galactic_ (not galactocentric!!) velocities to vl,vr,vb wrt. the observer (velocity units are unimportant)
__device__ void vel_xyz2lbr(float &vl, float &vb, float &vr, const float vx, const float vy, const float vz, const float l, const float b)
{
	float cl, sl, cb, sb;
	cl = cosf(l);
	sl = sinf(l);
	cb = cosf(b);
	sb = sinf(b);

	float tmp;
	vl  = -vx*sl  + vy*cl;
	tmp =  vx*cl  + vy*sl;
	vr  =  cb*tmp + vz*sb;
	vb  = -sb*tmp + vz*cb;
}

// convert vl,vr,vb velocities wrt. the observer to cartesian _Galactic_ (not galactocentric!!) velocities (velocity units are unimportant)
__device__ void vel_lbr2xyz(float &vx, float &vy, float &vz, const float vl, const float vb, const float vr, const float l, const float b)
{
	float cl, sl, cb, sb;
	cl = cosf(l);
	sl = sinf(l);
	cb = cosf(b);
	sb = sinf(b);

	float tmp = sb*vb - cb*vr;
	vx = -sl*vl - cl*tmp;
	vy =  cl*vl - sl*tmp;
	vz =  cb*vb + sb*vr;
}

// convert galactocentric cylindrical to galactocentric cartesian velocities
__device__ inline void vel_cyl2xyz(float &vx, float &vy, float &vz, const float vr, const float vphi, const float vz0, const float X, const float Y)
{
	float rho = sqrtf(X*X + Y*Y);
	float cphi = X / rho;
	float sphi = Y / rho;

	vx = -sphi * vphi + cphi * vr;
	vy =  cphi * vphi + sphi * vr;
	vz = vz0;
}

namespace galequ_constants
{
	static const double angp = ctn::d2r * 192.859508333; //  12h 51m 26.282s (J2000)
	static const double dngp = ctn::d2r * 27.128336111;  // +27d 07' 42.01" (J2000)
	static const double l0 = ctn::d2r * 32.932;	// galactic longitude of ascending node of galactic coordinate system (where b=0, dec=0)
	static const double ce = 0.88998740217659689; // cos(dngp)
	static const double se = 0.45598511375586859; // sin(dngp)
};

namespace float_galequ_constants
{
	static const float angp = ctn::d2r * 192.859508333; //  12h 51m 26.282s (J2000)
	static const float dngp = ctn::d2r * 27.128336111;  // +27d 07' 42.01" (J2000)
	static const float ce = (float)galequ_constants::ce;
	static const float se = (float)galequ_constants::se;
	static const float l0 = (float)galequ_constants::l0;
};

float inline __device__ sqrf(float x) { return x*x; }	

//static const float flt_d2r = (float)ctn::d2r;
static const float flt_r2d = (float)(1./ctn::d2r);
inline __device__ float deg(float radians)
{
	return flt_r2d * radians;
}

/*
	Convert proper motions in one coordinate system to de. The coordinate
	system orientations are given by specifying the pole of the source coordinate
	system in destination coordinates, and the longitude of ascending node
	of the source coordinate system, in source coordinates.

	transform_pmf: uses single-precision floating point arithmetic.

	Input:		l,b in radians; vl,vb (units not important)
	Output:		vlout, vbout in destination coordinate system
*/
__device__ void transform_pm(float &vlout, float &vbout, float l, float b, float vl, float vb,
	const double ce, const double se, const double l0)
{
	double cb = cos(b);
	double sb = sin(b);
	double cl = cos(l - l0);
	double sl = sin(l - l0);

	double tmp1 = cb * se - ce * sb * sl; // ??
	double tmp2 = ce * sb - cb * se * sl; // \propto cos(alpha)
	double tmp4 = cb * cl;                // \propto sin(alpha)
	double tmp3 = sb * se + cb * ce * sl; // sin(delta)

	if(fabs(tmp3) > 0.9999999)
	{
		// we're practically at the north/south pole.
		// better signal that computation here will be incorrect,
		// than give an incorrect result.
		vlout = vbout = 999.99;
	}
	else
	{
		double denom1 = 1. / ( tmp4*tmp4 + tmp2*tmp2 );
		double cd = sqrt(1. - tmp3*tmp3); // cos(delta);

		vl /= cb;
		vlout = (-ce * cl * vb + cb * tmp1 * vl) * denom1;
		vbout = (tmp1 * vb  + cb * ce * cl * vl) / cd;
		vlout *= cd;
	}
}

__device__ void transform_pmf(float &vlout, float &vbout, float l, float b, float vl, float vb,
	const float ce, const float se, const float l0)
{
	float cb = cosf(b);
	float sb = sinf(b);
	float cl = cosf(l - l0);
	float sl = sinf(l - l0);

	float tmp1 = cb * se - ce * sb * sl; // ??
	float tmp2 = ce * sb - cb * se * sl; // \propto cos(alpha)
	float tmp3 = sb * se + cb * ce * sl; // sin(delta)
	float cd = sqrtf(1 - sqrf(tmp3)); // cos(delta);

	vl /= cb;
	vlout = (-ce * cl * vb + cb * tmp1 * vl) / ( sqrf(cb * cl) + sqrf(tmp2) );
	vbout = (tmp1 * vb  + cb * ce * cl * vl) / cd;
	vlout *= cd;
}

/*
	Convert proper motions in galactic coordinate system to proper motions in equatorial
	coordinate system. Uses single-precision floating point arithmetic.

	Input:		l,b in radians; vl,vb (units not important)
	Output:		vra,vdec (in units of vl,vb)
*/
__device__ void pm_galequf(float &vra, float &vdec, float l, float b, float vl, float vb)
{
	using namespace float_galequ_constants;
	transform_pm(vra, vdec, l, b, vl, vb, ce, se, l0);
}

__device__ void pm_galequ(float &vra, float &vdec, float l, float b, float vl, float vb)
{
	using namespace galequ_constants;
	transform_pm(vra, vdec, l, b, vl, vb, ce, se, l0);
}

__device__ void pm_equgalf(float &vl, float &vb, float ra, float dec, float vra, float vdec)
{
	using namespace float_galequ_constants;
	transform_pm(vl, vb, ra, dec, vra, vdec, ce, se, angp - 0.5*ctn::pi);
}

__device__ void pm_equgal(float &vl, float &vb, float ra, float dec, float vra, float vdec)
{
	using namespace galequ_constants;
	transform_pm(vl, vb, ra, dec, vra, vdec, ce, se, angp - 0.5*ctn::pi);
}

__device__ float3 vcyl2pm(float l, float b, float vx, float vy, float vz, float X, float Y, float Z, os_vel2pm_data &par)
{
#if 0
	// These must produce (mu_l, mu_b, v_radial) = (12.72 mas/yr, -14.75 mas/yr, -27.34 km/s)
	X = 8100; Y = 100; Z = 3000;
	vx = 10; vy = 50; vz = -30;
	l = atan2(-Y, 8000-X);
	b = asin(Z / sqrt(sqrf(8000-X) + sqrf(Y) + sqrf(Z)));
	#if __DEVICE_EMULATION__
	printf("%f %f\n", deg(l), deg(b));
	#endif
#endif
	// convert the velocities from cylindrical to galactocentric cartesian system
	float3 pm;	// pm.x == mu_l, pm.y == mu_b, pm.z == radial velocity (km/s)
	vel_cyl2xyz(pm.x, pm.y, pm.z,   vx, vy, vz,   X, Y);

	// switch to Solar coordinate frame
	pm.x -= par.u0;
	pm.y -= par.v0 + par.vLSR;
	pm.z -= par.w0;

	// convert to velocities wrt. the observer
	vel_xyz2lbr(pm.x, pm.y, pm.z,   -pm.x, -pm.y, pm.z,  l, b);

	// convert to proper motions
	float D = sqrtf(sqrf(8000.f-X) + sqrf(Y) + sqrf(Z));
	pm.x /= 4.74 * D*1e-3;	// proper motion in mas/yr (4.74 km/s @ 1kpc is 1mas/yr)
	pm.y /= 4.74 * D*1e-3;

	return pm;
}
KERNEL(
	ks, 0, 
	os_vel2pm_kernel(
		otable_ks ks, os_vel2pm_data par, gpu_rng_t rng, 
		ct::cdouble::gpu_t lb0, 
		ct::cfloat::gpu_t XYZ,
		ct::cfloat::gpu_t vcyl,  
		ct::cfloat::gpu_t pmout),
	os_vel2pm_kernel,
	(ks, par, rng, lb0, XYZ, vcyl, pmout)
)
{

	for(uint32_t row=ks.row_begin(); row < ks.row_end(); row++)
	{
		// fetch prerequisites
		// NOTE: Single precision should be sufficient here for (l,b)
		float l = radf(lb0(row, 0));
		float b = radf(lb0(row, 1));
		float X = XYZ(row, 0);
		float Y = XYZ(row, 1);
		float Z = XYZ(row, 2);
		float vx = vcyl(row, 0);
		float vy = vcyl(row, 1);
		float vz = vcyl(row, 2);

		float3 pm = vcyl2pm(l, b, vx, vy, vz, X, Y, Z, par);

		// rotate to output coordinate system
		switch(par.coordsys)
		{
		case GAL:
			break;
		case EQU:
/*#if __DEVICE_EMULATION__
			fprintf(stderr, "before = %f %f %f %f\n", deg(l), deg(b), pm.x, pm.y);
#endif*/
			pm_galequ(pm.x, pm.y, l, b, pm.x, pm.y);
/*#if __DEVICE_EMULATION__
			fprintf(stderr, "after = %f %f\n", pm.x, pm.y);
			abort();
#endif*/
			break;
		default:
			//THROW(EAny, "Unknown coordinate system [id=" + str(coordsys) + "] requested");
			pm.x = pm.y = pm.z = -9.99;
			break;
		}

		pmout(row, 0) = pm.x;
		pmout(row, 1) = pm.y;
		pmout(row, 2) = pm.z;
	}
}


//======================================================================
//======================================================================
//    kinTMIII
//======================================================================
//======================================================================

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
	rng.load(threadID());

#define par os_kinTMIII_par
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
	for(int row=ks.row_begin(); row < ks.row_end(); row++)
	{
		// fetch prerequisites
		const int component = comp[row];
		float X = XYZ(row, 0);
		float Y = XYZ(row, 1);
		float Zpc = XYZ(row, 2);
		const float Rsquared = 1e-6 * (X*X + Y*Y);
		const float Z = 1e-3 * Zpc;

// #ifdef __DEVICE_EMULATION__
// 		printf("%d %d %d R=%f Z=%f\n", tid, (int)row, component, sqrtf(Rsquared), Z);
// #endif

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
/*#ifdef __DEVICE_EMULATION__
	printf("leaving.\n");
#endif*/
	rng.store(threadID());
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

inline __device__ double2 galequ(const double2 lb)
{
	using namespace galequ_constants;

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

	while(r.x < 0.)        { r.x += ctn::pi2; }
	while(r.x >= ctn::pi2) { r.x -= ctn::pi2; }

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

DEFINE_TEXTURE(secProb);
DEFINE_TEXTURE(cumLF);
DEFINE_TEXTURE(invCumLF);

__device__ bool draw_companion(float &M2, float M1, multiplesAlgorithms::algo algo, gpu_rng_t &rng)
{
#if 0 //__DEVICE_EMULATION__
	printf("%f %f\n", secProb.x0, secProb.inv_dx);
	{
		float tprob[] = {0., 2., 16.5, 18., 18.1};
		for(int i=0; i != sizeof(tprob)/sizeof(float); i++)
			printf("secProb(%f) = %f\n", tprob[i], (float)secProb.sample(tprob[i]));
	}
	{
		float tprob[] = {2.99, 3., 10., 17., 17.01};
		for(int i=0; i != sizeof(tprob)/sizeof(float); i++)
			printf("cumLF(%f) = %f\n", tprob[i], (float)cumLF.sample(tprob[i]));
	}
	{
		float tprob[] = {-0.01, 0., 0.364679, 1., 1.01};
		for(int i=0; i != sizeof(tprob)/sizeof(float); i++)
			printf("invCumLF(%f) = %f\n", tprob[i], (float)invCumLF.sample(tprob[i]));
	}
#endif
	// draw the probability that this star has a secondary
	float psec, u;

	psec = secProb.sample(M1);
	u = rng.uniform();
#if __DEVICE_EMULATION__
/*	for(float u=0; u <= 17; u += .1)
	{
		printf("secProb: %.3f %f\n", u, (float)secProb.sample(u));
	}
	abort();*/
#endif

#if __DEVICE_EMULATION__
//	assert(psec == 1.f);
#endif
	if(u > psec) { return false; }

	// draw the absolute magnitude of the secondary, subject to requested
	// algorithm
	using namespace multiplesAlgorithms;
	if(algo == EQUAL_MASS) { M2 = M1; return true; }

	float pprim = cumLF.sample(M1);
	u = rng.uniform();
	if(algo == LF_M2_GT_M1)
	{
		// draw subject to the requirement that it is fainter than the primary
		u = pprim + u * (1. - pprim);
	}
	M2 = invCumLF.sample(u);// + rng.gaussian(1.f);
//	M2 = 5.f + 12.f*u;
	if(algo == LF_M2_GT_M1 && M2 < M1)
	{
		// This can happen due to resampling of cumLF and invCumLF
		// (see the note in os_unresolvedMultiples::construct)
		M2 = M1;
	}

#if __DEVICE_EMULATION__
/*	for(float u=0; u <=1; u += 0.01)
 	{
 		printf("%.3f %f\n", u, (float)invCumLF.sample(u));
 	}
 	abort();*/
#endif

	return true;
}

KERNEL(
	ks, 3*4,
	os_unresolvedMultiples_kernel(otable_ks ks, gpu_rng_t rng, int nabsmag, ct::cfloat::gpu_t M, ct::cfloat::gpu_t Msys, ct::cint::gpu_t ncomp, multiplesAlgorithms::algo algo),
	os_unresolvedMultiples_kernel,
	(ks, rng, nabsmag, M, Msys, ncomp, algo)
)
{
	/*
		Input:	M -- absolute magnitude of the primary.
		Output:	Msys[] -- absolute magnitudes of system components.
			M      -- total absolute magnitude of the system
	*/
	rng.load(threadID());
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		float M1 = M[row];
		Msys(row, 0) = M1;
		float Ltot = exp10f(-0.4f*M1);
		int ncomps = 1;			/* Number of components of the multiple system */
		for(int i = 1; i < nabsmag; i++)
		{
			float M2;
			if(draw_companion(M2, M1, algo, rng))
			{
				Msys(row, i) = M2;
				Ltot += exp10f(-0.4f*M2);
				ncomps++;
			}
			else
			{
				Msys(row, i) = ABSMAG_NOT_PRESENT;
			}
		}
		float Mtot = -2.5*log10f(Ltot);
		    M[row] = Mtot;
		ncomp[row] = ncomps;
	}
	rng.store(threadID());
}

#include <vector>

#if BUILD_FOR_CPU && HAVE_CUDA
extern __TLS std::vector<xptrng::tptr<float> > *locuses;
extern __TLS std::vector<xptrng::tptr<uint> >   *flags;
#else
__TLS std::vector<xptrng::tptr<float> > *locuses;
__TLS std::vector<xptrng::tptr<uint> >   *flags;
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
	uint flags = 0;

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
		colors[ic] = (*locuses)[ic].elem(f, m);
		fl |= (*flags)[ic].elem(f, m);
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

void os_photometry_set_isochrones(const char *id, std::vector<xptrng::tptr<float> > *loc, std::vector<xptrng::tptr<uint> > *flgs)
{
	locuses = loc;
	flags = flgs;

#if HAVE_CUDA
	activeDevice dev(gpuExecutionEnabled("os_photometry_kernel")? 0 : -1);
	if(gpuGetActiveDevice() < 0) { return; }

	size_t width  = (*loc)[0].width();
	size_t height = (*loc)[0].height();

	xptrng::tptr<float4> texc;
	xptrng::tptr<uint4>  texf;
	int texid = 0;
	for(int i=0; i < loc->size(); i += 4)
	{
		// Get the pre-built arrays cached across kernel calls
		char idx[50];
		sprintf(idx, "%s%d", id, i);
		os_photometry_tex_get(idx, texc, texf);
		if(texc.isNull() || texf.isNull())
		{
			texc = xptrng::tptr<float4>(width, height);
			texf =  xptrng::tptr<uint4>(width, height);

			// Pack the lookups to float4
			for(int y=0; y != height; y++)
			{
				for(int x=0; x != width; x++)
				{
					texc.elem(x, y).x =                     (*loc)[i+0].elem(x, y)    ;
					texc.elem(x, y).y = i+1 < loc->size() ? (*loc)[i+1].elem(x, y) : 0;
					texc.elem(x, y).z = i+2 < loc->size() ? (*loc)[i+2].elem(x, y) : 0;
					texc.elem(x, y).w = i+3 < loc->size() ? (*loc)[i+3].elem(x, y) : 0;

					texf.elem(x, y).x =                      (*flgs)[i+0].elem(x, y)    ;
					texf.elem(x, y).y = i+1 < flgs->size() ? (*flgs)[i+1].elem(x, y) : 0;
					texf.elem(x, y).z = i+2 < flgs->size() ? (*flgs)[i+2].elem(x, y) : 0;
					texf.elem(x, y).w = i+3 < flgs->size() ? (*flgs)[i+3].elem(x, y) : 0;
				}
			}
			os_photometry_tex_set(idx, texc, texf);
		}
// 		printf("%f %f %f %f\n%f %f %f %f\n",
// 			texc(317, 28).x,
// 			texc(317, 28).y,
// 			texc(317, 28).z,
// 			texc(317, 28).w,
// 			(*loc)[i+0](317, 28),
// 			(*loc)[i+1](317, 28),
// 			(*loc)[i+2](317, 28),
// 			(*loc)[i+3](317, 28));

		// Bind isochrone array to texture reference
		cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindFloat);
		//cudaArray* cu_array = gpuMMU.mapToCUDAArray(texc, channelDesc);
		cudaArray* cu_array = texc.getCUDAArray(channelDesc);
		// set texture parameters & bind the array to the texture
		colorTextures[texid]->addressMode[0] = cudaAddressModeClamp;
		colorTextures[texid]->addressMode[1] = cudaAddressModeClamp;
		colorTextures[texid]->filterMode = cudaFilterModeLinear;
		colorTextures[texid]->normalized = false;    // access with normalized texture coordinates
		cuxErrCheck( cudaBindTextureToArray( *colorTextures[texid], cu_array, channelDesc) );

		// Bind flags array to texture reference
		channelDesc = cudaCreateChannelDesc(32, 32, 32, 32, cudaChannelFormatKindUnsigned);
		//cu_array = gpuMMU.mapToCUDAArray(texf, channelDesc);
		cu_array = texf.getCUDAArray(channelDesc);
		// set texture parameters & bind the array to the texture
		cflagsTextures[texid]->addressMode[0] = cudaAddressModeClamp;
		cflagsTextures[texid]->addressMode[1] = cudaAddressModeClamp;
		cflagsTextures[texid]->filterMode = cudaFilterModePoint;
		cflagsTextures[texid]->normalized = false;    // access with normalized texture coordinates
		cuxErrCheck( cudaBindTextureToArray( *cflagsTextures[texid], cu_array, channelDesc) );

		texid++;
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
	os_photometry_kernel(otable_ks ks, os_photometry_data lt, gcint flags, gcfloat DM, gcfloat Mr, int nabsmag, gcfloat mags, gcfloat FeH),
	os_photometry_kernel,
	(ks, lt, flags, DM, Mr, nabsmag, mags, FeH)
)
{
	float *c = ks.sharedMemory<float>();
	for(uint32_t row = ks.row_begin(); row < ks.row_end(); row++)
	{
		// construct colors given the absolute magnitude and metallicity
		float fFeH = FeH[row];
		float iFeH = (fFeH - lt.FeH0) / lt.dFeH;

		// generate system components, compute system luminosity
		for(int syscomp = 0; syscomp < nabsmag; syscomp++)
		{
			float fMr = Mr(row, syscomp);
			if(fMr >= ABSMAG_NOT_PRESENT) { break; }

			float iMr  = (fMr  -  lt.Mr0) / lt.dMr;
			int flag = sampleColors(c, iFeH, iMr, lt.ncolors);
			if(syscomp) { flag &= flags[row]; }
			flags[row] = flag;

			// compute absolute magnitudes in different bands, and store
			// as luminosity
			for(int b = 0; b <= lt.ncolors; b++)
			{
				float M = fMr;
				if(b < lt.bidx) { for(int i=b;       i != lt.bidx; i++) { M += c[i]; } }
				if(b > lt.bidx) { for(int i=lt.bidx; i != b;    i++)    { M -= c[i]; } }

				float L = exp10f(-0.4f*M);
				if(syscomp) { L += mags(row, b); }
				mags(row, b) = L;
			}
		}

		// convert luminosity to apparent magnitude of the system
		float dm = DM[row];
		for(int b = 0; b <= lt.ncolors; b++)
		{
			float Mtot = -2.5f * log10f(mags(row, b));
			float mtot = dm + Mtot;
			mags(row, b) = mtot;
		}

/*		if(row == 30)
		{
			for(int i=0; i != lt.ncolors; i++) { std::cerr << "c[" << i << "]=" << c[i] << "\n"; }
			for(int i=0; i != lt.ncolors+1; i++) { std::cerr << "mags[" << i << "]=" << mags(row, i) << "\n"; }
			exit(0);
		}*/
	}
}
