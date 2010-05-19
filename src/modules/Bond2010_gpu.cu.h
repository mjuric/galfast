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

#ifndef Bond2010_gpu_cu_h__
#define Bond2010_gpu_cu_h__

namespace Bond2010
{
	//
	// Module data passed to device kernel
	//
	struct farray5
	{
		float data[5];
		
		__device__ float& operator [] (int i) { return data[i]; }
		__device__ const float& operator [] (int i) const { return data[i]; }
	};
	
	struct vmoments
	{
		// first (means) and second (dispersion tensor) velocity moments
		farray5  m[3];	// means: vR,vPhi,vZ or vr,vTheta,vPhi
		farray5 ss[6];	// dispersion tensor components: s11, s12, s13, s22, s23, s33
	};

	struct os_Bond2010_data
	{
		bit_map comp_thin, comp_thick, comp_halo;

		// two-component disk
		float fk;
		vmoments disk_1, disk_2;

		// halo
		vmoments halo;
	
		// coordinate system definition
		float M[3][3];
		float3 T;
	};
}

//
// Device kernel implementation
//
#if !__CUDACC__ && !BUILD_FOR_CPU

	DECLARE_KERNEL(os_Bond2010_kernel(otable_ks ks, gpu_rng_t rng, cint_t::gpu_t comp, cint_t::gpu_t hidden, cfloat_t::gpu_t XYZ, cfloat_t::gpu_t vcyl));

#else // #if !__CUDACC__ && !BUILD_FOR_CPU

	#define K_IO
	#define K_OUT

namespace Bond2010
{
	/*
		Draws a 3D vector from a 3D gaussian defined by the symmetric tensor s11..s33.
	
		NOTE: WARNING: CAVEAT: Tensor components s11, s22, and s33 are SQUARE ROOTS
				OF VARIANCES (sigma) and not the variances (sigma^2)!!
	*/
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


	__device__ inline float modfun(float Rsquared, float Z, const farray5 &val)
	{
		// model function: result = a + b*|Z|^c + d*R^e
		return val[0] + val[1]*powf(fabs(Z), val[2]) + val[3]*powf(Rsquared, 0.5*val[4]);
	}

	__device__ void add_dispersion(K_IO float v[3], float Rsquared, float Z, const farray5 ellip[6], K_IO gpu_rng_t &rng)
	{
		// compute velocity dispersions at this position, and draw from trivariate gaussian
		// NOTE: ADDS THE RESULT TO v, DOES NOT ZERO v BEFORE !!
		float sigma[6];
		for (int i=0;i<6;i++)		
			sigma[i] = modfun(Rsquared, Z, ellip[i]);	
				
		trivar_gaussK_draw(v, sigma[0], sigma[1], sigma[2], sigma[3], sigma[4], sigma[5], rng);
	}

	__device__ void compute_means(K_OUT float v[3], float Rsquared, float Z, const farray5 means[3])
	{
		// returns means in v[3]
		for (int i=0;i<3;i++) 	
			v[i] = modfun(Rsquared, Z, means[i]);	
	}

	__device__ void get_disk_kinematics(K_OUT float v[3], float Rsquared, float Z, gpu_rng_t &rng, const os_Bond2010_data& par)
	{
		// select the gaussian
		float p = rng.uniform();
		const vmoments &disk = p < par.fk ? par.disk_1 : par.disk_2;

		// draw
		compute_means(v, Rsquared, Z, disk.m);
		if(v[1] > 0.) { v[1] = 0.; } // truncate v_phi > 0
		add_dispersion(v, Rsquared, Z, disk.ss, rng);
	}

	__device__ void get_halo_kinematics(K_OUT float v[3], float Rcyl_squared, float Z, K_IO gpu_rng_t &rng, const os_Bond2010_data& par)
	{
		// These are all in _SPHERICAL_ coordinate system, while the output velocity is expected in _CYLINDRICAL_ coordinates
		// v[0,1,2] = v_r, v_theta, v_phi
		compute_means(v, Rcyl_squared, Z, par.halo.m);
		add_dispersion(v, Rcyl_squared, Z, par.halo.ss, rng);

		// Rotate to cylindrical system (v_R, v_phi, v_Z)
		float invR = 1.f/sqrtf(Rcyl_squared + Z*Z);
		float cosa = sqrtf(Rcyl_squared) * invR;
		float sina =                   Z * invR;
		float vR = cosa * v[0] - sina * v[1];
		float vZ = sina * v[0] + cosa * v[1];

		v[0] = vR;
		v[1] = v[2]; // v_phi is unchanged
		v[2] = vZ;
	}
}

	__device__ __constant__ Bond2010::os_Bond2010_data os_Bond2010_par;

	KERNEL(
		ks, 3*4,
		os_Bond2010_kernel(
			otable_ks ks, gpu_rng_t rng, 
			cint_t::gpu_t comp,
			cint_t::gpu_t hidden,
			cfloat_t::gpu_t XYZ,
			cfloat_t::gpu_t vcyl),
		os_Bond2010_kernel,
		(ks, rng, comp, hidden, XYZ, vcyl)
	)
	{
		using namespace Bond2010;
		const os_Bond2010_data &par = os_Bond2010_par;

		rng.load(threadID());
		for(int row=ks.row_begin(); row < ks.row_end(); row++)
		{
			if(hidden(row)) { continue; }

			// fetch prerequisites
			const int cmp = comp(row);
			
#if 0
			float X = par.Rg - XYZ(row, 0);
			float Y = -XYZ(row, 1);
			float Zpc = XYZ(row, 2);
			const float Rsquared = 1e-6 * (X*X + Y*Y);
			const float Z = 1e-3 * Zpc;
#else
			float3 v = { XYZ(row, 0), XYZ(row, 1), XYZ(row, 2) };
			v = transform(v, par.T, par.M);
			const float Rsquared = 1e-6 * (v.x*v.x + v.y*v.y);
			const float Z = 1e-3 * v.z;
#endif

			if(par.comp_thin.isset(cmp) || par.comp_thick.isset(cmp))
			{
				float tmp[3]; 
				get_disk_kinematics(tmp, Rsquared, Z, rng, par);
				vcyl(row, 0) = tmp[0];
				vcyl(row, 1) = tmp[1];
				vcyl(row, 2) = tmp[2];
			}
			else if(par.comp_halo.isset(cmp))
			{
				float tmp[3]; 
				get_halo_kinematics(tmp, Rsquared, Z, rng, par);
				vcyl(row, 0) = tmp[0];
				vcyl(row, 1) = tmp[1];
				vcyl(row, 2) = tmp[2];
			}
		}
		rng.store(threadID());
	}


#endif // #else (!__CUDACC__ && !BUILD_FOR_CPU)

#endif // Bond2010_gpu_cu_h__
