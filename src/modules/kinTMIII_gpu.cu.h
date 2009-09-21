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

#ifndef kinTMIII_gpu_cu_h__
#define kinTMIII_gpu_cu_h__

//
// Module data passed to device kernel
//
struct farray5
{
	float data[5];
	
	__device__ float& operator [] (int i) { return data[i]; }
	__device__ const float& operator [] (int i) const { return data[i]; }
};

struct iarray5
{
	short int data[5];
	
	__device__ short int& operator [] (int i) { return data[i]; } 
	__device__ const short int& operator [] (int i) const { return data[i]; } 
};

struct i8array5
{
	char data[5];
	
	__device__ char& operator [] (int i) { return data[i]; } 
	__device__ const char& operator [] (int i) const { return data[i]; }
};

struct os_kinTMIII_data
{
	int comp_thin, comp_thick, comp_halo;
	float fk, DeltavPhi;
	farray5 	vPhi1, vPhi2, vR, vZ,
		sigmaPhiPhi1, sigmaPhiPhi2, sigmaRR, sigmaZZ, sigmaRPhi, sigmaZPhi, sigmaRZ,
		HvPhi, HvR, HvZ,
		HsigmaPhiPhi, HsigmaRR, HsigmaZZ, HsigmaRPhi, HsigmaZPhi, HsigmaRZ;
};

//
// Device kernel implementation
//
#if !__CUDACC__ && !BUILD_FOR_CPU

	DECLARE_KERNEL(os_kinTMIII_kernel(otable_ks ks, gpu_rng_t rng, cint_t::gpu_t comp, cfloat_t::gpu_t XYZ, cfloat_t::gpu_t vcyl));

#else // #if !__CUDACC__ && !BUILD_FOR_CPU

	#define K_IO
	#define K_OUT

	float inline __device__ sqrf(float x) { return x*x; }	

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
			cint_t::gpu_t comp,
			cfloat_t::gpu_t XYZ,
			cfloat_t::gpu_t vcyl),
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
			const int component = comp(row);
			float X = XYZ(row, 0);
			float Y = XYZ(row, 1);
			float Zpc = XYZ(row, 2);
			const float Rsquared = 1e-6 * (X*X + Y*Y);
			const float Z = 1e-3 * Zpc;

			if(component == par.comp_thin || component == par.comp_thick)
			{
				get_disk_kinematics(tmp, Rsquared, Z, rng, par, diskMeans, diskEllip);
			}
			else if(component == par.comp_halo)
			{
				get_halo_kinematics(tmp, Rsquared, Z, rng, haloMeans, haloEllip);
			}
			vcyl(row, 0) = tmp[0];
			vcyl(row, 1) = tmp[1];
			vcyl(row, 2) = tmp[2];
		}
		rng.store(threadID());
	#undef par
	}


#endif // #else (!__CUDACC__ && !BUILD_FOR_CPU)

#endif // kinTMIII_gpu_cu_h__
