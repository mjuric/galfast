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

#ifndef __simulate_base_h
#define __simulate_base_h

struct os_FeH_data
{
	float A[2], sigma[3], offs[3];
	float Hmu, muInf, DeltaMu;
};

struct os_vel2pm_data
{
	int coordsys;

	float vLSR,		// Local standard of rest velocity
 	      u0, v0, w0;	// Solar peculiar motion
};

static const int GAL = 0;
static const int EQU = 1;

struct trivar_gaussK
{
	float __host__ __device__ sqr(float x) { return x*x;}
#if 1

	double __host__ __device__ Mul(float(& a)[3][3],int row, float(& z)[3]) { return a[row][0]*z[0]+a[row][1]*z[1]+a[row][2]*z[2]; }
					

	void  draw(float y[3], float s11, float s12, float s13, float s22, float s23, float s33, rng_t &rng)
	{
		float a[3][3];
		double b[3][3];		

		a[0][1]=0.;	a[0][2]=0.;
					a[1][2]=0.;

		b[0][0]=sqr(s11); 
		b[1][0]=s12;      b[1][1]=sqr(s22); 
		b[2][0]=s13;      b[2][1]=s23;      b[2][2]=sqr(s33);
		a[0][0]=s11;
			a[1][0]=1./a[0][0] * ( b[1][0] );	
		
		a[1][1]=            sqrt( b[1][1] - a[1][0]*a[1][0] );
			a[2][0]=1./a[0][0] * (b[2][0] );
			a[2][1]=1./a[1][1] * (b[2][1] - a[2][0]*a[1][0]);		

		a[2][2]=            sqrt( b[2][2] - a[2][0]*a[2][0] - a[2][1]*a[2][1]);

		
		float z[3];
		
		z[0]=rng.gaussian(1.);
		z[1]=rng.gaussian(1.);
		z[2]=rng.gaussian(1.);

		y[0]+=Mul(a,0,z);
		y[1]+=Mul(a,1,z);
		y[2]+=Mul(a,2,z);
	}

	void __device__ draw(float y[3], float s11, float s12, float s13, float s22, float s23, float s33, gpu_rng_t &rng)
	{
		float a[3][3];
		double b[3][3];		

		a[0][1]=0.;	a[0][2]=0.;
					a[1][2]=0.;

		b[0][0]=sqr(s11); 
		b[1][0]=s12;      b[1][1]=sqr(s22); 
		b[2][0]=s13;      b[2][1]=s23;      b[2][2]=sqr(s33);
		a[0][0]=s11;
			a[1][0]=1./a[0][0] * ( b[1][0] );	
		
		a[1][1]=            sqrt( b[1][1] - a[1][0]*a[1][0] );
			a[2][0]=1./a[0][0] * (b[2][0] );
			a[2][1]=1./a[1][1] * (b[2][1] - a[2][0]*a[1][0]);		

		a[2][2]=            sqrt( b[2][2] - a[2][0]*a[2][0] - a[2][1]*a[2][1]);

		
		float z[3];
		
		z[0]=rng.gaussian(1.);
		z[1]=rng.gaussian(1.);
		z[2]=rng.gaussian(1.);

		y[0]+=Mul(a,0,z);
		y[1]+=Mul(a,1,z);
		y[2]+=Mul(a,2,z);
	}
#else
	gsl_matrix *A;
	gsl_vector *Z;

	trivar_gaussK()
	{
		A = gsl_matrix_alloc(3, 3);
		Z = gsl_vector_alloc(3);

		gsl_matrix_set_zero(A);
	}

	void set(float s11, float s12, float s13, float s22, float s23, float s33)
	{
		// populate A (assumes the upper triang is already 0), calculate Cholesky decomp
		gsl_matrix_set(A, 0, 0, sqr(s11));
		gsl_matrix_set(A, 1, 0,     s12 ); gsl_matrix_set(A, 1, 1, sqr(s22));
		gsl_matrix_set(A, 2, 0,     s13 ); gsl_matrix_set(A, 2, 1,     s23 ); gsl_matrix_set(A, 2, 2, sqr(s33));
		//print_matrix(A);

		int status = gsl_linalg_cholesky_decomp(A);
		ASSERT(status == 0);
		//print_matrix(A); std::cerr << "status=" << status << "\n";
	}

	//void __device__ draw(double y[3], rng_t &rng)	
	void draw(float ret[3], rng_t &rng, bool zero = false)
	{
		gsl_vector *y=gsl_vector_alloc(3);
		gsl_vector_set(y, 0, ret[0]);
		gsl_vector_set(y, 1, ret[1]);
		gsl_vector_set(y, 2, ret[2]);

		gsl_vector_set(Z, 0, rng.gaussian(1.));
		gsl_vector_set(Z, 1, rng.gaussian(1.));
		gsl_vector_set(Z, 2, rng.gaussian(1.));

		if(zero) { gsl_vector_set_zero(y); }
		//gsl_blas_dgemv(CblasNoTrans, 1., A, Z, 1., y);
		gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, A, Z);
		gsl_vector_add(y, Z);
		//std::cout << "XXXXX: " << y->data[0] << " " << y->data[1] << " " << y->data[2] << "\n";

		ret[0]=gsl_vector_get(y, 0 );
		ret[1]=gsl_vector_get(y, 1 );
		ret[2]=gsl_vector_get(y, 2 );
		
	}

	~trivar_gaussK()
	{
		gsl_vector_free(Z);
		gsl_matrix_free(A);
	}

#endif

};

struct farray5
{
	float data[5];
	
	__device__ float& operator [] (int i) { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
	__device__ const float& operator [] (int i) const { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
};

static const int BahcallSoneira_model_THIN = 0, BahcallSoneira_model_THICK = 1, BahcallSoneira_model_HALO = 2;

struct iarray5
{
	short int data[5];
	
	__device__ short int& operator [] (int i) { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
	__device__ const short int& operator [] (int i) const { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
};

struct i8array5
{
	char data[5];
	
	__device__ char& operator [] (int i) { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
	__device__ const char& operator [] (int i) const { 		
		//if (i<0 || i>4)
		//	THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
};

struct os_kinTMIII_data
{	
		float fk, DeltavPhi;
		farray5 	vPhi1, vPhi2, vR, vZ,
			sigmaPhiPhi1, sigmaPhiPhi2, sigmaRR, sigmaZZ, sigmaRPhi, sigmaZPhi, sigmaRZ,
			HvPhi, HvR, HvZ,
			HsigmaPhiPhi, HsigmaRR, HsigmaZZ, HsigmaRPhi, HsigmaZPhi, HsigmaRZ;
};

struct os_kinTMIII_data_int
{	
		float fk, DeltavPhi;
		iarray5 	vPhi1, vPhi2, vR, vZ,
			sigmaPhiPhi1, sigmaPhiPhi2, sigmaRR, sigmaZZ, sigmaRPhi, sigmaZPhi, sigmaRZ;
		
		i8array5 HvPhi, HvR, HvZ,
			HsigmaPhiPhi, HsigmaRR, HsigmaZZ, HsigmaRPhi, HsigmaZPhi, HsigmaRZ;		
};

struct os_kinTMIII_data_groupedA
{	
				
	farray5 *diskEllip[6], *haloEllip[6], *diskMeans[3], *haloMeans[3];
		
};





#endif
