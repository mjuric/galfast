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

// add Fe/H information
class os_FeH : public osink, os_FeH_data
{
public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool init(const Config &cfg, otable &t);
	virtual const std::string &name() const { static std::string s("FeH"); return s; }

	os_FeH() : osink()
	{
		prov.insert("FeH");
		req.insert("comp");
		req.insert("XYZ");
	}
};

// add Fe/H information
class os_fixedFeH : public osink
{
	protected:
		float fixedFeH;

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool init(const Config &cfg, otable &t);
		virtual const std::string &name() const { static std::string s("fixedFeH"); return s; }

		os_fixedFeH() : osink(), fixedFeH(0)
		{
			prov.insert("FeH");
		}
};

// convert velocities to proper motions
class os_vel2pm : public osink , public os_vel2pm_data
{
public:
	virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
	virtual bool init(const Config &cfg, otable &t);
	virtual const std::string &name() const { static std::string s("vel2pm"); return s; }

	os_vel2pm() : osink()
	{
		coordsys=GAL;
		req.insert("lb");
		req.insert("XYZ");
		req.insert("vcyl");
	}
};


struct trivar_gauss
{
	gsl_matrix *A;
	gsl_vector *Z;

	trivar_gauss()
	{
		A = gsl_matrix_alloc(3, 3);
		Z = gsl_vector_alloc(3);

		gsl_matrix_set_zero(A);
	}

	void set(double s11, double s12, double s13, double s22, double s23, double s33)
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

	void draw(gsl_vector *y, rng_t &rng, bool zero = false)
	{
		gsl_vector_set(Z, 0, rng.gaussian(1.));
		gsl_vector_set(Z, 1, rng.gaussian(1.));
		gsl_vector_set(Z, 2, rng.gaussian(1.));

		if(zero) { gsl_vector_set_zero(y); }
		//gsl_blas_dgemv(CblasNoTrans, 1., A, Z, 1., y);
		gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, A, Z);
		gsl_vector_add(y, Z);
		//std::cout << "XXXXX: " << y->data[0] << " " << y->data[1] << " " << y->data[2] << "\n";
	}

	~trivar_gauss()
	{
		gsl_vector_free(Z);
		gsl_matrix_free(A);
	}
};


struct trivar_gaussK
{
	double a[3][3];
	double z[3];

	void set(double s11, double s12, double s13, double s22, double s23, double s33)
	{
		double b[3][3];		

		a[0][1]=0.;	a[0][2]=0.;
					a[1][2]=0.;

		b[0][0]=sqr(s11); 
		b[1][0]=s12;      b[1][1]=sqr(s22); 
		b[2][0]=s13;      b[2][1]=s23;      b[2][2]=sqr(s33);
			
		// populate A (assumes the upper triang is already 0), calculate Cholesky decomp
		/*gsl_matrix_set(A, 0, 0, sqr(s11));
		gsl_matrix_set(A, 1, 0,     s12 ); gsl_matrix_set(A, 1, 1, sqr(s22));
		gsl_matrix_set(A, 2, 0,     s13 ); gsl_matrix_set(A, 2, 1,     s23 ); gsl_matrix_set(A, 2, 2, sqr(s33));*/
		//print_matrix(A);

		a[0][0]=s11;
			a[1][0]=1./a[0][0] * ( b[1][0] );	
		
		a[1][1]=sqrt( b[1][1] - a[1][0]*a[1][0] );
			a[2][0]=1./a[0][0] * (b[2][0] );
			a[2][1]=1./a[0][0] * (b[2][0] - b[2][0]*b[1][0]);		

		a[2][2]=sqrt( b[2][2] - a[2][0]*a[2][0] - a[2][1]*a[2][1]);
		
		
		//int status = gsl_linalg_cholesky_decomp(A);
		//ASSERT(status == 0);
		//print_matrix(A); std::cerr << "status=" << status << "\n";
	}

	double Mul(int row, double z[3]) { return a[row][0]*z[0]+a[row][1]*z[1]+a[row][2]*z[2]; }
		

	void draw(double y[3], rng_t &rng)
	{
		double z[3];
		
		z[0]=rng.gaussian(1.);
		z[1]=rng.gaussian(1.);
		z[2]=rng.gaussian(1.);

		y[0]+=Mul(0,z);
		y[1]+=Mul(1,z);
		y[2]+=Mul(2,z);

		//gsl_vector_set(Z, 0, rng.gaussian(1.));
		//gsl_vector_set(Z, 1, rng.gaussian(1.));
		//gsl_vector_set(Z, 2, rng.gaussian(1.));

				

		//gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, A, Z);
		//gsl_vector_add(y, Z);
		//std::cout << "XXXXX: " << y->data[0] << " " << y->data[1] << " " << y->data[2] << "\n";
	}

};


struct darray5
{
	double data[5];
	
	double& operator [] (int i) { 		
		if (i<0 || i>4)
			THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
	const double& operator [] (int i) const { 		
		if (i<0 || i>4)
			THROW(ENotImplemented, "We should have never gotten here");
		return data[i];
		} 
};

struct os_kinTMIII_data
{	
	trivar_gaussK tri_rnd;

		double fk, DeltavPhi;
		darray5 	vPhi1, vPhi2, vR, vZ,
			sigmaPhiPhi1, sigmaPhiPhi2, sigmaRR, sigmaZZ, sigmaRPhi, sigmaZPhi, sigmaRZ,
			HvPhi, HvR, HvZ,
			HsigmaPhiPhi, HsigmaRR, HsigmaZZ, HsigmaRPhi, HsigmaZPhi, HsigmaRZ;
		
		darray5 *diskEllip[6], *haloEllip[6], *diskMeans[3], *haloMeans[3];
		
};

	

// add kinematic information
class os_kinTMIII : public osink, os_kinTMIII_data
{	
	public:
		void add_dispersion(double v[3], double Rsquared, double Z, darray5 *ellip[6], rng_t &rng);
		void compute_means(double v[3], double Rsquared, double Z, darray5 *means[3]);

		void get_disk_kinematics(double v[3], double Rsquared, double Z, rng_t &rng, bool &firstGaussian);
		void get_halo_kinematics(double v[3], double Rsquared, double Z, rng_t &rng);

	public:
		virtual size_t process(otable &in, size_t begin, size_t end, rng_t &rng);
		virtual bool init(const Config &cfg, otable &t);
		virtual const std::string &name() const { static std::string s("kinTMIII"); return s; }

		os_kinTMIII() : osink()
		{
			prov.insert("vcyl");
			req.insert("comp");
			req.insert("XYZ");
		}
};

inline double modfun(double Rsquared, double Z, double a, double b, double c, double d, double e)
{
	return a + b*pow(fabs(Z), c) + d*pow(Rsquared, 0.5*e);
}

