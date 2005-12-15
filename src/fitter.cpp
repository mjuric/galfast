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

#include "analysis.h"
#include "model.h"

#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include <astro/system/options.h>
#include <astro/exceptions.h>
#include <astro/util.h>
#include <astro/math.h>
#include <astro/system/log.h>
#include <astro/io/format.h>
#include <astro/useall.h>

using namespace std;

extern "C"
int model_fdf (const gsl_vector * v, void *params, gsl_vector * f, gsl_matrix * J)
{
	model_fitter *m = (model_fitter *)params;

	m->set_parameters(v);
	m->fdf(f, J);
}

extern "C"
int model_f (const gsl_vector * v, void *params, gsl_vector * f)
{
	return model_fdf(v, params, f, NULL);
}

extern "C"
int model_df (const gsl_vector * v, void *params, gsl_matrix * J)
{
	return model_fdf(v, params, NULL, J);
}

int print_state (size_t iter, gsl_multifit_fdfsolver * s, int dof)
{
	fprintf (stderr, "iter: %3u x = ", iter);
	FOR(0, dof)
	{
	  	fprintf(stderr, "%15.8f ", gsl_vector_get(s->x, i));
	}
	fprintf(stderr, "|f(x)| = %g\n", gsl_blas_dnrm2 (s->f));
}

int model_fitter::fit(int cullIter, double nsigma)
{
	//
	// Generic fit driver
	//
	int n = ndata();	// number of data points

	// allocate arrays, copy non-constant parameters into the array
	gsl_vector *v = gsl_vector_alloc(ndof);
	get_parameters(v);

	// initialize and set up the solver
	gsl_multifit_fdfsolver *s = gsl_multifit_fdfsolver_alloc (gsl_multifit_fdfsolver_lmsder, n, ndof);

	gsl_multifit_function_fdf f;
	f.f = &model_f;
	f.df = &model_df;
	f.fdf = &model_fdf;
	f.n = n;
	f.p = ndof;
	f.params = (void *)this;
	gsl_multifit_fdfsolver_set (s, &f, v);

	// iterate
	int status = GSL_CONTINUE;
	for(int iter = 0; status == GSL_CONTINUE && iter < 500; iter++)
	{
		status = gsl_multifit_fdfsolver_iterate (s);

		print_state (iter, s, ndof);
		if (status && status != GSL_CONTINUE) { break; }; status = 0;
		//if (status) { break; }

		status = gsl_multifit_test_delta (s->dx, s->x, 1e-4, 1e-4);
	}
	fprintf (stderr, "status = %s\n", gsl_strerror (status));
	ASSERT(status == 0);

	// extract fitted parameters
	set_parameters(s->x);
	// expand the covariance matrix
	gsl_matrix *cov = gsl_matrix_alloc (ndof, ndof);
	gsl_multifit_covar (s->J, 0.0, cov);
	int np = nparams;
	covar.resize(np*np); fill(covar.begin(), covar.end(), 0.);
	int y = 0;
	FORj(r, 0, np) // rows
	{
		if(fixed[r]) continue;
		int x = 0;
		FORj(c, 0, np) // columns
		{
			if(fixed[c]) continue;
			covar[r*np + c] = gsl_matrix_get(cov, y, x);
			x++;
		}
		y++;
	}

	// chi^2/DOF
	double chi = gsl_blas_dnrm2(s->f);
	chi2_per_dof = sqr(chi)/ (n - ndof);

	// free the memory	
	gsl_matrix_free(cov);
	gsl_multifit_fdfsolver_free (s);

	if(cullIter--)
	{
		set_parameters(v);
		cull(nsigma);
		fit(cullIter, nsigma);
	}
	gsl_vector_free(v);
}

void load_disk(vector<rzpixel> *data, const std::string &filename)
{
	data->clear();
	text_input_or_die(in, filename);

	rzpixel p;
	//bind(in, p.r, 0, p.z, 1, p.N, 2, p.V, 3); // if loading from cleanedrz.txt
	bind(in, p.r, 0, p.z, 1, p.N, 5, p.V, 6); // if loading from output of median3d.pl
	while(in.next())
	{
		p.rho = p.N / p.V;
		p.sigma = sqrt(p.N) / p.V;
		data->push_back(p);
	}
};

void clean_disk(vector<rzpixel> *data, const std::string &how)
{
	double rbeam;
	if(how == "ngpbeam")
	{
		// find the r value closest to the Sun
		rbeam = 0;
		FOREACH(*data)
		{
			rzpixel &pix = *i;
			if(abs(pix.r-8000.) < abs(rbeam-8000))
				rbeam = pix.r;
		}
		cerr << "Beam radius: r = " << rbeam << "pc\n";
	}

	//
	// Remove the points near the plane of the galaxy
	// and close to the edges of the survey
	//
	vector<rzpixel> out;
	FOREACH(*data)
	{
		rzpixel &pix = *i;
		double phi = deg(atan2(pix.z,pix.r-8000));

 		if(-30 < phi && phi < 15) continue;
 		if(153 < phi && phi < 180) continue;
 		if(-180 < phi && phi < -150) continue;

		if(how == "thin")
		{
			if(!(abs(pix.z) >= 50 && abs(pix.z) <= 600)) continue;
		}
		else if(how == "thick")
		{
//			if(pix.z < 2500 && pix.z > 1000 && pix.r < 8000) continue;
//			if(pix.z < -2000 && pix.r > 8000) continue;

// /*			if(!(abs(pix.z) >= 1000 && abs(pix.z) <= 4000)) continue;
// 			if(!(abs(pix.z) >= 200 && abs(pix.z) <= 3000)) continue;*/
			
//			if(!(abs(pix.z) >= 200 && abs(pix.z) <= 1400)) continue;

			if(abs(pix.z) <= 500) continue;
			if(abs(pix.z) >= 2500) continue;
//			if(!(abs(pix.z) >= 200 && abs(pix.z) <= 3500)) continue;
			if(pix.z > 0 && deg(atan2(pix.z-500,pix.r-8000)) < 30. && pix.r < 10000) continue;
			if(pix.z > 0 && pix.r < 8000.) continue;
		}
		else if(how == "halo")
		{
			if(pix.z > 0 && deg(atan2(pix.z-500,pix.r-8000)) < 30. && pix.r < 10000) continue;
			if(pix.z > 0 && pix.z < 5000 && pix.r < 8000.) continue;

			if(pix.r < 5000) continue;
						
			if(!(abs(pix.z) >= 500)) continue;
		}
		else if(how == "ngpbeam")
		{
			if(pix.r != rbeam) continue;
			if(pix.z < 100) continue;
		}

		out.push_back(pix);
	}

	cerr << data->size() - out.size() << " pixels rejected.\n";
	cerr << out.size() << " pixels accepted.\n";
	*data = out;
}

int main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version);

	Options opts(
		"Program for fitting r-z dataplanes of the Galaxy",
		version,
		Authorship::majuric
	);

	//# add any arguments your program needs. eg:
	opts.argument("binsfile", "File from which to read a list of rz-dataplanes to fit.");

	// add any options your program might need. eg:
	// opts.option("meshFactor", "meshFactor", 0, "--", Option::required, "4", "Resolution decrease between radial steps");

	try {
		opts.parse(argc, argv);
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}

	/////// Start your application code here
	std::string binsfile = opts["binsfile"];
	text_input_or_die(in, binsfile);

//	std::string how = "halo";
//	std::string how = "thick";
//	std::string how = "thin";
	std::string how = "ngpbeam";

	model_fitter m;
	double r0, r1; std::string rzfile, modelname;
	if(how == "thin")
	{
		m.set_param("rho0", 0, false);
		m.set_param("l", 0, false);
		m.set_param("h", 0, false);
		m.set_param("z0", 0, true);
		
		m.set_param("f", 0.0, true);
		m.set_param("lt", 3500, true);
		m.set_param("ht", 400, true);
		m.set_param("fh", 0.0, true);
		m.set_param("q", 1.5, true);
		m.set_param("n", 3, true);

		bind(in, modelname, 0, r0, 1, r1, 2, rzfile, 3,
			m.param("rho0"), 4,
			m.param("l"), 5,
			m.param("h"), 6,
			m.param("z0"), 7
			);
	}
	else if(how == "thick")
	{
		m.set_param("rho0", 0, true);
		m.set_param("l", 0, false);
		m.set_param("h", 0, false);
		m.set_param("z0", 0, true);

//		m.set_param("f", 0.012, true);
		m.set_param("f", 0.03, false);
		m.set_param("lt", 4800, false);
		m.set_param("ht", 1500, false);

		m.set_param("fh", 0.0, true);
		m.set_param("q", 1.5, true);
		m.set_param("n", 3, true);

		bind(in, modelname, 0, r0, 1, r1, 2, rzfile, 3,
			m.param("rho0"), 4,
			m.param("l"), 5,
			m.param("h"), 6,
			m.param("z0"), 7
			);
	}
	else if(how == "halo")
	{
		m.set_param("rho0", 0, false);
		m.set_param("l", 0, true);
		m.set_param("h", 0, true);
		m.set_param("z0", 0, true);
		
		m.set_param("f", 0.0647, true);
		m.set_param("lt", 3300, true);
		m.set_param("ht", 1100, true);

		m.set_param("fh", 0.0015, false);
		m.set_param("q", 2.3, false);
		m.set_param("n", 1.16, true);

		bind(in, modelname, 0, r0, 1, r1, 2, rzfile, 3,
			m.param("rho0"), 4,
			m.param("l"), 5,
			m.param("h"), 6,
			m.param("z0"), 7
			);
	}
	else if(how == "ngpbeam")
	{
		m.set_param("rho0", 0.06, false);
		m.set_param("l", 2500, true);
		m.set_param("h", 270, false);
		m.set_param("z0", 24, true);

		m.set_param("f", 0.04, false);
		m.set_param("lt", 2500, true);
		m.set_param("ht", 1200, false);

		m.set_param("fh", 0.0, true);
		m.set_param("q", 1.5, true);
		m.set_param("n", 3, true);

		bind(in, modelname, 0, r0, 1, r1, 2, rzfile, 3);
	}

	m.print(cout, model_fitter::HEADING);
	cout << "\n";
//	cout << "# name ri0 ri1 chi2/dof rho0 l h z0 err(rho0 l h z0)\n";

	vector<rzpixel> data;
	gsl_vector *v = gsl_vector_alloc(m.ndof);
	while(in.next())
	{
		m.get_parameters(v);

		load_disk(&data, rzfile);
		clean_disk(&data, how);
		m.map = &data;

//		m.param("l") = 1.;

		cerr << "Fitting " << rzfile << "\n";
		m.fit(0, 5);
		m.print(cerr);
		cerr << "norm_thick = " << m.norm_at_Rg() << "\n";
		cerr << "rho(R=8kpc,Z=5kpc) = " << m.rho(8000,5000) << "\n";
		cerr << "rho1(R=8kpc,Z=5kpc) = " << m.rho_thin(8000,5000) << "\n";
		cerr << "rho2(R=8kpc,Z=5kpc) = " << m.rho_thick(8000,5000) << "\n";
		cerr << "rhoh(R=8kpc,Z=5kpc) = " << m.rho_halo(8000,5000) << "\n";

		cout << setw(10) << modelname << setw(7) << r0 << setw(7) << r1 << setw(40) << rzfile << setw(40) << rzfile;
		cout << " " << setw(10) << m.chi2_per_dof << " ";
		m.print(cout, model_fitter::LINE);
		cout << "\n";

		m.set_parameters(v);
	}
}
catch(EAny &e)
{
	e.print();
}
}
