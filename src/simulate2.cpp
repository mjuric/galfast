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

/**
	Basic algorithm outline, definitions:
	  * (x, y) - Lambert map coordinates (transformable to l, b)
	  * m - r band magnitude

	Integrations:
	  * split (x, y) lambert space into 1deg^2 grid
	  * For each nonempty (x,y) pixel (== a conical beam in 3D space):
	     * split r-i into 0.01mag bins
	     * For each r-i bin, calculate:
	        * N(x, y, ri) = Integrate[rho(m; x,y,ri) dm, m_min, m_max]
		  where m is r band magnitude, m_min and m_max are flux
		  limits of the survey
	  * Calculate the sum - this is the total number of stars expected 
	    in survey volume. Normalize by this number to convert the
	    density to probability.
	  * Calculate marginal distributions:
	     * N(ri|x,y), N(y|x), N(x)
	  * Calculate spline approximations to integrals of marginalized
	    distributions for use when transforming uniform variates
	    
	Sampling:
	  * Generating the sample of N stars, step 1
	     * 1: Sample p(x) to get x
	     * Sample p(y|x) to get y
	     * Reject and goto 1 if not in survey volume
	     * Sample p(ri|x,y) to get ri
	     * Store in array, goto 1 if number of stars != N
	  * Calculating magnitudes
	     * Sort the array of generated stars by x,y,ri
	     * For each star generated
	       * If not cached from previous star, calculate spline
	         approximation to p(m | x, y, ri)
	       * Draw a random variate from p(m | x, y, ri) - that is the
	         magnitude.
*/

namespace lapack
{
	extern "C"
	{
		#include <g2c.h>
		#include "common/clapack.h"
	}
}

extern "C"
{
	#include "gpc/gpc.h"
}

#include <nurbs++/nurbsS.h>

#include <boost/tuple/tuple.hpp>

#include "projections.h"
#include "model.h"

#include <vector>
#include <map>
#include <string>

#include <astro/math.h>
#include <astro/util.h>
#include <astro/system/log.h>
#include <astro/useall.h>

void poly_bounding_box(double &x0, double &x1, double &y0, double &y1, const gpc_polygon &p);

gpc_polygon poly_rect(double x0, double x1, double y0, double y1)
{
	static gpc_vertex v[4];
	static gpc_vertex_list vl = {4, v};
	static gpc_polygon rect = {1, NULL, &vl};

	v[0].x = x0; v[0].y = y0;
	v[1].x = x1; v[1].y = y0;
	v[2].x = x1; v[2].y = y1;
	v[3].x = x0; v[3].y = y1;

	return rect;
}

#if 0
typedef double (*fun_2d)(double x, double y);
void generic_nurbs_approximation(PlNurbsSurfacef &surf, fun_2d f,
	double x0, double x1, double dx,
	double y0, double y1, double dy
	)
{
	using namespace PLib;

	// calculate the number of matrix points needed
	int nx = (int)floor((x1 - x0) / dx) + 1; //if(nx * dx != (x1-x0)) { nx++; }
	int ny = (int)floor((y1 - y0) / dy) + 1; //if(ny * dy != (y1-y0)) { ny++; }

	Matrix_Point3Df pts(nx, ny);
	for(int i = 0; i < nx; i++)
	{
		double x = x0 + i*dx;
		if(i == nx-1) { x = x1; }
		for(int j = 0; j < ny; j++)
		{
			double y = y0 + j*dy;
			if(j == ny-1) { y = y1; }
			
			double fval = f(x, y);
			pts(i, j) = Point3Df(x, y, fval);
		}
	}

	surf.globalInterp(pts,3,3);
	surf.writeVRML97("surface.wrl");
}

double surffun(double x, double y)
{
	return 1+x*y;//*y*y;
}

void test_nurbs()
{
	using namespace PLib;

	PlNurbsSurfacef surf;
	generic_nurbs_approximation(surf, surffun,
		-1, 1, .1,
		-1, 1, .1);

	Point3Df pt(.66, .66, 1);
	float u = .66, v = .66;
	cerr << surf.projectOn(pt, u, v, 10000) << "\n";
	cerr << u << " " << v << "\n";
	cerr << surf(u, v) << "\n";
/*	for(int i = 0; i != 100; i++)
	{
		double x = math::rnd(-1, 1);
		double y = math::rnd(-1, 1);
		double v = surffun(x, y);
		double va = surf(x, y);
		cerr << "(x, y, v, vapprox, eps) = "
			<< x << ", " << y << ": " << v << " " << va << " => " <<
			(va - v)/v << "\n";
	}*/
}
#endif

lambert proj(rad(90), rad(90));
double dx = deg(1.); // model grid angular resolution
double dri = 0.01; // model CMD resolution
double ri0 = 0.0, ri1 = 1.5;	// color limits
double m0 = 14, m1 = 22;	// magnitude limits
double dm = 0.01;	// model CMD magnitude resolution
typedef std::map<float, model> models_t;
std::map<float, model> models; // a model for a given color ri

#if 0
distanceSampler::distanceSampler(const direction &d0, const double m0, const double m1)
{
	int n = ceil((m1-m0)/dm);
	std::valarray<double> m(n+1), y(n+1);
	FOR(0, xv.size()) { xv[i] = m0 + dm*i; }
	m[m.size()-1] = m1; // last entry is the flux limit

	// calculate N(m) splines for colors where the model has been given
	FOREACH(models_t::iterator, models)
	{
	}

	direction d = d0;
	sample_integral(xv, yv, func, (void *)&d);

	// normalization - convert the integrated density to probability
	N = yv[yv.size()-1];
	yv /= N;

	// construct spline approximation
	f.construct(&yv[0], &xv[0], yv.size());
}
#endif

void do_cmd_integrals(std::vector<double> &pdf, const double x, const double y)
{
	// deproject to (l,b)
	Radians l, b;
	proj.inverse(x, y, l, b);

	// TODO: continue here...
}

void integrations(const std::string &prefix)
{
	gpc_polygon sky;
	FILE *fp = fopen((prefix + ".gpc.txt").c_str(), "r");
	gpc_read_polygon(fp, 1, &sky);
	fclose(fp);

	std::map<std::pair<int, int>, gpc_polygon> skymap;	
	std::map<std::pair<int, int>, std::vector<double> > pdf; // N(x, y)(ri), marginalized over magnitude
	double x0, x1, y0, y1;
	poly_bounding_box(x0, x1, y0, y1, sky);
	int X = 0, Y = 0;
	for(double x = x0; x < x1; x += dx, X++)
	{
		double xa = x, xb = x+dx;
		for(double y = y0; y < y1; y += dx, Y++)
		{
			double ya = y, yb = y+dx;
			gpc_polygon r = poly_rect(xa, xb, ya, yb);

			gpc_polygon p;
			gpc_polygon_clip(GPC_INT, &sky, &r, &p);
			if(p.num_contours == 0) continue;

			skymap[std::make_pair(X, Y)] = p;

			// calculate CMD integrals
			do_cmd_integrals(pdf[std::make_pair(X, Y)], (xa + xb)/2., (ya + yb)/2.);
		}
	}
}

#if 0

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <valarray>

// template lapack_gb_matrix<typename T>
// 	class lapack_gb_matrix : public banded_matrix<T, boost::numeric::ublas::column_major>
// 	{
// 		
// 	}

namespace lapack {
// 	int dgbsv_(integer *n, integer *kl, integer *ku, integer *
// 	nrhs, doublereal *ab, integer *ldab, integer *ipiv, doublereal *b, 
// 	integer *ldb, integer *info);
	inline int solve(
		boost::numeric::ublas::vector<doublereal> &xy,
		boost::numeric::ublas::banded_matrix<doublereal> &m, 
		boost::numeric::ublas::vector<integer> &ipiv)
	{
		ASSERT(m.size1() == m.size2());
		ASSERT(xy.size() == m.size1());

		integer n = m.size1(), kl = m.lower(), ku = m.upper()-kl, nrhs = 1, ldab = 2*kl+ku+1;
/*		cout << n << " " << kl << " " << ku << " " << nrhs << " " << ldab << "\n";
		exit(-1);*/
		integer info;
		dgbsv_(&n, &kl, &ku, &nrhs, &m.data()[0], &ldab, &ipiv[0], &xy[0], &n, &info);
		return info;
	}

	inline int solve(
		boost::numeric::ublas::vector<doublereal> &xy,
		boost::numeric::ublas::banded_matrix<doublereal> &m)
	{
		boost::numeric::ublas::vector<integer> ipiv(m.size1());
		return solve(xy, m, ipiv);
	}
}

void prettyprint(double *v, int size)
{
	FOR(0, size)
	{
		if(i % 5 == 0) { cout << "\n"; }
		if(v[i] == 0) { cout << " * "; } else { cout << v[i] << " ";}
	}
	cout << "\n";
}

#include <iomanip>
void mprint(boost::numeric::ublas::banded_matrix<double> &m, boost::numeric::ublas::vector<double> *v = NULL)
{
	using namespace std;
	int n = m.size1(), kl = m.lower(), ku = m.upper();//-kl;
	FORj(r, 0, n)
	{
		FORj(c, 0, n)
		{
			r++; c++;
			if(std::max(1, c-ku) <= r && r <= std::min(n, c+kl))
			{
				cerr << setw(3) << m(r-1, c-1);
			} else {
				cerr << "  *";
			}
			r--; c--;
		}
		if(v != NULL)
		{
			cerr << setw(10) << (*v)[r];
		}
		cerr << "\n";
	}
	cerr << "\n";
}

void fit_int_spline();
void test_lapack()
{
 	fit_int_spline();
 	return;

	using namespace boost::numeric;
	using namespace lapack;
	using namespace std;

	integer n = 5, kl = 2, ku = 1, nrhs = 1;
	integer ldab = 2*kl+ku+1;
	ublas::banded_matrix<doublereal> m (n, n, kl, ku+kl);
	m(0, 0) = 11; m(0, 1) = 12;
	m(1, 0) = 21; m(1, 1) = 22; m(1, 2) = 23;
	m(2, 0) = 31; m(2, 1) = 32; m(2, 2) = 33; m(2, 3) = 34;
	              m(3, 1) = 42; m(3, 2) = 43; m(3, 3) = 44; m(3, 4) = 45;
	                            m(4, 2) = 53; m(4, 3) = 54; m(4, 4) = 55;
	double* v = &m.data()[0];
	prettyprint(v, m.data().size());
	mprint(m);
	ublas::vector<long int> ipiv(n);
	ublas::vector<double> rhs(n);
	FOR(0, 5) { rhs[i] = (double)i+1; }
	FOR(0, 5) { cerr << rhs[i] << "\t" << ipiv[i] << "\n"; }
	int info = lapack::solve(rhs, m, ipiv);

	cerr << "INFO = " << info << "\n";

	//cout << m << "\n";
	prettyprint(v, m.data().size()); cerr << "\n";
	mprint(m);
	FOR(0, 5) { cerr << rhs[i] << "\t" << ipiv[i] << "\n"; }
}

// double ib[] =
// {
// 	32, 27, 28, 15, 32, 8, 4, 1,
// 	1, 1, 1, 1, -1, 0, 0, 0,
// 	0, 1, 2, 3, 0, -1, 0, 0,
// 	0, 0, 2, 6, 0, 0, -2, 0
// };
// 
// double ib0[] =
// {
// 	32, 8, 4, 1,
// 	0, 1, 0, 0
// };
// 
// double ib1[] =
// {
// 	32, 27, 28, 15,
// 	0, 1, 2, 3
// };

void copy_block(boost::numeric::ublas::banded_matrix<double> &m, int r0, int c0, double *b, int w, int h)
{
	--b;
	FORj(i, 0, h)
		FORj(j, 0, w)
			if(*(++b) != 0.)
				m(r0+i, c0+j) = *b;
}

#if 0

double ib[] =
{
	1, 1, 1, 1, -1,  0,  0, 0,
	0, 1, 2, 3,  0, -1,  0, 0,
	0, 0, 1, 3,  0,  0, -1, 0,
	0, 0, 0, 0, 12,  6,  4, 3
};

double ib0[] =
{
	 0, 1, 0, 0,
	 0, 0, 1, 0,
	12, 6, 4, 3
};

double ib1[] =
{
	0, 1, 2, 3
};

#include <algorithm>

class ispline
{
public:
	double eval(double t, double *c)
	{
		return (((t*c[3]) + c[2])*t + c[1])*t + c[0];
	}
public:
	int nseg;
	std::valarray<double> x, h, c;
public:
	int find_x(double xv)
	{
		using namespace std;
		double *x = &this->x[0];
		int n = this->x.size();
		double *v = upper_bound(x, x+n, xv);
		--v;
		if(v-x >= n-1) { --v; } // if it's beyond the last segment, extrapolate
		return v-x;
	};

	double eval(double xv)
	{
		int i0 = find_x(xv);
		//cout << i0 << " ";
		double dx = x[i0+1] - x[i0];
		double t = (xv - x[i0])/dx;
		return eval(t, &c[4*i0]);
	};

	ispline(int np, double *vx, double *vy, double *coef)
		: nseg(np-1), x(vx, np), h(vy, np-1), c(coef, 4*(np-1))
	{}
};

void fit_int_spline()
{
	using namespace boost::numeric;
	using namespace lapack;
	using namespace std;

	valarray<double> x(6), h(x.size()-1);
	x[0] = 0; h[0] = 1;
	x[1] = 1; h[1] = 2;
	x[2] = 2; h[2] = 4;
	x[3] = 3; h[3] = 4;
	x[4] = 4; h[4] = 5;
	x[5] = 5;

	// construct and fill in the spline matrix and rhs vector
	int ns = h.size();	// number of segments	cout << "

	int n = 4*ns;
	ublas::vector<double> rhs(n); rhs *= 0.;
	ublas::banded_matrix<double> m(n, n, 3, 1 + 3);

	copy_block(m, 0, 0, ib0, 4, 3);
//	rhs[0] = 2.*(h[1]-h[0])/(x[2]-x[0]);	// D1
	rhs[2] = 12. * h[0] / (x[1] - x[0]);	// H1

	FOR(1, ns)
	{
		cerr << 4*i-1 << "\n";
		copy_block(m, 4*i-1, 4*(i-1), ib, 8, 4);
		rhs[4*i+2] = 12. * h[i] / (x[i] - x[i-1]);
	}

	copy_block(m, n-1, n-4, ib1, 4, 1);
//	rhs[n-1] = 2 * (h[n-1]-h[n-2]) / (x[n] - x[n-2]); // Dn
	mprint(m, &rhs);
	
	// solve spline matrix
 	int info = lapack::solve(rhs, m);
 	cerr << "LAPACK info = " << info << "\n\n";
 	mprint(m, &rhs);

	// repack to ispline data structure
	ispline isp(x.size(), &x[0], &h[0], &rhs[0]);
	for(double xx = 0; xx <= x[x.size()-1]+0.1; xx += 0.1)
	{
		cout << xx << "\t" << isp.eval(xx) << "\n";
	}
}

#else

double ib[] =
{
	1, 1, 1, -1, 0, 0,
	0, 1, 2, 0, -1, 0,
	0, 0, 0, 6,  3, 2
};

double ib0[] =
{
	 0, 1, 0,
	 6, 3, 1
};

double ib1[] =
{
	0, 1, 2
};

#include <algorithm>

class ispline
{
public:
	double eval(double t, double *c)
	{
		return ((c[2])*t + c[1])*t + c[0];
	}
public:
	int nseg;
	std::valarray<double> x, h, c;
public:
	int find_x(double xv)
	{
		using namespace std;
		double *x = &this->x[0];
		int n = this->x.size();
		double *v = upper_bound(x, x+n, xv);
		--v;
		if(v-x >= n-1) { --v; } // if it's beyond the last segment, extrapolate
		return v-x;
	};

	double eval(double xv)
	{
		int i0 = find_x(xv);
		//cout << i0 << " ";
		double dx = x[i0+1] - x[i0];
		double t = (xv - x[i0])/dx;
		return eval(t, &c[3*i0]);
	};

	ispline(int np, double *vx, double *vy, double *coef)
		: nseg(np-1), x(vx, np), h(vy, np-1), c(coef, 3*(np-1))
	{}
};

///////////////////////////////////

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <valarray>

class spline
{
private:
	gsl_interp *f;
	gsl_interp_accel *acc;
	std::valarray<double> xv, yv;
public:
	spline() : f(NULL), acc(NULL) {}
	spline(const double *x, const double *y, int n);
	void construct(const double *x, const double *y, int n);
	~spline();

 	double operator ()(double x)        { return gsl_interp_eval(f, &xv[0], &yv[0], x, acc); }
 	double deriv(double x)              { return gsl_interp_eval_deriv(f, &xv[0], &yv[0], x, acc); }
 	double deriv2(double x)             { return gsl_interp_eval_deriv2(f, &xv[0], &yv[0], x, acc); }
 	double integral(double a, double b) { return gsl_interp_eval_integ(f, &xv[0], &yv[0], a, b, acc); }

public:
	spline& operator= (const spline& a);
	spline(const spline& a) : f(NULL), acc(NULL) { *this = a; }
};

///////////////////////////////////

void fit_int_spline2(std::valarray<double> &x, std::valarray<double> &h)
{
	double sum = 0;
	FOR(0, h.size()-1)
	{
		h[i] = sum;
		sum += h[i+1];
		cerr << sum << "\n";
	}
	h[h.size()-1] = sum;
	cerr << x.size() << " " << h.size() << "\n";

	spline spl;
	spl.construct(&x[0], &h[0], h.size());
	for(double xx = 0; xx <= x[x.size()-1]+0.1; xx += 0.1)
	{
		cout << xx << "\t" << spl.deriv(xx) << "\n";
	}
	exit(0);
}

void fit_int_spline()
{
	using namespace boost::numeric;
	using namespace lapack;
	using namespace std;

	valarray<double> x(6), h(x.size()-1);
	x[0] = 0; h[0] = 1;
	x[1] = 1; h[1] = 2;
	x[2] = 2; h[2] = 4;
	x[3] = 3; h[3] = 4;
	x[4] = 4; h[4] = 5;
	x[5] = 5;

	fit_int_spline2(x, h);
	
	// construct and fill in the spline matrix and rhs vector
	int ns = h.size();	// number of segments	cout << "

	int n = 3*ns;
	ublas::vector<double> rhs(n); rhs *= 0.;
	ublas::banded_matrix<double> m(n, n, 2, 1 + 3);

	copy_block(m, 0, 0, ib0, 3, 2);
//	rhs[0] = 2.*(h[1]-h[0])/(x[2]-x[0]);	// D1
	rhs[1] = 6. * h[0] / (x[1] - x[0]);	// H1

	FOR(1, ns)
	{
		copy_block(m, 3*i-1, 3*(i-1), ib, 6, 3);
		rhs[3*i+1] = 6. * h[i] / (x[i] - x[i-1]);
	}

	copy_block(m, n-1, n-3, ib1, 3, 1);
//	rhs[n-1] = 2 * (h[n-1]-h[n-2]) / (x[n] - x[n-2]); // Dn
	mprint(m, &rhs);
	
	// solve spline matrix
 	int info = lapack::solve(rhs, m);
 	cerr << "LAPACK info = " << info << "\n\n";
 	mprint(m, &rhs);

	// repack to ispline data structure
	ispline isp(x.size(), &x[0], &h[0], &rhs[0]);
	for(double xx = 0; xx <= x[x.size()-1]+0.1; xx += 0.1)
	{
		cout << xx << "\t" << isp.eval(xx) << "\n";
	}
}

#endif
#endif
