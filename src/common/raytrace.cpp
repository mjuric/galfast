#include "raytrace.h"
#include "interval_arithmetic.h"
#include "dm.h"
#include "projections.h"

#include <astro/math/vector.h>
#include <astro/coordinates.h>
#include <astro/math.h>
#include <astro/util.h>
#include <astro/sdss/rungeometry.h>
#include <astro/io/fits.h>
#include <astro/system/options.h>

#include <astro/useall.h>

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>

#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sort.h>

#include "ximage.h"
#include "floodfill.h"
#include "analysis.h"
#include "binarystream.h"

using namespace std;

#define DBG 0

//
// Wrapper for printing out V3s in [abs(v), lon, lat] format
//

struct celestial
{
	const V3 &v;

	celestial(const V3 &v_) : v(v_) {}
};

OSTREAM(const celestial &c)
{
	out << "[r=" << abs(c.v)
	    << ", lon=" << deg(c.v.phi())
	    << ", lat=" << deg(ctn::pi/2 - c.v.theta())
	    << "]";
	return out;
}

OSTREAM(const M3 &m)
{
	FOR(0, 3) { out << m(i, 0) << " " << m(i, 1) << " " << m(i, 2) << "\n"; }
	return out;
}

M3 &M3::rotate(const V3 &axis, Radians angle)
{
	// TODO
	ASSERT(0);
}

V3 operator *(const M3 &a, const V3 &b)
{
	V3 ret(0.);
	FORj(r, 0, 3) FORj(c, 0, 3)
	{
		ret[r] += a(r, c) * b[c];
	}
	return ret;
}

M3 operator *(const M3 &a, const M3 &b)
{
	M3 ret(0.);
	FORj(i, 0, 3) FORj(j, 0, 3) FORj(k, 0, 3)
	{
		ret(i, k) += a(i, j) * b(j, k);
	}
	return ret;
}

M3 operator +(const M3 &a, const M3 &b)
{
	M3 ret;
	FORj(r, 0, 3) FORj(c, 0, 3) { ret(r, c) = a(r, c) + b(r, c); }
	return ret;
}

//
// structure describing a quadric form:
//	g(x, y) = a x^2 + b x y + c y^2 + d x + e y + f
//
struct Quadric
{
	double a, b, c, d, e, f;
	Quadric(double a_ = 0, double b_ = 0, double c_ = 0, double d_ = 0, double e_ = 0, double f_ = 0)
		: a(a_), b(b_), c(c_), d(d_), e(e_), f(f_)
	{ }

	double operator()(double x, double y) const { return g(x, y); }
	double g(double x, double y) const
	{
		return a*sqr(x) + b*x*y + c*sqr(y) + d*x + e*y + f;
	}
};

OSTREAM(const Quadric &q)
{
	return out << q.a << " " << q.b << " " << q.c << " " << q.d << " " << q.e << " " << q.f;
}

struct StandardQuadric : public Quadric
{
	// matrices for conversion from the coordinate system of the principal axes
	// quadric to the original coordinate system
	double cosf, sinf;	// rotation matrix from principal axes to zero-centered coordinate system
	double h, k;		// offset to add to move from zero-centered to original coordinate system

	StandardQuadric() {}
	StandardQuadric(const Quadric &q);

	void rotate(double x, double y, double &xr, double &yr, double dir = 1);	// rotate from standard to zero-centered coordinate system

	void toOriginal(double x, double y, double &xo, double &yo);
	void toPrincipal(double x, double y, double &xp, double &yp);
};

OSTREAM(const StandardQuadric &q)
{
	out << (Quadric &)q;
	out << " | " << q.h << " " << q.k << " | " << q.cosf << " " << q.sinf << " | " << deg(atan2(q.sinf, q.cosf));
	return out;
}

void StandardQuadric::rotate(double x, double y, double &xr, double &yr, double dir)
{
	// rotate x, y clockwise
	xr = cosf * x + dir * sinf * y;
	yr = -dir * sinf * x + cosf * y;
}

void StandardQuadric::toOriginal(double x, double y, double &xo, double &yo)
{
	rotate(x, y, xo, yo);
	xo += h; yo += k;
}

void StandardQuadric::toPrincipal(double x, double y, double &xp, double &yp)
{
	x -= h; y -= k;
	rotate(x, y, xp, yp, -1);
}

StandardQuadric::StandardQuadric(const Quadric &q)
{
	*((Quadric *)this) = q;

	// translate to zero-centered coordinate system, if not already there
	if(d != 0 || e != 0) {
		double D = sqr(b) - 4 * a * c;
		if(D == 0) { ASSERT(0); } // parabola or a line - TODO: handle it

		h = (2 * c * d - b * e) / D;
		k = (2 * a * e - b * d) / D;
		f = q(h, k);
		if(f == 0) { ASSERT(0); }  // degenerate case - point or lines - TODO: handle it
		d = e = 0;
	} else {
		h = k = 0;
	}

	// find eigenvalues
	if(b == 0) {
		// already diagonal
		cosf = 1; sinf = 0;
	} else {
/*		double D = sqrt(sqr(a - c) + sqr(b));
		double lc = ( a + c + D) / 2;
		double la = ( a + c - D) / 2;
		double a1 = (-a + c + D) / b;
		double a2 = (-a + c - D) / b;
		cosf = a1/sqrt(sqr(a1)+1);
		sinf = a2/sqrt(sqr(a2)+1);

		// set new values
		a = la;
		b = 0;
		c = lc;
*/
		double D = sqrt(sqr(a - c) + sqr(b));
		double l1 = (a + c + D) / 2,
		       l2 = (a + c - D) / 2;
		V2 e1(a - l1, b / 2), e2(a - l2, b / 2);
		e1 /= abs(e1); e2 /= abs(e2);
//		cout << "Det = " << e1.x * e2.y - e2.x * e1.y<< "\n";
		if(e1.x * e2.y - e2.x * e1.y < 0) { e2.x *= -1;  }
		cosf = e1.x; sinf = e2.x;
		a = l2;
		b = 0;
		c = l1;
#if 0
		gsl_eigen_symmv_workspace *w = gsl_eigen_symm_alloc(2);
		double m[] = { 
			a,   b/2,
			b/2, c     };
		gsl_matrix_view M = gsl_matrix_view_array(m, 2, 2);
		gsl_vector *eval = gsl_vector_alloc(2);
		gsl_matrix *evec = gsl_matrix_alloc(2, 2);

		gsl_eigen_symmv(A, eval, evec, w);

		gsl_eigen_symmv_free(w);
		gsl_vector_free(eval);
		gsl_matrix_free(evec);
#endif
	}
}

class UCP : public V2 // Unit Circle Point (has ordering)
{
public:
	Radians ph;
	enum Type { unset, cone, sphere, plane } type;
public:
	UCP() : V2(0, 0), ph(0), type(unset) {}
	UCP(double x_, double y_, double x0, double y0, Type type_) : type(type_), V2(x_ - x0, y_ - y0), ph(atan2(y, x)) { }
	UCP(const UCP &a) : V2(a), ph(a.ph), type(a.type) { }

	bool operator <(const UCP &a) const
	{
		return ph < a.ph;
	}

	Radians phi() const { return ph; }
};
OSTREAM(const UCP &u)
{
	out << (const V2 &)u << " " << deg(u.phi());

	switch(u.type) {
		case UCP::cone: out << " cone"; break;
		case UCP::sphere: out << " sphere"; break;
		case UCP::plane: out << " plane"; break;
		default: out << " unset"; break;
	}

	return out;
}

struct Intersector : public StandardQuadric
{
	gsl_poly_complex_workspace *w;

	// cone-circle intersection
	double bcoef[5];	// b4 and b3 must be preset
	double b20, b22, b10, b12, b00, b02, b04;	// these are the r coefficients
	double xp, yp;	// circle coordinates in principal axis coordinate system
	double xc, yc;	// circle coordinates in galactic coordinate system

	Intersector(const StandardQuadric &qq, double xc, double yc);
	~Intersector();

	int intersect(double r, vector<UCP> &roots);
};

Intersector::Intersector(const StandardQuadric &qq, double xc_, double yc_)
	: StandardQuadric(qq), xc(xc_), yc(yc_)
{
	toPrincipal(xc, yc, xp, yp);

	// precalculate various coefficients needed by intersect
	bcoef[4] = sqr(a - c);
	bcoef[3] = 4*a*(c-a)*yp;

	double x2 = sqr(xp);
	double y2 = sqr(yp);
	double a2 = sqr(a);
	double f2 = sqr(f);
	double rho2 = x2 + y2;
	double del2 = x2 - y2;

	b20 = 2*(a*(c*del2 - f) + c*f + a2*(x2 + 3*y2));
	b22 = 2*a*(-a + c);

	b10 = -4*a*yp*(-f + a*rho2);
	b12 = 4*a2*yp;

	b00 = 2*a*del2*f + f2 + a2*sqr(rho2);
	b02 = -2*a*(-f + a*rho2);
	b04 = a2;

	w = gsl_poly_complex_workspace_alloc(5);

#if 0
	cout << "--- xc --\n";
	cout << xc << " " << yc << "\n";
	cout << xp << " " << yp << "\n";
	cout << "--- b ---\n";
	cout << bcoef[4] << "\n";
	cout << bcoef[3] << "\n";
	cout << b20 << " " << b22 << "\n";
	cout << b10 << " " << b12 << "\n";
	cout << b00 << " " << b02 << " " << b04 << "\n";
	cout << "---   ---\n";
#endif

}

Intersector::~Intersector()
{
	gsl_poly_complex_workspace_free(w);
}

struct cxtmp {
	double x, y;
	bool operator <(const cxtmp &a) const { return x < a.x; }
};

#define DBG2 0

#if DBG2
static double zdebug;
#endif

// circle center must be in principal axis coordinate system
// roots are returned in original coordinate system
int Intersector::intersect(double r, vector<UCP> &roots)
{
	double r2 = sqr(r);

	bcoef[2] = b20 + r2 * b22;
	bcoef[1] = b10 + r2 * b12;
	bcoef[0] = b00 + r2 * (b02 + r2 * b04);

	// call GSL solve magic, store result in yy[] array
	cxtmp yy[4];
	int status = gsl_poly_complex_solve(bcoef, 5, w, (double *)yy);
	if(status != GSL_SUCCESS)
	{
		cout << "root solving qr method failed to converge\n";
		cout << "r = " << r << " ";
#if DBG2
		cout << "z = " << zdebug << "\n";
#endif
		cout << "\n";
		FOR(0, 6) { cout << bcoef[i] << " "; }
		cout << "\n";
		return 0;
	}
	sort(yy, yy + 4);

	//
	const double eps = 1e-5;
	double lasty;
	int n = 0;
	FOR(0, 4)
	{
//		cout << "Root: " << yy[i].x << " " << yy[i].y << "\n";
		if(abs(yy[i].y) > eps) continue; // complex root

		// calculate x
		double x, y = yy[i].x;
		double D = r2 - sqr(y - yp);
		if(D < 0) continue;
		D = sqrt(D);

		// check which is the right root and store it
		// also take care of roots with multiplicity > 2
		// assuming roots are sorted
//		cout << "lasty = " << lasty << ", y = " << y << " " << abs(lasty - y) << "\n";
		if(roots.size() && abs(lasty - y) < eps)
		{
			x = xp - D;
		} else {
			x = (abs(g(xp + D, y)) < eps) ? xp + D : xp - D;
		}
		lasty = y;

		// convert to original
//		cout << "root(sq coords) = " << UCP(x, y, xp, yp) << "\n";
		toOriginal(x, y, x, y);
//		cout << "root(orig coords) = " << UCP(x, y, xc, yc) << "\n";
		roots.push_back(UCP(x, y, xc, yc, UCP::cone));

		n++;
	}
//	cout << "N = " << n << "\n";
	return n;
}


///////////

struct PlaneIntersector : public Plane
{
	double coef[3];
	double c00, c02;
	double cy0, cy1;
	double x0;
	
	double xx(double y) const { return cy0 + cy1 * y; }
	int intersect(double r, vector<UCP> &roots);

	PlaneIntersector(const Plane &p, double xc, double yc, double z0);
};

PlaneIntersector::PlaneIntersector(const Plane &p, double x0, double y0, double z0)
{
	// assuming w == 0

	((Plane &)*this) = p;
	this->x0 = x0;

	double a2 = sqr(a), b2 = sqr(b), c2 = sqr(c), x02 = sqr(x0), y02 = sqr(y0), z02 = sqr(z0);

	c02 = -a2;
	c00 = a2*(x02 + y02) + 2*a*c*x0*z0 + c2*z02;
	coef[1] = 2*a*b*x0 - 2*a2*y0 + 2*b*c*z0;
	coef[2] = a2 + b2;
	
	cy0 = -c * z0/a;
	cy1 = -b/a;
}

int PlaneIntersector::intersect(double r, vector<UCP> &roots)
{
	double y0, y1;
	double r2 = sqr(r);
	coef[0] = c00 + c02 * r2;
	if(gsl_poly_solve_quadratic(coef[2], coef[1], coef[0], &y0, &y1) == 2)
	{
		if(a != 0) {
			roots.push_back(UCP(xx(y0), y0, x0, 0, UCP::plane));
			roots.push_back(UCP(xx(y1), y1, x0, 0, UCP::plane));
		} else {
			// Here we assume y0 == 0
			double D = sqrt(r2 - sqr(y0));
			roots.push_back(UCP(x0 + D, y0, x0, 0, UCP::plane));
			roots.push_back(UCP(x0 - D, y1, x0, 0, UCP::plane));
		}
		return 2;
	}
	return 0;
}

///////////

struct SphereIntersector
{
	double c00, c02;
	double d2;
	double xc;

	int intersect(double r, vector<UCP> &roots);

	SphereIntersector(double d, double xc, double z0);
};

SphereIntersector::SphereIntersector(double d, double x0, double z0)
	: xc(x0)
{
	d2 = sqr(d) - sqr(z0);
	c00 = (d2 + sqr(x0)) / (2 * x0);
	c02 = -1/(2 * x0);
}

int SphereIntersector::intersect(double r, vector<UCP> &roots)
{
	if(d2 < 0) return 0;

	double x = c00 + c02 * sqr(r);
	double x2 = sqr(x);
	double y2 = d2 - x2;
	if(y2 < 0) return 0;
	double y = sqrt(y2);

	roots.push_back(UCP(x, -y, xc, 0, UCP::sphere));
	roots.push_back(UCP(x, y, xc, 0, UCP::sphere));

	return 2;
}

///////////


typedef double Bangle;

inline Bangle barctan2(double y, double x)
{
	return y > 0 ? 1 - x : 3 + x;
}

inline void btan2(Bangle a, double &y, double &x)
{
	x = a < 2 ? 1 - a : a - 3;
	y = sqrt(1 - sqr(x));
}


// intersect z plane with a cone oriented according to matrix M and with an opening
// tangent squared t2
void intersect_cone(Quadric &q, const double t2, const M3 &M, const double z0)
{
	const double a = M(0, 0), b = M(0, 1), c = M(0, 2);
	const double d = M(1, 0), e = M(1, 1), f = M(1, 2);
	const double g = M(2, 0), h = M(2, 1), i = M(2, 2);

	double a2 = sqr(M(0, 0)), b2 = sqr(M(0, 1)), c2 = sqr(M(0, 2));
	double d2 = sqr(M(1, 0)), e2 = sqr(M(1, 1)), f2 = sqr(M(1, 2));
	double g2 = sqr(M(2, 0)), h2 = sqr(M(2, 1)), i2 = sqr(M(2, 2));

	q.a = g2 - (a2 + d2)*t2;
	q.b = -2*(-(g*h) + (a*b + d*e)*t2);
	q.c = h2 - (b2 + e2)*t2;
	q.d = 2*(g*i - (a*c + d*f)*t2)*z0;
	q.e = 2*(h*i - (b*c + e*f)*t2)*z0;
	q.f = (i2 - (c2 + f2)*t2)*sqr(z0);
}

//
// Rotation matrix for transforming vectors from geocentric galactic to run coordinates
//
void rotation_matrix(M3 &M, Radians node, Radians inc)
{
	// Coordinate system pointing _towards_ the galactic center
	//
	// x^ = (l =   0, b = 0)
	// y^ = (     90,     0)
	// z^ = (       ,    90)
	//

	double mu[3], nu[3];
// NOTE!! I changed _temporarily_ for testing galgcs->equgcs
	coordinates::galgcs(node, inc, rad(0),  rad(0), mu[0], nu[0]);
	coordinates::galgcs(node, inc, rad(90), rad(0), mu[1], nu[1]);
	coordinates::galgcs(node, inc, rad(0),  rad(90), mu[2], nu[2]);

	// convert to cartesian, in run coordinate system
	V3 x, y, z;
	x.spherical(1, ctn::pi/2 - nu[0], mu[0]);
	y.spherical(1, ctn::pi/2 - nu[1], mu[1]);
	z.spherical(1, ctn::pi/2 - nu[2], mu[2]);

	// set the matrix
	M(0) = x;
	M(1) = y;
	M(2) = z;
}

void rotation_matrix(const RunGeometry &geom, M3 &M)
{
	rotation_matrix(M, geom.node, geom.inc);
}

void matrix_transpose(M3 &Mt, const M3 &M)
{
	FORj(i, 0, 3) FORj(j, 0, 3) { Mt(j, i) = M(i, j); }
}

class ExtendedRunInfo
{
public:
	static auto_ptr<RunGeometryDB> geomDB;

public:
	sdss::RunGeometry geom;
	sdss::Mask mask;

	M3 M;				// rotation matrix from galactic to run coordinates
	M3 Minv;			// rotation matrix from run to galactic coordinates
	double t[6][2]; 		// tangents of limiting nu angles for each column
	double d0, d1;			// distance limits (parsecs)
	V3 p0, p1;			// normals of run bound planes (corresponding to mu0, mu1);
	V3 p0g, p1g;			// normals of run bound planes in galactic coordinates

	ExtendedRunInfo(int run, double d0_, double d1_)
		: d0(d0_), d1(d1_)
	{
		geomDB->getGeometry(run, geom);

// DEBUGGING
//geom.muStart = rad(0);
//geom.muEnd = rad(10);

		mask = sdss::Mask(geom);

		FOR(0, 6)
		{
			t[i][0] = tan(mask.lo(i));
			t[i][1] = tan(mask.hi(i));
		}

		rotation_matrix(geom, M);
		matrix_transpose(Minv, M);

		p0.cylindrical(1, geom.muStart, 0); p0.y = -p0.y; swap(p0.x, p0.y);
		p1.cylindrical(1, geom.muEnd,   0); p1.y = -p1.y; swap(p1.x, p1.y);

		p0g = Minv * p0;
		p1g = Minv * p1;
	}

	bool contains(const V3 &v, int col) const;
};
auto_ptr<RunGeometryDB> ExtendedRunInfo::geomDB;

OSTREAM(const ExtendedRunInfo &e)
{
	out << e.M << "\n" << e.Minv << "\n" << e.M * e.Minv << "\n";
	FOR(0, 6) { out << e.t[i][0] << " " << e.t[i][1] << " == " << deg(e.mask.lo(i)) << " " << deg(e.mask.hi(i)) << "\n"; }
	out << "\n";
	out << "muStart=" << deg(e.geom.muStart) << " muEnd=" << deg(e.geom.muEnd) << "\n";
	out << e.p0 << " " << e.p1 << "\n";
	out << e.p0g << " " << e.p1g << "\n";
	out << "\n";
	return out;
}

//
// v is in galactic cartesian coordinates
//
bool ExtendedRunInfo::contains(const V3 &x, int col) const
{
	// test distance limits
	double d = abs(x);
	if(d < d0 || d > d1) return false;

	V3 v(M*x);

	// test if it's within opening angle
	double tn = v.z / v.rho();
	if(tn < t[col][0] || tn > t[col][1]) return false;

	// test phi direction
	double phi = v.phi();
	if(phi < 0) { phi += ctn::pi2; }

	if(geom.muEnd < geom.muStart) { return phi < geom.muEnd || phi > geom.muStart; }
	return geom.muStart < phi && phi < geom.muEnd;
}

Radians sum_covered_angle(const ExtendedRunInfo &egeom, int col, vector<UCP> &roots, double r, double z)
{
	// sum the angle covered by intersections inside of camera column cones
	// roots must be sorted and roots.back() must be the same as .front(), but with phi shifted by 2pi
	//
	Radians total = 0;
	V3 v;
	FOR(1, roots.size())
	{
		v = (roots[i-1] + roots[i]) / 2.;	// midpoint on the circle, roots are in galactocentric coords.
		v *= r/abs(v);
		if(roots[i].phi() - roots[i-1].phi() > ctn::pi/2.) { v *= -1; }
#if DBG
		Radians phi = v.phi();
		V3 v0(v);
#endif
		
		v.z = z; v.x += Rg;			// switch to geocentric coordinates (g.c. is to the right)

#if DBG
		V3 p(egeom.M*v);
		cout << "col = " << col << " v= " << v0 << " phi=" << deg(phi) << " mu=" << deg(p.phi()) << " nu=" << deg(ctn::pi/2 - p.theta()) << "\n";
#endif
		if(egeom.contains(v, col))
		{
			total += roots[i].phi() - roots[i-1].phi();
#if DBG
				double xc = Rg, z0 = z;
				V3 r0(roots[i-1].x + xc, roots[i-1].y, z0), r1(roots[i].x + xc, roots[i].y, z0);
				r0 = egeom.M * r0; r1 = egeom.M*r1;
				double del_phi = roots[i].phi() - roots[i-1].phi();
				cout << "    [" << deg(r0.phi()) << ", " << deg(ctn::pi/2 - r0.theta()) << "] - ";
				cout << "[" << deg(r1.phi()) << ", " << deg(ctn::pi/2 - r1.theta()) << "] - del_phi = ";
				cout << deg(del_phi) << " dV = " << del_phi*r << "\n";
#endif
		}
	}
#if DBG
	cout << "pix = " << total*r << "\n";
#endif
	return total;
}

inline double cube(const double x) { return x*x*x; }

class raytrace_fill : public FloodFill<I2>
{
public:
	double dx, dx2;
	vector<UCP> roots;
	bool edge;

	double dv, vtotal, vexpected;
	ExtendedRunInfo &egeom;
	XImage &img; int xc, yc;
	int col;

	struct Intersectors
	{
		auto_ptr<Intersector> is0, is1;
		auto_ptr<PlaneIntersector> isp0, isp1;
		auto_ptr<SphereIntersector> iss0, iss1;

		ExtendedRunInfo egeom;
		int col;
		double z;

		Intersectors(ExtendedRunInfo &egeom_, int col_, double z_)
			: egeom(egeom_), col(col_), z(z_)
		{
			Quadric q0, q1;
			intersect_cone(q0, sqr(egeom.t[col][0]), egeom.M, z);	// find the intersecting quadrics, planes and spheres
			intersect_cone(q1, sqr(egeom.t[col][1]), egeom.M, z);
			StandardQuadric sq0(q0), sq1(q1);
			is0.reset(new Intersector(sq0, Rg, 0));
			is1.reset(new Intersector(sq1, Rg, 0));

			isp0.reset(new PlaneIntersector(egeom.p0g, Rg, 0, z));		// intersect on run boundaries
			isp1.reset(new PlaneIntersector(egeom.p1g, Rg, 0, z));

			iss0.reset(new SphereIntersector(egeom.d0, Rg, z));
			iss1.reset(new SphereIntersector(egeom.d1, Rg, z));		// intersect on limiting spheres
		}

		int intersect(double r, vector<UCP> &roots)
		{
			roots.clear();

			is0->intersect(r, roots);
			is1->intersect(r, roots);
			isp0->intersect(r, roots);
			isp1->intersect(r, roots);
			iss0->intersect(r, roots);
			iss1->intersect(r, roots);

			if(roots.empty()) return 0;

			sort(roots.begin(), roots.end());
			roots.push_back(roots.front());
			roots.back().ph += ctn::pi * 2;

			return roots.size() - 1;
		}

		Radians find_arcs(vector<UCP> &roots, double r, bool &edge);
	};

	typedef map<double, Intersectors *> isectmap;
	isectmap isect;
	Intersectors *getIntersectors(double z);

	raytrace_fill(ExtendedRunInfo &egeom_, int col, double dx, XImage &img_);
	I2 initial_flood_point();

	virtual bool test(const I2 &v);
	virtual void action(const I2 &v);
	void printroots(double r, double z, int col);

	virtual ~raytrace_fill();
};

void raytrace_fill::printroots(double r, double z, int col)
{
	r = r + Rg;
	FOR(0, roots.size())
	{
		UCP &v = roots[i];
		V3 gal = V3(v.x, v.y, z) + V3(Rg, 0, 0);
		cout << "ucp=" << v << ", gal=" << gal << ", run=" << celestial(egeom.M*gal);
		if(i != 0) {
			cout << " dphi=" << deg(roots[i].phi() - roots[i-1].phi());

			UCP &r0 = roots[i-1];
			UCP &r1 = roots[i];
			V3 vv((r0 + r1) / 2.);	// midpoint on the circle, roots are in galactocentric coords.
			vv *= r/abs(vv);
			if(r1.phi() - r0.phi() > ctn::pi) { vv *= -1; }

			vv.z = z; vv.x += Rg;			// switch to geocentric coordinates (g.c. is to the right)
//			cout << "\n   vvv=" << vv - V3(Rg, 0, z);
			if(egeom.contains(vv, col))
			{
				cout << " +";
			}
		}
		cout << "\n";
	}
}

raytrace_fill::raytrace_fill(ExtendedRunInfo &egeom_, int col_, double dx_, XImage &img_)
	: egeom(egeom_), col(col_), vtotal(0), img(img_), dx(dx_), dx2(sqr(dx_))
{
	vexpected = egeom.mask.width(col)*egeom.geom.length()*(cube(egeom.d1)-cube(egeom.d0))/3.;

	xc = img.x() / 2;
	yc = img.y() / 2;
}

raytrace_fill::~raytrace_fill()
{
	FOREACH(isectmap::iterator, isect)
	{
		delete (*i).second;
	}
}

//#define DBG 1

Radians raytrace_fill::Intersectors::find_arcs(vector<UCP> &roots, double r, bool &edge)
{
	// sum the angle covered by intersections inside of camera column cones
	// roots must be sorted and roots.back() must be the same as .front(), but with phi shifted by 2pi
	//
	edge = false;
	Radians total = 0;
	V3 v;
	FOR(1, roots.size())
	{
		UCP &r0 = roots[i-1];
		UCP &r1 = roots[i];
		v = (r0 + r1) / 2.;	// midpoint on the circle, roots are in galactocentric coords.
		v *= r/abs(v);
		if(r1.phi() - r0.phi() > ctn::pi) { v *= -1; }
#if DBG
		Radians phi = v.phi();
		V3 v0(v);
#endif

		v.z = z; v.x += Rg;			// switch to geocentric coordinates (g.c. is to the right)

#if DBG
		V3 p(egeom.M*v);
		cout << "col = " << col << " v= " << v0 << " phi=" << deg(phi) << " mu=" << deg(p.phi()) << " nu=" << deg(ctn::pi/2 - p.theta()) << "\n";
#endif
		if(egeom.contains(v, col))
		{
			total += r1.phi() - r0.phi();
			if(r0.type != UCP::cone || r1.type != UCP::cone) { edge = true; }
#if DBG
				double xc = Rg, z0 = z;
				V3 r0(roots[i-1].x + xc, roots[i-1].y, z0), r1(roots[i].x + xc, roots[i].y, z0);
				r0 = egeom.M * r0; r1 = egeom.M*r1;
				double del_phi = roots[i].phi() - roots[i-1].phi();
				cout << "    [" << deg(r0.phi()) << ", " << deg(ctn::pi/2 - r0.theta()) << "] - ";
				cout << "[" << deg(r1.phi()) << ", " << deg(ctn::pi/2 - r1.theta()) << "] - del_phi = ";
				cout << deg(del_phi) << " dV = " << del_phi*r << "\n";
#endif
		}
	}
#if DBG
	cout << "pix = " << total*r << "\n";
#endif
	return total;
}

raytrace_fill::Intersectors *raytrace_fill::getIntersectors(double z)
{
	// find the right intersector
	Intersectors *in;
	map<double, Intersectors *>::iterator isit = isect.find(z);
	if(isit == isect.end())
	{
		isect[z] = (in = new Intersectors(egeom, col, z));
	} else {
		in = (*isit).second;
	}
	return in;
}

bool raytrace_fill::test(const I2 &vpix)
{
	double r = vpix.x * dx + Rg;
	double z = vpix.y * dx;

#if DBG2
	zdebug = z;
#endif

	dv = 0;
	if(r > 0)
	{
		Intersectors *in = getIntersectors(z);

		// test for intersection and calculate volume covered if there was one
		if(in->intersect(r, roots))
		{
			Radians arc = in->find_arcs(roots, r, edge);
			dv = arc*r*dx2;
		}
	}

	// if there was no intersection, or if the edge flag is set
	// this pixel is near the edge of the volume
	// subsample it heavily to accurately account for the edge
	if(0 && (edge || dv == 0))
	{
		const int sx = 5;
		const double sdx = dx / sx;
		const double sdx2 = sqr(sdx);
		double z0 = z - 0.5*dx + 0.5*sdx,	// center of the first subpixel
		       r0 = r - 0.5*dx + 0.5*sdx;

		dv = 0;
		FORj(nz, 0, sx)
		{
			double z = z0 + nz * sdx;
#if DBG2
			zdebug = z;
#endif
			Intersectors *in = getIntersectors(z);
			FORj(nr, 0, sx)
			{
				double r = r0 + nr * sdx;

				if(r <= 0) { continue; }
				if(in->intersect(r, roots))
				{
					Radians arc = in->find_arcs(roots, r, edge);
					dv += arc*r*sdx2;
				}
			}
		}
	}

	return dv != 0;
}

void raytrace_fill::action(const I2 &vpix)
{
//	cout << vpix.x + xc << " " << vpix.y + yc << " " << dv <<"\n";
	img(vpix.x + xc, vpix.y + yc) += dv;
	vtotal += dv;
}

I2 raytrace_fill::initial_flood_point()
{
	// find (mu,nu) of stripe center
	double mu = egeom.mask.start + egeom.mask.length()/2.;
	if(mu > ctn::pi2) { mu -= ctn::pi2; }
	double nu = egeom.mask.lo(col) + egeom.mask.width(col)/2;

//	cout << egeom << "\n";
//	cout << "IFP: mu=" << deg(mu) << ", nu=" << deg(nu) << "\n";

	V3 v;
	v.spherical((egeom.d0 + egeom.d1) / 2., ctn::pi/2 - nu, mu);
	v = egeom.Minv * v;	// convert to galactic

	v.x = Rg - v.x;		// convert to galactocentric
	v.y *= -1;

	I2 vpix(int((v.rho() - Rg) / dx), int(v.z / dx));	// nearest pixel
//	cout << v << " " << vpix << " [" << vpix.x + xc << " " << vpix.y + yc << "]\n";
	if(test(vpix)) { return vpix; }

	// handle very thin runs
	FORj(i, -1, 2) FORj(j, -1, 2)
	{
		I2 vpix2(vpix.x + i, vpix.y + j);
		if(test(vpix2)) { return vpix2; }
	}

	// something very wrong just happened
	ASSERT(0);
}

///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////
///////////////////////////////////////

struct volume_map
{
	typedef map<S2, interval_set, less_S2> pixelmap;
	pixelmap pixels;

	float dx;
	int run;

	double Dmin, Dmax;	/// volume limits
};

BOSTREAM(const volume_map &vm)
{
	out << vm.run << vm.dx << vm.Dmin << vm.Dmax << vm.pixels;

	return out;
};

BISTREAM(volume_map &vm)
{
	in >> vm.run >> vm.dx >> vm.Dmin >> vm.Dmax >> vm.pixels;
	
	return in;
}

void test_intervals()
{
	float A[] = { 1, 2,   3, 4,   5, 6,   7, 8 };
	const int N = sizeof(A) / sizeof(int);
	
	interval_set is, is0;
	is0.resize(N);
	copy(A, A+N, is0.begin());
	cout << "Original: " << is0 << "\n";

	is = is0;
	interval k = make_pair(1.5, 5.5); // (in, in)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(2.5, 5.5); // (out, in)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(2.5, 6.5); // (out, out)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(1.5, 6.5); // (in, out)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(1.5, 1.7); // (in, in, same)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	
	// degenerates (just touching the boundaries)
	cout << "\n";
	is = is0;
	k = make_pair(1, 5); // (in, in)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(2, 5); // (out, in)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(2, 6); // (out, out)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(1, 6); // (in, out)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	// sides
	cout << "\n";
	is = is0;
	k = make_pair(0, 6);
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(1, 9); // (in, out)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	is = is0;
	k = make_pair(0, 10); // (in, out)
	add_interval(is, k);
	cout << k << " : " << is << " : " << length(is) << "\n";

	// intersections
	cout << "Intersections\n";
	k = make_pair(1.5, 5.5); // (in, in)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	k = make_pair(2.5, 5.5); // (out, in)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	k = make_pair(2.5, 6.5); // (out, out)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	k = make_pair(1.5, 6.5); // (in, out)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	k = make_pair(1.5, 1.7); // (in, in, same)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));
	
	k = make_pair(10, 11); // (empty)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	k = make_pair(0, 11); // (all)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));
	
	// degenerates (just touching the boundaries)
	cout << "\n";
	is = is0;
	k = make_pair(1, 5); // (in, in)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	is = is0;
	k = make_pair(2, 5); // (out, in)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	is = is0;
	k = make_pair(2, 6); // (out, out)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	is = is0;
	k = make_pair(1, 6); // (in, out)
	intersect_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0) && length(is) <= length(k));

	// intersections
	cout << "Subtractions\n";
	k = make_pair(1.5, 5.5); // (in, in)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));

	k = make_pair(2.5, 5.5); // (out, in)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));

	k = make_pair(2.5, 6.5); // (out, out)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));

	k = make_pair(1.5, 6.5); // (in, out)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));

	k = make_pair(1.5, 1.7); // (in, in, same)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));
	
	k = make_pair(10, 11); // (empty)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));

	k = make_pair(0, 11); // (all)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));
	
	// degenerates (just touching the boundaries)
	cout << "\n";
	is = is0;
	k = make_pair(1, 5); // (in, in)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));

	is = is0;
	k = make_pair(2, 5); // (out, in)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));

	is = is0;
	k = make_pair(2, 6); // (out, out)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));

	is = is0;
	k = make_pair(1, 6); // (in, out)
	subtract_interval(is, is0, k);
	cout << k << " : " << is << " : " << length(is) << "\n";
	ASSERT(length(is) <= length(is0));
}

//
// Z-ray intersectors
//
class ZIntersector
{
public:
	double z0y2, z0xy, z0x2;
	double z1y1, z1x1;
	double z2;
public:
	ZIntersector() {}

	void setup(double t2, const M3 &M);
	bool intersect(double x, double y, double &z0, double &z1);
};

OSTREAM(const ZIntersector &is)
{
	out << "z0(x2, xy, y2) = " << is.z0x2 << " " << is.z0xy << " " << is.z0y2 << "\n";
	out << "z1(x1, y1) = " << is.z1x1 << " " << is.z1y1 << "\n";
	out << "z2 = " << is.z2 << "\n";
	return out;
}

void ZIntersector::setup(double t2, const M3 &M)
{
	const double a = M(0, 0), b = M(0, 1), c = M(0, 2);
	const double d = M(1, 0), e = M(1, 1), f = M(1, 2);
	const double g = M(2, 0), h = M(2, 1), i = M(2, 2);

	double a2 = sqr(M(0, 0)), b2 = sqr(M(0, 1)), c2 = sqr(M(0, 2));
	double d2 = sqr(M(1, 0)), e2 = sqr(M(1, 1)), f2 = sqr(M(1, 2));
	double g2 = sqr(M(2, 0)), h2 = sqr(M(2, 1)), i2 = sqr(M(2, 2));

#define Power(x,n) x##n

	// Formulas calculated in 'doc/Z ray.nb'

	z0y2 = Power(h,2) - (Power(b,2) + Power(e,2))*Power(t,2);
	z0xy = 2*g*h - 2*(a*b + d*e)*Power(t,2);
	z0x2 = Power(g,2) - (Power(a,2) + Power(d,2))*Power(t,2);

	z1y1 = 2*h*i - 2*(b*c + e*f)*Power(t,2);
	z1x1 = 2*g*i - 2*(a*c + d*f)*Power(t,2);

	z2 = Power(i,2) - (Power(c,2) + Power(f,2))*Power(t,2);

#undef Power
}

bool ZIntersector::intersect(double x, double y, double &za, double &zb)
{
	double z0 = z0y2*sqr(y) + z0xy*x*y + z0x2*sqr(x);
	double z1 = z1y1*y + z1x1*x;
	
	return gsl_poly_solve_quadratic(z2, z1, z0, &za, &zb) == 2;
}

inline bool z_ray_sphere_intersect(double x, double y, double d2, double &za, double &zb)
{
	double D = d2 - (sqr(x) + sqr(y));

	if(D < 0) return false;

	za = -sqrt(D);
	zb = -za;

	return true;
}

inline bool z_ray_plane_intersect(double x, double y, const Plane &p, double &z)
{
	if(p.c == 0) return false;
	z = - (p.a*x + p.b*y) / p.c;
	return true;
}



class zray_fill : public FloodFill<I2>
{
public:
	double dx, dx2;

	vector<pair<double, bool> > roots;

	double vtotal;
	interval_set dv;
	ExtendedRunInfo &egeom;
	volume_map &vm;
	int col;
	double d02, d12;

	ZIntersector in0, in1;

	zray_fill(ExtendedRunInfo &egeom_, int col, double dx, volume_map &vm_);
	bool initial_flood_point(I2 &ifp);

	void intersect(double x, double y);
	void addroot(double z, bool edge = false);
	interval_set z_covered(double x, double y, bool &edge);
	void printroots(double x, double y);

	double flood();
	virtual bool test(const I2 &v);
	virtual void action(const I2 &v);
};

double zray_fill::flood()
{
	vtotal = 0;

	I2 ifp;
	if(initial_flood_point(ifp))
	{
		FloodFill<I2>::flood(ifp);
	}

	return vtotal;
}

zray_fill::zray_fill(ExtendedRunInfo &egeom_, int col_, double dx_, volume_map &vm_)
	: egeom(egeom_), col(col_), vtotal(0), vm(vm_), dx(dx_), dx2(sqr(dx_)),
	d02(sqr(egeom_.d0)), d12(sqr(egeom_.d1))
{
	in0.setup(sqr(egeom.t[col][0]), egeom.M);
	in1.setup(sqr(egeom.t[col][1]), egeom.M);
}

void zray_fill::addroot(double z, bool edge)
{
	roots.push_back(make_pair(z, edge)); 
}

void zray_fill::intersect(double x, double y)
{
	roots.clear();

	double za, zb;
	if(in0.intersect(x, y, za, zb)) { addroot(za); addroot(zb); }
	if(in1.intersect(x, y, za, zb)) { addroot(za); addroot(zb); }
	if(z_ray_sphere_intersect(x, y, d02, za, zb)) { addroot(za, true); addroot(zb, true); }
	if(z_ray_sphere_intersect(x, y, d12, za, zb)) { addroot(za, true); addroot(zb, true); }
	if(z_ray_plane_intersect(x, y, egeom.p0g, za)) { addroot(za, true); }
	if(z_ray_plane_intersect(x, y, egeom.p1g, za)) { addroot(za, true); }
	
	sort(roots.begin(), roots.end());
}

interval_set zray_fill::z_covered(double x, double y, bool &edge)
{
	interval_set dz;
	edge = false;
	FOR(1, roots.size())
	{
		double z0 = roots[i-1].first;
		double z1 = roots[i].first;

		double z = (z0 + z1) / 2.;
//		cout << "Testing: " << celestial(egeom.M*V3(x, y, z)) << "\n";
		if(egeom.contains(V3(x, y, z), col))
		{
			add_interval(dz, make_pair(z0, z1));
			edge = edge || roots[i-1].second || roots[i].second;
		}
	}
	return dz;
}

void zray_fill::printroots(double x, double y)
{
	FOR(0, roots.size())
	{
		cout << "(x,y,z) = " << x << ", " << y << ", " << roots[i].first;
		cout << " " << celestial(egeom.M * V3(x, y, roots[i].first)) << " edge=" << roots[i].second;
		if(i != 0) { cout << " dz=" << roots[i].first - roots[i-1].first; }
		cout << "\n";
	}
}

bool zray_fill::test(const I2 &vpix)
{
	double x = vpix.x * dx;
	double y = vpix.y * dx;

#define PIX_ORIG_CORNER
#ifdef PIX_ORIG_CORNER
	x += dx / 2.;
	y += dx / 2.;
#endif

	intersect(x, y);

	bool edge;
	dv = z_covered(x, y, edge);

	return dv.size();
}

void zray_fill::action(const I2 &vpix)
{
//	cout << vpix.x - xc << " " << vpix.y + yc << " " << dv <<"\n";
	// Note the -vpix.x and not +vpix.x - we invert this to have the galactic center to the left
	// and at the same time flip .y to keep the coordinatem righthanded
	// galactic longitude is still ccw, but starts at -x axis
#ifdef PIX_ORIG_CORNER
	// take into account that the coordinates are the coordinates of
	// bottom left corner of the pixel, but the volume of the pixel
	// is calculated in the center of the pixel
	add_interval(vm.pixels[S2(-vpix.x - 1, -vpix.y - 1)], dv);
#else
	add_interval(vm.pixels[S2(-vpix.x, -vpix.y)], dv);
#endif
	vtotal += length(dv) * dx2;
}

bool zray_fill::initial_flood_point(I2 &vpix)
{

	// find (mu,nu) of stripe center
	double mu = egeom.mask.start + egeom.mask.length()/2.;
	if(mu > ctn::pi2) { mu -= ctn::pi2; }
	double nu = egeom.mask.lo(col) + egeom.mask.width(col)/2;

	V3 v;
	v.spherical((egeom.d0 + egeom.d1) / 2., ctn::pi/2 - nu, mu);
	v = egeom.Minv * v;	// convert to galactic

	vpix = I2(int(v.x / dx), int(v.y / dx));	// nearest pixel
	if(test(vpix)) { return true; }

	// handle very thin runs
	FORj(i, -1, 2) FORj(j, -1, 2)
	{
		I2 vpix2(vpix.x + i, vpix.y + j);
		if(test(vpix2)) { vpix = vpix2; return true; }
	}

	// something very wrong just happened
	ASSERT(0);
}

int zrayng(const set<int> &runs, float dx, pair<float, float> r, pair<float, float> ri)
{
	double d0, d1;
	plx_gri_locus paralax;
	paralax.distance_limits(d0, d1, ri.first, ri.second, r.first, r.second);
	cout << "Distance limits: " << d0 << " " << d1 << "\n";
	ASSERT(d0 < d1);

 	ExtendedRunInfo::geomDB.reset(new RunGeometryDB);

	FOREACH(set<int>::const_iterator, runs)
	{
		int run = *i;

		volume_map vm;
		vm.dx = dx;
		vm.run = run;
		vm.Dmin = d0;
		vm.Dmax = d1;

		ExtendedRunInfo egeom(run, d0, d1);

		cout << "# run " << run << "\n";
		cout << "# col mapped predicted difference\n";
		double vtotexp = 0, vtot = 0;
		FORj(col, 0, 6)
		{
			zray_fill zf(egeom, col, dx, vm);
			zf.flood();
	
			double vexpected = egeom.mask.width(col)*egeom.geom.length()*(cube(egeom.d1)-cube(egeom.d0))/3.;
			vtotexp += vexpected;
			vtot += zf.vtotal;

			double pctdiff = 100 * (zf.vtotal / vexpected - 1);
			cout << "\t" << col << " " << zf.vtotal << " " << vexpected << " " << pctdiff << "%\n";
		}
		double pctdiff = 100 * (vtot / vtotexp - 1);
		cout << "\t" << "=" << " " << vtot << " " << vtotexp << " " << pctdiff << "%\n";

		string fn = io::format("volumes/vol.%05d.%5.3f-%5.3f.%5.3f-%5.3f.bin")
			<< run << r.first << r.second << ri.first << ri.second;
		binary_output_or_die(out, fn);
		out << vm;
	}
}

static bool bin3d = true;

// Dec10 debugging
double gx, gy, gz = 0, gdx, gn;

/*
	A "safe" version of floor() function, which takes into account the fact that
	machine representation of floating point numbers is not exact.

	If <x> is withing <eps> of ceil(x), <x> is assumed to be equal to ceil(x).
	
	Eg.: eps_floor(0.999999, 1e-4) == 1
*/
template<typename T>
T eps_floor(T x, T eps)
{
	T c = ceil(x);
	if(c - x < eps) { return c; }
	return c - T(1);
}

inline void store_dv(binned_run &br, const S3 &k, double dv, double dz)
{
	if(bin3d)
	{
		br.pixels[k].volume += dv;
		
		// Dec10 debugging
		double dik = ((gx - gdx/2.) / dz + .5);
		double djk = ((gy - gdx/2.) / dz + .5);
		double dkk = ((gz - gdx/2.) / dz + .5);

/*		short ik = (short)floor((gx - gdx/2.) / dz + .5);
		short jk = (short)floor((gy - gdx/2.) / dz + .5);
		short kk = (short)floor((gz - gdx/2.) / dz + .5);*/

		short ik = (short)eps_floor(dik, 1./(4.*gn));
		short jk = (short)eps_floor(djk, 1./(4.*gn));
		short kk = (short)eps_floor(dkk, 1./(4.*gn));
/*		short ik = (short)floor(dik);
		short jk = (short)floor(djk);
		short kk = (short)floor(dkk);*/
		if(ik != k.x || jk != k.y) {
		cerr << gx << " " << gy << " : " << k << " " << ik << " " << jk << " " << kk << " : ";
		cerr << setprecision(12) << dik << " " << setprecision(12) << djk << "\n";
		}
	}
	else
	{
#if 0
		// map pixel to meridianal plane
		static const double rg = Rg / dz;
		double	x = k.x + rg,
			y = k.y;
		double rho = sqrt(sqr(x) + sqr(y));
		// simple nearest-grid-point binning (note, the
		// coordinates are the coordinates of pixel centers)
		S3 rz((short)floor(rho + 0.5), k.z, 0);
		br.pixels[rz].volume += dv;
#else
		S3 rz(k.x, k.z, 0);	// because rho is stored in k.x
		br.pixels[rz].volume += dv;
#endif
	}
}

double pixelate_interval(binned_run &br, S2 idx, const interval_set &is0, float dz, float dv,
	double x, double y, double d02, double d12)
{
	// cut boundaries of the volume
	double za, zb;
	if(!z_ray_sphere_intersect(x, y, d12, za, zb)) return 0;	// outer sphere

	interval_set is;
	intersect_interval(is, is0, make_pair(float(za), float(zb)));

	if(!is.size()) return 0;

	// inner sphere
	if(z_ray_sphere_intersect(x, y, d02, za, zb))		// inner sphere
	{
		interval_set tmp;
		subtract_interval(tmp, is, make_pair(float(za), float(zb)));
		is = tmp;
		if(!is.size()) return 0;
	}

	gx = x; gy = y; // debugging

	double sum = 0;
	S3 k(idx);
	

	FOREACH(interval_set::const_iterator, is)
	{
		float x0 = *i / dz + 0.5; i++;
		float x1 = *i / dz + 0.5;

		float i0 = ceil(x0);
		float i1 = floor(x1);

		// last pixel
		k.z = (short)i1;
		const double ddv = (x1 - max(i1, x0)) * dv;
		store_dv(br, k, ddv, dz);
		sum += ddv;

		if(i0 <= i1)
		{
			// first pixel
			k.z = (short)i0 - 1;
			const double ddv = (i0 - x0) * dv;
			store_dv(br, k, ddv, dz);
			sum += ddv;

			// pixels in between
			for(k.z = (short)i0; k.z != i1; k.z++)
			{
				store_dv(br, k, dv, dz);
				sum += dv;
			}
		}
	}

#if 0 	
	// debugging code
	int bx = -1, by = -13;
	if(k.x == bx && k.y == by)
	{
		cerr << x << " " << y << " : " << " vol=" << sum << " " << is << "\n";
		FOREACH(br.pixels)
		{
			const S3 &k = (*i).first;
			binned_run::pixel &p = (*i).second;
			if(k.x == bx && k.y == by) { cerr << k << " " << p.volume << "\n"; }
		}
		exit(-1);
	}
#endif
	
	return sum;
}

double pixelate_volume(binned_run &br, const volume_map &vm, int n)
{
#ifdef PIX_ORIG_CORNER
	ASSERT(n % 2 == 0);
	int nhalf = n / 2;
#else
	ASSERT(n % 2 == 1);
#endif
	float dx = vm.dx * n;
	gdx = vm.dx;	// debugging
	gn = n;

	br.dx = dx;
	br.run = vm.run;
	
	double d02 = sqr(br.Dmin), d12 = sqr(br.Dmax);
	ASSERT(br.Dmin >= vm.Dmin);
	ASSERT(br.Dmax <= vm.Dmax);
	
	br.pixels.clear();

	float dv = sqr(vm.dx) * dx;

	double vvm = 0, vbr = 0;
	FOREACH(volume_map::pixelmap::const_iterator, vm.pixels)
	{
		S2 idx = (*i).first;
		const interval_set &is = (*i).second;

		// rescale idx
#ifdef PIX_ORIG_CORNER
		// physical coordinates of pixel center
		double x = (idx.x + 0.5) * vm.dx;
		double y = (idx.y + 0.5) * vm.dx;

		if(bin3d)
		{
			idx.x = (short)floor(float(idx.x + nhalf) / n);
			idx.y = (short)floor(float(idx.y + nhalf) / n);
		} else {
			double rho = sqrt(sqr(x + Rg) + sqr(y));
			// simple nearest-grid-point binning (note, the
			// coordinates are the coordinates of pixel centers)
			idx.x = (short)floor(rho / dx + 0.5);
		}

		vbr += pixelate_interval(br, idx, is, dx, dv, x, y, d02, d12);
#else
		ASSERT(0);
		// unimplemented - since n is odd, that has to be taken into account
		// (or not?)
#endif

		vvm += length(is) * sqr(vm.dx);
	}

	br.loaded = true;
	return vbr;
}

void rhoray(int argc, char **argv)
{
	VERSION_DATETIME(version);

	Options opts(
		"\
	Program for 'flattening' volume maps of SDSS runs. Used for mstars project. \n\
	Takes into consideration limiting apparent and absolute magnitudes and produces \n\
	a volume map which equally covers all objects in a given ri magnitude bin. \n\
	The program uses the output of volume.x",
		version,
		Authorship::majuric
	);

	opts.argument("outputFITSFile", "Filename of the output FITS file");
	opts.argument("run", "Run number or a file which contains runs to flatten");
	opts.argument("x", "Width of the image [pixels]");
	opts.argument("y", "Height of the image [pixels]");
	opts.argument("dx", "Scale [pc/pixel]");
	opts.argument("d0", "Minimum geocentric distance [pc]");
	opts.argument("d1", "Maximum geocentric distance [pc]");

	try {
		opts.parse(argc, argv);
	} catch(EOptions &e) {
		cout << opts.usage(argv);
		e.print();
		exit(-1);
	}

	/////// Now start with the real busyness //////

	// parse the command line
	string outputFile = opts["outputFITSFile"];
	set<int> runs;
	int run = opts["run"];
	if(run == 0)
	{
		text_input_or_die (in, (std::string)opts["run"]);
		load(in, runs, 0);
	} else {
		runs.insert(run);
	}

	double dx = opts["dx"];
	double d0 = opts["d0"];
	double d1 = opts["d1"]; ASSERT(d0 < d1);
	int x = opts["x"];
	int y = opts["y"];

	XImage img(x, y);

	gsl_set_error_handler_off ();

	cout << "# col mapped predicted difference\n";
	FOREACH(set<int>::iterator, runs)
	{
		int run = *i;
		ExtendedRunInfo egeom(run, d0, d1);

		cout << "# run " << run << "\n";
		FORj(col, 0, 6)
		{
			raytrace_fill rf(egeom, col, dx, img);
			rf.flood(rf.initial_flood_point());

			double pctdiff = 100 * (rf.vtotal / rf.vexpected - 1);
			cout << "\t" << col << " " << rf.vtotal << " " << rf.vexpected << " " << pctdiff << "%";
			if(abs(pctdiff) > 0.1) { cout << " WARNING - LARGE DIFFERENCE!"; }
			cout << "\n";
		}
#if 0
		FOREACHj(k, img)
		{
			if(run == 2650 || run == 2766) break;
			if(*k == 100000)
			{
				cout << "run === " << run << " " << k.x << " " << k.y << "\n";
				ASSERT(0);
			}
		}
#endif
	}

#ifdef HAVE_LIBCCFITS
	io::fits::write(img, outputFile);
#else
	ASSERT(0) { cerr << "libCCfits support not compiled in.\n"; }
#endif
}

int bin_volumes(const std::set<int> &runs, double dx, int ndx, pair<float, float> r, pair<float, float> ri)
{
	double Dmin, Dmax;
	plx_gri_locus paralax;
	paralax.distance_limits(Dmin, Dmax, ri.first, ri.second, r.first, r.second);
	cout << "Distance limits: " << Dmin << " " << Dmax << "\n";

	ExtendedRunInfo::geomDB.reset(new RunGeometryDB);

	FOREACH(std::set<int>::const_iterator, runs)
	{
		int run = *i;
		cout << run << " ... "; cout.flush();

		volume_map vm;
		binned_run br;

		br.Dmin = Dmin;
		br.Dmax = Dmax;

//		string fn = io::format("volumes/vol.%05d.%5.3f-%5.3f.%5.3f-%5.3f.bin")
//			<< run << r.first << r.second << ri.first << ri.second;
		string fn = io::format("volumes/vol.%05d.14.000-21.800.0.100-0.150.bin")
			<< run;
//fn = "merged.bin";
		gz_binary_input_or_die(in, fn);
		in >> vm;

		// adjust the dx and D in volume
		double fact = dx / vm.dx;
		vm.Dmax *= fact;
		vm.Dmin *= fact;
		vm.dx = dx;
		{
		FOREACH(volume_map::pixelmap::iterator, vm.pixels) {
			FOREACHj(interval_set::iterator, j, (*i).second) {
				*j *= fact;
			}
		}
		}

		if(i == runs.begin())
		{
			cout << "[Available volume: " << vm.Dmin << " " << vm.Dmax << "] "; cout.flush();
		}

		// pixelize volume
		double vbr = pixelate_volume(br, vm, ndx);

		// sanity check
		ExtendedRunInfo egeom(run, Dmin, Dmax);
		double vexpected = 0;
		FORj(col, 0, 6)
		{
			vexpected += egeom.mask.width(col)*egeom.geom.length()*(cube(egeom.d1)-cube(egeom.d0))/3.;
		}
		double pctdiff = 100 * (vbr/vexpected - 1);
		cout << vbr << " " << vexpected << " " << pctdiff << "%\n";
		if(abs(pctdiff) >= 0.2) { cout << "	WARNING: large pctdiff!\n"; }

 		fn = io::format("maps/map.%s%05d.%5.3f-%5.3f.%5.3f-%5.3f.bin")
 			<< (bin3d ? "" : "rz.") << run << r.first << r.second << ri.first << ri.second;
 		binary_output_or_die(out, fn);
 		out << br;

output_or_die(tout, "vol.txt");
tout << br;
	}
	return 0;
}


/* 
      ///////////////////////////////////////////////////////////////////////////////
	Merge a set of interval volume maps for each run, into a merged volume map
	which traces uniquely covered volume
      ///////////////////////////////////////////////////////////////////////////////
*/

struct volstream
{
	ibinarystream in;

	size_t nleft;

	S2 k;
	interval_set is;

	volstream(const char *name) : in(*new gzstream::ifstream(name)), nleft(0) {}
	~volstream() { delete (gzstream::ifstream *)&in.f; }
};

struct priority_volstream : public less_S2 {
	bool operator()(const volstream *v1, const volstream *v2) const
	{
		return !less_S2::operator ( )(v1->k, v2->k);
	}
};

class unistream
{
protected:
	priority_queue<volstream *, vector<volstream *>, priority_volstream> q;
public:
	volume_map vm;
	size_t size;
	void add(const std::string &fn)
	{
		auto_ptr<volstream> vs(new volstream((fn+".gz").c_str()));

		vs->in  >> vm.run >> vm.dx >> vm.Dmin >> vm.Dmax
			>> vs->nleft;
		ASSERT(vs->nleft);
		cout << "\t" << vs->nleft << "\n"; cout.flush();

		vs->in	>> vs->k >> vs->is;
		size += vs->nleft;

		q.push(vs.release());
	}

	bool next(S2 &k, interval_set &is)
	{
		if(q.empty()) return false;

		volstream *vs = q.top();
		q.pop();

		k = vs->k;
		is = vs->is;

		--vs->nleft;
		if(vs->nleft)
		{
			vs->in >> vs->k >> vs->is;
			q.push(vs);
		} else {
			delete vs;
		}

		return true;
	}

	unistream() : size(0) {}

	~unistream()
	{
		while(!q.empty()) { delete q.top(); q.pop(); }
	}
};

int merge_volumes(const std::string &outfn, const std::set<int> &runs)
{
	volume_map merged;
	unistream ss;
	binary_output_or_die(out, outfn);

	int n = 1;
	FOREACH(std::set<int>::const_iterator, runs)
	{
		int run = *i;
		cout << setw(4) << n << " : " << run << "\n"; cout.flush();

		string fn = io::format("volumes/vol.%05d.14.000-21.800.0.100-0.150.bin")
			<< run;

		ss.add(fn);
		n++;
	}

	// write out volume header
	ss.vm.run = 0;
	out << ss.vm.run << ss.vm.dx << ss.vm.Dmin << ss.vm.Dmax
	    << ss.size;

	cerr << "Total intervals to process: " << ss.size << "\n";
	ss.size = 0;

 	// merge intervals and write them out
	double volume = 0, vtot = 0;
	interval_set is, tmp;
	S2 k, kprev;
	bool first = true;
	ticker tick(10000);
	while(ss.next(k, tmp))
	{
		volume += length(tmp);

		if(first)	// bootstrap
		{
			first = false;
			kprev = k;
		}
		
		if(k == kprev)
		{
			add_interval(is, tmp);
		}
		else
		{
			out << kprev << is;
			ss.size++;
			vtot += length(is);

			kprev = k;
			is = tmp;
		}
		tick.tick();
	}

	// write correct header
	out_stream.seekp(0);
	out << ss.vm.run << ss.vm.dx << ss.vm.Dmin << ss.vm.Dmax
	    << ss.size;

	double pctdiff = 100 * (vtot/volume - 1);
	cout << vtot << " " << volume << " " << pctdiff << "%\n";

	return 0;
}

/*


*/
#if 1
class pixelize_volume
{
public:
	unistream ss;

	// current interval and its index (in dx units, the coordinate of _bottom_left_ of the pixel)
	S2 idx; interval_set is;

	// current pixel volume, linear dimension
	float dv, dx;
	
	// physical coordinates of the _center_ of current pixel
	double x, y, z;

public:
	double dxfactor;
	bool next()
	{
		if(!ss.next(idx, is)) return false;
		FOREACH(interval_set::iterator, is) { *i *= dxfactor; }
		return true;
	}

	void init(const std::string &volfn, double dx)
	{
		ss.add(volfn);

		// adjust the dx and D in volume
		// it's basically a bugfix (intervals should have been dimensionless)...
		volume_map &vm = ss.vm;
		dxfactor = dx / vm.dx;
		vm.Dmax *= dxfactor;
		vm.Dmin *= dxfactor;
		this->dx = vm.dx = dx;
		cout << "[Available volume: " << vm.Dmin << " " << vm.Dmax << "] "; cout.flush();

		dv = cube(vm.dx);
	}

	virtual bool prepareInterval() = 0;
	virtual double action(float ddv) = 0;

	/*
		Bins the interval into pixels. For each pixel, sets (x, y, z) to coordinates
		of its _center_ and calls action(ddv), where ddv is the volume covered within
		that pixel.
	*/
	virtual double pixelate_interval()
	{
		if(!prepareInterval()) return 0;

		double sum = 0;
		FOREACH(interval_set::iterator, is)
		{
			float x0 = *i / dx; i++;
			float x1 = *i / dx;

			float i0 = ceil(x0);
			float i1 = floor(x1);

			// last pixel
			z = (i1  + 0.5)*dx;
			const double ddv = (x1 - max(i1, x0)) * dv;
			sum += action(ddv);

			if(i0 <= i1)
			{
				// first pixel
				z = ((i0 - 1) + 0.5)*dx;
				const double ddv = (i0 - x0) * dv;
				sum += action(ddv);
	
				// pixels in between
				int kz;
				for(kz = (short)i0; kz != i1; kz++)
				{
					z = (kz + 0.5)*dx;
					sum += action(dv);
				}
			}
		}

		return sum;
	}

	double run()
	{

		double vvm = 0, vbr = 0;
		ticker tick(10000);
		while(next())
		{
			// set (x, y) physical coordinates of pixel center
			x = (idx.x + 0.5) * dx;
			y = (idx.y + 0.5) * dx;

			vbr += pixelate_interval();
			vvm += length(is) * sqr(dx);
			
			tick.tick();
		}

		return vbr;
	}
};

class bin3d_pixelize_volume : public pixelize_volume
{
public:
	double d02, d12; // square of inner and outer boundary of the runs
	binned_run br;
	float ndx;
	double flooreps;

	bin3d_pixelize_volume(const std::string &volfn, float dx, int n, pair<float, float> r, pair<float, float> ri)
	{
		ASSERT(n % 2 == 0);

		init(volfn, dx);

		flooreps = 1./(4.*n);

		ndx = dx * n;
		br.dx = ndx;
		br.run = ss.vm.run;

		double Dmin, Dmax;
		plx_gri_locus paralax;
		paralax.distance_limits(Dmin, Dmax, ri.first, ri.second, r.first, r.second);
		cout << "Distance limits: " << Dmin << " " << Dmax << "\n";
		br.Dmin = Dmin;
		br.Dmax = Dmax;
		ASSERT(br.Dmin >= ss.vm.Dmin);
		ASSERT(br.Dmax <= ss.vm.Dmax);

		d02 = sqr(br.Dmin), d12 = sqr(br.Dmax);
		br.pixels.clear();
	}

	virtual bool prepareInterval()
	{
		// cut boundaries of the volume
		double za, zb;
		if(!z_ray_sphere_intersect(x, y, d12, za, zb)) return 0;	// outer sphere

		interval_set is0 = is;
		intersect_interval(is, is0, make_pair(float(za), float(zb)));
	
		if(!is.size()) return false;
	
		// inner sphere
		if(z_ray_sphere_intersect(x, y, d02, za, zb))		// inner sphere
		{
			interval_set tmp;
			subtract_interval(tmp, is, make_pair(float(za), float(zb)));
			is = tmp;
			if(!is.size()) return false;
		}

		return true;
	}
#if 0
	// Dec10 debugging
	double pixelate_interval()
	{
		double sum = pixelize_volume::pixelate_interval();

		short i = (short)eps_floor((x - dx/2.) / ndx + .5, flooreps);
		short j = (short)eps_floor((y - dx/2.) / ndx + .5, flooreps);
		int bx = -1, by = -13;
		if(i == bx && j == by)
		{
			cerr << x << " " << y << " : " << " vol=" << sum << " " << is << "\n";
			FOREACH(br.pixels)
			{
				const S3 &k = (*i).first;
				binned_run::pixel &p = (*i).second;
				if(k.x == bx && k.y == by) { cerr << k << " " << p.volume << "\n"; }
			}
			exit(-1);
		}
	}
#endif
	S3 binnedPixelCoords(double x, double y, double z)
	{
/*		short i = (short)floor((x - dx/2.) / ndx + .5);
		short j = (short)floor((y - dx/2.) / ndx + .5);
		short k = (short)floor((z - dx/2.) / ndx + .5);*/
		short i = (short)eps_floor((x - dx/2.) / ndx + .5, flooreps);
		short j = (short)eps_floor((y - dx/2.) / ndx + .5, flooreps);
		short k = (short)eps_floor((z - dx/2.) / ndx + .5, flooreps);
		return S3(i, j, k);
	}

	virtual double action(float ddv)
	{
		br.pixels[binnedPixelCoords(x, y, z)].volume += ddv;
		return ddv;
	}
};

//
//	Bins a merged volume, stored in file <volfn>, taking the pixel size to be <dx>
//	into pixels of <dx>*<ndx>, and stores the result into <outfile>
//
int bin_volumes2(const std::string &outfile, const std::string &volfn, double dx, int ndx, pair<float, float> r, pair<float, float> ri)
{
	bin3d_pixelize_volume binner(volfn, dx, ndx, r, ri);
	cout << binner.run() << "\n";
	binner.br.loaded = true;

	binary_output_or_die(out, outfile);
	out << binner.br;

	output_or_die(tout, "vol.txt");
	tout << binner.br;

	return 0;
}

//
//	Traverses <uniqMapFn> which contains a merged unique volume as created
//	by bin_volumes2() and adds the unique volume information to each pixel
//	of <brs>
//
void mergeUniqVolume(binned_runset &brs, const std::string &uniqMapFn, pair<float, float> r, pair<float, float> ri)
{
	// load uniq volume
	binned_run br;
	binary_input_or_die(in, uniqMapFn);
	in >> br;
	
	// Calculate the solid angle covered by uniq volume
	// this is mostly for sanity check/debugging purposes
	double Dmin, Dmax;
	plx_gri_locus paralax;
	paralax.distance_limits(Dmin, Dmax, ri.first, ri.second, r.first, r.second);
	
	double vsouth = 0, vnorth = 0;
	FOREACH(binned_run::pixelmap::iterator, br.pixels)
	{
		const S3 &v = (*i).first;
		binned_run::pixel &p = (*i).second;
		if(v.z >= 0) { vnorth += p.volume; }
		else { vsouth += p.volume; }
	}
	double volume = vnorth + vsouth;
	double Omega = volume / (1./3.*(cube(Dmax) - cube(Dmin)));
	Omega /= sqr(ctn::d2r);
	double OmegaSouth = (vsouth / (1./3.*(cube(Dmax) - cube(Dmin)))) / sqr(ctn::d2r);
	double OmegaNorth = (vnorth / (1./3.*(cube(Dmax) - cube(Dmin)))) / sqr(ctn::d2r);
	cerr << "Total volume       : " << volume << "\n";
	cerr << "Distance boundaries: " << Dmin << " " << Dmax << "\n";
	cerr << "Solid angle covered: " << Omega << "\n";
	cerr << "Solid angle north  : " << OmegaSouth << "\n";
	cerr << "Solid angle south  : " << OmegaNorth << "\n";

	FOREACH(binned_runset::pixelmap::iterator, brs.pixels)
	{
		const S3 &k = (*i).first;
		binned_runset::pixel &p = (*i).second;

		ASSERT(br.pixels.count(k))
		{
			cerr << k << " does not exist in unique map\n";
			cerr << "Pixel details: N=" << p.uniqueN << ", V=" << p.volume << "\n";
		}

		p.uniqueVolume = br.pixels[k].volume;

		double pct;
		if(p.uniqueVolume > p.volume && abs(pct = p.uniqueVolume / p.volume - 1) > 1e-3)
//		if(p.uniqueVolume != p.volume && abs(pct = p.uniqueVolume / p.volume - 1) > 1e-3)
		{
			cerr << k << " " << p.uniqueVolume << " > " << p.volume << "(" << pct*100 <<")\n";
		}
	}
	
	// Dec 8th debugging: dump unique map
	text_output_or_die(to, "uniqmap.txt");
	to << br;
}

class gc_pixelize_volume : public bin3d_pixelize_volume
{
public:
	Radians nu, tnu;
	Radians node, inc;
	ZIntersector in;
	M3 M;

	gc_pixelize_volume(const std::string &volfn, float dx, int n, pair<float, float> r, pair<float, float> ri,
		Radians node_, Radians inc_, Radians nu_)
		: bin3d_pixelize_volume(volfn, dx, n, r, ri), nu(nu_), node(node_), inc(inc_)
		{
			// setup intersectors
			rotation_matrix(M, node, inc);
			tnu = tan(nu);
			in.setup(sqr(tnu), M);

			cerr << "node = " << deg(node) << "\n";
			cerr << "inc  = " << deg(inc) << "\n";
			cerr << "nu   = [" << deg(nu) << ", " << deg(-nu) << ")\n";
		}
	
	//
	// v is in galactic cartesian coordinates
	//
	bool contains(const V3 &x) const
	{
		// test distance limits
//		double d2 = dot(x, x);
//		if(d2 < d02 || d2 > d12) return false;
	
		V3 v(M*x);

		// test if it's within opening angle
		double t = v.z / v.rho();
		return abs(t) < tnu;
	}

#if 1
	interval_set is0;
	virtual bool prepareInterval()
	{
		if(!bin3d_pixelize_volume::prepareInterval()) { return false; }

		// coordinates are stored with gc in -x direction
		double x = -this->x;
		double y = -this->y;
//		x = 10; y = 10;

		double za, zb;
		// not intersecting the cone at all - everything on this
		// ray is within the cone
		if(!in.intersect(x, y, za, zb)) return true;

		ASSERT(za <= zb);
		
		double z = (za + zb) / 2.;
		is0 = is;
		if(!contains(V3(x, y, z)))
		{
			// za and zb are intersecting the same half-cone - things
			// with z < za and z > zb are within the cone
			subtract_interval(is, is0, make_pair(float(za), float(zb)));
		} else {
			// za and zb are intersecting different half-cones
			// they define the interval on the ray which is within the cone
			intersect_interval(is, is0, make_pair(float(za), float(zb)));
		}

		return is.size() != 0;
	}
#endif
	
	virtual double action(float ddv)
	{
		ASSERT(is.size());
		
		// transform to GC coordinate system (note: GC is in -x direction)
		V3 v(-x, -y, z);
		Radians lon = v.phi(), lat = ctn::pi/2. - v.theta();
		coordinates::galgcs(node, inc, lon, lat, lon, lat);
		if(abs(lat) > nu) {
			if(abs(lat)/nu > 1.05) {
		 		cerr << deg(lat) << "[x=" << x << " y=" << y << "]\n";
			}
			return 0;
		}

		// bin
		v.celestial(abs(v), lon, lat);
		br.pixels[binnedPixelCoords(v.x, v.y, v.z)].volume += ddv;
		//cerr << "+";
		return ddv;
	}
};

//
//
int bin_gc_cut(
	const std::string &outfile, const std::string &volfn, 
	double dx, int ndx,
	pair<float, float> r, pair<float, float> ri,
	Radians node, Radians inc, pair<Radians, Radians> nu)
{
	gc_pixelize_volume binner(volfn, dx, ndx, r, ri, node, inc, nu.second);
	cout << binner.run() << "\n";
	binner.br.loaded = true;

	binary_output_or_die(out, outfile);
	out << binner.br;

	output_or_die(tout, "gccut.txt");
	tout << binner.br;

	return 0;
}

void plane_transformer::setup(const V3 &x1, const V3 &x2, const V3 &x3, const V3 &origin, double delta_, bool earth_on_x_axis)
{
	delta = delta_;

	x[0] = x1; x[1] = x2; x[2] = x3;

	// plane normal and d
	V3 n = cross(x2 - x1, x3 - x1);
	n /= abs(n);		// normal
	double d = dot(n,origin);	// distance from the origin (fourth plane param)
	p  = Plane(n, -d);	// the plane
	
	cerr << "\tPlane = " << p.a << " " << p.b << " " << p.c << " " << p.w << "\n";

	// setup coordinate system and transformation matrix
	V3 a0(-8000, 0, 0);	// galactic center
	V3 a1(0, 0, 0);		// earth

	// basis vectors
	V3 e3 = n;		// z axis
	V3 e1;
		if(earth_on_x_axis)
		{
			e1 = (a1 - a0);
			e1 = e1 - n*dot(n,e1);	// x axis (projection of a1-a0 to p plane)
		} else {
			e1 = x2 - x1;
			cerr << "Earth NOT on x axis.";
		}
		e1 /= abs(e1);
	V3 e2 = cross(e3, e1);	// z cross x = y for righthanded c.s.
	matrix_transpose(M, M3(e1, e2, e3));	// galactic -> plane CS transform. matrix

	cerr << "Basis vectors: " << e1 << e2 << e3 << "\n";	
	cerr << "\tMatrix =\n" << M << " [earthcentric -> plane]\n";

	// zero point of the coordinate system, in earthcentric coordinates
	//
	// CAVEAT: origin can move the origin of the z axis, to above or below the plane,
	// if origin point itself is not in the plane.
	//
	t0 = origin;
	cerr << "\tOrigin = " << t0 << " [earthcentric]\n";
}

V3 plane_transformer::toPlane(const V3 &v)
{
	V3 tmp = v - t0;
	tmp = M*tmp;
	return tmp;
}

class plane_pixelize_volume : public bin3d_pixelize_volume
{
public:
	plane_transformer pt;
public:
	plane_pixelize_volume(const std::string &volfn, float dx, int n, pair<float, float> r, pair<float, float> ri)
		: bin3d_pixelize_volume(volfn, dx, n, r, ri)
		{}

	virtual double action(float ddv)
	{
		// transform to plane coordinate system
//		x = 0; y = -1000; z = 0;
//		cerr << "v0 = " << V3(x, y, z) << "\n";
		V3 v = pt.toPlane(V3(x, y, z));
//		cerr << "vt = " << v << "\n";
//		exit(-1);
		if(abs(v.z) > pt.delta) return 0;

		br.pixels[binnedPixelCoords(v.x, v.y, v.z)].volume += ddv;
		return ddv;
	}
};

V3 equecen(const double d, const Radians ra, const Radians dec)
{
	Radians l, b;
	coordinates::equgal(ra, dec, l, b);
//l = ra; b = dec;
	V3 v; v.celestial(d, l, b);
	
	// convert standard cartesian galactic to earthcentric galactic coordinates
	v.x = -v.x;
	v.y = -v.y;
	
	return v;
}

void ecenequ(const V3 &v_, double &d, Radians ra, Radians dec)
{
	V3 v(-v_.x, -v_.y, v_.z);

	d = abs(v);

	coordinates::galequ(v.phi(), v.lat(), ra, dec);
}

void setupPlaneTransformer(
	plane_transformer &pt,
	const std::string &coordsys,
	double d1, pair<double, double> p1,
	double d2, pair<double, double> p2,
	double d3, pair<double, double> p3,
	double d0, pair<double, double> p0,
	double delta, bool earth_on_x_axis)
{
	// calculate earthcentric coordinates
	double d[] = {d0, d1, d2, d3};
	pair<double, double> x[] = {p0, p1, p2, p3};
	V3 v[4];
	
	FOR(0, 4)
	{
		Radians l, b;
		Radians ra = rad(x[i].first);
		Radians dec = rad(x[i].second);

		cerr << coordsys << "\n";
		if(coordsys == "galcart")
		{
			v[i].set(d[i] - Rg, x[i].first, x[i].second);
		}
		else
		{
			if(coordsys == "equ")
			{
				coordinates::equgal(ra, dec, l, b);
			}
			else if(coordsys == "gal")
			{
				l = ra; b = dec;
				cerr << deg(l) << " " << deg(b) << "\n";
			}
	
			v[i].celestial(d[i], l, b);
			
			// convert standard cartesian galactic to earthcentric galactic coordinates
			v[i].x = -v[i].x;
			v[i].y = -v[i].y;
		}

		cerr << "\tx" << i << " = " << d[i] << "pc, ra=" << deg(ra) << ", dec=" << deg(dec) << " " << v[i] << "\n";
	}
	cerr << "\tdelta = +-" << delta << "pc\n";

	pt.setup(v[1], v[2], v[3], v[0], delta, earth_on_x_axis);
}

int bin_plane_cut(
	const std::string &outfile, const std::string &volfn, 
	double dx, int ndx,
	pair<float, float> r, pair<float, float> ri,
	const std::string &coordsys,
	double d1, pair<double, double> p1,
	double d2, pair<double, double> p2,
	double d3, pair<double, double> p3,
	double d0, pair<double, double> p0,
	double delta, bool earth_on_x_axis)
{
	plane_pixelize_volume binner(volfn, dx, ndx, r, ri);
	setupPlaneTransformer(binner.pt, coordsys.size() ? coordsys : "equ", d1, p1, d2, p2, d3, p3, d0, p0, delta, earth_on_x_axis);
	cout << binner.run() << "\n";
	binner.br.loaded = true;

	binary_output_or_die(out, outfile);
	out << binner.br;

	//output_or_die(tout, "gccut.txt");
	//tout << binner.br;

	return 0;
}
#endif


class cylindrical_pixelize_volume : public bin3d_pixelize_volume
{
public:
	Radians phi0;
public:
	cylindrical_pixelize_volume(const std::string &volfn, float dx, int n, pair<float, float> r, pair<float, float> ri, Radians phi0_ = 0)
		: bin3d_pixelize_volume(volfn, dx, n, r, ri), phi0(phi0_)
		{}

	virtual double action(float ddv)
	{
		// bin in galactocentric cylindrical (rho, phi, z) coordinates
		V3 v(x + Rg, y, z);

		// find the median rho for this pixel, use it to determine phi
		// this is to minimize individual pixel distortions
		double rho = binnedPixelCoords(v.rho(), 0, 0).x * ndx;

		Radians phi = modulo(v.phi() - phi0, 2*ctn::pi);
		double rphi = phi*rho;
		double z = v.z;

		//cerr << v << " == " << V3(rho, rphi, z) << " == " << binnedPixelCoords(rho, rphi, z) << " " << ndx << " " << dx << "\n";

		br.pixels[binnedPixelCoords(rho, rphi, z)].volume += ddv;
		return ddv;
	}
};

int bin_cylindrical(
	const std::string &outfile, const std::string &volfn, 
	double dx, int ndx,
	pair<float, float> r, pair<float, float> ri,
	Radians phi0
)
{
	cylindrical_pixelize_volume binner(volfn, dx, ndx, r, ri, phi0);
	cerr << "phi0 = " << deg(phi0) << "\n";
	cout << binner.run() << "\n";
	binner.br.loaded = true;

	binary_output_or_die(out, outfile);
	out << binner.br;

	output_or_die(tout, "gccut.txt");
	tout << binner.br;

	return 0;
}
