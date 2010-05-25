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

#include "galfast_config.h"

#include "transform.h"

#include <astro/coordinates.h>
#include <astro/system/config.h>
#include <astro/exceptions.h>
#include <astro/util.h>
#include <astro/math.h>
#include <cux.h>

#include <iostream>
#include <iomanip>
#include <cstdlib>

#include "../modules/module_lib.h"

#include <astro/useall.h>

//
// Multiply two 3D matrices
//
void matmul3d(double C[3][3], double A[3][3], double B[3][3])
{
	for(int i = 0; i != 3; i++)
	{
		for(int j = 0; j != 3; j++)
		{
			C[i][j] = 0;
			for(int k = 0; k != 3; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void lbd2xyz(double v[3])
{
	// Convert (l, b, D) -> (x, y, z). lb must be in radians.
	double cl = cos(v[0]), cb = cos(v[1]);
	double sl = sin(v[0]), sb = sin(v[1]);
	double d = v[2];

	v[0] = d*cl*cb;
	v[1] = d*sl*cb;
	v[2] =    d*sb;
}

void get_gal2xxx_matrix(double go[3][3], const std::string &coordsys, const peyton::system::Config &cfg)
{
	using namespace peyton::coordinates;

	// get a matrix that rotates a point xyz in Galactic coordinate
	// system to coordinate system coordsys
	//
	// Valid values for coordsys are 'gal', 'equ' and 'galplane'
	//
	if(coordsys == "gal" || coordsys == "galcentric")
	{
		// nothing to do -- return a unit matrix
		memset(go, 0, sizeof(go[0][0])*9);
		go[0][0] = go[1][1] = go[2][2] = 1.f;
	}
	else if(coordsys == "equ")
	{
		// Galactic -> equatorial transformation
		double x[3], y[3], z[3];
		galequ(0., 0., 		x[0], x[1]); x[2] = 1;	// Galactic x axis in equatorial system
		galequ(rad(90.), 0.,	y[0], y[1]); y[2] = 1;	// Galactic y axis in equatorial system
		galequ(0., rad(90.),	z[0], z[1]); z[2] = 1;	// Galactic z axis in equatorial system

#if 0
		std::cout << "x(lbd)=" << std::setprecision(10) << deg(x[0]+(x[0]<0)*ctn::twopi) << " " << deg(x[1]) << " " << x[2] << "\n";
		std::cout << "y(lbd)=" << std::setprecision(10) << deg(y[0]+(y[0]<0)*ctn::twopi) << " " << deg(y[1]) << " " << y[2] << "\n";
		std::cout << "z(lbd)=" << std::setprecision(10) << deg(z[0]+(z[0]<0)*ctn::twopi) << " " << deg(z[1]) << " " << z[2] << "\n";
#endif

		// transform to cartesian
		lbd2xyz(x);
		lbd2xyz(y);
		lbd2xyz(z);
		
		// pack into the rotation matrix M = (x, y, z)
		go[0][0] = x[0]; go[0][1] = y[0]; go[0][2] = z[0];
		go[1][0] = x[1]; go[1][1] = y[1]; go[1][2] = z[1];
		go[2][0] = x[2]; go[2][1] = y[2]; go[2][2] = z[2];

#if 0
		// test
		print_matrix(go);

		float3 v = { 0.f, 0.f, 1.f };
		v = matmul3d(go, v);
		std::cout << std::setprecision(10) << v.x << " " << v.y << " " << v.z << "\n";
		Radians ra = atan2(v.y, v.x); ra += (ra < 0)*ctn::twopi;
		Radians dec = asin(v.z / sqrt(sqr(v.x) + sqr(v.y) + sqr(v.z)));
		std::cout << "ra=" << deg(ra) << " dec=" << deg(dec) << "\n";
		abort();
#endif
	}
	else if(coordsys == "galplane")
	{
		// parameter-defined system!
		double Rg = cfg.get("Rg");
		double z0 = cfg.get("z0");

		// transform _from_ Galactic _to_ Galactic plane system.
		// This is a rotation by -alpha around y axis, followed
		// by a 180deg rotation around the z axis.
		double sina = -z0/Rg, cosa = sqrt(1 - sqr(sina));
		go[0][0] = -cosa; go[0][1] =  0; go[0][2] = sina;
		go[1][0] =     0; go[1][1] = -1; go[1][2] =    0;
		go[2][0] =  sina; go[2][1] =  0; go[2][2] = cosa;
	}

#if 0
	//std::cerr << deg(phi) << " " << deg(theta) << " " << deg(psi) << "\n";
	print_matrix(go);

	// Do some testing here...
	//float3 v = { -8000.f, 0.f, 0.f };
	float3 v = { 0.f, 0.f, 1.f };
	v = matmul3d(go, v);
	std::cout << std::setprecision(10) << v.x << " " << v.y << " " << v.z << "\n";
	exit(0);
#endif
}

void xxx2galxyz(const peyton::system::Config &cfg, const std::string &coordsys, const std::string &type, double v[3])
{
	using namespace peyton::coordinates;
	bool unkt = false, unkcs = false;

	if(type.empty())
	{
		v[0] = v[1] = v[2] = 0.;
	}

	// Convert a vector in coordinate system coordsys to Galactic cartesian
	// coordinate system
	if(coordsys == "gal")
	{
		if(type == "sph")
		{
			v[0] = rad(v[0]);
			v[1] = rad(v[1]);
			lbd2xyz(v);
		}
		else if(type == "xyz" || type.empty())
		{
			// no need to do anything
		}
		else unkt = true;
	}
	else if(coordsys == "equ")
	{
		if(type == "sph" || type.empty())
		{
			v[0] = rad(v[0]);
			v[1] = rad(v[1]);
			equgal(v[0], v[1], v[0], v[1]);
			lbd2xyz(v);
		}
		else unkt = true;
	}
	else if(coordsys == "galcentric")
	{
		if(type == "xyz" || type.empty())
		{
			// parameter-defined system!
			double Rg = cfg.get("Rg");
			
			// transform _from_ Galactocentric _to_ Galactic cartesian
			v[0] = Rg - v[0];
			v[1] = -v[1];
		}
		else unkt = true;
	}
	else if(coordsys == "galplane")
	{
		if(type == "xyz" || type.empty())
		{
			// parameter-defined system!
			double Rg = cfg.get("Rg");
			double z0 = cfg.get("z0");
			
			// transform _from_ Galactic plane system _to_ Galactic cartesian
			// This is a 2D rotation by -alpha, where alpha = asin(z/Rg)
			// followed by the usual Galactocentric->Galactic conversion
			double sina = -z0/Rg, cosa = sqrt(1 - sqr(sina));
			double x = cosa * v[0] - sina * v[2];
			double z = sina * v[0] + cosa * v[2];
			v[0] = Rg - x;
			v[1] = -v[1];
			v[2] = z;
		}
		else unkt = true;
	}
	else unkcs = true;
	
	if(unkcs) { THROW(EAny, "Unknown coordinate system " + coordsys + "."); }
	if(unkt)  { THROW(EAny, "Unknown rotation type " + type + "."); }
}

void get_rotation_matrix(double M[3][3], const std::string &spec, const peyton::system::Config &cfg)
{
	//
	// Returns a matrix transforming a point in the Galactic
	// frame to the xyz frame specified by the orientation specification <spec>.
	//
	// The general format of <spec> is that of '<coordsys> <kind> [params...]' 
	// where <kind> specifies the way in which the rotation is given (e.g., 
	// Euler angles, axis+angle, etc.) while <coordsys> specifies the starting 
	// coordinate system in which the rotation happens
	//

	// no rotation at all
	if(spec.empty())
	{
		memset(M, 0, sizeof(M[0][0])*9);
		M[0][0] = M[1][1] = M[2][2] = 1.;
		return;
	}

	// starting coordinate system (starting orientation)
	std::istringstream in(spec.c_str());
	std::string coordsys;
	in >> coordsys;
	double R[3][3];
	get_gal2xxx_matrix(R, coordsys, cfg); // galactic->"origin" system transformation

	std::string kind;
	if(!(in >> kind))
	{
		// no extra rotation beyond coordinate system rotation
		memcpy(M, R, sizeof(R[0][0])*9);
		return;
	}

	double rot[3][3];	
	if(kind == "axisangle")
	{
		// computes: a rotation matrix transforming a point in "origin" frame to the xyz frame of the target
		//
		// axisangle <x> <y> <z> <angle>
		//
		// The (xyz,angle) quad specifies the orientation of the target coordinate system wrt. the origin.
		// Imagine it as having the target system rotated by an angle alpha ccw around the axis (x,y,z).
		//
		double x, y, z, alpha;
		in >> x >> y >> z >> alpha;
		if(!in) { THROW(EAny, "Four numbers (x, y, z, angle) required for 'axisangle' rotation."); }

		// normalize the vector
		double norm = sqrt(x*x + y*y + z*z);
		x /= norm; y /= norm; z /= norm;

		// flip the angle (we want to go from origin->target frame)
		alpha = -rad(alpha);

		// construct the rotation matrix
		// picked up from http://en.wikipedia.org/wiki/Rotation_matrix#Axis_and_angle
		double c = cos(alpha), s = sin(alpha), C = 1.-c;
		double	xs  = x*s,   ys = y*s,   zs = z*s,
			xC  = x*C,   yC = y*C,   zC = z*C,
			xyC = x*yC, yzC = y*zC, zxC = z*xC;
		rot[0][0] = x*xC+c;   rot[0][1] = xyC-zs;   rot[0][2] = zxC+ys;
		rot[1][0] = xyC+zs;   rot[1][1] = y*yC+c;   rot[1][2] = yzC-xs;
		rot[2][0] = zxC-ys;   rot[2][1] = yzC+xs;   rot[2][2] = z*zC+c;
	}
	else if(kind == "euler")
	{
		// computes: a rotation matrix transforming a point in "origin" frame to the xyz frame of the target
		//
		// If imagining this as the rotation of the coordinate system from the origin to target, the following
		// rotation sequence occurs:
		//
		//	1) ccw around the z axis by an angle phi
		//	2) ccw around the x axis by theta
		//	3) ccw around the z axis by psi
		//
		// See http://mathworld.wolfram.com/EulerAngles.html, eqs 6-14
		//
		double phi, theta, psi;
		in >> phi >> theta >> psi;
		if(!in) { THROW(EAny, "Three angles (phi, theta, psi) required for 'euler' rotation."); }
		phi = rad(phi); theta = rad(theta); psi = rad(psi);
	
		rot[0][0] =  cos(psi)   * cos(phi)   - cos(theta) * sin(phi) * sin(psi);
		rot[0][1] =  cos(psi)   * sin(phi)   + cos(theta) * cos(phi) * sin(psi);
		rot[0][2] =  sin(psi)   * sin(theta);
		rot[1][0] = -sin(psi)   * cos(phi)   - cos(theta) * sin(phi) * cos(psi);
		rot[1][1] = -sin(psi)   * sin(phi)   + cos(theta) * cos(phi) * cos(psi);
		rot[1][2] =  cos(psi)   * sin(theta);
		rot[2][0] =  sin(theta) * sin(phi);
		rot[2][1] = -sin(theta) * cos(phi);
		rot[2][2] =  cos(theta);
	}
	else
	{
		THROW(EAny, "Unknown rotation type '" + kind + "'");
	}

	// compose the rotations from Galactic->origin and origin->target
	// This is just 3D matrix multiplication: M = rot*R
	matmul3d(M, rot, R);

#if 0
	//std::cerr << deg(phi) << " " << deg(theta) << " " << deg(psi) << "\n";
	print_matrix(rot);
	std::cout << "     x\n";
	print_matrix(R);
	std::cout << "     =\n";
	print_matrix(M);

	// Do some testing here...
	float3 v = { -8000.f, 1000.f, 0.f };
	v = matmul3d(M, v);
	std::cout << std::setprecision(10) << v.x << " " << v.y << " " << v.z << "\n";
	abort();
#endif
}

void load_transform(float tran[3], float rot[3][3], const std::string &center_str, const std::string &orientation_str, const peyton::system::Config &cfg)
{
	// Returns the transformation vector and matrix that transform 
	// a point xyz in Galactic cartesian system to internal module
	// coordinate system xyz'
	//
	// The returned transformations are to be applied as follows:
	//
	//     xyz' = rot * (xyz + tran)
	//
	// that is, translation by tran, followed by multiplication by
	// matrix rot (usually the rotation matrix)
	//
	// Currently, the two keyword are read to construct the transformation:
	//	center = <coordsys> <vector>
	//	rotation = <coordsys> <rotationtype> <parameters>
	//
	// Possible values of coordsys are:
	//	equ, gal, galxyz, galcentric_xyz, galplane_xyz
	// see docs/coordinatesystems.txt for a description

	// load translation vector
	{
		double ctr[3] = { 0 };
		std::string coordsys, type;
		std::istringstream ss(center_str);
		ss >> coordsys >> type >> ctr[0] >> ctr[1] >> ctr[2];
		xxx2galxyz(cfg, coordsys, type, ctr);
		tran[0] = -ctr[0];
		tran[1] = -ctr[1];
		tran[2] = -ctr[2];
	}

	// load rotation vector
	{
		double M[3][3];
		get_rotation_matrix(M, orientation_str, cfg);

		// copy the output matrix
		for(int i = 0; i != 3; i++)
			for(int j = 0; j != 3; j++)
				rot[i][j] = M[i][j];
	}
}
