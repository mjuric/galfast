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

#ifndef projections__h
#define projections__h

#include <astro/math.h>
#include <cmath>

namespace peyton {
namespace math {

	inline bool between(const double x, const double x0, const double x1) { return x0 <= x && x < x1; }

	template<typename T>
	inline bool between(const T x, const std::pair<T, T> cmp) { return cmp.first <= x && x < cmp.second; }
	
	class lambert
	{
		Radians l0, phi1, cosphi1, sinphi1;
	
	public:
		lambert(Radians l0_ = peyton::math::rad(90), Radians phi1_ = peyton::math::rad(90))
			: phi1(phi1_), cosphi1(cos(phi1_)), sinphi1(sin(phi1_)), l0(l0_)
		{ }

		template<typename T> // usually, T = Radians (double), but it can also be valarray<double>
		void convert(const T l, const T phi, T &x, T &y) const
		{
			T kp = sqrt(2./(1+sinphi1*sin(phi)+cosphi1*cos(phi)*cos(l - l0)));
			x = kp*cos(phi)*sin(l - l0);
			y = kp*(cosphi1*sin(phi)-sinphi1*cos(phi)*cos(l-l0));
		}

		template<typename T> // usually, T = Radians (double), but it can also be valarray<double>
		void inverse(const T x, const T y, T &l, T &phi) const
		{
			T r = sqrt(x*x + y*y);
			T c = 2*asin(0.5*r);

			phi =     r != 0. ? asin(cos(c)*sinphi1 + y*sin(c)*cosphi1/r) : phi1;
			l = l0 + (r != 0. ? atan2(x*sin(c), r*cosphi1*cos(c) - y*sinphi1*sin(c)) : 0.);
		}
	};



	class gnomonic
	{
		double l0, phi1, cosphi1, sinphi1;
	
	public:
		gnomonic(Radians l0_ = 0., Radians phi1_ = 0.)
			: phi1(phi1_), cosphi1(cos(phi1_)), sinphi1(sin(phi1_)), l0(l0_)
		{ }

		template<typename T> // usually, T = Radians (double), but it can also be valarray<double>
		void convert(const T l, const T phi, T &x, T &y) const
		{
			const T cosc = sinphi1*sin(phi)+cosphi1*cos(phi)*cos(l-l0);
			x = cos(phi)*sin(l - l0) / cosc;
			y = (cosphi1*sin(phi)-sinphi1*cos(phi)*cos(l-l0)) / cosc;
		}

		template<typename T> // usually, T = Radians (double), but it can also be valarray<double>
		void inverse(const T x, const T y, T &l, T &phi) const
		{
			T r = sqrt(x*x + y*y);
			T c = atan(r);

			phi =     r != 0. ? asin(cos(c)*sinphi1 + y*sin(c)*cosphi1/r) : phi1;
			l = l0 + (r != 0. ? atan2(x*sin(c), r*cosphi1*cos(c) - y*sinphi1*sin(c)) : 0.);
		}
	};
} // math
} // peyton
	
#endif
