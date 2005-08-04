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
		double l0, cosphi1, sinphi1;
	
	public:
		lambert(Radians l0_ = peyton::math::rad(90), Radians phi1 = peyton::math::rad(90))
			: cosphi1(cos(phi1)), sinphi1(sin(phi1)), l0(l0_)
		{ }

		void convert(const Radians l, const Radians phi, double &x, double &y) const
		{
			double kp = sqrt(2./(1+sinphi1*sin(phi)+cosphi1*cos(phi)*cos(l - l0)));
			x = kp*cos(phi)*sin(l - l0);
			y = kp*(cosphi1*sin(phi)-sinphi1*cos(phi)*cos(l-l0));
		}
	};

} // math
} // peyton
	
#endif
