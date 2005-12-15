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

#ifndef simulate_h__
#define simulate_h__

#include "projections.h"
#include "model.h"
#include <string>


class cumulative_dist
{
public:
	typedef std::pair<double, double> value_type;
	std::vector<value_type> hist;
public:
	double operator()(double prob);
	void construct(const std::vector<double> &x, const std::vector<double> &y);
};

class sim
{
public:
	std::string prefix;	// prefix for input/output datafiles
	galactic_model *model;	// models for this simulation

	double dx; 		// model grid angular resolution
	double dri;		// model CMD resolution
	double ri0, ri1;	// color limits
	double m0, m1;		// magnitude limits
	double dm;		// model CMD magnitude resolution

public: // internal variables - usually you don't need to touch these (for northern sky)
	peyton::math::lambert proj;	// lambert projector object (by default, centered at north pole)
	double x0, x1, y0, y1;		// lambert survey footprint bounding box

public:
	sim(const std::string &pfx, galactic_model *model);

	struct star // storage structure for Monte Carlo generated stars
	{
		static sim *gsim;

		double x, y, ri, m;
		int X() const { return (int)((x - gsim->x0)/gsim->dx); }
		int Y() const { return (int)((y - gsim->y0)/gsim->dx); }
		int RI() const { return (int)((ri - gsim->ri0)/gsim->dri); }
	};

	void precalculate_mpdf();
	void magnitude_mpdf(cumulative_dist &mspl, double x, double y, double ri);
	void montecarlo(unsigned int K);

protected:
	double ri_mpdf(std::vector<double> &pdf, const double x, const double y);
};

#endif // simulate_h__
