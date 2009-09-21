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

#include "gpu.h"
#include <astro/constants.h>

/**
	Things "C"-ish enough for CUDA to swallow (and to upload to GPU) go into this header.
*/

// returns numeric_limits::epsilon() for the type of x
#define EPSILON_OF(x) std::numeric_limits<typeof(x)>::epsilon()

static const float ABSMAG_NOT_PRESENT = 99.999f;

static const int GAL = 0;
static const int EQU = 1;

	namespace galequ_constants
	{
		static const double angp = peyton::ctn::d2r * 192.859508333; //  12h 51m 26.282s (J2000)
		static const double dngp = peyton::ctn::d2r * 27.128336111;  // +27d 07' 42.01" (J2000)
		static const double l0   = peyton::ctn::d2r * 32.932;	// galactic longitude of ascending node of galactic coordinate system (where b=0, dec=0)
		static const double ce   = 0.88998740217659689; // cos(dngp)
		static const double se   = 0.45598511375586859; // sin(dngp)

		static const double halfpi = peyton::ctn::halfpi;
	};

	namespace float_galequ_constants
	{
		static const float angp = (float)galequ_constants::angp; //  12h 51m 26.282s (J2000)
		static const float dngp = (float)galequ_constants::dngp;  // +27d 07' 42.01" (J2000)
		static const float ce   = (float)galequ_constants::ce;
		static const float se   = (float)galequ_constants::se;
		static const float l0   = (float)galequ_constants::l0;

		static const float halfpi = (float)galequ_constants::halfpi;
	};

namespace peyton { namespace system { class Config; }};
class otable;
class osink;
class opipeline;

struct skygenConfig;
struct skypixel;
struct skyConfigInterface
{
	virtual bool init(
		otable &t,
		const peyton::system::Config &cfg,	// model cfg file
		const skygenConfig &sc,
		const skypixel *pixels) = 0;
	virtual void initRNG(rng_t &rng) = 0;		// initialize the random number generator from CPU RNG
	virtual double integrateCounts() = 0;		// return the expected starcounts contributed by this model
	virtual void setDensityNorm(float norm) = 0;
	virtual size_t run(otable &in, osink *nextlink) = 0;
	virtual ~skyConfigInterface() {};
};

#endif
