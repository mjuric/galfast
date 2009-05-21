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

#include "config.h"

#include "paralax.h"
#include <astro/system/log.h>
#include <astro/util.h>
#include <astro/useall.h>

double Rg = 8000.;

inline OSTREAM(const std::vector<double> &a)
{
	if(a.empty()) { return out; }
	out << a[0];
	FOR(1, a.size()) { out << " " << a[i]; }
	return out;
}

void plx_gri_locus_ng::setParalaxCoefficients(const std::vector<double> &Mcoef)
{
	Mrc = Mcoef;

	DLOG(verb1) << "color-absmag coef. = " << Mcoef;
}
