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

