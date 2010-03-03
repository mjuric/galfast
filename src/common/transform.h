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

#ifndef transform_h__
#define transform_h__

#include <astro/system/config.h>
#include <string>
#include <cstdio>

// loads a coordinate system transform (translation+rotation) given the keywords and config object
// see the comments in transform.cpp for more details.
void load_transform(float tran[3], float rot[3][3], const std::string &center_str, const std::string &orientation_str, const peyton::system::Config &cfg);

template<typename T>
void print_matrix(T rot[3][3])
{
	printf("% 11.9f % 11.9f % 11.9f\n", rot[0][0], rot[0][1], rot[0][2]);
	printf("% 11.9f % 11.9f % 11.9f\n", rot[1][0], rot[1][1], rot[1][2]);
	printf("% 11.9f % 11.9f % 11.9f\n", rot[2][0], rot[2][1], rot[2][2]);
}

#endif // transform_h__
