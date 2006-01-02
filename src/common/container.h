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

#ifndef _container_h__
#define _container_h__

#include <sstream>

template<typename C>
class container_aux
{
	C &c;
	typedef typename C::value_type value_type;

	class aux
	{
		adapt_c<C> c;
	public:
		aux(C &c_) : c(c_) { c.clear(); }
		aux insert(const value_type &v) { c.push_back(v); return *this; }
		aux operator ,(const value_type &v) { return insert(v); }
	};

public:
	container_aux(C &c_) : c(c_) {}
	aux operator=(const value_type &v) { return aux(c).insert(v); }
};
template<typename C>
container_aux<C> container(C &c) { return container_aux<C>(c); }



template<typename C>
std::string join(const std::string &separator, const C& c)
{
	std::ostringstream ss;
	FOREACH(c)
	{
		if(i != c.begin()) { ss << separator; }
		ss << *i;
	}
	return ss.str();
}

#endif
