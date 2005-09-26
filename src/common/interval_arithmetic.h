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

#ifndef interval_arithmetic_h
#define interval_arithmetic_h

#include <vector>
#include <utility>
#include <iostream>

#include <astro/util.h>
#include <astro/system/log.h>

#include "binarystream.h"

// Interval and interval set classes
typedef std::pair<float, float> interval;
typedef std::vector<float> interval_set;

inline OSTREAM(const interval &k)
{
	return out << "(" << k.first << ", " << k.second << ")";
}

inline OSTREAM(const interval_set &is)
{
	ASSERT(is.size() % 2 == 0);
	
	out << "[";
	FOREACH(is)
	{
		out << " (" << *i;
		out << ", " << *(++i) << ")";
	}
	out << " ]";
	
	return out;
}

inline BOSTREAM(const interval_set &is)
{
	ASSERT(is.size() % 2 == 0);
	
	out << is.size();
	FOREACH(is)
	{
		out << *i;
		++i;
		out << *i;
	}
	
	return out;
}

inline BISTREAM(interval_set &is)
{
	size_t size;
	in >> size;
	ASSERT(is.size() % 2 == 0);

	is.reserve(size);
	is.clear();
	interval_set::value_type v;
	while(size--)
	{
		in >> v;
		is.push_back(v);
	}

	return in;
}

inline void add_interval(interval_set &is, const interval &k)
{
	ASSERT(k.first <= k.second);

	if(k.first == k.second) return; // do not add 0-length intervals

	interval_set::iterator	x0 = lower_bound(is.begin(), is.end(), k.first),
				x1 = upper_bound(is.begin(), is.end(), k.second);
	int idx0 = x0 - is.begin();

	bool	odd0 = (idx0) % 2,
		odd1 = (x1 - is.begin()) % 2;

	is.erase(x0, x1);

	if(!odd0)
	{
		is.insert(is.begin() + idx0, k.first);
		idx0++;
	}
	
	if(!odd1)
	{
		is.insert(is.begin() + idx0, k.second);
	}
}

inline void add_interval(interval_set &is, const interval_set &is0)
{
	interval_set::value_type tmp;
	FOREACH(is0)
	{
		tmp = *i;
		++i;
		add_interval(is, std::make_pair(tmp, *i));
	}
}

inline void intersect_interval(interval_set &out, const interval_set &is, const interval &k)
{
	ASSERT(k.first <= k.second);

	interval_set::const_iterator	x0 = upper_bound(is.begin(), is.end(), k.first),
					x1 = lower_bound(is.begin(), is.end(), k.second);

	bool	odd0 = (x0 - is.begin()) % 2,
		odd1 = (x1 - is.begin()) % 2;

	out.clear();
	out.reserve(x1 - x0 + 2);

	if(odd0)
	{
		out.push_back(k.first);
	}

	std::insert_iterator<interval_set> ii(out, out.end());
	copy(x0, x1, ii);

	if(odd1)
	{
		out.push_back(k.second);
	}
}

inline void subtract_interval(interval_set &out, const interval_set &is, const interval &k)
{
	ASSERT(k.first <= k.second);

	interval_set::const_iterator	x0 = lower_bound(is.begin(), is.end(), k.first),
					x1 = upper_bound(is.begin(), is.end(), k.second);

	bool	odd0 = (x0 - is.begin()) % 2,
		odd1 = (x1 - is.begin()) % 2;

	out.clear();
	out.reserve(x1 - x0 + 2);

	copy(is.begin(), x0, std::insert_iterator<interval_set>(out, out.end()));
	if(odd0) { out.push_back(k.first); }
	if(odd1) { out.push_back(k.second); }
	copy(x1, is.end(), std::insert_iterator<interval_set>(out, out.end()));
}

inline interval_set::value_type length(const interval_set &is)
{
	ASSERT(is.size() % 2 == 0);

	interval_set::value_type sum = 0;
	FOREACH(is)
	{
		sum -= *i;
		++i;
		sum += *i;
	}
	return sum;
}

inline float length(const interval &k)
{
	return k.second - k.first;
}

#endif
