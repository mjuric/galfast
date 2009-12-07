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

#ifndef compmgr_h__
#define compmgr_h__

#include <vector>
#include <map>
#include <astro/macros.h>

/**
	componentMap_t -- map between sparse compIDs and internally used sequential indices

	Sequential indices are used internally to refer to individual models.
	They're easier to index for fast lookup (via bitmaps), than the sparse
	compIDs.
*/
struct componentMap_t
{
public:
	std::map<uint32_t, uint32_t> comp2seq;	// map compID -> compIdx
	std::vector<uint32_t> seq2comp;		// vector of compIDs

public:
	uint32_t seqIdx(uint32_t compID);	// return an index corresponding to compID, creating one if it doesn't exist
	uint32_t compID(uint32_t seqIdx);	// return the componentID corresponding to component index

	size_t size() const { return seq2comp.size(); }
};

// globals
extern componentMap_t componentMap;	// global map of compIDs<->indices (defined in pipeline.cpp)

#endif // compmgr_h__
