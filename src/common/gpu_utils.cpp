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

#include <iostream>

///////////////////////////////////////////////////////////
// CUDA helpers

//
// Find the dimensions (bx,by) of a 2D grid of blocks that 
// has as close to nblocks blocks as possible
//
void find_best_factorization(int &bx, int &by, int nblocks)
{
	bx = -1;
	int best_r = 100000;
	for(int bytmp = 1; bytmp != 65536; bytmp++)
	{
		int r  = nblocks % bytmp;
		if(r < best_r && nblocks / bytmp < 65535)
		{
			by = bytmp;
			bx = nblocks / bytmp;
			best_r = r;
			
			if(r == 0) { break; }
			bx++;
		}
	}
	if(bx == -1) { std::cerr << "Unfactorizable?!\n"; exit(-1); }
}

//
// Given a total number of threads, their memory requirements, and the
// number of threadsPerBlock, compute the optimal allowable grid dimensions.
// Returns false if the requested number of threads are impossible to fit to
// shared memory.
//
bool calculate_grid_parameters(int gridDim[3], int threadsPerBlock, int neededthreads, int dynShmemPerThread, int staticShmemPerBlock)
{
	const int shmemPerMP =  16384;

	int dyn_shared_mem_required = dynShmemPerThread*threadsPerBlock;
	int shared_mem_required = staticShmemPerBlock + dyn_shared_mem_required;
	if(shared_mem_required > shmemPerMP) { return false; }

	// calculate the total number of threads
	int nthreads = neededthreads;
	int over = neededthreads % threadsPerBlock;
	if(over) { nthreads += threadsPerBlock - over; } // round up to multiple of threadsPerBlock

	// calculate the number of blocks
	int nblocks = nthreads / threadsPerBlock;
	if(nthreads % threadsPerBlock) { nblocks++; }

	// calculate block dimensions so that there are as close to nblocks blocks as possible
	find_best_factorization(gridDim[0], gridDim[1], nblocks);
	gridDim[2] = 1;

	return true;
}

