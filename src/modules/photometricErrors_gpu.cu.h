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

#ifndef photometricErrors_gpu_cu_h__
#define photometricErrors_gpu_cu_h__

//
// Module data passed to device kernel
//

// -- none for now (this is a CPU module; TODO: port this module to GPU)

//
// Device kernel implementation
//
#if (__CUDACC__ || BUILD_FOR_CPU)

// -- no kernel yet (this is a CPU module; TODO: port this module to GPU)

#endif // (__CUDACC__ || BUILD_FOR_CPU)

#endif // photometricErrors_gpu_cu_h__
