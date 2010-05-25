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

#include "galfast_config.h"

#include "column.h"
#include "modules/module_lib.h"

#include "modules/FeH_gpu.cu.h"
#include "modules/GaussianFeH_gpu.cu.h"
#include "modules/fixedFeH_gpu.cu.h"
#include "modules/unresolvedMultiples_gpu.cu.h"
#include "modules/photometry_gpu.cu.h"
#include "modules/kinTMIII_gpu.cu.h"
#include "modules/Bond2010_gpu.cu.h"
#include "modules/vel2pm_gpu.cu.h"
#include "modules/gal2other_gpu.cu.h"
