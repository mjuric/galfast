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

#if 1
#include "dm.h"

void reprocess_driver();

main(int argc, char **argv)
{
try
{
	std::set<int> runs;
	loadRuns(runs, "/home/mjuric/projects/galaxy/workspace/catalogs/runs.txt");

	gsl_set_error_handler_off ();
#if 0
	//makelookup(runs, "dm_tmpcat.dmm", "dm_tmpcat_index.dmm", "dm_run_index.dmm");
	make_run_index_offset_map();
#endif

#if 0
	average_magnitudes();
	//starinfo(58932874);
#endif

#if 1
	reprocess_driver();
#endif

} catch (peyton::exceptions::EAny &e)
{
	e.print();
} catch (...)
{
	std::cout << "Uncaught exception!\n";
}

}
#else
#include <iostream>

int main(int argc, char **argv)
{
	std::cerr << "This exe has not been compiled because of the lack of CCfits library.\n";
	return -1;
}
#endif // HAVE_LIBCCFITS
