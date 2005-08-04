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

#include <sstream>
#include <astro/system/options.h>

#include "dm.h"

#include <sqlplus.hh>

#include <astro/useall.h>
using namespace std;

int main(int argc, char **argv)
{
try
{
	VERSION_DATETIME(version);
	DMMArray<mobject> arr("dm_unique_stars.dmm");

	try {
		Connection con("galaxy");
		Query query = con.query();
		
		query << "select * from stars";
		Result res = query.store();
	
		cout << "Existing Records Found: " << res.size() << endl << endl;

		cout << "Clearing table.\n";
		query.execute("delete from stars");

		cout << "Importing data...\n";
		ticker tick(10000);
		FOR(0, arr.size())
		{
			mobject m = arr[i];

			FORj(j, 0, 3) { if(!isfinite(m.magErr[j])) { m.magErr[j] = 0; } }

			query << "insert into stars values (" <<
				i << ',' << setprecision(15) << m.obs_offset << ',' << setprecision(15) << m.n << ',' << setprecision(15) << m.Ar << ',' << setprecision(15) <<
				m.mag[0] << ',' << setprecision(15) << m.magErr[0] << ',' << setprecision(15) << m.N[0] << ',' << setprecision(15) <<
				m.mag[1] << ',' << setprecision(15) << m.magErr[1] << ',' << setprecision(15) << m.N[1] << ',' << setprecision(15) <<
				m.mag[2] << ',' << setprecision(15) << m.magErr[2] << ',' << setprecision(15) << m.N[2] << ',' << setprecision(15) <<
				m.flags << ',' << setprecision(15) <<
				m.ml_mag[0] << ',' << setprecision(15) << m.ml_mag[1] << ',' << setprecision(15) << m.ml_mag[2] << ',' << setprecision(15) <<
				m.D << ',' << setprecision(15) <<
				m.ra << ',' << setprecision(15) << m.dec
				<< ")";
//			cout << "Query : " << query.preview() << endl; 
			query.execute(RESET_QUERY);
			tick.tick();
		}
		
	} catch (BadQuery er) { // handle any connection or
				// query errors that may come up
		cerr << "Error: " << er.error << endl;
		return -1;
	} catch (BadConversion er) { // handle bad conversions
		cerr << "Error: Tried to convert \"" << er.data << "\" to a \""
			<< er.type_name << "\"." << endl;
		return -1;
	} 

}
catch(EAny &e)
{
	e.print();
}
}
