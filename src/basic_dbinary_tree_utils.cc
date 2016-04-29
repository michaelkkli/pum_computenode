/*
	Copyright (C) 2009 Michael Li
	This file is part of the Computenode Library.

	The Computenode Library is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "basic_dbinary_tree_utils.hh"
#include "box.hh"
#include "box_utils.hh"

template <int dim>
void
gnuplot_output ( vector<string> &             keys,
                 basic_dbinary_tree<dim> &    dtree,
                 ostream &                   out,
                 bool                         labels )
{
	int num_keys = keys.size();
	
	vector<box<dim> >    boxes;
	dtree.get_box ( keys, boxes );
	
	out << "#!/usr/bin/gnuplot -persist\n";
	out << "set size square\n";
	
	double centre[dim];
	if ( labels ) {
		for ( int i=0; i<num_keys; ++i ) {
			out << "set label \""
			    << keys[i] << "\" at ";
			boxes[i].get_centre_point ( centre );
			for ( int d=0; d<dim; ++d ) {
				if ( d!=0 ) {
					out << ",";
				}
				out << centre[d];
			}
			out << "\n";
		}
	}
	
	out << "plot '-' w l\n";
	
	gp_draw ( boxes, out );
}

//

template void gnuplot_output<2> ( vector<string>&, basic_dbinary_tree<2>&, ostream&, bool );
template void gnuplot_output<3> ( vector<string>&, basic_dbinary_tree<3>&, ostream&, bool );
