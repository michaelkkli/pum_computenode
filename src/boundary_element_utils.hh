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

#ifndef _BOUNDARY_ELEMENT_UTILS_HH_
#define _BOUNDARY_ELEMENT_UTILS_HH_

#include "boundary_element.hh"
#include "box.hh"

#include <iostream>
#include <vector>

using std::ostream;
using std::vector;

template<int dim>
void line_segment_get_point ( const double* line, double t, double* );

template<int dim>
bool intersect_line_with_plane ( const double* line,
				 const double* outward_normal,
				 const double* point_on_plane,
			         double* intersection_parameter );

void gp_draw_single( const boundary_element<1>&, ostream&, double level=0.0 );
void gp_draw_single( const boundary_element<2>&, ostream&, double level=0.0 );
void gp_draw_single( const boundary_element<3>&, ostream&, double ignored=0.0 );

template<int dim>
void gp_draw( const vector<boundary_element<dim> >&, ostream& );

#endif // _BOUNDARY_ELEMENT_UTILS_HH_
