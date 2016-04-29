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

#ifndef _COMPUTENODE_UTILS_HH_
#define _COMPUTENODE_UTILS_HH_  

#include "computenode.hh"

#include <algorithms>
using std::max;

template<int dim>
void create_grid_on_bounding_box ( const box<dim>&,
				   vector<double>&,
				   int points_in_x_direction,
				   int points_in_y_direction=1,
				   int points_in_z_direction=1 );

template<int dim>  
void create_boxes_for_points ( const vector<double>&,
			       vector<box<dim> >&,
			       double x_size,
			       double y_size = 0.0,
			       double z_size = 0.0 );

#endif // _COMPUTENODE_UTILS_HH_
