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

#include "computenode_utils.hh"

template<int dim>
void create_grid_on_bounding_box ( const box<dim>& bx,
				   vector<double>& vec,
				   int num_x,
				   int num_y=1,
				   int num_z=1 )
{
  if ( num_x <= 0 || num_y <= 0 || num_z<=0 || bx.empty() ) {
    vec.clear();
    return;
  }
  int num_points;
  int max_num;
  if ( dim == 1 ) {
    num_points = num_x;
    max_num = num_x;
  } else if ( dim == 2 ) {
    num_points = num_x * num_y;
    max_num = std::max( max_num, num_y );
  } else if ( dim == 3 ) {
    num_points = num_x * num_y * num_z;
    max_num = std::max( max_num, num_z );
  }
  
  std::vector<double> temp( max_num );

  const double* vals = bx.get();

  double nums[3] = { num_x, num_y, num_z };

  for ( int i=0; i<dim; ++i ){
    // Separation between points.
    double d = (vals[i+1]-vals[i]) / (nums[i]+1);
    double half_d = 0.5 * d;

    for ( int j=0; j<nums[i]; ++j ) {
      temp[j] = half_d + d * j;
    }
  }
}
