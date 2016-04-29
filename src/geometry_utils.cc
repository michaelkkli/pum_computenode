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

#include "geometry_utils.hh"

#include <cassert>
#include <cmath>

using std::cos;
using std::sin;
using std::slice;


void generate_circular_boundary ( int                   num_pts,
                                  double                radius,
                                  const double *        centre,
                                  valarray<double> &    out,
                                  bool                  anticlockwise )
{
	// 2D only.
	out.resize ( 2 * num_pts );
	
	double angle_step = 2.*M_PI/num_pts;
	
	for ( int i=0; i<num_pts; ++i ) {
		out[2*i]   = cos ( angle_step*i );
		out[2*i+1] = sin ( angle_step*i );
	}
	
	out *= radius;

	if ( !anticlockwise ) {
		out[ slice(1,num_pts,2) ] *= valarray<double>(-1.0, num_pts);
	}
	
	if ( centre ) {
		double centre_x = centre[0];
		double centre_y = centre[1];
		
		out[ slice(0,num_pts,2) ] += valarray<double>(centre_x, num_pts);
		out[ slice(1,num_pts,2) ] += valarray<double>(centre_y, num_pts);
	}
	
}


