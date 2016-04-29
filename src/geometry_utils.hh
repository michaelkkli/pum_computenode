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

#ifndef _GEOMETRY_UTILS_HH_
#define _GEOMETRY_UTILS_HH_

#include <valarray>

using std::valarray;

void generate_circular_boundary ( int                   num_pts,
                                  double                radius,
                                  const double *        centre,
                                  valarray<double> &    out,
                                  bool                  anticlockwise=true ); 

#endif // _GEOMETRY_UTILS_HH_
