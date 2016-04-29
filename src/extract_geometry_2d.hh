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

#ifndef _EXTRACT_GEOMETRY_2D_HH_
#define _EXTRACT_GEOMETRY_2D_HH_

#include <string>
#include <valarray>
using std::string;
using std::valarray;


template <typename Container>
void extract_geometry_2d ( const string& filename,
                           Container&,
#if 0
                           int gaussian_kernel_size = 5,
#endif
                           double binarize_threshold = 0.17*255.0,
                           int close_kernel_size = 20,
                           int island_keep_size = 512*512/10 );

void apply_gaussian ( const string& filename, int gaussian_kernel_size=5 );

template <typename value_t>
void simple_off_input_geometry_2d ( const string &    filename, valarray<value_t> &    geom );

#endif // _EXTRACT_GEOMETRY_2D_HH_
