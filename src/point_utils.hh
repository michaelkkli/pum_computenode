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

#ifndef _POINT_UTILS_HH_
#define _POINT_UTILS_HH_

#include "box.hh"

#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <valarray>

#include <cassert>

using std::binder1st;
using std::ostream;
using std::string;
using std::unary_function;
using std::valarray;
using std::vector;




void make_subdivisions ( double a, double b, int n, valarray<double>& out );

void make_point_vec ( const valarray<double>&,
		      const valarray<double>&,
		      valarray<double>& );

void make_point_vec ( const valarray<double>&,
		      const valarray<double>&,
		      const valarray<double>&,
		      valarray<double>& );

template<int dim>
void make_point_vec_on_box ( const box<dim>&,
			     valarray<double>&,
			     int points_in_x_direction,
			     int points_in_y_direction=1,
			     int points_in_z_direction=1 );

template<int dim>
void make_point_vec_interior_box ( const box<dim>&,
				   valarray<double>&,
				   int points_in_x_direction,
				   int points_in_y_direction=1,
				   int points_in_z_direction=1 );

void output_values ( const vector<valarray<double> >&, ostream& );

template<int dim>
void gp_draw_points ( const valarray<double>&, ostream&, double z_offset=0.0 );

template<int dim>
void gp_draw_points ( const valarray<double>& pts,
	const valarray<double>& vals,
	ostream&,
	int block=-1,
	double z_offset=0.0 );

template<int dim>
void gp_draw_labels ( string* prefix,
		      const valarray<double>& pts,
		      const vector<string>* labels,
		      std::ostream& out,
		      double offset_x=0.0,
		      double offset_y=0.0,
		      double offset_z=0.0 );

#if 0
// Implicit template instantiation will allow Pred to be deduced.
template<int dim>
void not_predicate_point_zero_values ( const interior_predicate<dim>& pred, const valarray<double>& points, valarray<double>& values );
#endif

#endif // _POINT_UTILS_HH_
