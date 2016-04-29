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

#ifndef _BOUNDARY_HH_
#define _BOUNDARY_HH_

#include "box.hh"
#include "boundary_element.hh"

#include "hacker_shimrat_algorithm112.hh"

#include <algorithm>
#include <iostream>
#include <valarray>
#include <vector>

using std::ostream;
using std::valarray;
using std::vector;

template<int dim>
class boundary {
public:
  boundary();
  ~boundary();
  boundary( const boundary<dim>& other );
  boundary<dim>& operator= ( const boundary<dim>& other );
public:
  void deep_copy( const boundary<dim>& other );
  void simple_set( valarray<double>& );
//  void set_simplicies ( const valarray<double>& );
  void update();
  void get_bounding_box ( box<dim>& ) const;
  bool closed_intersect_box ( const box<dim>& ) const;
  bool open_intersect_box ( const box<dim>& ) const;
  bool open_intersect_point ( const double* ) const;
  void gp_draw_boundary( ostream&, double level=0.0) const;
private:
  bool     updated;
  box<dim> bounding_box;
  vector<boundary_element<dim> > elems;
  valarray<double>               set_vals; // Mostly for the benefit of the 1D case.
};

#endif // _BOUNDARY_HH_
