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

#ifndef _BOUNDARY_ELEMENT_HH_
#define _BOUNDARY_ELEMENT_HH_

#include <string>
#include <vector>

using std::string;
using std::vector;

// Only references used to box<dim> and boundary_element<dim>
// so we do not need to introduce unnecessary dependencies by
// including unneeded class declarations (forward declaration is sufficient).
template<int dim> class box;
template<int dim> class boundary_element;

class box_intersect_line_scratch;

template<int dim>
class boundary_element {
public:
  boundary_element();
  ~boundary_element();
  void deep_copy( const boundary_element<dim>& );

// Unneeded.

  boundary_element( const boundary_element<dim>& );
  boundary_element<dim>& operator= ( const boundary_element<dim>& );


public:
  /**
  	1D - take one pair of double (two doubles in total)  to give location and direction.
	2D - take two pairs of double (four doubles in total) to define line and deduce outward normal.
  	3D - take three triples of double (nine doubles in total) to define triangle.
  */
  void set_oriented_simplex ( const double* );
  const double* get_point( int ) const;
  void get_centre_point( double* ) const;
  double distance_centre_to_point ( const double* ) const;
  void finalize();
	bool is_finalized() const;
  void includes_point ( const double* ) const;
  void get_bounding_box ( box<dim>& bx ) const;
	
  bool closed_intersect_box ( const box<dim>&, box_intersect_line_scratch & ) const;
  bool open_intersect_box ( const box<dim>&, box_intersect_line_scratch & ) const;
	
private: // Member variables to be copied.
  bool   finalized; // Normalized plus other details.
  double endpoints[dim*dim]; // Point in 1D, line in 2D, triangle in 3D.
  double outward[dim];
  bool   axis_aligned;
  // "+0" or "-0" along x-axis; "+1" or "-1" along y-axis; "+2" or "-2" along z-axis.
  string preferred_axis;
#if 0
private:
  friend bool closed_intersect<> ( const box<dim>&, const boundary_element<dim>& );
  friend bool   open_intersect<> ( const box<dim>&, const boundary_element<dim>& );
#endif // 0
};

#endif // _BOUNDARY_ELEMENT_HH_
