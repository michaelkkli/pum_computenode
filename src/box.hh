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

#ifndef _BOX_HH_
#define _BOX_HH_

#include <boost/shared_ptr.hpp>

#include <valarray>
#include <vector>
using boost::shared_ptr;
using std::valarray;
using std::vector;

class box_intersect_line_scratch;

/**
   The class box is a concrete class (see ``The C++ Programming Language''
   by Bjarne Stroustrup). It is to be light-weight and testable; without
   dependency on anything else.
 */
template<int dim>
class box {
public:
	typedef shared_ptr<box<dim> > ptr;
	typedef vector<ptr>           vec_ptr;
public:
  // These are necessary as we sometime want to modify the base
  // box but sometimes we want to keep the base box unchanged.
  static void expand_by_factor ( const box<dim>&, double, box<dim>& );
  static void translate_centre ( const box<dim>&, const double*, box<dim>& );
public:
  box();
  ~box();
  // Default copy constructor, and assignment operator are correct.

public: // Non-const member functions.
  void set( const double* );
  void scale( const double* );
  void scale( double );
  void translate( const double* );
  void recentre( const double* );
  void scale_and_translate ( double a0, double a1, double a2,
			     double d0, double d1, double d2 );
  void clip_against( const box<dim>& );
public: // Const member functions.
  inline
  const double* get() const
  {
    return extents;
  }
  
  void get_centre_point ( double* ) const;

  bool empty() const;
  double measure() const;
  
  bool closed_intersect_box( const box<dim>& ) const;
  bool open_intersect_box( const box<dim>& ) const;
  
  bool closed_contains_box( const box<dim>& ) const;
  bool open_contains_box( const box<dim>& ) const;

  bool closed_intersect_point ( const double* ) const;
  bool open_intersect_point ( const double* ) const;
	
	void closed_intersect_point ( valarray<double>&, valarray<bool>& ) const;
  
  void intersect_line ( const double*,
                        box_intersect_line_scratch &,
                        bool* intersects_closed,
                        bool* intersects_open ) const;
	
  bool closed_intersect_line ( const double*, box_intersect_line_scratch & ) const;
  bool open_intersect_line ( const double*, box_intersect_line_scratch & ) const;
  
  // For change of coordinates.
  void global_partial_local ( double* ) const;
  
  // For change of coordinates.
  void local_partial_global ( double* ) const;


	/**
		The pointer is to the first of 2*dim doubles which represent a line
		lying within the closed box. This member function works out the length
		of the corresponding line in local-coordinates.
	*/
	double local_line_length ( double * ) const;


public: // Const member functions related to coordinate transformations of single points.
  /**
     Map local coordinates to global coordinates.
     x(e) = (a+b)/2 + e*(b-a)/2
	@param const double* Local coordinates.
	@param double* Global coordinates.
   */
  void map_local_to_global ( const double*, double* ) const;

  /**
     Map global coordinates to local coordinates.
     e(x) = ( x - (a+b)/2 ) * 2/(b-a)
     @param const double* Global coordinates.
     @param double* Local coordinates.
   */
  void map_global_to_local ( const double*, double* ) const;

  /**
     Map restricted local coordinates to local coordinates.
     For a \< c \< d \< b where (c,d) is contained in (a,b).
     x(r) = (c+d)/2 + r*(d-c)/2
     x(e) = (a+b)/2 + e*(b-a)/2
     so
     e(r) = ( (c+d-a-b) + r*(d-c) )/(b-a)
     @param const double* Restricted local coordinates.
     @param double* Local coordinates.
  */
  void map_restricted_local_to_local ( const box<dim>& restricted_local,
				       const double*,
				       double* ) const;
  
  // The functionality of a ``map_restricted_local_to_global'' is the same
  // as applying ``map_local_to_global'' on the restricted local box.

public: // Const member functions related to coordinate transformations of multiple points.
  void map_local_to_global ( valarray<double>&,
			     valarray<double>& ) const;
	
	void map_local_to_global ( valarray<double>& in, valarray<bool>& pred, valarray<double>& out ) const;

  void map_global_to_local ( valarray<double>&,
			     valarray<double>& ) const;
	
	void map_global_to_local ( valarray<double>& in, valarray<bool>&, valarray<double>& out ) const;

  void map_restricted_local_to_local ( const box<dim>& restricted_local,
				       valarray<double>&,
				       valarray<double>& ) const;
  
  // The functionality of a ``map_restricted_local_to_global'' is the same
  // as applying ``map_local_to_global'' on the restricted local box.

private:
  inline
  double* get()
  {
    return &extents[0];
  }

private:
  double extents[ 2*dim ];
};

struct box_intersect_line_scratch {
	vector<box<1> >    box_vector;
	box<1>             all_intersected;
};

#endif // _BOX_HH_
