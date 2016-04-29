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

#include "boundary_element.hh"

#include "box.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>

using std::copy;
using std::cout;
using std::min;
using std::max;
using std::abs;
using std::sqrt;
using std::ostringstream;

template<int dim>
boundary_element<dim>::boundary_element () : finalized( false ), axis_aligned(false), preferred_axis("")
{
}

template<int dim>
boundary_element<dim>::~boundary_element()
{
}
template<int dim>
void
boundary_element<dim>::deep_copy( const boundary_element<dim>& other )
{
	finalized      = other.finalized;
	
	copy( other.endpoints, other.endpoints + dim*dim, this->endpoints  );
	copy( other.outward, other.outward + dim, this->outward );
	
	axis_aligned   = other.axis_aligned;
	preferred_axis = other.preferred_axis;
}

// Unneeded.


template<int dim> 
boundary_element<dim>::boundary_element ( const boundary_element<dim>& other )
{
	this->deep_copy( other );
}


template<int dim>
boundary_element<dim>&
boundary_element<dim>::operator= ( const boundary_element<dim>& other )
{
	this->deep_copy( other );
	return *this;
}


template <>
void
boundary_element<1>::set_oriented_simplex ( const double* sim )
{
	assert( sim );
	endpoints[0] = sim[0];
	outward[0]   = sim[1];
}

/**
	To get outward pointing vector, work out vector from first
	point to the second point. Reflect in y-axis (so multiply first coordinate
	by -1) and then reflect in the line y=x (swap the coordinates).

	Vector is (sim[2]-sim[0], sim[3]-sim[1]) so outward normal (assuming
	coordinates given in an anticlockwise winding sense and the inside is on
	the left)

		( sim[3] - sim[1], sim[0] - sim[2] )
*/
template <>
void
boundary_element<2>::set_oriented_simplex ( const double* sim )
{
	assert( sim );
	endpoints[0] = sim[0];
	endpoints[1] = sim[1];
	endpoints[2] = sim[2];
	endpoints[3] = sim[3];
	
	outward[0] = sim[3] - sim[1];
	outward[1] = sim[0] - sim[2];
}

/**
	Assume anticlockwise winding of the boundary of the simplex. We calculate the
	vector of the first edge out from the first point.
		( sim[3]-sim[0], sim[4]-sim[1], sim[5]-sim[2] )
	From the first point to the last point is
		( sim[6]-sim[0], sim[7]-sim[1], sim[8]-sim[2] )
	Let the above be (a,b,c)
	             and (d,e,f).
	The cross product gives outward vector
		( bf-ce, cd-af, ae-bd )
*/
template <>
void
boundary_element<3>::set_oriented_simplex ( const double* sim )
{
	assert( sim );
	endpoints[0] = sim[0];
	endpoints[1] = sim[1];
	endpoints[2] = sim[2];
	endpoints[3] = sim[3];
	endpoints[4] = sim[4];
	endpoints[5] = sim[5];
	endpoints[6] = sim[6];
	endpoints[7] = sim[7];
	endpoints[8] = sim[8];
	
	double a = sim[3]-sim[0];
	double b = sim[4]-sim[1];
	double c = sim[5]-sim[2];
	
	double d = sim[6]-sim[0];
	double e = sim[7]-sim[1];
	double f = sim[8]-sim[2];
	
	outward[0] = b*f-c*e;
	outward[1] = c*d-a*f;
	outward[2] = a*e-b*d;
}

template<int dim>
const double*
boundary_element<dim>::get_point( int d ) const
{
	assert( 0<=d && d<dim );
	if ( 0>d || d>=dim ) {
		return 0;
	} else {
		return &endpoints[ d*dim ];
	}
}

template <>
void
boundary_element<1>::get_centre_point( double* co ) const
{
	assert( co );
	co[0] = endpoints[0];
}

template<int dim>
void
boundary_element<dim>::get_centre_point( double* co ) const
{
	assert( co );
	double tmp;
	for ( int d=0; d<dim; ++d ) {
		tmp = 0.0;
		for ( int i=0; i<dim; ++i ) {
			tmp += endpoints[ dim*i + d ];
		}
		co[d] = tmp;
	}
}

template<int dim>
double
boundary_element<dim>::distance_centre_to_point ( const double* co ) const
{
	assert( co );
	
	double cen[dim];
	get_centre_point( cen );
	
	double tmp = 0.0;
	double term;
	for ( int d=0; d<dim; ++d ) {
		term = ( co[d] - cen[d] );
		tmp += term*term;
	}
	return sqrt( tmp );
}

template  <int dim>
void
boundary_element<dim>::finalize()
{
	double tmp = 0.0;
	if ( dim == 1 ) {
	    if ( outward[0] > 0.0 ) {
		outward[0] = 1.0;
	    } else {
		outward[0] = -1.0;
	    }
		axis_aligned = true;
		if ( outward[0] > 0.0 ) {
			preferred_axis = "+0";
		} else {
			preferred_axis = "-0";
		}
	} else {
		for ( int d=0; d<dim; ++d ) {
			tmp += outward[d]*outward[d];
		}
		tmp = sqrt( tmp );
		assert ( tmp > 0.0 );
		for ( int d=0; d<dim; ++d ) {
			outward[d] /= tmp; // Normalize outward pointing normal.
		}
		// This lower loop relies on the upper loop
		// having performed the normalization.
		axis_aligned = false;
		int max_d = 0;
		double largest_projection = abs( static_cast<long double>(outward[0]) );
		for ( int d=0; d<dim; ++d ) {
			if( abs( static_cast<long double>(outward[d]) ) > largest_projection ) {
				max_d = d;
				largest_projection = abs( static_cast<long double>(outward[d]) );
			}
		}
		if ( outward[max_d] > 0.0 ) {
			ostringstream oss;
			oss << "+" << max_d;
			preferred_axis = oss.str();
			
			// Preferred axis now correct.
			// 2008-09-01 Mike Li.
			// cout << " preferred_axis is " << preferred_axis <<"\n";
		} else {
			ostringstream oss;
			oss << "-" << max_d;
			preferred_axis = oss.str();
		}
		// We take 2*pi/10000 to be a very small angle and if
		// projection of the normal to the axis 'max_d' results in
		// angle less than 2*pi/10000 (which is equivalent to the
		// projected value being larger than
		// 	cos( 2*pi/10000 ) = 0.999999802607918431)
		// then we regard the normal as being aligned to the particular axis.
		if ( abs( static_cast<long double>(outward[max_d])) > 0.999999802607918431 ) {
			axis_aligned = true;
		} else {
			// Valgrind memcheck helped to find this was missing.
			// Doesn't solve problem of uninitialized values though.
			// 2008-09-01 Mike Li.
			axis_aligned = false;
		}
	} // dim != 1
	finalized = true;
}
template <int dim>
bool
boundary_element<dim>::is_finalized() const
{
	return finalized;
}

template <>
bool
boundary_element<1>::closed_intersect_box ( const box<1>& bx, box_intersect_line_scratch & scratch ) const
{
	return bx.closed_intersect_point( this->endpoints );
}

template <>
bool
boundary_element<2>::closed_intersect_box ( const box<2>& bx, box_intersect_line_scratch & scratch ) const
{
	return bx.closed_intersect_line( endpoints, scratch );
}

template <>
bool
boundary_element<1>::open_intersect_box ( const box<1>& bx, box_intersect_line_scratch & scratch ) const
{
	return bx.open_intersect_point( this->endpoints );
}

template <>
bool
boundary_element<2>::open_intersect_box ( const box<2>& bx, box_intersect_line_scratch & scratch ) const
{
	return bx.open_intersect_line( endpoints, scratch );
}


template <>
bool
boundary_element<3>::closed_intersect_box ( const box<3>& bx, box_intersect_line_scratch & scratch ) const
{
	assert( !"Not coded." ); // TODO: code this. Placeholder to allow compilation for 3D.
	return false;
}

template <>
bool
boundary_element<3>::open_intersect_box ( const box<3>& bx, box_intersect_line_scratch & scratch ) const
{
	assert( !"Not coded." ); // TODO: code this. Placeholder to allow compilation for 3D.
	return false;
}


template <>
void
boundary_element<1>::get_bounding_box ( box<1>& bx ) const
{
	assert( finalized );
	bx.set( &endpoints[0] );
}

template <>
void
boundary_element<2>::get_bounding_box ( box<2>& bx ) const
{
	assert( finalized );
	double tmp[4];
	
	// x-direction
	tmp[0] = min( endpoints[0], endpoints[2] );
	tmp[1] = max( endpoints[0], endpoints[2] );
	
	// y-direction
	tmp[2] = min( endpoints[1], endpoints[3] );
	tmp[3] = max( endpoints[1], endpoints[3] );
	
	bx.set( &tmp[0] );
}

template <>
void
boundary_element<3>::get_bounding_box ( box<3>& bx ) const
{
	assert( finalized );
	double tmp[6];
	// x-direction
	tmp[0] = min( endpoints[0], min( endpoints[3], endpoints[6] ) );
	tmp[1] = max( endpoints[0], max( endpoints[3], endpoints[6] ) );
	// y-direction
	tmp[2] = min( endpoints[1], min( endpoints[4], endpoints[7] ) );
	tmp[3] = max( endpoints[1], max( endpoints[4], endpoints[7] ) );
	// z-direction
	tmp[4] = min( endpoints[2], min( endpoints[5], endpoints[8] ) );
	tmp[5] = max( endpoints[2], max( endpoints[5], endpoints[8] ) );
	
	bx.set( &tmp[0] );
}

// Explicit template initialization.

//template class boundary_element<1>;
template class boundary_element<2>;
//template class boundary_element<3>;

