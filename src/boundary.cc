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

#include "boundary.hh"
#include "boundary_element_utils.hh"

#include "box_utils.hh"

#include <iostream>
#include <cassert>

using std::cout;

template<int dim>
boundary<dim>::boundary() : updated( false ), set_vals(0)
{
}

template<int dim>
boundary<dim>::~boundary()
{
}

template<int dim>
boundary<dim>::boundary( const boundary<dim>& other )
{
	deep_copy( other );
}

template<int dim>
boundary<dim>&
boundary<dim>::operator= ( const boundary<dim>& other )
{
	deep_copy( other );
	return *this;
}

template<int dim>
void
boundary<dim>::deep_copy( const boundary<dim>& other )
{
	updated      = other.updated;
	bounding_box = other.bounding_box;
	elems        = other.elems;
	
	set_vals.resize( other.set_vals.size() );
	set_vals     = other.set_vals;
}

template<int dim>
void
boundary<dim>::simple_set( valarray<double>& sim )
{
	set_vals.resize( sim.size() );
	set_vals = sim;
	assert( set_vals.size() == sim.size() );
	
	// TODO: remove. Debug point in poly.
#if 0
	cout << "Saved points \n";
	for ( int i=0; i<set_vals.size(); ++i ) {
		cout << set_vals[i] << ", ";
	}
	cout << "\b\b\n";
#endif
	
	if ( dim==1 ) {
		assert( sim.size() >= 4  );
		elems.resize(2);

		elems[0].set_oriented_simplex( &sim[0] );
		elems[1].set_oriented_simplex( &sim[2] );

	} else if ( dim==2 ) {
		assert( sim.size() >= 6 );
		double pts[4];
		assert( sim.size()%dim==0 );
		int num_pts = sim.size()/dim;
		elems.resize( num_pts );
		for( int i=0; i<num_pts; ++i ) {
			pts[0] = sim[2*i];
			pts[1] = sim[2*i+1];
			if ( i != num_pts-1 ) {
				assert ( 2*(i+1)+1 < sim.size() );
				
				pts[2] = sim[2*(i+1)];
				pts[3] = sim[2*(i+1)+1];
			} else {
				pts[2] = sim[0];
				pts[3] = sim[1];
			}
			elems[i].set_oriented_simplex( &pts[0] );
		}
		assert( elems.size() == num_pts );
	} else {
		assert( dim == 3 );
		assert( !"Not implemented yet." );
	}
}

#if 0
template<int dim>
void
boundary<dim>::set_simplicies ( const valarray<double>& sim )
{
	assert( sim.size() > 0 );
	
	int doubles_per_sim = ( dim == 1 ? 2 : dim*dim ); // Number of doubles needed to define a simplex.
	
	int size = sim.size();
	
	assert( size%doubles_per_sim == 0 );
	
	int num_simplices = size / doubles_per_sim;
	
	elems.resize( num_simplices );
	
	for ( int i=0; i<num_simplices; ++i ) {
		elems[i].set_oriented_simplex( &sim[ i*doubles_per_sim ] );
	}
	
	set_vals = sim;
}
#endif

template<int dim>
void
boundary<dim>::update()
{
	assert( elems.size() > 0 );
	
	vector<box<dim> > tmp_vec( elems.size() );
	typename vector<box<dim> >::iterator tmp_it( tmp_vec.begin() );
	
	typename vector<boundary_element<dim> >::iterator it( elems.begin() );
	typename vector<boundary_element<dim> >::iterator end( elems.end() );
	for( ; it!=end; ++it, ++tmp_it ) {
		it->finalize();
		it->get_bounding_box( *tmp_it );
	}
	
	make_bounding_box( tmp_vec, bounding_box );
	
	updated = true;
}

// The bounding box is formed in update() so that
// this get member function can be const. This is
// logical as access should not require modification.
template<int dim>
void
boundary<dim>::get_bounding_box ( box<dim>& bx ) const
{
	assert( updated );
	bx = bounding_box;
}


template<int dim>
bool
boundary<dim>::closed_intersect_box ( const box<dim>& bx ) const
{

	size_t elems_size = elems.size();
	assert ( elems_size > 0 );
	
	assert ( updated );
	
	if ( !bounding_box.closed_intersect_box( bx ) ) {
		return false;
	}
	
	box_intersect_line_scratch    scratch;
	
	for ( int i=0; i<elems_size; ++i ) {
		assert ( elems[i].is_finalized() );
		if ( elems[i].closed_intersect_box( bx, scratch ) ) {
			return true;
		}
	}
	return false;
#if 0

	/*	cout << "Box is ";
	output_description( bx, cout );
	cout << "There are " << elems.size() << " boundary elements : ";
	gp_draw ( elems, cout );
	*/
	
	assert( elems.size() > 0 );
	typename vector<boundary_element<dim> >::const_iterator it( elems.begin() );
	typename vector<boundary_element<dim> >::const_iterator end( elems.end() );
	
	assert ( it != end );
	
	int count = 0;
	
	while ( it!=end ) {
		assert ( it->is_finalized() );
		//cout << " = = = = = = = " << count << "\n";
		assert ( &elems[count] == &(*it) );
		
		assert ( elems[count].closed_intersect_box(bx) == it->closed_intersect_box( bx ) );
		if ( it->closed_intersect_box( bx ) ) {
			return true;
		}
		++it;
		++count;
	}
	return false;
#endif
}

template<int dim>
bool
boundary<dim>::open_intersect_box ( const box<dim>& bx ) const
{
	assert( elems.size() > 0 );
	
	assert ( updated );
	
	if ( !bounding_box.open_intersect_box( bx ) ) {
		return false;
	}
	
	box_intersect_line_scratch    scratch;
	
	typename vector<boundary_element<dim> >::const_iterator it( elems.begin() );
	typename vector<boundary_element<dim> >::const_iterator end( elems.end() );
	while ( it!=end ) {
		if ( it->open_intersect_box( bx, scratch ) ) {
			return true;
		}
		++it;
	}
	return false;
}

template <>
bool
boundary<1>::open_intersect_point ( const double* co ) const
{
	assert( elems.size() > 0 );
	// Assumes the boundary is made up of two ends of
	// an interval.
	// TODO: allow for boundary made up of multiple ends of
	// intervals.
	return ( set_vals[0] < co[0] && co[0] < set_vals[2] );
}

template <>
bool
boundary<2>::open_intersect_point ( const double* co ) const
{
	assert( co );
	assert( elems.size() > 0 );
	
	assert ( updated );

	if ( !bounding_box.open_intersect_point( co ) ) {
		return false;
	}

	
	int num_pts = set_vals.size()/2;
	valarray<double> in_x( num_pts );
	valarray<double> in_y( num_pts );
	
	assert( set_vals.size() > 0 );
	
	for ( int i=0; i<num_pts; ++i ) {
		in_x[i] = set_vals[ 2*i   ];
		in_y[i] = set_vals[ 2*i+1 ];
	}
	return point_in_polygon( num_pts, &in_x[0], &in_y[0], co[0], co[1] );
}

template <>
bool
boundary<3>::open_intersect_point ( const double* co ) const
{
	assert( elems.size() > 0 );
	assert( !"Not implemented." );
	// TODO: this could be hard to code.
	return false;
}

template<int dim>
void boundary<dim>::gp_draw_boundary( ostream& out, double level ) const
{
	
	int size = elems.size();
	for ( int i=0; i<size; ++i ) {
		gp_draw_single( elems[i], out, level );
	}
}

//

// template class boundary<1>;
template class boundary<2>;
// template class boundary<3>;

