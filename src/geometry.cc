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

#include "geometry.hh"
#include "box_utils.hh"

#include <cassert>
#include <iostream>

using std::cout;

template<int dim>
geometry<dim>::geometry()
{
}

template<int dim>
geometry<dim>::~geometry()
{
}

template<int dim>
void
geometry<dim>::create_boundary ( const string& name, valarray<double>& sim )
{
	assert( created.find( name ) == created.end() );
	created.insert( name );
	// boundaries[ name ].set_simplices( sim );
	boundaries[name].simple_set( sim );
}

template<int dim>
void
geometry<dim>::create_region ( const string& name, const string& oriented_boundary )
{
	assert( created.find( name ) == created.end() );
	created.insert( name );
	if ( oriented_boundary == "" ) {
		// No boundary has been specified.
		regions[name]; // Add to regions map (creates an empty vector).
	} else {
		// A boundary has be specified already.
		
		// The first character must be '+' or '-' to denote
		// the inside and the outside respectively.
		assert( oriented_boundary[0] == '+' || oriented_boundary[0] == '-' );
		
		// Take the name of the boundary.
		string boundary_name( oriented_boundary, 1, string::npos );
		
		// The boundary name specified must be an existing boundary.
		assert( boundaries.find( boundary_name ) != boundaries.end() );

		// Make sure to push back the oriented boundary and not just
		// the name of the boundary.
		regions[name].push_back( oriented_boundary );
	}
}

template<int dim>
const map<string, boundary<dim> >&
geometry<dim>::access_boundaries () const
{
	assert ( boundaries.size() > 0 );
	return this->boundaries;
}

template<int dim>
const map<string, vector<string> >&
geometry<dim>::access_regions () const
{
	return this->regions;
}

template<int dim>
void
geometry<dim>::update()
{
	typename map<string, boundary<dim> >::iterator it( boundaries.begin() );
	typename map<string, boundary<dim> >::iterator end( boundaries.end() );
	
	for ( ; it!=end; ++it ) {
		it->second.update();
	}
}

template<int dim>
bool
geometry<dim>::is_boundary( const string& name ) const
{
	return boundaries.find( name ) != boundaries.end();
}

template<int dim>
bool
geometry<dim>::is_region( const string& name ) const
{
	return regions.find( name ) != regions.end();
}

template<int dim>
bool 
geometry<dim>::is_support ( const string& in ) const
{
	return is_boundary ( in ) || is_region ( in );
}

template<int dim>
bool geometry<dim>::gp_draw ( const string& name, ostream& out, double level ) const
{
	if ( is_boundary( name ) ) {
		assert( boundaries.find(name) != boundaries.end() );
		const boundary<dim>& bdry = (boundaries.find(name))->second;
		bdry.gp_draw_boundary( out, level );
		return true;
	} else {
		assert( is_region( name ) );
		string boundary_name;
		assert( regions.find(name) != regions.end() );
		const vector<string>& oriented_bdry = (regions.find(name))->second;
		
		int size = oriented_bdry.size();
		for ( int i=0; i<size; ++i ) {
			boundary_name.assign( oriented_bdry[i], 1, string::npos );
			assert( is_boundary( boundary_name ) );
			
			assert( boundaries.find(boundary_name) != boundaries.end() );
			const boundary<dim>& bdry = (boundaries.find(boundary_name))->second;
			bdry.gp_draw_boundary( out, level );
			out << "\n";	
		}
	}
	return false;
}

/**
 * Determination is based on an 'OR' relationship as written.
 */
template<int dim>
bool
geometry<dim>::inside_region( const string& region, const double* co ) const
{
	assert( regions.find( region ) != regions.end() );
	const vector<string>& oriented_boundaries = regions.find( region )->second;
	
	typename vector<string>::const_iterator it( oriented_boundaries.begin() );
	typename vector<string>::const_iterator end( oriented_boundaries.end() );
	string boundary_name;
	for ( ; it!=end; ++it ) {
		boundary_name.assign( *it, 1, string::npos );
		if ( (*it)[0] == '+' ) {
			assert( boundaries.find( boundary_name ) != boundaries.end() );
			if ( (boundaries.find( boundary_name )->second).open_intersect_point( co ) )
			{
				return true;
			}
		} else {
			assert( (*it)[0] == '-' );
			assert( boundaries.find( boundary_name ) != boundaries.end() );
			if ( !(boundaries.find( boundary_name )->second).open_intersect_point( co ) )
			{
				return true;
			}
		}
	}
	return false;
}

/**
	For viz_refinement of eight and greater, this member function becomes
	a highly significant cost.
*/
template <int dim>
void
geometry<dim>::intersect_region ( const string &              region,
                                  const vector<box<dim> > &   in,
                                  valarray<bool>&              out ) const
{
	assert ( region != "" );
	assert( regions.find( region ) != regions.end() );
	const vector<string>& oriented_boundaries = regions.find( region )->second;
	assert ( oriented_boundaries.size() > 0 );
	
	int in_size = in.size();
	
	out.resize ( in_size );
	
	typename vector<string>::const_iterator it;
	typename vector<string>::const_iterator	end ( oriented_boundaries.end() );
	
	double centre[dim];
	bool cont_outer = false;
	for ( int b=0; b<in_size; ++b ) {
		it  = oriented_boundaries.begin();
		assert ( it != end );
		
		cont_outer = false;
		
		in[b].get_centre_point ( centre );
		out[b] = inside_region( region, centre );
		
		if ( out[b] ) {
			continue;
		}
		
		// Even if the centre of the box is not inside the region, the box may still
		// intersect the boundary.
		
		string boundary_name;
		for ( ; it!=end; ++it ) {
			boundary_name.assign( *it, 1, string::npos );
			assert ( boundary_name.size() > 0 );
			
			if ( (*it)[0] == '+' ) {
				assert( boundaries.find( boundary_name ) != boundaries.end() );
				
				if ( (boundaries.find( boundary_name )->second).open_intersect_box ( in[b] ) ) {
					out[b] = true;
					cont_outer = true;
					break;
				}
				
			} else {
				assert( (*it)[0] == '-' );
				assert( boundaries.find( boundary_name ) != boundaries.end() );
				if ( !(boundaries.find( boundary_name )->second).open_intersect_box ( in[b] ) )
				{
					out[b] = true;
					cont_outer = true;
					break;
				};
				
			}
		}
		if ( cont_outer ) {
			continue;
		}
	}
	
}

template<int dim>
void
geometry<dim>::get_bounding_box ( box<dim>& bounding_box )
{
	assert( boundaries.size() > 0 );
	
	vector<box<dim> > tmp_vec( boundaries.size() );
	typename vector<box<dim> >::iterator tmp_it( tmp_vec.begin() );
	
	typename map<string, boundary<dim> >::iterator it( boundaries.begin() );
	typename map<string, boundary<dim> >::iterator end( boundaries.end() );
	for( ; it!=end; ++it, ++tmp_it ) {
		(it->second).get_bounding_box( *tmp_it );
	}
	
	make_bounding_box( tmp_vec, bounding_box );
}

//

// template class geometry<1>;
template class geometry<2>;
// template class geometry<3>;
