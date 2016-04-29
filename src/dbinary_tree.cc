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

#include "dbinary_tree.hh"
#include "dbinary_tree_utils.hh"
#include "box_utils.hh"

#include <petscksp.h>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <cassert>
#include <cmath>
#include <cstring>
#include <valarray>

using std::inserter;
using std::back_inserter;
using std::copy;
using std::cout;
using std::insert_iterator;
using std::merge;
using std::pow;
using std::set_difference;
using std::set_intersection;
using std::set_union;
using std::strcpy;
using std::valarray;

#if 0
template<int dim>
int key_level( const string& k )
{
	return (k.size()-1)/2;
}
#endif

template<int dim>
dbinary_tree<dim>::dbinary_tree ()
  : updated_neighbours(false),
    updated( false ),
    updated_global_boxes( false ),
    have_parent_box( false ),
    cover_factor( 1.2 )
{
	set_max_level( 100 );
	present_tree.insert( "0" );
}

template<int dim>
dbinary_tree<dim>::~dbinary_tree() { }

template<int dim>
void
dbinary_tree<dim>::set_geometry ( shared_ptr<geometry<dim> >& g, double expand )
{
	assert( g );
	geom = g;
	
	box<dim> bbox;
	geom->get_bounding_box( bbox );
	
	assert( expand > 0.0 );
	
	double scale[dim];
	for ( int d=0; d<dim; ++d ) {
		scale[d] = expand;
	}
	
	bbox.scale( scale );
	
	set_parent_box( bbox );
	
	updated = false;
}

template<int dim>
const geometry<dim>&
dbinary_tree<dim>::access_geometry () const
{
	assert( geom );
	return *geom;
}

template<int dim>
const set<string>&
dbinary_tree<dim>::access_region_boundary_keys( const string& region ) const
{
	assert( updated );
	assert( present_region_boundary_keys.find( region ) != present_region_boundary_keys.end() );
	return ( present_region_boundary_keys.find( region ) )->second;
}

template<int dim>
const set<string>&
dbinary_tree<dim>::access_region_interior_keys( const string& region ) const
{
	assert( updated );
	assert( present_region_interior_keys.find( region ) != present_region_interior_keys.end() );
	return ( present_region_interior_keys.find( region ) )->second;
}
template<int dim>
const set<string>&
dbinary_tree<dim>::access_region_all_keys( const string& region ) const
{
	assert( updated );
	assert( present_region_all_keys.find( region ) != present_region_all_keys.end() );
	return ( present_region_all_keys.find( region ) )->second;
}

template<int dim>
const set<string>&
dbinary_tree<dim>::access_boundary_keys ( const string& boundary ) const
{
	assert( updated );
	assert( geom->is_boundary( boundary ) );
	assert( present_boundary_keys.find( boundary ) != present_boundary_keys.end() );
	return ( present_boundary_keys.find( boundary ) )->second;
}

template<int dim>
const box<dim>&
dbinary_tree<dim>::access_parent_box() const
{
	return parent_box;
}

template<int dim>
string
dbinary_tree<dim>::find_patch_containing_point( const string& region, const double* co ) const
{
	assert( updated );
	if( !geom->inside_region( region, co ) ) {
		// Point is not inside region.
		return "";
	} else {
		box<dim> tmp_box;
		{
			assert( present_region_boundary_keys.find( region ) != present_region_boundary_keys.end() );
			const set<string>& boundary_keys = present_region_boundary_keys.find( region )->second;
			typename set<string>::const_iterator it( boundary_keys.begin() );
			typename set<string>::const_iterator end( boundary_keys.end() );
			while ( it!=end ) {
				get_global_box( *it, tmp_box );
				if ( tmp_box.open_intersect_point( co ) ) {
					return *it;
				}
				++it;
			}
		}
		{
			assert( present_region_interior_keys.find( region ) != present_region_interior_keys.end() );

			const set<string>& interior_keys = present_region_interior_keys.find( region )->second;
			
			// TODO: remove. Debug no interior patches problem.
			// cout << "Interior_keys size is " << interior_keys.size() << "\n";
			
			typename set<string>::const_iterator it( interior_keys.begin() );
			typename set<string>::const_iterator end( interior_keys.end() );
			while ( it!=end ) {
				get_global_box( *it, tmp_box );
				if ( tmp_box.open_intersect_point( co ) ) {
					return *it;
				}
				++it;
			}
		}
		assert( !"The cover must be incomplete." );
		return "";
	}
}

// Regardless of support.
template <int dim>
void
dbinary_tree<dim>::get_patch_containing_point ( const double* co, vector<string>& out ) const
{
	assert ( parent_box.closed_intersect_point( co ) );
	
	basic_dtree->get_key_expanded_intersect_point ( initial_level, co, cover_factor, out );
}

template<int dim>
void
dbinary_tree<dim>::get_boundary_neighbours( const string& boundary, const string& key, set<string>& out ) const
{
	assert(updated_neighbours);
	assert(geom->is_boundary(boundary));
	
	typename map<string, set<string> >::const_iterator tmp_it(present_neighbours.find( key ));
	assert( tmp_it!=present_neighbours.end() );
	const set<string>& all_neighbours = tmp_it->second;
	
	tmp_it = present_boundary_keys.find( boundary );
	assert( tmp_it!=present_boundary_keys.end() );
	
	const set<string>& boundary_keys=tmp_it->second;
	out.clear();
	set_intersection( boundary_keys.begin(), boundary_keys.end(),
			  all_neighbours.begin(), all_neighbours.end(),
			  inserter( out, out.end() )
			);
}

template<int dim>
void
dbinary_tree<dim>::get_region_neighbours( const string& region, const string& key, set<string>& res ) const
{
	assert( updated );
	assert( updated_neighbours );
	assert( geom->is_region( region ) );
	
	assert( present_neighbours.find( key ) != present_neighbours.end() );
	const set<string>& all_neighbours = (present_neighbours.find( key ))->second;
	
	// gdb helped find that `key' was incorrectly used as `region' here.
	assert( present_region_boundary_keys.find(region) != present_region_boundary_keys.end() );
	const set<string>& boundary_keys = ( present_region_boundary_keys.find(region) )->second;
	
	// gdb helped find that `key' was incorrectly used as `region' here.
	assert( present_region_interior_keys.find(region) != present_region_interior_keys.end() );
	const set<string>& interior_keys = ( present_region_interior_keys.find(region) )->second;
	
	set<string> keys_covering_region;
	merge( boundary_keys.begin(), boundary_keys.end(),
	       interior_keys.begin(), interior_keys.end(),
	       inserter( keys_covering_region, keys_covering_region.end() )
	     );
	
	res.clear();
	set_intersection( all_neighbours.begin(), all_neighbours.end(),
			  keys_covering_region.begin(), keys_covering_region.end(),
			  inserter( res, res.end() )
			);
}


template<int dim>
void
dbinary_tree<dim>::get_region_neighbours_incident_point ( const string& region, const string& key, const double* co, set<string>& out ) const
{
	set<string> region_neighbours;
	get_region_neighbours( region, key, region_neighbours );
	out.clear();
	
	box<dim> tmp_box;
	
	typename set<string>::iterator it( region_neighbours.begin() );
	typename set<string>::iterator end( region_neighbours.end() );
	for ( ; it!=end; ++it ) {
		
		// GDB helped find that key was incorrectly substituted
		// for *it so that boxes were returned that were not
		// incident on the point.
		get_global_box( *it, tmp_box );
		if ( tmp_box.open_intersect_point( co ) ) {
			out.insert( *it );
		}
	}
	
}

template<int dim>
void
dbinary_tree<dim>::get_patch_keys( const string& name, set<string>& res ) const
{
	// Name is region or boundary name.
	res.clear();
	if ( geom->is_boundary( name ) ) {
		res = (present_boundary_keys.find( name ))->second;
	} else {
		assert( geom->is_region( name ) );
		res = (present_region_all_keys.find(name))->second;
	}
}

template <int dim>
void
dbinary_tree<dim>::get_neighbour_keys( const string& name, map<string, set<string> >& out ) const
{
	const set<string>& keys = access_patch_keys( name );
	
	out.clear();
	
	set<string>::const_iterator begin( keys.begin() );
	set<string>::const_iterator it( keys.begin() ), end( keys.end() );
	
	/**
		2009-01-22 The key neighbours will no longer take the domain
		cover into account at this part of the code. It was found using
		valgrind/callgrind profiling that the use of set_intersection
		was a massive problem here around levels 7 and 8.
		
		The index-neighbours will be correct for the cover when
		the cover_structure is created.
		
		At level 7 with -g running through valgrind/callgrind,
		get_cover_structure drops from 58.93% to 16.57% of cost which
		is a lot to do with wiping out the 51.20% of cost of set_intersection
		at this level.
		
		generate_index_information must be used to place the key-neighbours
		back in step and this is required before some of the visualization
		routines.
	*/
	for ( ; it!=end; ++it ) {
		assert ( present_neighbours.find( *it ) != present_neighbours.end() );
		const set<string>& present_all_neigh = (present_neighbours.find( *it ))->second;
		
		set<string> &    out_set = out[*it];
		
		out_set.clear();
		out_set = present_all_neigh;
	}
#if 0
	for ( ; it!=end; ++it ) {
		assert ( present_neighbours.find( *it ) != present_neighbours.end() );
		const set<string>& all_neigh = (present_neighbours.find( *it ))->second;

		// This command does not scale beyond level 8.
		// 2009-01-21. ML.
		// Test remove.
#if 1
		set_intersection( begin,
		                  end,
		                  all_neigh.begin(),
		                  all_neigh.end(),
		                  inserter( out[*it], out[*it].end() )
		                );
#else
		copy ( all_neigh.begin(), all_neigh.end(),
		       inserter(out[*it], out[*it].end() ) );
		//out[*it] = (present_neighbours.find( *it ))->second;
#endif
	}
#endif
}


template<int dim>
const set<string>&
dbinary_tree<dim>::access_patch_keys( const string& region_or_boundary ) const
{
	if ( geom->is_boundary(region_or_boundary) ) {
		assert( present_boundary_keys.find(region_or_boundary)!=present_boundary_keys.end() );
		return (present_boundary_keys.find(region_or_boundary))->second;
	} else {
		assert( geom->is_region(region_or_boundary) );
		assert( present_region_all_keys.find(region_or_boundary)!=present_region_all_keys.end() );
		return (present_region_all_keys.find(region_or_boundary))->second;
	}
}

template<int dim>
int
dbinary_tree<dim>::get_num_patches( const string& region_or_boundary ) const
{
	if ( geom->is_boundary( region_or_boundary ) ) {
		return ((present_boundary_keys.find( region_or_boundary ))->second).size();
	} else {
		assert( geom->is_region( region_or_boundary ) );
		return ((present_region_all_keys.find(region_or_boundary))->second).size();
	}	
}

template<int dim>
bool
dbinary_tree<dim>::support_has_key( const string& name, const string& key ) const
{
	if ( geom->is_boundary( name ) ) {
		const set<string>& tmp_set = (present_boundary_keys.find( name ))->second;
		return tmp_set.find( key ) != tmp_set.end();
	} else {
		assert( geom->is_region( name ) );
		
		const set<string>& tmp_set = (present_region_all_keys.find(name))->second;

		return tmp_set.find( key ) != tmp_set.end();
	}
}

template<int dim>
void
dbinary_tree<dim>::get_patches( const set<string>& in, vector<box<dim> >& res ) const
{
	// TODO: more asserts here.
	assert( updated );
	// assert( in.size() > 0 );
	
	res.resize( in.size() );
	
	typename vector<box<dim> >::iterator res_it( res.begin() );
	
	typename set<string>::const_iterator it( in.begin() );
	typename set<string>::const_iterator end( in.end() );
	for ( ; it!=end; ++it, ++res_it ) {
		get_global_box( *it, *res_it );
	}
}

template<int dim>
void
dbinary_tree<dim>::get_region_patches( const string& region, vector<box<dim> >& res ) const
{
	assert( updated );
	assert( geom->is_region( region ) );
	assert( present_region_all_keys.find(region) != present_region_all_keys.end() );
	get_patches( (present_region_all_keys.find(region))->second, res );
}

template<int dim>
void
dbinary_tree<dim>::get_region_boundary_patches ( const string& region, vector<box<dim> >& out ) const
{
	assert( updated );
	assert( geom->is_region( region ) );
	get_patches( access_region_boundary_keys(region), out );
}
  
template<int dim>
void
dbinary_tree<dim>::get_region_interior_patches ( const string& region, vector<box<dim> >& out ) const
{
	assert( updated );
	assert( geom->is_region( region ) );
	get_patches( access_region_interior_keys(region), out );
}
  
template<int dim>
void
dbinary_tree<dim>::get_boundary_patches( const string& boundary, vector<box<dim> >& res ) const
{
	assert( updated );
	assert( geom->is_boundary( boundary ) );
	get_patches( access_boundary_keys(boundary), res );
}

template<int dim>
void
dbinary_tree<dim>::set_parent_box( const box<dim>& par_box )
{
	// global_boxes.clear();
	// updated_global_boxes = false;
	updated = false;
	// global_boxes["0"] = par_box;
	parent_box = par_box;
	have_parent_box = true;
}

template<int dim>
bool
dbinary_tree<dim>::is_updated() const
{
	return updated;
}

template<int dim>
void
dbinary_tree<dim>::update_neighbours()
{
	// assert( updated_global_boxes );
	
	typename set<string>::const_iterator it(present_tree.begin());
	typename set<string>::const_iterator end(present_tree.end());
	typename set<string>::const_iterator inner_it;

	box<dim> outer_box, inner_box;
	
	for ( ; it!=end; ++it ) {
		
		
		get_global_box( *it, outer_box );
		
		// access is expensive for large covers.
		// const box<dim>& outer_box = access_global_box(*it);
		
		// Based on profiling have changed how the inner loop works to avoid
		// dereferencing two iterators and comparison between strings.
		//for ( inner_it = present_tree.begin(); *inner_it < *it; ++inner_it ) {
		
		for ( inner_it=it; inner_it!=end; ++inner_it ) {
			
			if ( inner_it == it ) {
				// Each patch is a neighbour of itself for simplicity
				// of handling with assembly.
				present_neighbours[*it].insert(*it);
				continue;
			}
			
			
			get_global_box( *inner_it, inner_box );
			
			// access is expensive for large covers.
			// const box<dim>& inner_box = access_global_box(*inner_it);
			
			if ( outer_box.open_intersect_box( inner_box ) ) {
				present_neighbours[ *it ].insert( *inner_it );
				present_neighbours[ *inner_it ].insert( *it );
			}
		}
	}
	updated_neighbours = true;
}

template<int dim>
void
dbinary_tree<dim>::update()
{
	assert( have_parent_box );
	assert( present_tree.size() > 0 );
	
	if ( !basic_dtree ) {
		basic_dtree.reset ( new basic_dbinary_tree<dim> );
	}
	
	basic_dtree->initialize ( parent_box, max_level );
	basic_dtree->get_neighbours ( present_tree, present_neighbours );
	updated_neighbours = true;
	
// Caching of created boxes is inefficient for large covers.

	PetscLogDouble log_t0, log_t1;
	
	PetscGetTime ( &log_t0 );
	// Create boundary keys.
	{
		present_boundary_keys.clear();
		
		const map<string, boundary<dim> >& boundaries = geom->access_boundaries();
		
		assert ( boundaries.size() > 0 );
		
		typename map<string, boundary<dim> >::const_iterator it( boundaries.begin() );
		typename map<string, boundary<dim> >::const_iterator end( boundaries.end() );

		box<dim> tmp_box;
		
		for ( ; it!=end; ++it ) {
			set<string>& keys_for_boundary = present_boundary_keys[ it->first ];
			
			typename set<string>::const_iterator tree_it( present_tree.begin() );
			typename set<string>::const_iterator tree_end( present_tree.end() );
			
			assert ( present_tree.size() > 0 );
			
			for ( ; tree_it!=tree_end; ++tree_it ) {
				get_global_box( *tree_it, tmp_box );
				assert ( !tmp_box.empty() );
				if ( (it->second).closed_intersect_box( tmp_box ) ) {
					keys_for_boundary.insert( *tree_it );
				}
			}
		}
	}
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Updating of boundary keys took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscGetTime ( &log_t0 );
	// Create interior keys.
	{
		present_interior_keys.clear();
		
		set<string> all_interior_keys;
		
		typename map<string, set<string> >::const_iterator it( present_boundary_keys.begin() );
		typename map<string, set<string> >::const_iterator end( present_boundary_keys.end() );
		
		box<dim> test_box;
		double centre[dim];
		
		for ( ; it!=end; ++it ) {
			const set<string>& boundary_keys = it->second;
			
			// Make a new set of keys that represent patches
			// that do not intersect the boundary of interest.
			// We get a set that we know to be interior or
			// exterior to the boundary of interest.
			set_difference( present_tree.begin(), present_tree.end(),
					boundary_keys.begin(), boundary_keys.end(),
					inserter( all_interior_keys, all_interior_keys.end() )
				      );
			// Allow.
			// assert( all_interior_keys.size() > 0 && "Every patch is a boundary patch." );
			
			// Get access to the boundaries so we can check if a particular
			// point is inside or outside a particular boundary.
			const map<string, boundary<dim> >& boundaries = geom->access_boundaries();
			
			typename set<string>::const_iterator sep_it( all_interior_keys.begin() );
			typename set<string>::const_iterator sep_end( all_interior_keys.end() );
			
			// Do each non-boundary key one at a time.
			for ( ; sep_it!=sep_end; ++sep_it ) {
				
				get_global_box( *sep_it, test_box );
				test_box.get_centre_point( centre );
				
				// TODO: think this can be removed.
				//cout << "Test box " << *sep_it << " is ";
				//gp_draw_single( test_box, cout );
				
				// Separate here.
				assert( boundaries.find( it->first ) != boundaries.end() );
				if ( (boundaries.find( it->first )->second).open_intersect_point( centre ) ) {
					string interior_name( "+" );
					interior_name.append( it->first );
					
					// TODO: remove. Debugging no interior patches at all.
					// cout << "Appending something to " << interior_name << "\n";
					
					present_interior_keys[ interior_name ].insert( *sep_it );
				} else {
					string interior_name( "-" );
					interior_name.append( it->first );
					present_interior_keys[ interior_name ].insert( *sep_it );
				}
			}
		}
		// Allow.
		// assert( present_interior_keys.size() > 0 );
		
	}
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Updating of interior keys took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscGetTime ( &log_t0 );
	// Divide by regions.
	{
		present_region_boundary_keys.clear();
		present_region_interior_keys.clear();
		present_region_all_keys.clear();
		
		const map<string, vector<string> >& region = geom->access_regions();
		
		typename map<string, vector<string> >::const_iterator reg_it( region.begin() );
		typename map<string, vector<string> >::const_iterator reg_end( region.end() );
		
		string region_name, interior_name, boundary_name;
		set<string> boundary_set;
		set<string> interior_set;
		
		for ( ; reg_it!=reg_end; ++reg_it ) {
			region_name = reg_it->first;
			typename vector<string>::const_iterator bound_it( reg_it->second.begin() );
			typename vector<string>::const_iterator bound_end( reg_it->second.end() );
			
			// gdb helped identify missing for enclosing statement.
			for ( ; bound_it!=bound_end; ++bound_it ) {
				
				interior_name = *bound_it;
				
				// Make sure we are getting the +/- in front of the boundary name.
				assert( interior_name[0]=='-' || interior_name[0]=='+' );
				
				// Don't take the '+' or '-' from the front
				// and we get the boundary name.
				boundary_name.assign( interior_name, 1, string::npos );
				
				// Won't exist on first run through.assert( 
				// present_boundary_keys.find( boundary_name ) != present_boundary_keys.end() );
				// Boundaries bring all their covering patches with them hence
				// we take the union of the relevant keys.
				{
					set<string>& union_set = boundary_set;
					union_set.clear();
					set_union( present_boundary_keys[ boundary_name ].begin(),
						present_boundary_keys[ boundary_name ].end(),
						present_region_boundary_keys[ region_name ].begin(),
						present_region_boundary_keys[ region_name ].end(),
						inserter( union_set, union_set.end() )
						);
					present_region_boundary_keys[ region_name ] = union_set;
					assert( union_set.size() > 0 );
				}
				
				// Won't exist on first run through.
				//assert( present_interior_keys.find( interior_name ) != present_interior_keys.end() );
				
				// Again the union of interior keys unless we're doing something silly
				// with our boundary inclusions/exclusions.
				if ( present_interior_keys.find( interior_name ) == present_interior_keys.end() ) {
					// There are no interior patches for this boundary.
					present_region_interior_keys[ region_name ];
				} else {
					set<string>& union_set = interior_set;
					union_set.clear();
					
					
					assert( present_interior_keys.find( interior_name ) != present_interior_keys.end() );
					
					set_union( present_interior_keys[ interior_name ].begin(),
						present_interior_keys[ interior_name ].end(),
						present_region_interior_keys[ region_name ].begin(),
						present_region_interior_keys[ region_name ].end(),
						inserter( union_set, union_set.end() )
						);
					present_region_interior_keys[ region_name ] = union_set;
					assert( union_set.size() > 0 && "No interior patches for a particular boundary. Either the patch size is too small or something strange is happening." );
				}
				
				merge( boundary_set.begin(), boundary_set.end(),
				interior_set.begin(), interior_set.end(),
				inserter( present_region_all_keys[ region_name ],
						present_region_all_keys[ region_name ].begin() )
				);
			     }
		}
	}
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Updating of region keys took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	if ( !updated_neighbours ) {
		PetscGetTime ( &log_t0 );
		update_neighbours();
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "Updating of neighbours took %f minutes.\n", (log_t1-log_t0)/60.0 );
	}
	
	updated = true;
}

template<int dim>
void
dbinary_tree<dim>::set_initial_refinement( int level )
{
	assert( level >= 0 );
	
	initial_level = level;
	
	present_tree.clear();
#if 0
	present_tree.insert( "0" );
	
	stack<string> working_stack;
	
	set<string>::iterator it( present_tree.begin() ), end( present_tree.end() );
	for ( ; it!=end; ++it ) {
		working_stack.push( *it );
	}
	
	present_tree.clear();
	string tmp;
	while ( !working_stack.empty() ) {
		tmp = working_stack.top();
		working_stack.pop();
		if ( key_level<dim>( tmp ) == level ) {
			present_tree.insert( tmp );
		} else {
			split_dbinary_node_and_push( tmp, working_stack );
		}
	}
#endif

	vector<string> tmp;

	generate_keys<dim>( level, tmp );

	copy( tmp.begin(), tmp.end(), inserter( present_tree, present_tree.begin() ) );
	
	updated_global_boxes = false;
	updated = false;
}

template<int dim>
void
dbinary_tree<dim>::set_cover_factor( double cf )
{
	assert( cf > 0 && cf < 5.0);
	
	// global_boxes.clear(); // Global boxes will be invalidated.
	// updated_global_boxes = false;
	
	cover_factor = cf;
}

template<int dim>
void
dbinary_tree<dim>::dump_cout( int num_lines ) const
{
	assert( updated );
	cout << "dbinary_tree:dump_cout()\n";
	int c=0;
	{
		set<string>::const_iterator it( present_tree.begin() );
		set<string>::const_iterator end( present_tree.end() );
		for ( ; it!=end && c<num_lines; ++c, ++it ) {
			cout << *it << ", ";
		}
		if ( present_tree.size() <= num_lines ) {
			// Output complete.
			cout << "\b\b. OUTPUT COMPLETE\n";
		} else {
			// Output incomplete.
			cout << "...\n";
		}
		cout << "Number of patches " << present_tree.size() << "\n\n";
	}
	{
		cout << "present_region_all_keys by region\n";
		typename map<string, set<string> >::const_iterator r_it( present_region_all_keys.begin() );
		typename map<string, set<string> >::const_iterator r_end( present_region_all_keys.end() );
		for ( ; r_it!=r_end; ++r_it ) {
			cout << "Species " << r_it->first << ".\n";
			typename set<string>::const_iterator k_it( (r_it->second).begin() );
			typename set<string>::const_iterator k_end( (r_it->second).end() );
			c=0;
			for ( ; k_it!=k_end && c < num_lines; ++k_it, ++c ) {
				cout << *k_it << ", ";
			}
			if ( (r_it->second).size()  <= num_lines ) {
				// Output complete.
				cout << "\b\b. OUTPUT COMPLETE\n";
			} else {
				// Output incomplete.
					cout << "...\n";
			}
			cout << "Number of patches in cover " << (r_it->second).size() << ".\n";
		}
		cout << "\n";
	}
	{
		// Dump of neighbours.
		cout << "dbinary_tree dump of present_neighbours.\n";
		typename map<string, set<string> >::const_iterator k_it( present_neighbours.begin() );
		typename map<string, set<string> >::const_iterator k_end( present_neighbours.end() );
		for ( ; k_it!=k_end; ++k_it ) {
			cout << k_it->first << " has " << (k_it->second).size() << " neighbours : ";
			typename set<string>::const_iterator it( (k_it->second).begin() );
			typename set<string>::const_iterator end( (k_it->second).end() );
			for ( ; it!=end; ++it ) {
				cout << *it << ", ";
			}
			cout << "\b\b\n";
		}
		cout << "\n";
	}
	cout << "End of dump_cout.\n";
}

template<int dim>
void
dbinary_tree<dim>::dump_gp_template_boxes( ostream& out )
{
	gp_draw( template_boxes, out );
}

template<int dim>
void
dbinary_tree<dim>::dump_gp( ostream& out ) const
{
	assert( updated );
	vector<box<dim> > vec_box( present_tree.size() );
	typename vector<box<dim> >::iterator vit( vec_box.begin() );
	
	set<string>::const_iterator it( present_tree.begin() );
	set<string>::const_iterator end( present_tree.end() );
	
	for ( ; it!=end; ++vit, ++it ) {
		get_global_box( *it, *vit );
	}
	gp_draw( vec_box, out );
}

template<int dim>
void
dbinary_tree<dim>::gp_draw_parent_box( ostream& out, double level ) const
{
	assert( have_parent_box );
	gp_draw_single( parent_box, out, level );
}

template<int dim>
void
dbinary_tree<dim>::gp_draw_boundary_cover ( const string& boundary, ostream& out, double level ) const
{
	assert( updated );
	assert( geom->is_boundary( boundary ) );
	
	vector<box<dim> > boundary_patches;
	get_boundary_patches( boundary, boundary_patches );
	gp_draw( boundary_patches, out, level );
}
 
template<int dim>
void
dbinary_tree<dim>::gp_draw_region_interior_cover ( const string& region, ostream& out, double level ) const
{
	assert( updated );
	assert( geom->is_region( region ) );
	
	vector<box<dim> > region_interior_patches;
	get_region_interior_patches( region, region_interior_patches );
	gp_draw( region_interior_patches, out, level );
}

template<int dim>
void
dbinary_tree<dim>::gp_draw_region_boundary_cover ( const string& region, ostream& out, double level ) const
{
	assert( updated );
	assert( geom->is_region( region ) );
	
	vector<box<dim> > region_boundary_patches;
	get_region_boundary_patches( region, region_boundary_patches );
	gp_draw( region_boundary_patches, out, level );
}
template<int dim>
bool
dbinary_tree<dim>::region_key_is_boundary_key ( const string& region, const string& key ) const
{
	map<string, set<string> >::const_iterator it( present_region_boundary_keys.find(region) );
	assert( it!=present_region_boundary_keys.end() );
	
	set<string>::const_iterator key_it( (it->second).find( key ) );
	return key_it!=(it->second).end();
}

template<int dim>
void
dbinary_tree<dim>::set_max_level ( int new_max_level )
{
  assert( 0<new_max_level && new_max_level < 1000 );
  
  // Checking that new_max_level > max_level is premature
  // optimization and requires that max_level be given a default
  // value in the dbinary_tree default constructor.

    powers_of_one_half.resize( new_max_level + 1);
    for ( int i= 0; i<new_max_level + 1; ++i ) {
      powers_of_one_half[i] = std::pow( 0.5, i );
    }
    powers_of_two.resize( new_max_level + 1 );
    for ( int i=0; i<new_max_level + 1; ++i ) {
      powers_of_two[i] = static_cast<int>( 1<<i ); // Gives 2^i. pow won't work here.
    }
    template_boxes.resize( new_max_level+1 );
    double ext[ 2*dim ];
    for ( int d=0; d<dim; ++d ) {
	ext[ 2*d ] = 0.0; // Lower bound in each dimension is zero for the template box.
    }
    for ( int i=0; i<new_max_level+1; ++i ) {
	double upper;
	if ( i==0 ) {
		upper = 2.0;
	} else if ( i==1 ) {
		upper = 1.0;
	} else {
		upper = powers_of_one_half[ i - 1 ];
	}
	for ( int d=0; d<dim; ++d ) {
		ext[ 2*d + 1 ] = upper;
	}
	template_boxes[i].set( &ext[0] );
    }
    max_level = new_max_level;

}

/**
	We work backwards from the end of the key (e.g. 0-1-1-1-0-0-0).
	The level `l' determines the finest increment/decrement 0.5^l.
	The value for level zero does not cause an increment/decrement.
	For l-1>0, the increment/decrement is twice that of the level `l'.
	A '0' causes a decrement whereas a '1' causes an increment of
	the corresponding multiple of 0.5^l.
*/
template<int dim>
double
dbinary_tree<dim>::centre_for_projected_key ( const string& pkey ) const
{
	int multiple = 0;
	int lev = key_level<dim>( pkey );
	for ( int i=lev; i>0; --i ) {
		if ( pkey[2*i] == '0' ) {
			multiple -= powers_of_two[lev-i];
		} else {
			multiple += powers_of_two[lev-i];
		}
//		cout << "powers_of_two "<<powers_of_two[lev-i]<<"\n";
	}
	double tmp = multiple*powers_of_one_half[lev];
	assert( -1.0<= tmp && tmp <= 1.0 );
//	cout << "claimed centre is "<<tmp<<" (multiple "<<multiple<<", powers_of_one_half "
//			<< powers_of_one_half[lev] <<")\n";
	return tmp;
}

#if 0
template<int dim>
void
dbinary_tree<dim>::create_global_box ( const string& key, box<dim>& out ) const
{
	basic_dtree->get_box ( key, out );
	out.scale ( cover_factor );
	

	out = template_boxes[ key_level<dim>( key ) ];
	string projected_keys[dim];
	for ( int d=0; d<dim; ++d ) {
		project_key( key, d, projected_keys[d] );
//		cout << "Projected key is " << projected_keys[d] << "\n";
	}
	double centre[dim];
	for ( int d=0; d<dim; ++d ) {
		centre[d] = centre_for_projected_key( projected_keys[d] );
	}
	
	double scale[dim];
	
	// Partial derivative of global coordinates written as
	// functions of local coordinates. As the coordinate transformation
	// from a box to a box is affine, these partial derivatives are constants
	// of the form (b_i - a_i)/2 for coordinate i.
	parent_box.global_partial_local( &scale[0] );
	
	for ( int d=0; d<dim; ++d ) {
		scale[d] *= cover_factor;
	}
	out.scale( &scale[0] );
	
	double global_centre[dim];
	parent_box.map_local_to_global( &centre[0], &global_centre[0] );
	
	out.recentre( &global_centre[0] );

}
#endif

#if 0
template <>
void
dbinary_tree<1>::project_key ( const string& in, int coord, string& out ) const
{
	assert( coord == 0 );
	
	out = in;
}

template <>
void
dbinary_tree<2>::project_key ( const string& in, int coord, string& out ) const
{
	assert( 0<=coord && coord<2 );
	
	int size = in.size();
	
	valarray<char> tmp(size);
	{
		const char* c_str = in.c_str();
		strcpy ( &tmp[0], c_str );
	}
	
	// Regression caused by two things. Last entry of key was not projected.
	// Key character equal to '0' was unhandled. 2008-08-20 Mike Li.
	
	// Make sure we reach the last entry of key.
	int end = (size-1)/2 + 1;
	
	if ( coord == 0 ) {
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			
			// Regression caused by missing c == '0' check. 2008-08-20 Mike Li.
			
			if ( c == '0' || c == '2' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}
	} else {
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			
			// Regression caused by missing c == '0' check. 2008-08-20 Mike Li.
			
			if ( c == '0' || c == '1' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}	
	}
	out.assign ( &tmp[0], size );
	
#if 0 //ndef NDEBUG // Keep for future debug.
	std::clog << "\t\t\t" << in << " projected to coordinate " << coord << " is " << out << "\n";
#endif
}

template <>
void
dbinary_tree<3>::project_key ( const string& in, int coord, string& out ) const
{
	assert( 0<=coord && coord<3 );
	
	int size = in.size();
	
	char tmp[size];
	{
		const char* c_str = in.c_str();
		strcpy ( tmp, c_str );
	}
	
	// Make sure we reach last entry of key.
	int end = (size-1)/2 + 1;
	
	if ( coord == 0 ) {
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			if ( c == '0' || c == '2' || c == '4' || c== '6' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}
	} else if ( coord == 1 ) {
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			if ( c == '0' || c == '1' || c == '4' || c == '5' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}	
	} else {
		assert ( coord == 2 );
		
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			if ( c == '0' || c == '1' || c == '2' || c == '3' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}	
	}
	out.assign ( tmp, size );
}
#endif

// Too much redundancy from code called many times.
#if 0
template<int dim>
void
dbinary_tree<dim>::project_key( const string& s, int coord, string& ret ) const
{
  assert( 0<=coord && coord<dim );
  ret = s;
//    cout << "ret is "<<ret<<"\n";
  if ( coord == 0 ) {
    for ( int i=0; i<ret.length(); ++i ) {
	if ( ret[i] == '-' ) {
		continue;
	}
      if ( ret[i] == '0' || ret[i] == '2' || ret[i] == '4' || ret[i] == '6' ) {
	ret[i] = '0';
      } else {
	ret[i] = '1';
      }
    }
    return;
  } else if ( coord == 1 ) {
    for ( int i=0; i<ret.length(); ++i ) {
	if ( ret[i] == '-' ) {
		continue;
	}
      if ( ret[i] == '0' || ret[i] == '1' || ret[i] == '4' || ret[i] == '5' ) {
	ret[i] = '0';
      } else {
	ret[i] = '1';
      }
    }
    return;
  } else if ( coord == 2 ) {
    for ( int i=0; i<ret.length(); ++i ) {
	if ( ret[i] == '-' ) {
		continue;
	}
      if ( ret[i] == '0' || ret[i] == '1' || ret[i] == '2' || ret[i] == '3' ) {
	ret[i] = '0';
      } else {
	ret[i] = '1';
      }
    }
    return;
  }
}
#endif

#if 0
template<int dim>
void
dbinary_tree<dim>::get_endpoints_for_projected_key ( string s,
						     double& e0,
						     double& e1 ) const
{
  e0 = 0.0;
  for ( int i=2; i<s.length(); i+=2 ) {
    if ( s[i] == '1' ) {
      e0 += powers_of_one_half[i/2];
    }
  }
  e1 = e0 + powers_of_one_half[ this->level( s ) ];
}
#endif

#if 0
template<int dim>
void
dbinary_tree<dim>::get_box ( string s, box<dim>& bx ) const
{
  double ext[6];
  string tmp_str;
  for ( int d=0; d<dim; ++d ) {
    project_key( s, d, tmp_str );
    get_endpoints_for_projected_key( tmp_str, ext[2*d], ext[2*d+1] );
  }
  bx.set( &ext[0] );
  const double* lzbx = level_zero_box.get();

  if ( dim == 1 ) {
    bx.scale_and_translate( 0.5*(lzbx[1]-lzbx[0]), 0.0, 0.0,
			    0.5*(lzbx[0]+lzbx[1] ), 0.0, 0.0 );
  } else if ( dim == 2 ) {
    bx.scale_and_translate( 0.5*(lzbx[1]-lzbx[0]), 0.5*(lzbx[3]-lzbx[2]), 0.0,
			    0.5*(lzbx[0]+lzbx[1]), 0.5*(lzbx[2]+lzbx[3]), 0.0 );
  } else if ( dim == 3 ) {
    bx.scale_and_translate( 0.5*(lzbx[1]-lzbx[0]), 0.5*(lzbx[3]-lzbx[2]), 0.5*(lzbx[5]-lzbx[4]),
			    0.5*(lzbx[0]+lzbx[1]), 0.5*(lzbx[2]+lzbx[3]), 0.5*(lzbx[4]+lzbx[5]) );
  }
}
#endif

template<int dim>
void
dbinary_tree<dim>::split_dbinary_node_and_push ( string s, stack<string>& keep )
{
  string tmp_s = s + "-0";
  int i = tmp_s.size()-1;
  
  keep.push( tmp_s );
  
  tmp_s[i] = '1';
  keep.push( tmp_s );
  
  if ( dim >= 2 ) {
    tmp_s[i]='2';
    keep.push( tmp_s );
    tmp_s[i]='3';
    keep.push( tmp_s );
  }

  if ( dim >= 3 ) {
    tmp_s[i]='7';
    keep.push( tmp_s );
    
    tmp_s[i]='6';
    keep.push( tmp_s );
    
    tmp_s[i]='5';
    keep.push( tmp_s );
    
    tmp_s[i]='4';
    keep.push( tmp_s );
  }
}

template<int dim>
void
dbinary_tree<dim>::get_global_box( const string& key, box<dim>& out ) const
{
	basic_dtree->get_box ( key, out );
	out.scale ( cover_factor );
}

template<int dim>
void
dbinary_tree<dim>::get_global_box( const vector<string> & in, vector<box<dim> >& out ) const
{
	basic_dtree->get_box ( in, out );
	
	int in_size = in.size();
	for ( int i=0; i<in_size; ++i ) {
		out[i].scale ( cover_factor );
	}
}

#if 0 // Should not be necessary yet.
template<int dim>
void
dbinary_tree<dim>::get_patch_keys_incident_point ( const string&      region,
                                                   const double*      co,
                                                   vector<string>&    out ) const
{
	vector<string> keys;
	basic_dtree->get_key_expanded_intersect_point ( initial_level,
	                                                co,
	                                                cover_factor,
	                                                keys );
	vector<box<dim> > boxes;
	get_global_box ( keys, boxes );
	
	valarray<bool> intersects;
	geom->intersect_region ( region, boxes, intersects );
	
	out.clear();
	
	int keys_size = keys.size ();
	for ( int i=0; i<keys_size; ++i ) {
		if ( intersects[i] ) {
			out.push_back ( keys[i] );
		}
	}
	
#ifndef NDEBUG
	string certain;
	basic_dtree->get_key_intersect_point ( initial_level, co, certain );
	vector<string>::iterator it = find ( out.begin(), out.end(), certain );
	assert ( it != out.end() );
#endif
}
#endif

template <int dim>
void
dbinary_tree<dim>::get_intersect_region ( const string &      region,
                                          vector<string> &    in,
                                          vector<string> &    out )
{
#if 1 // Takes boundary patches also.
	vector<box<dim> > boxes;
	get_global_box ( in, boxes );
	
	valarray<bool> intersects;
	
	assert ( boxes.size() > 0 );
	geom->intersect_region ( region, boxes, intersects );
	
	out.clear();
	
	int in_size = in.size ();
	for ( int i=0; i<in_size; ++i ) {
		if ( intersects[i] ) {
			out.push_back ( in[i] );
		}
	}

	
#else // does not take into account boundary intersecting box.
	int in_size = in.size();
	
	out.clear();
	assert ( out.size() == 0 );
	
	assert ( basic_dtree && geom );
	
	box<dim> unexpanded_box;
	double centre[dim];
	
	for ( int i=0; i<in_size; ++i ) {
		basic_dtree->get_box ( in[i], unexpanded_box );
		unexpanded_box.get_centre_point( centre );
		if ( geom->inside_region ( region, centre ) ) {
			out.push_back ( in[i] );
		}
	}
#endif
}


template <int dim>
void
dbinary_tree<dim>::get_cover_structure ( const string &                          name,
                                         typename cover_structure<dim>::ptr &    out_pcs )
{
	PetscLogDouble log_t0, log_t1;
	
	if ( cover_struct_ ) {
		out_pcs = cover_struct_;
	} else {
		cover_struct_.reset ( new cover_structure<dim> );
	}
	
	cover_struct_->support_name = name;
	
	PetscGetTime ( &log_t0 );
	get_patch_keys ( name, cover_struct_->patch_keys );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "get_patch_keys took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscGetTime ( &log_t0 );
	get_neighbour_keys( name, cover_struct_->patch_neighbour_keys );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "get_neighbour_keys took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	// generate_index_information is straight after get_neighbour_keys because the
	// neigbhour keys need to be adjusted to take into account the cover as defined by
	// the boundary. 2009-01-22 ML.
	PetscGetTime ( &log_t1 );
	generate_index_information ( *cover_struct_ );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "generate_index_information took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscGetTime ( &log_t0 );
	get_patches( cover_struct_->patch_keys, cover_struct_->patches );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "get_patches took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	vector<vector<box<dim> > >& neighs = cover_struct_->neighbour_patches;
	neighs.resize( (cover_struct_->patch_neighbour_keys).size() );
	
	map<string, set<string> >::const_iterator it( cover_struct_->patch_neighbour_keys.begin() );
	map<string, set<string> >::const_iterator end( cover_struct_->patch_neighbour_keys.end() );
	
	PetscGetTime ( &log_t0 );
	for ( int count=0; it!=end; ++count, ++it ) {
		get_patches( it->second, neighs[count] );
	}
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "get_patches for-loop took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	// Fix up an interface mis-match with partition_of_unity_function.
	// FIXME: settle on keeping copies or keeping shared pointers.
	{	
		PetscGetTime ( &log_t0 );
		copy_vector_to_vec_ptr ( cover_struct_->patches,
		                         cover_struct_->ptr_patches );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "copy_vector_to_vec_ptr took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		PetscGetTime ( &log_t0 );
		int num = cover_struct_->neighbour_patches.size();
		(cover_struct_->ptr_neighbour_patches).resize( num );
		for ( int i=0; i<num; ++i ) {
			copy_vector_to_vec_ptr( cover_struct_->neighbour_patches[i], cover_struct_->ptr_neighbour_patches[i] );
		}
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "neighbours copy_vector_to_vec_ptr took %f minutes.\n", (log_t1-log_t0)/60.0 );
	}
	
	out_pcs = cover_struct_;
	
#if 0 // ndef NDEBUG // Keep for future debug.
	// Give information about the generated cover structure.
	int num_patches = pcs->patch_keys.size();
	
	set<string>::iterator pk_it( pcs->patch_keys.begin() );
	map<string, set<string> >::iterator pnk_map_it ( (pcs->patch_neighbour_keys).begin() );
	set<string>::iterator pnk_it, pnk_end;
	
	for ( int cc=0; cc<num_patches; ++cc, ++pk_it, ++pnk_map_it ) {
		std::clog << *pk_it << " has neighbours ";
		pnk_end = (pnk_map_it->second).end();
		for ( pnk_it = (pnk_map_it->second).begin(); pnk_it!=pnk_end; ++pnk_it ) {
			if ( pnk_it != (pnk_map_it->second).begin() ) {
				std::clog << ", ";
			}
			std::clog << *pnk_it;
		}
		std::clog << "\n"
			<< "\tThis corresponds to index " << cc << " and neighbours ";
		int size = pcs->index_to_neighbour_indices[cc].size();
		for ( int i=0; i<size; ++i ) {
			std::clog << pcs->index_to_neighbour_indices[cc][i];
			if ( i!=size-1 ) {
				std::clog << ", ";
			}
		}
		std::clog << "\n";
	}
#endif
}

#if 0
template<int dim>
int
dbinary_tree<dim>::level ( string& s ) const
{
  assert( s.length()%2 == 1 );
  return (s.length()-1)/2;
}
#endif

//

#if 0
template int key_level<1>( const string& );
template int key_level<2>( const string& );
template int key_level<3>( const string& );
#endif

template class dbinary_tree<2>;
