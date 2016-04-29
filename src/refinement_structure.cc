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

#include "refinement_structure.hh"
#include <set>
using std::set;

template <int dim>
void generate_keys ( refinement_structure<dim> & in, vector<string> & out )
{
	typedef map<string, int, levelwise_less>::iterator   mit_t;
	mit_t    it(in.levels.begin());
	mit_t    end(in.levels.end());

	set<string>    keys;
	
	assert ( it != end );
	assert ( it->first == "0" );
	generate_keys<dim> ( it->second, out );
	keys.insert ( out.begin(), out.end() );
	++it;
	
	for ( ; it!=end; ++it ) {
		generate_keys<dim> ( it->first, it->second, out );
#ifndef NDEBUG
		int tmp = keys.erase ( it->first );
		assert ( tmp == 1 );
#else
		keys.erase ( it->first );
#endif

		// This will not take into account the ordered
		// nature of the vector.
		keys.insert ( out.begin(), out.end() );
	}

	out.reserve ( keys.size() );
	out.assign ( keys.begin(), keys.end() );
}


template <int dim>
void generate_keys_and_neighbours ( refinement_structure<dim> &     rs,
                                    vector<string> &                out,
                                    map<string, set<string> >&      nbrs,
                                    line_segments*                  boundary_segs,
                                    map<string, pair<int,int> >*    boundary_detail )
{
	// Template version of this does not support boundary things yet. See the
	// 2D version for efficient finding of boundary patches.
	assert( !boundary_segs && !boundary_detail );

	double expansion = rs.cover_factor;
	
	// Refinement structure map iterator type.
	typedef map<string, int, levelwise_less>::iterator   rsmit_t;
	rsmit_t    rsm_it(rs.levels.begin());
	rsmit_t    rsm_end(rs.levels.end());
	
	set<string>    keys;
	nbrs.clear();

	// Make sure we have been given a non-empty refinement request and
	// that we always being the refinement from the parent box with key "0".
	assert ( (rsm_it!=rsm_end) && (rsm_it->first=="0") );
	
	// Generate the initial set of keys. The return is
	// by reference into a vector to allow proper ordering of
	// the keys during recursive creation.
	generate_keys<dim> ( rsm_it->second, out );
	
	// Transfer to a set to allow insertion and removal as we proceed.
	keys.insert ( out.begin(), out.end() );
	++rsm_it;
	
	assert ( rs.dtree );
	basic_dbinary_tree<dim> &    tree = *(rs.dtree);
	
	// Use the basic_dbinary_tree to get complete neighbours.
	// This can only be done efficiently for some fixed level.
	// We use a step-wise approach later to modify this
	// neighbour map as we further refine the keys.
	tree.get_neighbours ( keys, nbrs );
	
	typedef set<string>::iterator       sit_t; // Set iterator type.
	typedef vector<string>::iterator    vit_t; // Vector iterator type.
	
	for ( ; rsm_it!=rsm_end; ++rsm_it ) {
		// Starting from rsm_it->first, carry out rsm_it->second
		// complete refinements and place the result temporarily
		// in `out'.
		generate_keys<dim> ( rsm_it->first, rsm_it->second, out );
		
		keys.erase ( rsm_it->first );
		
		// Insert the new keys.
		keys.insert ( out.begin(), out.end() );
		
		set<string>& known_nbrs = nbrs[rsm_it->first];
		
		// Erase the self-neighbour of the refined/replaced key.
		known_nbrs.erase(rsm_it->first);
		
		sit_t    s_it  = known_nbrs.begin();
		sit_t    s_end = known_nbrs.end();
		
		// Establish neighbour relationship for the new keys in 'out'
		// by doing an exhaustive local check.
		// This makes sure that the neighbour information for existing keys
		// is correct.
		for ( ; s_it != s_end; ++s_it ) {
			vit_t    v_it = out.begin();
			vit_t    v_end = out.end();
			for ( ; v_it!=v_end; ++v_it ) {
				if ( tree.are_neighbours ( *s_it, *v_it, expansion ) ) {
					nbrs[*s_it].insert ( *v_it );
					nbrs[*v_it].insert ( *s_it );
				}
			}
		}
		
		// Final record of the refined/replaced key is erased.
		nbrs.erase ( rsm_it->first );
		
		// Make sure the neighbour information for the new child keys
		// added is also correct.
		if ( rsm_it->second == 1 ) {
			// If we only refine by one level, we are guaranteed
			// that the new children are all neighbours of each other.
			vit_t    v_it = out.begin();
			vit_t    v_beg = out.begin();
			vit_t    v_end = out.end();
			for ( ; v_it!=v_end; ++v_it ) {
				nbrs[*v_it].insert ( v_beg, v_end );
			}
		} else {
			// If we refine by more than one level, we conduct an
			// exhaustive search for neighbours.
			
			vit_t    v_outer_it;
			vit_t    v_end;
			
			// Exhaustive comparison.
			for ( ; v_outer_it!=v_end; ++v_outer_it ) {
				vit_t    v_inner_it = v_outer_it;
				
				nbrs[*v_outer_it].insert ( *v_inner_it ); // Self-neighbour.
				++v_inner_it;
				
				for ( ; v_inner_it!=v_end; ++v_inner_it ) {
					if ( tree.are_neighbours ( *v_outer_it, *v_inner_it, expansion ) ) {
						nbrs[*v_outer_it].insert ( *v_inner_it );
						nbrs[*v_inner_it].insert ( *v_outer_it );
					}
				}
			}
		}
	}
	
	// Return the keys by reference.
	
	// Changed from reserve to resize.
	// 2009-02-17 ML.
	out.resize ( keys.size() );
	out.assign ( keys.begin(), keys.end() );
}

template <>
void generate_keys_and_neighbours<2> ( refinement_structure<2> &     rs,
                                    vector<string> &                out,
                                    map<string, set<string> >&      nbrs,
                                    line_segments*                  boundary_segs,
                                    map<string, pair<int,int> >*    boundary_detail )
{
	double expansion = rs.cover_factor;
	
	// Refinement structure map iterator type.
	typedef map<string, int, levelwise_less>::iterator   rsmit_t;
	rsmit_t    rsm_it(rs.levels.begin());
	rsmit_t    rsm_end(rs.levels.end());
	
	set<string>    keys;
	nbrs.clear();

	// Make sure we have been given a non-empty refinement request and
	// that we always being the refinement from the parent box with key "0".
	assert ( (rsm_it!=rsm_end) && (rsm_it->first=="0") );
	
	// Generate the initial set of keys. The return is
	// by reference into a vector to allow proper ordering of
	// the keys during recursive creation.
	generate_keys<2> ( rsm_it->second, out );
	
	// Transfer to a set to allow insertion and removal as we proceed.
	keys.insert ( out.begin(), out.end() );
	++rsm_it;
	
	assert ( rs.dtree );
	basic_dbinary_tree<2> &    tree = *(rs.dtree);
	
	// Use the basic_dbinary_tree to get complete neighbours.
	// This can only be done efficiently for some fixed level.
	// We use a step-wise approach later to modify this
	// neighbour map as we further refine the keys.
	tree.get_neighbours ( keys, nbrs );
	
	if ( boundary_segs ) {
		// If boundary segments are supplied, the boundary_detail
		// data structure must be available for output.
		assert ( boundary_detail );
	}
	
	int num_segs = 0; // Unused unless requested with boundary_detail != 0.
	vector<string>    boundary_keys;
	vector<int>       entry_segments;
	
	if ( boundary_detail ) {
		// We cannot populate boundary_detail without boundary segments
		// information.
		assert ( boundary_segs );
		num_segs = (boundary_segs->segments.size()/2)-1;
		
		boundary_detail->clear();
		
		get_intersect_keys_entry_indices ( rs,
                                                   *boundary_segs,
                                                   boundary_keys,
                                                   entry_segments );

		int num_bkeys = boundary_keys.size();
		
		// Some keys are boundary keys regardless of the expansion factor
		// whereas others become boundary keys after expansion.
		// Keep also first entry.
		map<string,int>    neighbours_of_permanent_boundary_keys;
		
		// Store details for permanent boundary keys and make a list of boundary
		// key neighbours with a guess at the first entry.
		for ( int i=0; i<num_bkeys; ++i ) {
			// Store the details that we know so far.
			pair<int, int> &    num_entry_num_inside = (*boundary_detail)[boundary_keys[i]];
			num_entry_num_inside.first               = entry_segments[i];
			
			// We do not know the number of boundary points inside yet.
			num_entry_num_inside.second              = -1;
			
			set<string> &    nb = nbrs[boundary_keys[i]];
			
			set<string>::iterator   nb_it  (nb.begin());
			set<string>::iterator   nb_end (nb.end());
			// Store neighbours in a way that gives us a head start on the first
			// entry.
			for ( ; nb_it!=nb_end; ++nb_it ) {
				// Subtract 10 to give a lead-in to the intersection segment.
				// It may cause needless extra box-line intersection checks
				// but this is likely better than going all around the boundary
				// conducting these checks.
			
				int guess_entry = num_entry_num_inside.first-10;
				if ( guess_entry > num_segs-1 ) {
					guess_entry -= num_segs;
					assert ( guess_entry<num_segs );
				}
			
				neighbours_of_permanent_boundary_keys[*nb_it] = guess_entry;
			}
		}
		
		
		map<string,int>::iterator msi_it  ( neighbours_of_permanent_boundary_keys.begin() );
		map<string,int>::iterator msi_end ( neighbours_of_permanent_boundary_keys.end() );
		
		box<2>    bx;
	
		// This bit is clever as we only check neighbours of permanent boundary
		// patches when looking for patches that become boundary patches
		// due to expansion.	
		for ( ; msi_it!=msi_end; ++msi_it ) {
			map<string, pair<int,int> >::iterator   bd_it  = boundary_detail->find(msi_it->first);
			map<string, pair<int,int> >::iterator   bd_end = boundary_detail->end();
			
			// Only handle keys that we have not seen already.
			// These are the keys that intersect the boundary
			// only when expanded.
			if ( bd_it==bd_end ) {
				
				tree.get_box( msi_it->first, bx, expansion );
				
				pair<int,int> tmp;
				
				get_first_entry_num_inside ( *boundary_segs,
                                                             bx,
                                                             msi_it->second,
                                                             -1,
                                                             tmp.first,
                                                             tmp.second );

				if ( tmp.first != -1 ) {
					// Actual first entry was not -1 hence we have found a boundary patch.
					// A new entry is made in *boundary_detail.
					pair<int,int> &    entry_and_inside = (*boundary_detail)[msi_it->first];
					entry_and_inside.first  = tmp.first;
					entry_and_inside.second = tmp.second;
				}
			}
		}
	}
	
	// Everything above has been purely to handle the initial refinement of
	// the parent "0".
	
	typedef set<string>::iterator       sit_t; // Set iterator type.
	typedef vector<string>::iterator    vit_t; // Vector iterator type.
	
	// TODO: boundary handling for general refinement.
	// Boundary handling for general refinement has not been written yet
	// in order to get full boundary handling in the complete refinement
	// case first.
	for ( ; rsm_it!=rsm_end; ++rsm_it ) {
		// Starting from rsm_it->first, carry out rsm_it->second
		// complete refinements and place the result temporarily
		// in `out'.
		generate_keys<2> ( rsm_it->first, rsm_it->second, out );
		
		keys.erase ( rsm_it->first );
		
		// Insert the new keys.
		keys.insert ( out.begin(), out.end() );
		
		set<string>& known_nbrs = nbrs[rsm_it->first];
		
		// Erase the self-neighbour of the refined/replaced key.
		known_nbrs.erase(rsm_it->first);
		
		sit_t    s_it  = known_nbrs.begin();
		sit_t    s_end = known_nbrs.end();

		// This section is not entered for a single complete refinement
		// but more detail refinement requires removal of deleted key
		// from all neighoubrs.
		abort();
		
		
		// Establish neighbour relationship for the new keys in 'out'
		// by doing an exhaustive local check.
		// This makes sure that the neighbour information for existing keys
		// is correct.
		for ( ; s_it != s_end; ++s_it ) {
			vit_t    v_it = out.begin();
			vit_t    v_end = out.end();
			for ( ; v_it!=v_end; ++v_it ) {
				if ( tree.are_neighbours ( *s_it, *v_it, expansion ) ) {
					nbrs[*s_it].insert ( *v_it );
					nbrs[*v_it].insert ( *s_it );
				}
			}
		}
		
		// Final record of the refined/replaced key is erased.
		nbrs.erase ( rsm_it->first );
		
		// Make sure the neighbour information for the new child keys
		// added is also correct.
		if ( rsm_it->second == 1 ) {
			// If we only refine by one level, we are guaranteed
			// that the new children are all neighbours of each other.
			vit_t    v_it = out.begin();
			vit_t    v_beg = out.begin();
			vit_t    v_end = out.end();
			for ( ; v_it!=v_end; ++v_it ) {
				nbrs[*v_it].insert ( v_beg, v_end );
			}
		} else {
			// If we refine by more than one level, we conduct an
			// exhaustive search for neighbours.
			
			vit_t    v_outer_it;
			vit_t    v_end;
			
			// Exhaustive comparison.
			for ( ; v_outer_it!=v_end; ++v_outer_it ) {
				vit_t    v_inner_it = v_outer_it;
				
				nbrs[*v_outer_it].insert ( *v_inner_it ); // Self-neighbour.
				++v_inner_it;
				
				for ( ; v_inner_it!=v_end; ++v_inner_it ) {
					if ( tree.are_neighbours ( *v_outer_it, *v_inner_it, expansion ) ) {
						nbrs[*v_outer_it].insert ( *v_inner_it );
						nbrs[*v_inner_it].insert ( *v_outer_it );
					}
				}
			}
		}
	}
	
	// Return the keys by reference.
	
	// Changed from reserve to resize.
	// 2009-02-17 ML.
	out.resize ( keys.size() );
	out.assign ( keys.begin(), keys.end() );
}


template <int dim>
void get_key_intersect_point ( refinement_structure<dim> &    in,
                               const double*                  co,
                               string &                       out )
{
	assert ( (in.dtree != 0) && (co != 0));
	
	typedef map<string, int, levelwise_less>::iterator   mit_t;
	mit_t    it(in.levels.begin());
	mit_t    end(in.levels.end());

	// This may not be efficient.
	box<dim>    bx;
	basic_dbinary_tree<dim> &    dtree = *(in.dtree);
	
	assert ( it != end );
	int level = 0;
	assert ( it->first == "0" );
	level += it->second;
	++it;

	int skip_size = 1;
	
	for ( ; it!=end; ++it ) {
		const string& key = it->first;
		if ( key.size() == skip_size ) {
			continue;
		}
		
		dtree.get_box ( key, bx );
		if ( bx.open_intersect_point( co ) ) {
			level += it->second;
		}
	}
	dtree.get_key_intersect_point( level, co, out );
}

//

template struct refinement_structure<2>;
template struct refinement_structure<3>;

template void generate_keys<2> ( refinement_structure<2>&, vector<string>& );
template void generate_keys<3> ( refinement_structure<3>&, vector<string>& );

// template void generate_keys_and_neighbours<2> ( refinement_structure<2>&, vector<string>&, map<string, set<string> >& );


template void get_key_intersect_point<2> ( refinement_structure<2>&, const double*, string& );
template void get_key_intersect_point<3> ( refinement_structure<3>&, const double*, string& );
