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

#include "box_utils.hh"
#include "cover_structure.hh"

#include <algorithm>
#include <iostream>
using std::clog;
using std::copy;
using std::find;

template <int dim>
void generate_index_information ( cover_structure<dim>& cs )
{
	int num_patches = cs.patch_keys.size();
	
	// Allow determination of the index given the patch key.
	cs.index_to_patch_key.resize( num_patches );
	copy( cs.patch_keys.begin(), cs.patch_keys.end(), cs.index_to_patch_key.begin() );
	
	cs.patch_key_to_index.clear();
	// Allow determination of the patch key given the index.
	for ( int i=0; i<num_patches; ++i ) {
		cs.patch_key_to_index[ cs.index_to_patch_key[i] ] = i;
	}
	
	cs.index_to_neighbour_indices.resize( num_patches );
	
	vector<string>    actual_key_neighbours;
	actual_key_neighbours.reserve( static_cast<size_t>( pow(3,dim) ) );

	// GDB helped find that the neighbour information was being kept about
	// patch keys that did not appear in the cover. These where used to
	// create the patch boxes and the result was approximation functions
	// that were zero everywhere.
	// Keep only keys that appear in the cover.
	// 2009-02-12 ML.
	{
		map<string, set<string> >    wanted_neigh_info;
		
		set<string>::iterator    s_it ( cs.patch_keys.begin() );
		set<string>::iterator    s_end ( cs.patch_keys.end() );
		for ( ; s_it!=s_end; ++s_it ) {
			wanted_neigh_info[*s_it] = cs.patch_neighbour_keys[*s_it];
		}
		assert ( wanted_neigh_info.size() <= cs.patch_neighbour_keys.size() );
		swap ( cs.patch_neighbour_keys, wanted_neigh_info );
	}
	
	/**
		2009-01-22 Optimization that will put the key-neighbours
		and index-neighbours out of step. It has been found that it
		is too expensive to take into account the boundary when creating
		key-neighbours. Key-neighbours will be key-neighbours in the tree
		whereas index-neighbours will be neighbours that are in the cover.
		
		For a given key/index in the cover, the number of index-neighbours
		may be fewer than the number of key-neighbours.
	*/
	set<string>::iterator    s_it  = cs.patch_keys.begin();
	set<string>::iterator    s_end = cs.patch_keys.end();
	for ( int ptch=0; s_it!=s_end; ++ptch, ++s_it ) {
		// GDB helped to find that the entire neighbour map
		// patch_neighbour_keys was being iterated over so ptch
		// was indexing outside the valid range.
		// 2009-02-12 ML.
		assert ( ptch<cs.index_to_neighbour_indices.size() );
	
		vector<int>& neigh_ind = cs.index_to_neighbour_indices[ptch];
		neigh_ind.clear();
		
		map<string, set<string> >::iterator m_it = cs.patch_neighbour_keys.find ( *s_it );
		
		assert ( m_it != cs.patch_neighbour_keys.end() );
		
		int max_num_neigh = (m_it->second).size();
		
		neigh_ind.reserve( max_num_neigh );
		
		set<string> &    key_neighbours = m_it->second;
		
		set<string>::const_iterator n_it ( key_neighbours.begin() );
		set<string>::const_iterator n_end ( key_neighbours.end() );
		
		map<string,int>::const_iterator    pkti_end ( cs.patch_key_to_index.end() );
		
		actual_key_neighbours.clear();
		
		// Add index-neighbour if the key-neighbour lies within the cover.
		for ( ; n_it!=n_end; ++n_it ) {
			map<string,int>::const_iterator    pkti ( cs.patch_key_to_index.find( *n_it ) );
			
			if ( pkti != pkti_end ) {
				neigh_ind.push_back ( pkti->second );
				actual_key_neighbours.push_back ( *n_it );
			}
		}
		
		// Correct the key-neighbours taking into account the cover.
		key_neighbours.clear();
		copy ( actual_key_neighbours.begin(), actual_key_neighbours.end(),
		       inserter ( key_neighbours, key_neighbours.begin() ) );
	}
}

template <int dim>
void generate_patch_information ( cover_structure<dim> & cs )
{
	assert ( cs.ref_structure );
	
	basic_dbinary_tree<dim> &    dtree = *((cs.ref_structure)->dtree);
	
	double expansion = (cs.ref_structure)->cover_factor;
	
	dtree.get_box( cs.patch_keys, cs.patches, expansion );
	
	vector<vector<box<dim> > >& neighs = cs.neighbour_patches;
	neighs.resize( (cs.patch_neighbour_keys).size() );
	
	map<string, set<string> >::const_iterator it( cs.patch_neighbour_keys.begin() );
	map<string, set<string> >::const_iterator end( cs.patch_neighbour_keys.end() );
	
	for ( int count=0; it!=end; ++count, ++it ) {
		dtree.get_box ( it->second, neighs[count], expansion );
	}
	
	// Fix up an interface mis-match with partition_of_unity_function.
	// FIXME: settle on keeping copies or keeping shared pointers.
	{	

		copy_vector_to_vec_ptr ( cs.patches,
		                         cs.ptr_patches );

		int num = cs.neighbour_patches.size();
		(cs.ptr_neighbour_patches).resize( num );
		for ( int i=0; i<num; ++i ) {
			copy_vector_to_vec_ptr( cs.neighbour_patches[i], cs.ptr_neighbour_patches[i] );
		}
	}
}

// This must be called AFTER generate_index_information in order to ensure
// the neighbours all lie within the cover. 2009-07-06 ML.
template <int dim>
void generate_patch_boundary_information ( cover_structure<dim> & inout,
                                           set<string> & bdry_keys,
                                           line_segments& domain )
{
	
	set<string>::iterator    bd_it ( bdry_keys.begin() );
	set<string>::iterator    bd_end ( bdry_keys.end() );
	
	{
		// We must find which strictly interior boxes expand to
		// become boundary patches. These are important for
		// restricting the integration to within the domain.
		// 2009-07-06 ML.
		
		assert ( inout.ref_structure );
		
		basic_dbinary_tree<dim> &    dtree = *((inout.ref_structure)->dtree);
		
		double expansion = (inout.ref_structure)->cover_factor;
		
		set<string>   int_bdry_keys;
		
		assert ( inout.patch_neighbour_keys.size() > 0 );
		
		
		box<dim>    expanded_box;
		
		for ( ; bd_it!=bd_end; ++bd_it ) {
			set<string> &    key_neighbours = inout.patch_neighbour_keys[*bd_it];
			
			set<string>::const_iterator n_it ( key_neighbours.begin() );
			set<string>::const_iterator n_end ( key_neighbours.end() );
			
			for ( ; n_it!=n_end; ++n_it ) {
				dtree.get_box( *n_it, expanded_box, expansion );
				
				int    intersection_first_entry  = -1;
				int    intersection_num_inside   =  0;
			
				get_first_entry_num_inside ( domain,
				                             expanded_box,
				                             0, -1, intersection_first_entry,
				                             intersection_num_inside );
			
				if ( intersection_first_entry != -1 ) {
					int_bdry_keys.insert ( *n_it );
				}
			}
		}
		
		set<string> perm_bdry_keys;
		swap ( perm_bdry_keys, bdry_keys );
		
		set_union ( perm_bdry_keys.begin(), perm_bdry_keys.end(),
		            int_bdry_keys.begin(),  int_bdry_keys.end(),
		            inserter ( bdry_keys, bdry_keys.begin() ) );
	}
	
	// Form the vector<bool>.
	{
		inout.is_boundary_patch.reserve ( inout.patch_keys.size() );
		set<string>::iterator    it ( inout.patch_keys.begin() );
		set<string>::iterator    end ( inout.patch_keys.end() );
		
		// Update the iterators.
		bd_it = bdry_keys.begin();
		bd_end = bdry_keys.end();
		
		inout.is_boundary_patch.clear();
		for ( ; it!=end; ++it ) {
			if ( (bd_it!=bd_end) && (*it == *bd_it) ) {
				++bd_it;
				inout.is_boundary_patch.push_back ( true );
			} else {
				inout.is_boundary_patch.push_back ( false );
			}
		}
	}
}

template <int dim>
void get_patch_keys_incident_point ( cover_structure<dim>& cs, const double* co, vector<string> & out )
{
	int patches_size = cs.patches.size();
	assert ( cs.index_to_patch_key.size() == patches_size );
	
	out.clear();
	
	string    main_key;
	
	assert ( cs.ref_structure );
	get_key_intersect_point ( *(cs.ref_structure), co, main_key );
	
	typedef set<string>::const_iterator    s_cit;
	
	// GDB helped find that evaluation outside the cover was not being allowed
	// as an assert was present to guarantee the incident main_key was found
	// in the neighbour map. 2009-02-12 ML.
	if ( cs.patch_neighbour_keys.find(main_key)==cs.patch_neighbour_keys.end() ) {
		// Point does not intersect the cover.
		return;
	}
	
	const set<string> &    nbrs = cs.patch_neighbour_keys[main_key];
	
	s_cit    it  = nbrs.begin();
	s_cit    end = nbrs.end();
	
	assert ( (cs.ref_structure)->dtree );
	basic_dbinary_tree<dim> &    dtree = *((cs.ref_structure)->dtree);
	
	double expansion = (cs.ref_structure)->cover_factor;
	
	for ( ; it!=end; ++it ) {
		if ( dtree.box_intersect_point ( *it, co, expansion ) ) {
			out.push_back ( *it );
		}
	}
	

// cover_structure now has access to a basic_dbinary_tree via refinement structure.
// It appears this is the most logical approach now. 2009-02-11 ML.

	

#if 0
	// Turns out this is more efficient than using valarray<bool> to
	// carry out all the intersections in one go in a free function
	// that avoids calling closed_intersect_point on each box.
	// 2008-09-10 Mike Li.
	for ( int i=0; i<patches_size; ++i ) {
		if ( cs.patches[i].closed_intersect_point ( co ) ) {
			out.push_back ( cs.index_to_patch_key[i] );
		}
	}
#endif
}

//

template class cover_structure<2>;

template void generate_index_information<2> ( cover_structure<2>& cs );
template void generate_patch_information<2> ( cover_structure<2> & cs );
template void generate_patch_boundary_information<2> ( cover_structure<2>&, set<string>&, line_segments& );

template void get_patch_keys_incident_point<2> ( cover_structure<2>&, const double*, vector<string> & );


