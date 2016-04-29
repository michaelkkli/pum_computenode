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

#include "global_approximation_space.hh"
#include "cubic_spline_weight.hh"
#include "polynomial.hh"

#include <algorithm>
#include <iostream>
#include <set>

using std::clog;
using std::copy;
using std::cout;
using std::find;
using std::set;

template<int dim>
global_approximation_space<dim>::global_approximation_space( local_basis_type lbt )
{
	default_weight_.reset( new cubic_spline_weight<dim> );
	
	default_local_basis_.resize( dim + 1 );
	
	if ( dim==1 ) {
		clog << "USING LINEAR LOCAL BASIS.\n";
		default_local_basis_[0].reset( new polynomial<dim>( 1.0 ) );
		default_local_basis_[1].reset( new polynomial<dim>( 1.0, 1 ) );
	} else if ( dim==2 ) {
		switch ( lbt ) {
		case monomial1:
			default_local_basis_[0].reset( new polynomial<dim>( 1.0, 0, 0 ) );
			default_local_basis_[1].reset( new polynomial<dim>( 1.0, 1, 0 ) );
			default_local_basis_[2].reset( new polynomial<dim>( 1.0, 0, 1 ) );
			break;
		case monomial2:
			default_local_basis_.resize ( 6 );
			default_local_basis_[0].reset( new polynomial<dim>( 1.0, 0, 0 ) );
			default_local_basis_[1].reset( new polynomial<dim>( 1.0, 1, 0 ) );
			default_local_basis_[2].reset( new polynomial<dim>( 1.0, 1, 1 ) );
			default_local_basis_[3].reset( new polynomial<dim>( 1.0, 2, 0 ) );
			default_local_basis_[4].reset( new polynomial<dim>( 1.0, 0, 1 ) );
			default_local_basis_[5].reset( new polynomial<dim>( 1.0, 0, 2 ) );
			break;
		default:
			default_local_basis_[0].reset( new polynomial<dim>( 1.0, 0, 0 ) );
			default_local_basis_[1].reset( new polynomial<dim>( 1.0, 1, 0 ) );
			default_local_basis_[2].reset( new polynomial<dim>( 1.0, 0, 1 ) );
			break;
		}
	} else if ( dim==3 ) {
		clog << "USING LINEAR LOCAL BASIS.\n";
		default_local_basis_[0].reset( new polynomial<dim>( 1.0, 0, 0, 0 ) );
		default_local_basis_[1].reset( new polynomial<dim>( 1.0, 1, 0, 0 ) );
		default_local_basis_[2].reset( new polynomial<dim>( 1.0, 0, 1, 0 ) );
		default_local_basis_[3].reset( new polynomial<dim>( 1.0, 0, 0, 1 ) );
	}
	
	// Very important as global functions are allowed in the local basis.
	// Found that monomials were not being set as local functions when
	// high numbers seen returned from them.
	// 2008-08-20 Mike Li.
	int size = default_local_basis_.size();
	for ( int i=0; i<size; ++i ) {
		default_local_basis_[i]->set_local_function();
	}
}

template<int dim>
global_approximation_space<dim>::~global_approximation_space()
{
}

template <int dim>
void
global_approximation_space<dim>::give_dbinary_tree ( typename dbinary_tree<dim>::ptr & dt )
{
	assert ( dt );
	dtree_ = dt;
}

template <int dim>
void
global_approximation_space<dim>::set_support ( const string & in )
{
	// 2009-02-11 Trying to eliminate use of dbinary_tree::ptr dtree_ in preference of basic_dbinary_tree.
	// assert ( dtree_ && ( dtree_->access_geometry() ).is_support(in) );
	support_ = in;
}

// DEPRECATED.
template <int dim>
void
global_approximation_space<dim>::give_cover_structure ( typename cover_structure<dim>::ptr& pcs )
{
	assert( pcs );
	
	cover_struct_ = pcs;
//	assert ( pcs->ref_structure );
//	assert ( cover_struct_->ref_structure );

	// The boundary keys and entry segment data must be used somewhere in here to
	// allow restriction of integration. Currently investigating partition_of_unity_function.
	// 2009-06-24 ML.
	
	{
		set<string>::const_iterator it  ( cover_struct_->patch_keys.begin() );
		set<string>::const_iterator end ( cover_struct_->patch_keys.end()   );
		
		for ( int c=0 ; it!=end; ++c, ++it ) {
			cached_pu_functions_[*it].reset( new partition_of_unity_function<dim> );
			cached_pu_functions_[*it]->set ( cover_struct_->ptr_patches[c],
						cover_struct_->ptr_neighbour_patches[c],
						default_weight_ );
		}
	}
	
	{
		local_approx_functions_.clear();
		
		set<string>::const_iterator it ( cover_struct_->patch_keys.begin() );
		set<string>::const_iterator end( cover_struct_->patch_keys.end()   );
		
		for ( ; it!=end; ++it ) {
			local_approx_functions_[*it] = default_local_basis_;
		}
	}
}

template <int dim>
void
global_approximation_space<dim>::get_patch_to_num_dof ( valarray<int>& out ) const
{
	
	out.resize( local_approx_functions_.size() );
	
	typename map<string, typename differentiable_function<dim>::vec_ptr>::const_iterator
		it ( local_approx_functions_.begin() );
	
	typename map<string, typename differentiable_function<dim>::vec_ptr>::const_iterator
		end( local_approx_functions_.end()   );
	
	for ( int c=0 ; it!=end; ++c, ++it ) {
		out[c] = (it->second).size();
	}
}

template <int dim>
int
global_approximation_space<dim>::est_patch_num_dof () const
{
	return default_local_basis_.size();
}

template <int dim>
void
global_approximation_space<dim>::get_global_basis_by_patch ( int in, vector<global_basis_function<dim> >& out ) const
{
	string key = cover_struct_->index_to_patch_key[in];
	
	assert( local_approx_functions_.find( key )
	                                      != local_approx_functions_.end() );
	
	const typename differentiable_function<dim>::vec_ptr& patch_local_basis
		= (local_approx_functions_.find( key ))->second;
	
	int size = patch_local_basis.size();
	
	assert ( size > 0 );
	
	out.resize( size );
	
	for ( int i=0; i<size; ++i ) {
		assert ( cached_pu_functions_.find(key) != cached_pu_functions_.end() );
			
		out[i].set( (cached_pu_functions_.find(key))->second,
		            patch_local_basis[i] );
		
		assert ( patch_local_basis[i] );
	}
}

template <int dim>
void 
global_approximation_space<dim>::get_incident_patch_keys_and_indices ( const double* co,
                                                                       vector<string>& keys,
                                                                       vector<int>& indices ) const
{
	keys.clear();
	indices.clear();

	// We can do better than this by eliminating a search
	// that grows quickly.
#if 0
	 get_patch_keys_incident_point( *cover_struct_, co, keys );
	
	int keys_size = keys.size();
	indices.resize ( keys_size );
		
	
	typename map<string, int>::const_iterator it;
	for ( int i=0; i<keys_size; ++i ) {
		it = cover_struct_->patch_key_to_index.find ( keys[i] );
		
		assert ( it != cover_struct_->patch_key_to_index.end() );
		indices[i] = it->second;
	}
#else
	
	string    main_key ( "-1" );
	
	
	// Find a key containing the point that is also in the cover.
	
		vector<string>    many_keys;
		
		assert ( cover_struct_ );
		cover_structure<dim> &    cs = *cover_struct_;
		
		// dtree_->get_patch_containing_point ( co, many_keys );
		// Removing all use of binary_tree.
		// Take into account the cover_structure.
		// 2009-02-12 ML.
		get_patch_keys_incident_point ( cs, co, many_keys );
		if ( many_keys.size()>0 ) {
			// Previously we had to ensure one of the keys contained the point
			// but now we are assured that all keys contain the point.
			// GDB helped to find that the indices were not being calculated at all.
			// 2009-02-13 ML.
			main_key = many_keys[0];
		}

#if 0
		const set<string> &    cover_keys = cs.patch_keys;
		
		set<string>::const_iterator it, end_it ( cover_keys.end() );
		
		int many_keys_size = many_keys.size();
		
		for ( int i=0; i<many_keys_size; ++i ) {
			it = cover_keys.find( many_keys[i] );
			
			if ( it != end_it ) {
				main_key = *it;
				
#ifndef NDEBUG
			{
				typename map<string, int>::const_iterator m_it ( cs.patch_key_to_index.find ( main_key ) );
				assert ( m_it != cs.patch_key_to_index.end() );
			}
#endif
				
				break;
			}
		}
#endif	
		// An attempt may be made to evaluate outside the support and we allow this.
		if ( main_key == "-1" ) {
#if 0//ndef NDEBUG
			{
				vector<string> tmp_keys;
				get_patch_keys_incident_point( *cover_struct_, co, tmp_keys );
				assert ( tmp_keys.size() == 0 );
			}
#endif
			keys.clear ();
			indices.clear ();
			return;
		}
	
	
	keys.clear ();
	indices.clear ();
	
#if 0 // ndef NDEBUG
	{
		typename set<string>::const_iterator s_it ( cs.patch_keys.begin() );
		typename set<string>::const_iterator s_end ( cs.patch_keys.end() );
		cout << "patch_keys check\n";
		for ( ; s_it!=s_end; ++s_it ) {
			cout << *s_it << "\n";
		}
		cout << "\n";
		
		typename map<string,int>::const_iterator b_it ( cs.patch_key_to_index.begin() );
		typename map<string,int>::const_iterator b_end ( cs.patch_key_to_index.end() );
		
		cout << "patch_key_to_index check\n";
		for ( ; b_it!=b_end; ++b_it ) {
			cout << b_it->first << "\t" << b_it->second << "\n";
		}
		cout << "\n";
	}
#endif
	
	typename map<string, int>::const_iterator m_it ( cs.patch_key_to_index.find ( main_key ) );
	assert ( m_it != cs.patch_key_to_index.end() );
	
	const vector<int>&       potential_neigh_indices = cs.index_to_neighbour_indices[ m_it->second ];
	const vector<box<dim> >& potential_neigh_patches = cs.neighbour_patches[ m_it->second ];
	
	assert ( potential_neigh_indices.size () == potential_neigh_patches.size () );
	assert ( potential_neigh_indices.size () > 0 );
	
	int pot_size = potential_neigh_indices.size ();
	
	for ( int i=0; i<pot_size; ++i ) {
		if ( potential_neigh_patches[i].open_intersect_point ( co ) ) {
			keys.push_back ( cs.index_to_patch_key[ potential_neigh_indices[i] ] );
			indices.push_back ( potential_neigh_indices[i] );
		}
	}
	
#if 0//ndef NDEBUG
	{
			vector<string> tmp_keys;
			get_patch_keys_incident_point( *cover_struct_, co, tmp_keys );
			vector<string>::iterator tmp_it ( find ( tmp_keys.begin(), tmp_keys.end(), main_key ) );
			
			assert ( tmp_keys == keys );
			
			typename map<string, int>::const_iterator it;
			
			int tmp_size = tmp_keys.size();
			
			vector<int> tmp_indices ( tmp_size );
			
			for ( int i=0; i<tmp_size; ++i ) {
				it = cover_struct_->patch_key_to_index.find ( tmp_keys[i] );
				
				assert ( it != cover_struct_->patch_key_to_index.end() );
				tmp_indices[i] = it->second;
			}
			assert ( indices == tmp_indices );
	}
#endif
	
#endif
}

template <int dim>
void
global_approximation_space<dim>::get_solution_evaluation_structure ( valarray<double> &                 in,
                                                                     solution_evaluation_structure &    out ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size / dim;
	
	out.eval_keys.resize ( num_pts );
	out.eval_indices.resize ( num_pts );
	
	string    main_key ( "-1" );
	
	assert ( cover_struct_ );
	cover_structure<dim>& cs = *cover_struct_;
	
	vector<string>    many_keys;
	
	for ( int pt=0; pt<num_pts; ++pt ) {
		// The first job is to find a key containing the point that is also in the cover.
		
		// Set a known value so we know when a suitable key is not found.
		main_key = "-1";
		
		// Removal of use of dbinary_tree.
		// 2009-02-12 ML.
		// dtree_->get_patch_containing_point ( &in[dim*pt], many_keys );
		
		get_patch_keys_incident_point ( cs, &in[dim*pt], many_keys );
		if (many_keys.size() > 0) {
			// Any key will do now.
			main_key = many_keys[0];
		}
#if 0
		const set<string> &    cover_keys = cs.patch_keys;
		
		set<string>::const_iterator it, end_it ( cover_keys.end() );
		
		int many_keys_size = many_keys.size();
		
		for ( int i=0; i<many_keys_size; ++i ) {
			it = cover_keys.find( many_keys[i] );
			
			if ( it != end_it ) {
				main_key = *it;
				break;
			}
		}
#endif
		
		// An attempt may be made to evaluate outside the support and we allow this.
		if ( main_key == "-1" ) {
#if 0//ndef NDEBUG
			{
				vector<string> tmp_keys;
				get_patch_keys_incident_point( *cover_struct_, &in[dim*pt], tmp_keys );
				assert ( tmp_keys.size() == 0 );
			}
#endif
			
			out.eval_keys[pt].clear();
			out.eval_indices[pt].clear();
			continue;
		}
		
		out.eval_keys[pt].clear();
		out.eval_indices[pt].clear ();
		
		typename map<string, int>::const_iterator m_it ( cs.patch_key_to_index.find ( main_key ) );
		assert ( m_it != cs.patch_key_to_index.end() );
		
		const vector<int> &       potential_neigh_indices = cs.index_to_neighbour_indices[ m_it->second ];
		const vector<box<dim> > & potential_neigh_patches = cs.neighbour_patches[ m_it->second ];
		
		assert ( potential_neigh_indices.size () == potential_neigh_patches.size () );
		assert ( potential_neigh_indices.size () > 0 );
		
		int pot_size = potential_neigh_indices.size ();
		
		for ( int i=0; i<pot_size; ++i ) {
			if ( potential_neigh_patches[i].open_intersect_point ( &in[dim*pt] ) ) {
				out.eval_keys[pt].push_back ( cs.index_to_patch_key[ potential_neigh_indices[i] ] );
				out.eval_indices[pt].push_back ( potential_neigh_indices[i] );
			}
		}
		
#if 0//ndef NDEBUG
	{
		vector<string> tmp_keys;
		get_patch_keys_incident_point( *cover_struct_, &in[dim*pt], tmp_keys );
		vector<string>::iterator tmp_it ( find ( tmp_keys.begin(), tmp_keys.end(), main_key ) );
		
		assert ( tmp_keys == out.eval_keys[pt] );
		
		typename map<string, int>::const_iterator it;
		
		int tmp_size = tmp_keys.size();
		
		vector<int> tmp_indices ( tmp_size );
		
		for ( int i=0; i<tmp_size; ++i ) {
			it = cover_struct_->patch_key_to_index.find ( tmp_keys[i] );
			
			assert ( it != cover_struct_->patch_key_to_index.end() );
			tmp_indices[i] = it->second;
		}
		assert ( out.eval_indices[pt] == tmp_indices );
	}
#endif
		
	}
}


template <int dim>
double
global_approximation_space<dim>::global_evaluate ( const double*                    co,
                                                   const vector<string> &           incident_patch_keys,
                                                   const vector<const double*> &    coeff_starts ) const
{
	assert ( co );
	
	int size = incident_patch_keys.size();
	
	assert ( coeff_starts.size() == size );
	
	double res = 0.0;
	
	double present_weight = 0.0;
	double weight_sum     = 0.0;
	
	double tmp_res = 0.0;
	
	box<dim> local_box;

	typename map<string, typename differentiable_function<dim>::vec_ptr>::const_iterator map_it;
	
	for ( int i=0; i<size; ++i ) {
		local_box = cover_struct_->patches[ cover_struct_->patch_key_to_index[ incident_patch_keys[i] ] ];
		present_weight = default_weight_->global_evaluate ( local_box, co );
		
		weight_sum += present_weight;
		
		map_it = local_approx_functions_.find ( incident_patch_keys[i] );
		
		assert ( map_it != local_approx_functions_.end() );
		
		const typename differentiable_function<dim>::vec_ptr & local_basis = map_it->second;
		
		int local_basis_size = local_basis.size();
		
		tmp_res = 0.0;
		
		for ( int j=0; j<local_basis_size; ++j ) {
			assert ( local_basis[j] );
			assert ( coeff_starts[i] );
			
			tmp_res += local_basis[j]->global_evaluate ( local_box, co ) * coeff_starts[i][j];
		}
		
		res += present_weight * tmp_res;
	}
	
	return res/weight_sum;
}

#if 0
template<int dim>
bool
global_approximation_space<dim>::is_updated() const
{
	return updated;
}
#endif

#if 0
template<int dim>
int
global_approximation_space<dim>::size() const
{
	
}
#endif

#if 0

template<int dim>
void
global_approximation_space<dim>::set_dbinary_tree ( typename dbinary_tree<dim>::ptr & t )
{
	assert( t );
	dtree = t;
}

template<int dim>
void
global_approximation_space<dim>::set_support ( const string & n )
{
	assert( dtree );
	assert( (dtree->access_geometry()).is_region( n ) ); // Support only regions for now.
	support = n;
}

template<int dim>
void
global_approximation_space<dim>::update()
{
	assert( dtree );
	assert( support != "" );
	{
		const set<string>& all_keys = dtree->access_region_all_keys( support );
		patch_keys.resize( all_keys.size() );
		copy( all_keys.begin(), all_keys.end(), patch_keys.begin() );
	}
	{
		const set<string>& interior_keys = dtree->access_region_interior_keys( support );
		patch_keys_interior.resize( interior_keys.size() );
		copy( interior_keys.begin(), interior_keys.end(), patch_keys_interior.begin() );
	}
	{
		const set<string>& boundary_keys = dtree->access_region_boundary_keys( support );
		patch_keys_boundary.resize( boundary_keys.size() );
		copy( boundary_keys.begin(), boundary_keys.end(), patch_keys_boundary.begin() );
	}
	
	patch_weight = default_weight;
	
	{
		patch_local_basis.clear();
		
		typename vector<string>::iterator it(patch_keys.begin()), end(patch_keys.end());
		for ( ; it!=end ; ++it ) {
			patch_local_basis[*it] = default_local_basis;
		}
	}
	
	{
		int num_patches = patch_keys.size();
		int size_local_basis = default_local_basis.size();
		
		vector<box<dim> > patches;
		
		dtree->get_patches( dtree->access_region_all_keys(support), patches );
		
		global_basis.resize( num_patches * size_local_basis );
		
		for ( int ptch=0; ptch<num_patches; ++ptch ) {
			for ( int i=0; i<size_local_basis; ++i ) {
				typename global_basis_function<dim>::ptr& gb_fn = global_basis[ ptch*size_local_basis + i ];
				if ( !gb_fn ) {
					gb_fn.reset( new global_basis_function<dim> );
				}
				//gb_fn->set_all( patches[ptch], i, 
			}
		}
	}

	updated = true;
}
#endif

//

template class global_approximation_space<2>;
