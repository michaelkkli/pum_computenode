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

#include "pum_discretization.hh"
#include "polynomial.hh"
#include "point_utils.hh"

#include "cubic_spline_weight.hh"

#include "interior_predicate.hh"

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>

using std::bind1st;
using std::cout;
using std::distance;
using std::max;
using std::mem_fun;
using std::unary_function;

template<int dim>
pum_discretization<dim>::pum_discretization() :
	have_coefficients(false),
	have_weights(true),
	have_local_basis(true),
	have_integration_scheme(false),
	updated(false)
{
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
	MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );
	
#if 0
	// have_weights and have_local_basis default to true as default functions are
	// set on creation of species. One must work to invalidate this.
	shared_ptr<polynomial<dim> > tmp_weight( new polynomial<dim>( 1.0 ) );
	tmp_weight->set_local_function();
#endif
	shared_ptr<function<dim> > tmp_weight( new cubic_spline_weight<dim> );
	
	default_weight = tmp_weight;
	assert( default_weight );
	
	vector<shared_ptr<polynomial<dim> > > tmp_local_basis;
	tmp_local_basis.resize( dim + 1 );
	if ( dim==1 ) {
		tmp_local_basis[0].reset( new polynomial<dim>( 1.0 ) );
		tmp_local_basis[1].reset( new polynomial<dim>( 1.0, 1 ) );
	} else if ( dim==2 ) {
		tmp_local_basis[0].reset( new polynomial<dim>( 1.0, 0, 0 ) );
		tmp_local_basis[1].reset( new polynomial<dim>( 1.0, 1, 0 ) );
		tmp_local_basis[2].reset( new polynomial<dim>( 1.0, 0, 1 ) );
	} else if ( dim==3 ) {
		tmp_local_basis[0].reset( new polynomial<dim>( 1.0, 0, 0, 0 ) );
		tmp_local_basis[1].reset( new polynomial<dim>( 1.0, 1, 0, 0 ) );
		tmp_local_basis[2].reset( new polynomial<dim>( 1.0, 0, 1, 0 ) );
		tmp_local_basis[3].reset( new polynomial<dim>( 1.0, 0, 0, 1 ) );
	}
	
	default_local_basis.clear();
	
	typename vector<shared_ptr<polynomial<dim> > >::iterator it( tmp_local_basis.begin() );
	typename vector<shared_ptr<polynomial<dim> > >::iterator end( tmp_local_basis.end() );
	for ( ; it!=end; ++it ) {
		(*it)->set_local_function();
		assert( (*it)->is_local_function() );
		default_local_basis.push_back( *it );
	}
}

template<int dim> pum_discretization<dim>::~pum_discretization() { }

template<int dim>
void
pum_discretization<dim>::set_dbinary_tree ( shared_ptr<dbinary_tree<dim> >& dt )
{
	assert( dt );
	dtree = dt;
	
	updated = false;
}

template<int dim>
void
pum_discretization<dim>::create_species ( const string& name, const string& region_or_boundary )
{
	assert( dtree );
	assert( dtree->is_updated() );
	
	assert( species_support.find( name ) == species_support.end() );
	const geometry<dim>& rg = dtree->access_geometry();
	assert( rg.is_boundary( region_or_boundary ) || rg.is_region( region_or_boundary ) );
	species_support[name] = region_or_boundary;
	
	// Set up patch weight.
	assert( present_patch_weight.find( name ) == present_patch_weight.end() );
	present_patch_weight[name] = default_weight;
			
	// Set up the local basis.
	assert( present_local_basis.find( name ) == present_local_basis.end() );
	map<string, vector<const_pfunction> >& species_local_basis = present_local_basis[name];;
	
	set<string> patch_keys;
	dtree->get_patch_keys( region_or_boundary, patch_keys );
	assert( patch_keys.size() > 0 );
	typename set<string>::const_iterator k_it( patch_keys.begin() );
	typename set<string>::const_iterator k_end( patch_keys.end() );
	for ( ; k_it!=k_end; ++k_it ) {
		// GDB found that the rhs had been omitted.
		species_local_basis[ *k_it ] = default_local_basis;
	}
	
	updated          = false;
}

#if 0 // TODO: remove.
template<int dim>
void
pum_discretization<dim>::species_set_function( const string& species, const pfunction& ic )
{
	assert( is_species(species) );
	assert(ic);
	species_initial_conditions[species] = ic;
	assert( species_initial_conditions[species]->is_global_function() );
}
#endif

// TODO: remove
#if 0
template<int dim>
const typename pum_discretization<dim>::const_pfunction&
pum_discretization<dim>::get_species_function( const string& species ) const
{
	assert( is_species( species ) );
	typename map<string, const_pfunction>::const_iterator it( species_initial_conditions.find( species ) );
	assert( it != species_initial_conditions.end() );
	return it->second;
}
#endif


template<int dim>
void
pum_discretization<dim>::set_integration_scheme( const pintegration_scheme int_sch )
{
	assert( int_sch );
	integrator = int_sch;
	have_integration_scheme = true;
}


template<int dim>
bool
pum_discretization<dim>::is_species ( const string& species_name ) const
{
	return species_support.find( species_name ) != species_support.end();
}
template<int dim>
bool
pum_discretization<dim>::is_boundary_species( const string& species ) const
{
	typename map<string, string>::const_iterator tmp_it(species_support.find(species));
	assert( tmp_it!=species_support.end() );
	return (dtree->access_geometry()).is_boundary(tmp_it->second);
}

template<int dim>	
bool
pum_discretization<dim>::is_interior_species( const string& species ) const
{
	typename map<string, string>::const_iterator tmp_it(species_support.find(species));
	assert( tmp_it!=species_support.end() );
	return (dtree->access_geometry()).is_region(tmp_it->second);
}

#if 0
template<int dim>
int pum_discretization<dim>::get_num_patches( const string& species_name ) const
{
	// TODO: boundary species.
	int boundary_patches = (dtree->access_region_boundary_keys.find( species_name ))->size();
	int interior_patches = (dtree->access_region_interior_keys.find( species_name ))->size();
	assert( boundary_patches + interior_patches > 0 );
	return boundary_patches + interior_patches;
}
#endif // 0

template<int dim>
int
pum_discretization<dim>::get_num_dof ( const string& species ) const
{
	assert( is_species( species ) );
	assert( present_local_basis.find( species ) != present_local_basis.end() );
	const map<string, vector<const_pfunction> >& local_approx =
			 present_local_basis.find( species )->second;
	
	int total = 0;
	
	typename map<string, vector<const_pfunction> >::const_iterator it( local_approx.begin() );
	typename map<string, vector<const_pfunction> >::const_iterator end( local_approx.end() );
	for( ; it!=end; ++it ) {
		total += (it->second).size();
	}
	return total;
}

template<int dim>
int
pum_discretization<dim>::get_num_patches ( const string& species ) const
{
	assert( is_species(species) );
	
	return dtree->get_num_patches( get_species_support(species) );
}

// TODO: remove.
#if 0
template<int dim>
void
pum_discretization<dim>::set_coefficients( const string& species, const valarray<double>& coeff )
{
	assert( updated );
	assert( coeff.size() == get_num_dof( species ) );
	assert( is_species( species ) );
	
	valarray<double>& valarr = present_coefficients[species];
	
	valarr.resize ( coeff.size() );
	valarr = coeff;
	
	// For some reason simple assignment like this doesn't work.
	// Resize and then assign is ok.
	// present_coefficients[ species ] = coeff;
	
	// Simple assignment of valarray doesn't work. Keep this assert.
	assert( present_coefficients[ species ].size() == coeff.size() );
	
	have_coefficients = true;
}
#endif


// TODO: remove.
#if 0
template<int dim>
const valarray<double>&
pum_discretization<dim>::access_species_coefficients ( const string& species ) const
{
	assert( is_species( species ) );
	assert( updated );
	assert( have_coefficients );
	
	assert( present_coefficients.find( species ) != present_coefficients.end() );
	
	return (present_coefficients.find( species ))->second;
}
#endif

template<int dim>
void
pum_discretization<dim>::set_default_weight_and_local_basis()
{
	typename map<string, string>::const_iterator sp_it( species_support.begin() );
	typename map<string, string>::const_iterator sp_end( species_support.end() );
	for ( ; sp_it!=sp_end; ++sp_it ) {
		present_patch_weight[ sp_it->first ] = default_weight;
		map<string, vector<const_pfunction> >& species_local_basis = present_local_basis[ sp_it->first ];
		typename map<string, vector<const_pfunction> >::iterator k_it( species_local_basis.begin() );
		typename map<string, vector<const_pfunction> >::iterator k_end( species_local_basis.end() );
		for ( ; k_it!=k_end; ++k_it ) {
			k_it->second = default_local_basis;
		}
	}
	have_weights     = true;
	have_local_basis = true;
}

template<int dim>
void
pum_discretization<dim>::update_global_index_start ()
{
	assert( dtree );
	assert( have_weights && have_local_basis );
	
	present_global_index_start.clear();
	
	// Iterate over species to make the code easier to read.
	typename map<string, string>::const_iterator sp_it( species_support.begin() );
	typename map<string, string>::const_iterator sp_end( species_support.end() );
	set<string> patch_keys;
	int previous_index;
	for ( ; sp_it!=sp_end; ++sp_it ) {
		map<string, vector<const_pfunction> >& species_present_local_basis =
				present_local_basis[sp_it->first];
				
		map<string, int>& species_global_index_start =
				present_global_index_start[sp_it->first];
		
		typename map<string, vector<const_pfunction> >::iterator it( species_present_local_basis.begin() );
		typename map<string, vector<const_pfunction> >::iterator end( species_present_local_basis.end() );
		
		previous_index = 0;
		for ( ; it!=end; ++it ) {
			species_global_index_start[it->first] = previous_index;
			previous_index += (it->second).size();
		}
	}
}

template<int dim>
void
pum_discretization<dim>::update_parallel_partitioning()
{
	// Decide on how to allocate the species patches.
	{
		valarray<int> species_num_patches( species_support.size() );
		int* v_it = &species_num_patches[0];
		
		map<string, string>::const_iterator sp_it( species_support.begin() );
		map<string, string>::const_iterator sp_end( species_support.end() );
		for ( ; sp_it!=sp_end; ++sp_it ) {
			*(v_it++) = get_num_patches( sp_it->first );
		}
		
		assert( comm_size > 0 );
		int n = species_num_patches.sum();
		
		assert( n >= comm_size ); // At least patch per compute node.
		
		int spare = n % comm_size;
		int each  = (n-spare)/comm_size;
		
		present_patch_allocation.resize(comm_size);
		
		for ( int i=0; i<comm_size; ++i ) {
			present_patch_allocation[i] = each;
			if ( spare > 0 ) {
				// Distribute the extra dofs.
				++present_patch_allocation[i];
				--spare;
			}
		}
		// Sanity check to make sure we have not mis-allocated.
		assert( present_patch_allocation.sum() == n );
	}
	
	// Allocate the patches. The species are kept separate and
	// it is species-patches that are allocated so a patch may
	// appear as many times as there are species using that patch.
	{
		present_parallel_partition.clear();
		
		int check_num_allocated = 0;
		int r=0;
		int allocate = present_patch_allocation[r];
		typename map<string, string>::const_iterator sp_it( species_support.begin() );
		typename map<string, string>::const_iterator sp_end( species_support.end() );
		for ( ; sp_it!=sp_end; ++sp_it ) {
			const set<string>& species_patch_keys = access_patch_keys(sp_it->first);
			typename set<string>::const_iterator k_it( species_patch_keys.begin() );
			typename set<string>::const_iterator k_end( species_patch_keys.end() );
			for ( ; k_it!=k_end; ++k_it ) {
				if ( allocate > 0 ) {
					present_parallel_partition[r][sp_it->first].insert( *k_it );
					--allocate;
					++check_num_allocated;
				} else {
					++r;
					allocate = present_patch_allocation[r];
					present_parallel_partition[r][sp_it->first].insert( *k_it );
					--allocate;
					++check_num_allocated;
				}
			}
		}
		
		assert( check_num_allocated == present_patch_allocation.sum() );
	}
	
	// Allocate dof.
	{
		int n;
		present_allocation_of_dof.resize(comm_size, 0);
		assert( present_allocation_of_dof.size() == comm_size );
		for ( int r=0; r<comm_size; ++r ) {
			n = 0;
			assert( present_parallel_partition.find(r)!=present_parallel_partition.end() );
			const map<string, set<string> >& species_patch_keys_on_r = (present_parallel_partition.find(r))->second;
			typename map<string, set<string> >::const_iterator sp_it( species_patch_keys_on_r.begin() );
			typename map<string, set<string> >::const_iterator sp_end( species_patch_keys_on_r.end() );
			for ( ; sp_it!=sp_end; ++sp_it ) {
				typename set<string>::const_iterator it( (sp_it->second).begin() );
				typename set<string>::const_iterator end( (sp_it->second).end() );
				for ( ; it!=end; ++it ) {
					assert( present_local_basis.find(sp_it->first)!=present_local_basis.end() );
					assert( present_local_basis[sp_it->first].find(*it)!=present_local_basis[sp_it->first].end() );
					
					// Count the size of the local basis in each species-patch.
					n += present_local_basis[sp_it->first][*it].size();
				}
			}
			present_allocation_of_dof[r] = n;
		}
		assert( present_allocation_of_dof.sum() == get_num_dof() );
	}
}

template<int dim>
void
pum_discretization<dim>::update_local_basis_sizes ()
{
	present_local_basis_size.clear();
	
	// GDB showed loop was written incorrectly as was assuming present_local_basis_size
	// was partially populated with patch keys when in fact this member function must
	// do the entire population. The correct version gives a greatly simplified inner loop!
	
	typename map<string, string>::const_iterator sp_it( species_support.begin() );
	typename map<string, string>::const_iterator sp_end( species_support.end() );
	set<string> patch_keys;
	for ( ; sp_it!=sp_end; ++sp_it ) {
		map<string, int>& species_local_basis_size =
				present_local_basis_size[sp_it->first];
		map<string, vector<const_pfunction> >& species_present_local_basis =
				present_local_basis[sp_it->first];
		
		typename map<string, vector<const_pfunction> >::iterator it( species_present_local_basis.begin() );
		typename map<string, vector<const_pfunction> >::iterator end( species_present_local_basis.end() );
		
		for ( ; it!=end; ++it ) {
			species_local_basis_size[ it->first ] = (it->second).size();
		}
	}
}

template<int dim>
int
pum_discretization<dim>::get_local_basis_size( const string& species, const string& key ) const
{
	assert( is_species_patch_key( species, key ) );
	
	assert( present_local_basis_size.find( species ) != present_local_basis_size.end() );
	
	const map<string, int>& species_local_basis_size =
			(present_local_basis_size.find( species ))->second;
			
	assert( species_local_basis_size.find( key ) != species_local_basis_size.end() );
			
	return (species_local_basis_size.find( key ))->second;
}

template<int dim>
void
pum_discretization<dim>::update()
{
	assert( dtree );
	assert( dtree->is_updated() );
	
	update_global_index_start();
	
	update_local_basis_sizes();
	
	update_parallel_partitioning();
	
	// TODO: remove
	// Initial condition set for all species.
	// assert( species_initial_conditions.size() == species_support.size() );
	
	updated = true;
}

template<int dim>
bool
pum_discretization<dim>::is_updated() const
{
	return updated;
}

template<int dim>
int
pum_discretization<dim>::get_global_index_start( const string& species, const string& key ) const
{
	assert( updated );
	assert( is_species( species ) );
	
	const map<string, int>& species_global_index_start = (present_global_index_start.find( species ))->second;
	
	return (species_global_index_start.find( key ))->second;
}

template<int dim>
const typename pum_discretization<dim>::const_pfunction&
pum_discretization<dim>::access_patch_weight( const string& species ) const
{
	assert( updated );
	assert( species_support.find( species ) != species_support.end() );
	assert( present_patch_weight.find( species ) != present_patch_weight.end() );
	return (present_patch_weight.find( species ))->second;
}

template<int dim>
const map<string, vector<typename pum_discretization<dim>::const_pfunction> >&
pum_discretization<dim>::access_local_basis( const string& species ) const
{
	assert( updated );
	assert( species_support.find( species ) != species_support.end() );
	assert( present_local_basis.find( species ) != present_local_basis.end() );
	return (present_local_basis.find( species ))->second;
}
template<int dim>
bool
pum_discretization<dim>::species_open_intersect_point( const string& species, const double* co ) const
{
	assert( is_species(species) );
	assert( co );
	if (!co) { return false; }
	
	return (dtree->access_geometry()).inside_region( get_species_support(species), co );
}

/**
	This was written before pum_discretization and allows different weights for each patch whereas
	pum_discretization keeps one weight for the whole species.
*/
template<int dim>
void
pum_discretization<dim>::get_patches_and_weights ( const string& species,
				    const vector<string>& vec_keys,
				       vector<box<dim> >& patches,
				       vector<const_pfunction>& patch_weights ) const
{
	int size = vec_keys.size();
	patches.resize(size);
	patch_weights.resize(size);
	
	typename vector<string>::const_iterator it(vec_keys.begin());
	typename vector<string>::const_iterator end(vec_keys.end());
	typename vector<box<dim> >::iterator box_it(patches.begin());
	typename vector<const_pfunction>::iterator f_it(patch_weights.begin());
	for ( ; it!=end; ++it, ++box_it, ++f_it ) {
		get_patch( species, *it, *box_it );
		get_weight( species, *f_it );
	}
}

template<int dim>
string
pum_discretization<dim>::find_patch_containing_point( const string& species, const double* co ) const
{
	assert( updated );
	
	string support_name = (species_support.find( species ))->second;
	
	// Check we have been give a species that is defined on a region.
	assert( (dtree->access_geometry()).is_region( support_name ) );
	
	return dtree->find_patch_containing_point( support_name, co );
}

template<int dim>
string
pum_discretization<dim>::get_species_support( const string& species ) const
{
	assert( is_species( species ) );
	return (species_support.find( species ))->second;
}

template<int dim>
void
pum_discretization<dim>::get_patch_keys( const string& species, set<string>& out ) const
{
	assert( is_species( species ) );
	string support = get_species_support( species );
	dtree->get_patch_keys( support, out );
}

template<int dim>
const set<string>&
pum_discretization<dim>::access_patch_keys( const string& species ) const
{
	assert( is_species( species ) );
	return dtree->access_patch_keys( get_species_support(species) );
}

template<int dim>
bool
pum_discretization<dim>::is_species_patch_key( const string& species, const string& key ) const
{
	// TODO: remove. Present for assert debug.
#if 0
	cout << "Support for species " << species << " is " << get_species_support( species ) << ".\n";
	if ( (dtree->access_geometry()).is_region( get_species_support( species ) ) ) {
		cout << get_species_support( species ) << " is a region.\n";
	} else {
		assert( (dtree->access_geometry()).is_boundary( get_species_support( species ) ) );
		cout << get_species_support( species ) << " is a boundary.\n";
	}
#endif
	// Debugging found that single keyword "return" was missing.
	return dtree->support_has_key( get_species_support( species ), key );
}
template<int dim>
void
pum_discretization<dim>::get_species_patch_neighbours( const string& species, const string& key, set<string>& out ) const
{
	if ( is_interior_species(species) ) {
		dtree->get_region_neighbours(get_species_support(species),key,out);
	} else {
		assert( is_boundary_species(species) );
		dtree->get_boundary_neighbours(get_species_support(species),key,out);
	}
}

template<int dim>
void
pum_discretization<dim>::get_patch( const string& species, const string& key, box<dim>& out ) const
{
	assert( is_species_patch_key( species, key ) );

	dtree->get_global_box( key, out );
}

template<int dim>
void
pum_discretization<dim>::get_patches( const string& species, vector<box<dim> >& cover ) const
{
	assert( updated );
	assert( is_species( species ) );
	string support_name = get_species_support( species );
	if ( (dtree->access_geometry()).is_boundary( support_name ) ) {
		// We only need to return boundary patches.
		dtree->get_boundary_patches( support_name, cover );
	} else {
		dtree->get_region_patches( support_name, cover );
	}
}


template<int dim>
void
pum_discretization<dim>::get_weight( const string& species, const_pfunction& weight ) const
{
	assert( updated );
	assert( is_species( species ) );
	assert( present_patch_weight.find( species ) != present_patch_weight.end() );
	
	weight = (present_patch_weight.find( species ))->second;
}

template<int dim>
void
pum_discretization<dim>::get_local_basis( const string& species,
			   const string& key,
			   vector<const_pfunction>& res ) const
{
	assert( is_species_patch_key( species, key ) );
	
	assert( present_local_basis.find( species ) != present_local_basis.end() );
	
	const map<string, vector<const_pfunction> >& species_local_basis =
			(present_local_basis.find( species ))->second;
			
	assert( species_local_basis.find( key ) != species_local_basis.end() );
			
	res = (species_local_basis.find( key ))->second;
}

template<int dim>
const vector<typename pum_discretization<dim>::const_pfunction>&
pum_discretization<dim>::access_local_basis( const string& species, const string& key ) const
{
	assert( is_species_patch_key( species, key ) );
	typename map<string, map<string, vector<const_pfunction> > >::const_iterator
			sp_it( present_local_basis.find( species ) );
	assert( sp_it!=present_local_basis.end());
	
	const map<string, vector<const_pfunction> >& species_local_basis = sp_it->second;
	
	typename map<string, vector<const_pfunction> >::const_iterator
			k_it( species_local_basis.find(key) );
	assert( k_it!=species_local_basis.end());
	return k_it->second;
}

template<int dim>
int
pum_discretization<dim>::get_num_dof() const
{
	int n=0;
	map<string, string>::const_iterator it( species_support.begin() );
	map<string, string>::const_iterator end( species_support.end() );
	for ( ; it!=end; ++it ) {
		n += get_num_dof( it->first );
	}
	return n;
}

template<int dim>
const valarray<int>&
pum_discretization<dim>::access_dof_allocation( int comm_size ) const
{
	assert( updated );
	return present_allocation_of_dof;
}

template<int dim>
void
pum_discretization<dim>::assemble_mass_matrix( Mat mass_matrix ) const
{
	// Get iterator that represents species and patches this computenode is
	// responsible for.
	typename map<int, map<string, set<string> > >::const_iterator tmp_it( present_parallel_partition.find(comm_rank));
	
	assert( tmp_it!=present_parallel_partition.end() );
	
	// Map representing species and patches this computenode is responsible for.
	const map<string, set<string> >& my_parallel_partition = tmp_it->second;
	assert( my_parallel_partition.size() > 0 );
	
	// Helper variables for use with PETSc.
	int m, n;
	valarray<int> idxm, idxn;
	valarray<double> values;
	
	set<string> species_patch_neighbours;
	string species;
	string first_key, second_key;
	
	// Iterate over species.
	typename map<string, set<string> >::const_iterator sp_it(my_parallel_partition.begin());
	typename map<string, set<string> >::const_iterator sp_end(my_parallel_partition.end());
	for ( ; sp_it!=sp_end; ++sp_it ) {
		// Iterate over patches of a species - outer patch key loop.
		typename set<string>::const_iterator k_it((sp_it->second).begin());
		typename set<string>::const_iterator k_end((sp_it->second).end());
		for ( ; k_it!=k_end; ++k_it ) {
			
			species = sp_it->first;
			first_key = *k_it;
			
			// Prepare the weight function for species here
			// in the outer loop.
			const_pfunction weight;
			get_weight( species, weight );
			
			m = access_local_basis( species, first_key ).size();
			assert( m>0 );
			idxm.resize(m);
			idxm[0] = get_global_index_start( species, first_key );
			// Index from 1.
			for ( int i=1; i<m; ++i ) {
				// Notice this is +i because we add to idxm[0]
				// rather than adding to idxm[i-1].
				// idxm[0]+i == idxm[i-1]+1;
				idxm[i] = idxm[0] + i;
			}
			assert( idxm[m-1] == idxm[0]+ m-1 );
			
			// Work through neighbours and fill in the rows for this patch *k_it.
			
			get_species_patch_neighbours(species, first_key, species_patch_neighbours);
			
			// To get testing against itself.
			species_patch_neighbours.insert( first_key );
			
			// Iterate over patches of a species - inner patch key loop.
			typename set<string>::const_iterator other_it( species_patch_neighbours.begin() );
			typename set<string>::const_iterator other_end( species_patch_neighbours.end() );
			for ( ; other_it!=other_end; ++other_it ) {
				second_key = *other_it;
				
				n = access_local_basis( species, second_key ).size();
				idxn.resize(n);
				idxn[0] = get_global_index_start( species, second_key );
				// Index from 1.
				for ( int i=1; i<n; ++i ) {
					idxn[i] = idxn[0] + i;
				}
				
				// We now know where to put the values. All that is left
				// is to calculate the values we'll put in.
				
				values.resize( m*n );
				
				vector<box<dim> > patches;
				dtree->get_patches( species_patch_neighbours, patches );
				
				int id_one = distance(species_patch_neighbours.begin(),species_patch_neighbours.find(first_key));
				
				int id_two = distance( species_patch_neighbours.begin(), species_patch_neighbours.find(second_key)  );
				
				vector<const_pfunction> local_approx_functions_one;
				vector<const_pfunction> local_approx_functions_two;
				
				get_local_basis(species, first_key, local_approx_functions_one);
				get_local_basis(species, second_key, local_approx_functions_two);
				
				assert( integrator );
				
				if ( dtree->region_key_is_boundary_key( get_species_support(species), first_key ) ) {
					interior_predicate<dim> pred( *this, species );
					
					integrator->compute_local_mass_matrix(  id_one,	id_two, patches, weight,
						local_approx_functions_one,
						local_approx_functions_two,
						&values[0], &pred );
				} else {
				
					integrator->compute_local_mass_matrix(  id_one,	id_two, patches, weight,
						local_approx_functions_one,
						local_approx_functions_two,
						&values[0], 0 );
				}
				
				// TODO: remove
#ifndef NDEBUG
				cout << "Mike idxm\n";
				for ( int i=0; i<m; ++i ) {
					cout << idxm[i] << ", ";
				}
				cout << "\nMike idxn\n";
				for ( int i=0; i<n; ++i ) {
					cout << idxn[i] << ", ";
				}
				cout << "\n";
#endif
				
				MatSetValues( mass_matrix, m, &(idxm[0]), n, &(idxn[0]), &values[0], INSERT_VALUES );
				
			}
		}
	}
	MatAssemblyBegin( mass_matrix, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( mass_matrix, MAT_FINAL_ASSEMBLY );

// TODO: evaluate whether removing is good or not.

#if 0 // ndef NDEBUG
	PetscInt rows;
	MatGetSize( mass_matrix, &rows, PETSC_NULL );
	if ( rows <= 12 ) {
		MatView( mass_matrix, PETSC_VIEWER_STDOUT_WORLD );
//		MatView( mass_matrix, PETSC_VIEWER_DRAW_WORLD );
	}
#endif

}

template<int dim>
void
pum_discretization<dim>::assemble_rhs_vector( Vec rhs_vec ) const
{	
	// Get iterator that represents species and patches this computenode is
	// responsible for.
	typename map<int, map<string, set<string> > >::const_iterator tmp_it( present_parallel_partition.find(comm_rank));
	
	assert( tmp_it!=present_parallel_partition.end() );
	
	// Map representing species and patches this computenode is responsible for.
	const map<string, set<string> >& my_parallel_partition = tmp_it->second;
	assert( my_parallel_partition.size() > 0 );
	
	// Helper variables for use with PETSc.
	int ni;
	valarray<int> ix;
	valarray<double> values;
	
	set<string> species_patch_neighbours;
	string species;
	string first_key;
	
	// Iterate over species.
	typename map<string, set<string> >::const_iterator sp_it(my_parallel_partition.begin());
	typename map<string, set<string> >::const_iterator sp_end(my_parallel_partition.end());
	for ( ; sp_it!=sp_end; ++sp_it ) {
		
		species = sp_it->first;
		
		const const_pfunction& species_ic = get_species_function( species );
		
		assert( species_ic->is_global_function() );
		
		// Iterate over patches of a species - outer patch key loop.
		typename set<string>::const_iterator k_it((sp_it->second).begin());
		typename set<string>::const_iterator k_end((sp_it->second).end());
		for ( ; k_it!=k_end; ++k_it ) {
			first_key = *k_it;
			
			// Prepare the weight function for species here
			// in the outer loop.
			const_pfunction weight;
			get_weight( species, weight );
			
			ni = access_local_basis( species, first_key ).size();
			assert( ni>0 );
			ix.resize(ni);
			ix[0] = get_global_index_start( species, first_key );
			// Index from 1.
			for ( int i=1; i<ni; ++i ) {
				ix[i] = ix[0] + i;
			}
			assert( ix[ni-1] == ix[0]+ ni-1 );

			values.resize( ni );

			// Work through neighbours and fill in the rows for this patch *k_it.
			
			get_species_patch_neighbours(species, first_key, species_patch_neighbours);
			
			// To get testing against itself.
			species_patch_neighbours.insert( first_key );
			
			vector<box<dim> > patches;
			dtree->get_patches( species_patch_neighbours, patches );
				
			int id_one = distance(species_patch_neighbours.begin(),species_patch_neighbours.find(first_key));
			
			vector<const_pfunction> local_approx_functions_one;
				
			get_local_basis(species, first_key, local_approx_functions_one);
			
			assert( integrator );
			
			if ( dtree->region_key_is_boundary_key( get_species_support(species), first_key ) ) {
				interior_predicate<dim> pred( *this, species );
				integrator->compute_local_projection_vector( species_ic,
								id_one,
								patches,
								weight,
								local_approx_functions_one,
								&values[0], &pred );
			} else {
				integrator->compute_local_projection_vector( species_ic,
								id_one,
								patches,
								weight,
								local_approx_functions_one,
								&values[0], 0 );
			}
			
			VecSetValues( rhs_vec, ni, &ix[0], &values[0], INSERT_VALUES );
		}
	}
	VecAssemblyBegin( rhs_vec );
	VecAssemblyEnd( rhs_vec );
	
	// TODO: evaluate whether removing is good or not.

#ifndef NDEBUG
	PetscInt rows;
	VecGetSize( rhs_vec, &rows );
	if ( rows <= 12 ) {
		cout << "Hello?";
		VecView( rhs_vec, PETSC_VIEWER_STDOUT_WORLD );
//		MatView( rhs_vec, PETSC_VIEWER_DRAW_WORLD );
	}
#endif

}

template<int dim>
void
pum_discretization<dim>::get_local_coefficients( const string& species,
					const valarray<double>& coeff,
				  const string& key,
      				  valarray<double>& out ) const
{
	assert( is_species( species ) );
	assert( have_coefficients );
	
	int size = get_local_basis_size( species, key );
	out.resize( size );
	int start = get_global_index_start( species, key );
	
	assert( start + size -1 < get_num_dof( species ) );
	
	for ( int i=0; i<size; ++i ) {
		out[i] = coeff[ start + i ];
	}
}

template<int dim>
void
pum_discretization<dim>::evaluate_partition_of_unity( const string& species,
					    const double* co,
					    valarray<double>& res ) const
{
	assert( updated );
	assert( is_species( species ) );
	
	assert( species_open_intersect_point(species,co) );
	if ( !species_open_intersect_point(species,co) ) {
		res.resize(0);
		assert( res.size() == 0 );
		return;
	}
	
	const_pfunction weight;
	get_weight( species, weight );
	
	assert( weight );

	assert ( co );
	
	string key = find_patch_containing_point( species, co );
	
	if ( key == "" ) {
		res.resize(0);
		assert( res.size() == 0 );
		return;
	}
	
	assert( key != "" );
	
	set<string> incident_keys;
	dtree->get_region_neighbours_incident_point( get_species_support(species), key, co, incident_keys );
	incident_keys.insert(key);
	
	assert( incident_keys.size() > 0 );
	
	vector<box<dim> > patches;
	dtree->get_patches( incident_keys, patches );
	assert( incident_keys.size() == patches.size() );
	
	// TODO: remove after finding cause of all points incident only on a single patch.
	// cout << "patches size is " << patches.size() << "\n";
	
	assert( patches.size() > 0 );
	
	valarray<double> values( patches.size() );
	
	typename vector<box<dim> >::const_iterator it( patches.begin() );
	typename vector<box<dim> >::const_iterator end( patches.end() );
	
	int counter = 0;
	
	assert( weight->is_local_function() );
	
	for ( ; it!=end; ++it, ++counter ) {
		assert( !it->empty() );
		if ( !it->closed_intersect_point( co ) ) {
			values[counter] = 0.0;
		} else {
			values[counter] = weight->global_evaluate( *it, co );
		}
	}
	
	double sum_all_weights = values.sum();
	
	// TODO: remove. negative weights returned as absolute value missing
	// in weight calculation.
#if 0
	cout << "values :\n";
	for ( int i=0; i<values.size(); ++i ) {
		cout << values[i] << ", ";
	}
	cout << "\nsum is " << sum_all_weights <<"\n";
#endif
	assert( sum_all_weights > 0.0 );

	int size = values.size();
	res.resize( size );
	for ( int i=0; i<size; ++i ) {
		res[i] = values[i]/sum_all_weights;
	}
}
template<int dim>
double
pum_discretization<dim>::evaluate_local_approximation( const string& species,
				const valarray<double>& coeff,
			      const string& key,
			      const double* co ) const
{
	assert( updated );
	
	assert( is_species_patch_key( species, key ) );
	
	assert( species_open_intersect_point(species,co) );
	if ( !species_open_intersect_point(species,co) ) {
		return 0.0;
	}
	
	box<dim> patch;
	get_patch( species, key, patch );
	assert( !patch.empty() );
	
	assert( patch.open_intersect_point(co) );
	
	vector<const_pfunction> single_local_basis;
	get_local_basis( species, key, single_local_basis );
	
	valarray<double> local_coefficients;
	get_local_coefficients( species, coeff, key, local_coefficients );
	
	assert( single_local_basis.size() == local_coefficients.size() );
	
	double res = 0.0;
	
	typename vector<const_pfunction>::const_iterator f_it( single_local_basis.begin() );
	typename vector<const_pfunction>::const_iterator f_end( single_local_basis.end() );
	const double* coeff_ptr = &local_coefficients[0];
	for ( ; f_it!=f_end; ++f_it, ++coeff_ptr ) {
		assert( (*f_it)->is_local_function() );
		
		// TODO: remove after small local approx values bug sorted.
		// cout << "believed coefficient is " << *coeff << "\n";
		// cout << "basis function eval value is " << (*f_it)->global_evaluate( patch, co ) << "\n";
		
		if ( patch.closed_intersect_point( co ) ) {
			res += (*coeff_ptr) * ( (*f_it)->global_evaluate( patch, co ) );
		}
			// Else do nothing.
			// res+= 0;
	}

	return res;
}

template<int dim>
void
pum_discretization<dim>::evaluate_local_approximation( const string& species, const valarray<double>& coeff, const double* co, valarray<double>& res ) const
{
	assert( updated );
	assert( is_species( species ) );
	
	assert( species_open_intersect_point(species,co) );
	if ( !species_open_intersect_point(species,co) ) {
		res.resize( 0 );
		return;
	}
	
	string key = find_patch_containing_point( species, co );
	
	if ( key == "" ) {
		res.resize( 0 );
		return;
	}
	
	assert( key != "" );
	
	set<string> incident_keys;
	dtree->get_region_neighbours_incident_point( get_species_support(species), key, co, incident_keys );
	incident_keys.insert(key);
	
	res.resize( incident_keys.size() );
	double* val = &res[0];
	
	typename set<string>::const_iterator sp_it( incident_keys.begin() );
	typename set<string>::const_iterator sp_end( incident_keys.end() );
	
	for ( ; sp_it!=sp_end; ++sp_it, ++val ) {
		*val = evaluate_local_approximation( species, coeff, *sp_it, co );
	}
}

template<int dim>
double
pum_discretization<dim>::evaluate_pum_discretization_approximation( const string& species,
const valarray<double>& coeff,
			    const double* co ) const
{
	assert( updated );
	
	assert( species_open_intersect_point(species,co) );
	if ( !species_open_intersect_point(species,co) ) {
		return 0.0;
	}
	
	valarray<double> tmp;
	valarray<double> local_approx;
	
	evaluate_partition_of_unity( species, co, tmp );

	evaluate_local_approximation( species, coeff, co, local_approx );
	
	assert( tmp.size() == local_approx.size() );

		
	if ( tmp.size() == 0 ) {
		// Point is not in the region of definition of the species.
		return 0.0;
	}
	
	
	evaluate_local_approximation( species, co, local_approx );
	
	assert ( tmp.size() > 0 );
	assert( tmp.size() == local_approx.size() );
	
	tmp *= local_approx;
	
	return tmp.sum();
}

template<int dim>
void
pum_discretization<dim>::evaluate_pum_discretization_approximation( const string& species,
					const valarray<double>& coeff,
					 const valarray<double>& in,
					 valarray<double>& res ) const
{
	assert( updated );
	assert( in.size() % dim == 0 );
	int num_pts = in.size()/dim;
	assert( num_pts > 0 );
	res.resize( num_pts );
	
	for ( int i=0; i<num_pts; ++i ) {
		assert( dim*i < in.size() );

		if ( species_open_intersect_point( species, &in[dim*i] ) ) {
			res[i] = evaluate_pum_discretization_approximation( species, &in[dim*i] );
		} else {
			res[i] = 0.0;
		}

	}
}

template<int dim>
void
pum_discretization<dim>::gp_draw_approximation( const string& species, const valarray<double>& coeff, ofstream& out, int n, double level )
{
	assert( updated );
	assert( n > 0 );
	valarray<double> eval_pts;
	valarray<double> results;
	make_point_vec_on_box<dim>( dtree->access_parent_box(), eval_pts, n, n, n );
	evaluate_pum_discretization_approximation( species, coeff, eval_pts, results );
	
	gp_draw_points<dim>( eval_pts, results, out, n, level );
}

//

template class pum_discretization<2>;
template class pum_discretization<3>;

