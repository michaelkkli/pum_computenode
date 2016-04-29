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

#include "pum_approx.hh"

#include "point_utils.hh"

#include <algorithm>
using std::min;
using std::copy;

#include <fstream>
using std::ofstream;

#include <cmath>
using std::pow;


template <int dim, typename key_t>
pum_approx<dim, key_t>::pum_approx()
{
}

template <int dim, typename key_t>
pum_approx<dim, key_t>::~pum_approx()
{
}

template <int dim, typename key_t>
pum_approx<dim, key_t>&
pum_approx<dim, key_t>::add_patch( const key_type& key, const box<dim>& bx )
{
  global_keys.insert( key );
  patches[key] = bx;
  register_change();
  return *this;
}

template <int dim, typename key_t>
int
pum_approx<dim, key_t>::num_local_coefficients_needed() const
{
  return local_matrix_size;
}

template <int dim, typename key_t>
void
pum_approx<dim, key_t>::give_local_coefficients ( const valarray<double>& c )
{
//  assert( c.size() == num_local_coefficients_needed );
  local_coefficients = c;
}

template <int dim, typename key_t>
double
pum_approx<dim, key_t>::get_coefficient( const key_type& key, int i ) const
{
	assert( local_keys.find( key ) != local_keys.end() );
	assert( i >= 0 );
	assert( i < neighbours[ key ].size() );
	return local_coefficients[ local_index_start[ key ] + i  ];
}

template <int dim, typename key_t>
double
pum_approx<dim, key_t>::seq_global_evaluate ( const double* co ) const
{
	// Implementation is not correct for parallel execution.
	
	assert( co );
	
	set<key_type>     incident;
	vector<box<dim> > patches;
	vector<const_pfunction> weights;
	
	key_type key = find_patch_containing_point( co );
	if ( global_keys.find( key ) == global_keys.end() ) {
		return 0.0;
	}
	get_neighbours( key, incident );
	incident.insert( key );
	
	vector<key_type> keys_v( incident.size() );
	
	copy( incident.begin(), incident.end(), keys_v.begin() );
	
	get_patches_and_weights( incident, patches, weights );
	
	valarray<double> weight_vals( incident.size() );
	
	// Precalculate weights.
	{
		int end = patches.size();
		for ( int i=0; i<end; ++i ) {
			weight_vals[i] = weights[i]->local_evaluate( patches[i], co );
		}
	}
	
	double sum_all_weights = weight_vals.sum();
	double result = 0.0;
	
	{
		vector<const_pfunction> local_approx;
		int global_start;
		double tmp;
		for ( int i=0; i<patches.size(); ++i ) {
			tmp = 0.0;
			get_local_approx_functions( keys_v[i], local_approx );
			if ( local_approx.size() == 0 ) {
				tmp = get_coefficient( keys_v[i], 0 );
			} else {
				global_start = global_index_start[ keys_v[i] ];
				for ( int j=0; j<local_approx.size(); ++j ) {
					tmp +=
						get_coefficient( keys_v[i],j)
							* local_approx[j]->local_evaluate( patches[i], co );
				}
			}
			result += tmp * weight_vals[i]; // Division by all weights left to the end.
		}
	}
	return result / sum_all_weights;
}

template <int dim, typename key_t>
void
pum_approx<dim, key_t>::output_gnuplot ( const char* filename,
				  int num_pts,
				  bool on_boundary ) const
{
  assert( filename );
  assert( num_pts > 0 );
  assert( coefficients_provided );
  int n_to_dim = dim;
  for ( int i=1; i<dim; ++i ) {
    n_to_dim *= dim;
  }
  valarray<double> pts( n_to_dim * dim  );
  if ( on_boundary ) {
    make_point_vec_on_box ( bounding_box, pts, num_pts, num_pts, num_pts );
  } else {
    make_point_vec_interior_box ( bounding_box, pts, num_pts, num_pts, num_pts );
  }

  valarray<double> vals( n_to_dim );
//must be replaced
 // evaluate_approximation( pts, vals );

  ofstream file ( filename );
  assert( file );

  for ( int p=0; p<n_to_dim; ++p ) {
    for ( int d=0; d<dim; ++d ) {
      file << pts[ p*dim +d ] << " ";
    }
    file << vals[p] << "\n";
  }
}



template <int dim, typename key_t>
key_t
pum_approx<dim, key_t>::find_patch_containing_point ( const double* co,
					       bool lowest_index ) const
{
	typename map<key_type, box<dim> >::const_iterator  it( patches.begin() );
	typename map<key_type, box<dim> >::const_iterator end( patches.end()   );
	for ( ; it != end; ++it ) {
		if ( (it->second).open_intersect_point( co ) ) {
			return it->first;
		}
	}
	return end;
}

// Perhaps useful for greater decomposition of the method.
#if 0
template <int dim, typename key_t>
double
pum_approx<dim, key_t>::evaluate_local_approx ( int p,
					 const double* co ) const
{
  assert( 0 <= p && p < patches.size() );
  assert( co );
  assert( patches[p].does_intersect( co ) );
  double tmp = 0.0;
  vector<const_pfunction>& local_approx = local_approx_functions[ p ];
  for ( int i=0; i<local_approx.size(); ++i ) {
    tmp += local_approx[i]->evaluate( co );
  }
  return tmp;
}
#endif

template <int dim, typename key_t>
void
pum_approx<dim, key_t>::register_change()
{
  updated               = false;
  coefficients_provided = false;
}

template <int dim, typename key_t>
void
pum_approx<dim, key_t>::update_neighbours()
{
  vector<key_t> global_keys_vec( global_keys.size() );
  {
	copy( global_keys.begin(), global_keys.end(), global_keys_vec.begin() );
	// Make sure ``neighbours'' is up to date.
	box<dim> tmp_box;
	for ( int i=0; i<global_keys_vec.size(); ++i ) {
		tmp_box = patches[ global_keys_vec[i] ];
		neighbours[global_keys_vec[i]].clear();
		for ( int j=0; j<i; ++j ) {
			if ( tmp_box.open_intersect_box( patches[global_keys_vec[i]] ) ) {
				neighbours[ global_keys_vec[i] ].insert( global_keys_vec[i] );
				neighbours[ global_keys_vec[j] ].insert( global_keys_vec[i] );
			}
		}
	}
  }
}

template <int dim, typename key_t>
void
pum_approx<dim, key_t>::update() {

  update_neighbours();
  
  int previous_index = 0;
  typename vector<key_t>::iterator it, end;
  it  = global_keys.begin();
  end = global_keys.end();
  int tmp;
  
// Make sure ``global_index_start'' is up to date.
  for ( ; it!=end; ++it ) {
	global_index_start[ it ] = previous_index;
	tmp = local_approx_functions[ it ].size();
	// The PU function is used for approximation if no local approx
	// functions are assigned.
	previous_index += tmp > 0 ? tmp : 1;
  }



#if 0
  int num_skipable_entries = global_index_start[ first_local_patch ];
  local_index_start.resize( num_local_patches );
  for ( int i=0; i<num_local_patches; ++i ) {
    local_index_start[i] =
      global_index_start[ first_local_patch + i ] - num_skipable_entries;
  }
#endif


// Load balancing not done yet. Local_keys is unpopulated.
// Should be populated depending on load balancing.
  num_local_patches = local_keys.size();

  updated = true;
}



//

template class pum_approx<1, string>;
template class pum_approx<2, string>;
template class pum_approx<3, string>;
