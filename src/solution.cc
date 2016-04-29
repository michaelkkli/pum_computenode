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

#include "solution.hh"

template<int dim>
solution<dim>::solution ()
	: initialized_(false)
{
}

template <int dim>
int
solution<dim>::set_global_approximation_space ( const typename global_approximation_space<dim>::ptr & gas, int size_hint )
{
	if ( gas ) {
		global_approx_space_ = gas;
		if ( size_hint>0 ) {
			global_coefficients_.resize ( size_hint );
		} else {
			valarray<int> tmp;
			global_approx_space_->get_patch_to_num_dof ( tmp );
			
			int num_dof = tmp.sum();
			
			global_coefficients_.resize ( num_dof );
			
			return num_dof;
		}
	}
}

template <int dim>
void
solution<dim>::set_dof_structure ( const typename dof_structure::ptr & ds )
{
	assert ( ds );
	if ( ds ) {
		dof_struct_ = ds;
	}
}

template <int dim>
void
solution<dim>::set_coefficients ( valarray<double>& in )
{
	if ( in.size() != global_coefficients_.size() ) {
		abort();
	}
	global_coefficients_ = in;
}

#if 0
template <int dim>
void 
solution<dim>::global_evaluate ( valarray<double> & in, valarray<double>& out )
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	
	int num_pts = in_size/dim;
	
	out.resize ( num_pts );
	
	vector<string>            incident_patch_keys;
	vector<int>               incident_patch_indices;
	vector<const double*>     coeff_starts;
	
	assert ( global_approx_space_ );
	
	// g++-4.2.3
	// src/solution.cc:79: internal compiler error: in lower_stmt, at gimple-low.c:282
	//#pragma omp parallel for private(incident_patch_keys,incident_patch_indices,coeff_starts)
	for ( int i=0; i<num_pts; ++i ) {
		double* curr = &in[ dim*i ];
		
		global_approx_space_->get_incident_patch_keys_and_indices ( curr, incident_patch_keys, incident_patch_indices );
		
		int num_incident = incident_patch_keys.size();
		
		assert ( incident_patch_indices.size() == num_incident );
		
		if ( num_incident == 0 ) {
			// Point outside cover.
			out[i] = 0.0;
			continue;
		}
		
		coeff_starts.clear();
		coeff_starts.resize( num_incident );
		
		assert ( dof_struct_ );
		
		for ( int j=0; j<num_incident; ++j ) {
			int start_index = dof_struct_->patch_to_start_index[ incident_patch_indices[j] ];
			coeff_starts[j] = &global_coefficients_[ start_index ];
		}
		
		out[i] = global_approx_space_->global_evaluate ( curr, incident_patch_keys, coeff_starts );
	}
}
#endif

template<int dim>
void solution<dim>::global_evaluate ( valarray<double> &                 pts,
                                      solution_evaluation_structure &    ses,
                                      valarray<double> &                 out )
{
	int pts_size = pts.size();
	assert ( pts_size % dim == 0 );
	
	int num_pts = pts_size/dim;
	
	out.resize ( num_pts );
	
	
	
	assert ( global_approx_space_ );
	
	
	
	// g++-4.2.3
	// src/solution.cc:79: internal compiler error: in lower_stmt, at gimple-low.c:282
#pragma omp parallel for
	for ( int i=0; i<num_pts; ++i ) {
		// In here purely for benefit of g++4.3
		vector<const double*>     coeff_starts;
		
		
		const vector<string> &    incident_patch_keys    = ses.eval_keys[i];
		const vector<int>    &    incident_patch_indices = ses.eval_indices[i];
		
		int num_incident = incident_patch_keys.size();
		
		assert ( incident_patch_indices.size() == num_incident );
		
		if ( num_incident == 0 ) {
			// Point outside cover.
			out[i] = 0.0;
			continue;
		}
		
		coeff_starts.clear();
		coeff_starts.resize( num_incident );
		
		assert ( dof_struct_ );
		
		for ( int j=0; j<num_incident; ++j ) {
			int start_index = dof_struct_->patch_to_start_index[ incident_patch_indices[j] ];
			coeff_starts[j] = &global_coefficients_[ start_index ];
		}
		out[i] = global_approx_space_->global_evaluate ( &pts[dim*i], incident_patch_keys, coeff_starts );
	}
}

#if 0
template<int dim>
void
solution<dim>::set_coefficients ( valarray<double>& rhs )
{
	global_coefficients_.resize( rhs.size() );
	global_coefficients_ = rhs;
	assert( rhs.size() > 0 ) ;
	have_coefficients_ = true;
}
#endif

//

template class solution<2>;
