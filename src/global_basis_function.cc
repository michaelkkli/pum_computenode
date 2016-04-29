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

#include "box_utils.hh" // Presently only for debug.
#include "global_basis_function.hh"
#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>

using std::abs;
using std::cout;
using std::fabs;
using std::fill;
using std::slice;


template <int dim>
void
global_basis_function<dim>::set ( const typename box<dim>::vec_ptr&                 patches,
                                  int                                               main_patch_index,
                                  const typename differentiable_function<dim>::ptr& patch_weight,
                                  const typename differentiable_function<dim>::ptr& local_approx )
{
	int num = patches.size();
	assert( 0<=main_patch_index && main_patch_index<num );

	typename box<dim>::vec_ptr neighbours;
	for ( int i=0; i<num; ++i ) {
		if ( i!= main_patch_index ) {
			neighbours.push_back( patches[i] );
		}
	}

	typename partition_of_unity_function<dim>::ptr tmp_pu_function( new partition_of_unity_function<dim> );
	tmp_pu_function->set( patches[main_patch_index], neighbours, patch_weight );

	this->set( tmp_pu_function, local_approx );
}

template <int dim>
void
global_basis_function<dim>::set( const typename partition_of_unity_function<dim>::ptr& pu_fn,
                                 const typename differentiable_function<dim>::ptr&     local_approx )
{
	assert( pu_fn && local_approx );
	pu_function_           = pu_fn;
	local_approx_function_ = local_approx;
}

template <int dim>
const box<dim>&
global_basis_function<dim>::access_box () const
{
	assert ( pu_function_ );
	return pu_function_->access_box();
}

// Currently used by petsc_solver-3 with double global basis function version used by pum_convergence-1.
template <int dim>
bool
global_basis_function<dim>::get_integration_domain_decomposition( integration_domain_decomposition<dim>& out,
                                                                  int                                    decomposition_level ) const
{	
	switch ( decomposition_level ) {
	case 0:
		out.box_vector.resize(1);
		out.box_vector[0] = pu_function_->access_box();
		
		return true;
	case 1:
		// Only uses current patch information. 2009-05-24 ML.
		pu_function_->get_patch_decomposition_by_neighbours( out.box_vector );
		
		return true;
	default:
		assert( "No such decomposition_level" && 0 );
		return false;
	}
}



template <int dim>
bool
global_basis_function<dim>::get_integration_domain_decomposition( const global_basis_function<dim>& first,
                                                          const global_basis_function<dim>& second,
                                                          integration_domain_decomposition<dim>& out,
                                                          int level )
{
	switch( level ) {
	case 0:
		out.box_vector.resize(1);
		out.box_vector[0] = first.access_box();
		out.box_vector[0].clip_against ( second.access_box() );
#if 0 //ndef NDEBUG // Keep for future debug.
		std::clog << "First box is ";
		output_description ( first.access_box(), std::clog );
		std::clog << ". Second box is ";
		output_description ( second.access_box(), std::clog );
		std::clog << "\n";
		assert ( !out.box_vector[0].empty() );
#endif
		return true;
	case 1:
	{
		const partition_of_unity_function<dim>& first_pu = *(first.pu_function_);
		const partition_of_unity_function<dim>& second_pu = *(second.pu_function_);

		partition_of_unity_function<dim>::get_patch_decomposition_by_neighbours( first_pu, second_pu, out.box_vector );
	}
		return true;
	default:
		assert( "Invalid level." && 0 );
		return false;
	}
}

template <int dim>
double
global_basis_function<dim>::global_evaluate( const double* co ) const
{
	assert( pu_function_ && local_approx_function_ );
	double tmp = pu_function_->global_evaluate( co );
	return tmp * local_approx_function_->global_evaluate( pu_function_->access_box(), co );
}

template <int dim>
double
global_basis_function<dim>::local_evaluate( const double* co ) const
{
	assert( pu_function_ && local_approx_function_ );
	double tmp = pu_function_->local_evaluate( co );
	
	double local_approx_value = local_approx_function_->local_evaluate( pu_function_->access_box(), co );
	
#if 0 //ndef NDEBUG // Keep for future debug.
	std::clog << "Partition of unity function evaluation is " << tmp
		<< " and local approx is " << local_approx_value << ".\n";
	
	// Monomials should be small.
	assert ( fabs ( static_cast<long double>( local_approx_value ) )<3.0 );
#endif
	
	return tmp * local_approx_value;
}

template <int dim>
void
global_basis_function<dim>::global_evaluate_grad( const double * co, double* out ) const
{
	pu_function_->global_evaluate_grad( co, out );

	valarray<double> tmp(dim);

	local_approx_function_->global_evaluate_grad( pu_function_->access_box(), co, &tmp[0] );	

	for ( int d=0; d<dim; ++d ) {
		out[d] += tmp[d];
	}
}

template <int dim>
void
global_basis_function<dim>::local_evaluate_grad ( const double * co, double* grad ) const
{
	assert ( co && grad );
	valarray<double> pu_grad ( dim );
	pu_function_->local_evaluate_grad( co, &pu_grad[0] );

	valarray<double> local_grad (dim);
	const box<dim>& patch = pu_function_->access_box();
	local_approx_function_->local_evaluate_grad( patch,
	                                             co,
	                                             &local_grad[0] );
	
	pu_grad    *= local_approx_function_->local_evaluate ( patch, co );
	local_grad *= pu_function_->local_evaluate ( co );
	
	assert ( pu_grad.size() == local_grad.size() );
	
	valarray<double>& product_rule = pu_grad;
	product_rule += local_grad;

	for ( int d=0; d<dim; ++d ) {
		grad[d] = product_rule[d];
	}
}

template <int dim>
void
global_basis_function<dim>::local_evaluate_grad ( valarray<double> &                 in,
                                                  global_basis_function_scratch &    scratch,
                                                  valarray<double> &                 out_grad ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size / dim;
	
	out_grad.resize ( in_size );
	
	assert ( pu_function_ );
	const box<dim>& patch = pu_function_->access_box();
	
	// Evaluate gradients of partition of unity function and local approximation function.
	valarray<double>& pu_grad = out_grad; // Make use of passed-in object.
	
	valarray<double> &    local_approx_grad = scratch.local_approx_grad;
	local_approx_grad.resize ( in_size );

	pu_function_->local_evaluate_grad ( in, scratch.pu_fun_scratch, pu_grad );
	local_approx_function_->local_evaluate_grad( patch, in, local_approx_grad );
	
	// Evaluate values of partition of unity function and local approximation function.
	valarray<double> &    pu_vals = scratch.pu_vals;
	pu_vals.resize ( num_pts );
	valarray<double> &    local_approx_vals = scratch.local_approx_vals;
	local_approx_vals.resize ( num_pts );
	
	
	pu_function_->local_evaluate ( in, scratch.pu_fun_scratch, pu_vals );
	local_approx_function_->local_evaluate ( patch, in, local_approx_vals );
	
	for ( int i=0; i<num_pts; ++i ) {
		for ( int d=0; d<dim; ++d ) {
			pu_grad[ dim*i+d ] *= local_approx_vals[i];
		}
		// pu_grad[ slice( dim*i, dim, 1 ) ] *= valarray<double>( local_approx_vals[i], dim );
	}
	
	for ( int i=0; i<num_pts; ++i ) {
		for ( int d=0; d<dim; ++d ) {
			local_approx_grad[dim*i+d] *= pu_vals[i];
		}
		// local_approx_grad[ slice( dim*i, dim, 1 ) ] *= valarray<double>( pu_vals[i], dim);
	}
	
	// Remember pu_grad is the same as out_grad at this point.
	assert ( out_grad.size() == local_approx_grad.size() );
	out_grad += local_approx_grad;
	
	assert ( out_grad.size() == num_pts*dim );
	
	// Internal consistency testing.
#ifndef NDEBUG
	{
		double grad[dim];
		for ( int i=0; i<num_pts; ++i ) {
			local_evaluate_grad ( &in[dim*i], grad );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs ( out_grad[dim*i+d] - grad[d] ) < 1e-10 );
			}
		}
	}
#endif
}

template <int dim>
void
global_basis_function<dim>::local_evaluate( valarray<double> &                 in,
                                            global_basis_function_scratch &    scratch,
                                            valarray<double> &                 out ) const
{	
	// Profiling has shown this is the perfect place to attempt optimization.

	assert( pu_function_ && local_approx_function_ );
	
	int in_size = in.size();
	assert( in_size%dim == 0 );

	int num_pts = in_size/dim;
	out.resize( num_pts );

	pu_function_->local_evaluate( in, scratch.pu_fun_scratch, out );
	assert ( out.max() < 1.0 + 1e-10 );
	
	valarray<double> &    local_vals = scratch.local_approx_vals;
	local_vals.resize ( num_pts );

	local_approx_function_->local_evaluate( pu_function_->access_box(), in, local_vals );
	
	assert ( local_vals.size() == out.size() );
	out *= local_vals;
	
	// Internal consistency testing.
#ifndef NDEBUG
	for ( int i=0; i<num_pts; ++i ) {
		assert ( fabs ( out[i] - local_evaluate ( &in[dim*i] ) ) < 1e-10 );
	}
#endif
}


// Highly inefficient and should be removed.
// Making sure not used anywhere.
// Appears not to be. 2009-02-20 ML.
#if 0
template <int dim>
void
global_basis_function<dim>::local_evaluate_and_grad ( valarray<double> &                in,
                                                      global_basis_function_scratch &   scratch,
                                                      valarray<double> &                out_vals,
                                                      valarray<double> &                out_grad ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size / dim;
	
	out_grad.resize ( in_size );
	
	assert ( pu_function_ );
	const box<dim>& patch = pu_function_->access_box();

	// Evaluate values of partition of unity function and local approximation function.
	valarray<double> &    pu_vals = scratch.pu_vals;
	pu_vals.resize ( num_pts );
	valarray<double> &    local_approx_vals = scratch.local_approx_vals;
	local_approx_vals.resize ( num_pts );
	
	// Evaluate gradients of partition of unity function and local approximation function.
	valarray<double>& pu_grad = scratch.pu_grad;
	
	valarray<double> &    local_approx_grad = scratch.local_approx_grad;
	local_approx_grad.resize ( in_size );

	// pu_function_->local_evaluate ( in, scratch.pu_fun_scratch, pu_vals );
	// pu_function_->local_evaluate_grad ( in, scratch.pu_fun_scratch, pu_grad );
	pu_function_->local_evaluate_and_grad ( in, scratch.pu_fun_scratch, pu_vals, pu_grad );
	
	local_approx_function_->local_evaluate ( patch, in, local_approx_vals );
	local_approx_function_->local_evaluate_grad( patch, in, local_approx_grad );
	
// Stepwise replacement of code.
#if 1
	local_evaluate_and_grad ( in,
	                          pu_vals,
	                          pu_grad,
	                          local_approx_vals,
                                  local_approx_grad,
                                  scratch,
                                  out_vals,
                                  out_grad );


#else
	out_vals.resize ( pu_vals.size() );
	out_vals = pu_vals;
	out_vals *= local_approx_vals;
	
	out_grad.resize ( pu_grad.size() );
	out_grad = pu_grad;
	
	for ( int i=0; i<num_pts; ++i ) {
		for ( int d=0; d<dim; ++d ) {
			out_grad[ dim*i+d ] *= local_approx_vals[i];
		}
		// pu_grad[ slice( dim*i, dim, 1 ) ] *= valarray<double>( local_approx_vals[i], dim );
	}
	
	for ( int i=0; i<num_pts; ++i ) {
		for ( int d=0; d<dim; ++d ) {
			local_approx_grad[dim*i+d] *= pu_vals[i];
		}
		// local_approx_grad[ slice( dim*i, dim, 1 ) ] *= valarray<double>( pu_vals[i], dim);
	}
	
	assert ( out_grad.size() == local_approx_grad.size() );
	out_grad += local_approx_grad;
	
	assert ( out_grad.size() == num_pts*dim );
#endif
	
	// Internal consistency testing.
#ifndef NDEBUG
	{
		double val;
		double grad[dim];
		for ( int i=0; i<num_pts; ++i ) {
			val = local_evaluate ( &in[dim*i] );
			assert ( fabs ( out_vals[i] - val ) < 1e-10 );
			local_evaluate_grad ( &in[dim*i], grad );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs ( out_grad[dim*i+d] - grad[d] ) < 1e-10 );
			}
		}
	}
#endif
}
#endif

template <int dim>
void
global_basis_function<dim>::local_evaluate_and_grad ( valarray<double> &                in,
                                                      valarray<double> &                pu_vals,
                                                      valarray<double> &                pu_grad,
                                                      valarray<double> &                local_approx_vals,
                                                      valarray<double> &                local_approx_grad,
                                                      global_basis_function_scratch &   scratch,
                                                      valarray<double> &                out_vals,
                                                      valarray<double> &                out_grad ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size / dim;

	out_vals.resize ( pu_vals.size() );
	out_vals = pu_vals;
	out_vals *= local_approx_vals;

	out_grad.resize ( pu_grad.size() );
	out_grad = pu_grad;
	
	for ( int i=0; i<num_pts; ++i ) {
		for ( int d=0; d<dim; ++d ) {
			out_grad[ dim*i+d ] *= local_approx_vals[i];
		}
		// pu_grad[ slice( dim*i, dim, 1 ) ] *= valarray<double>( local_approx_vals[i], dim );
	}
	
	for ( int i=0; i<num_pts; ++i ) {
		for ( int d=0; d<dim; ++d ) {
			out_grad[dim*i+d] += local_approx_grad[dim*i+d] * pu_vals[i];
		}
		// local_approx_grad[ slice( dim*i, dim, 1 ) ] *= valarray<double>( pu_vals[i], dim);
	}
	
	// Internal consistency testing.
#if 0// ndef NDEBUG // Accepted. 2009-07-06 ML.
	{
		double val;
		double grad[dim];
		for ( int i=0; i<num_pts; ++i ) {
			val = local_evaluate ( &in[dim*i] );
			assert ( fabs ( out_vals[i] - val ) < 1e-10 );
			local_evaluate_grad ( &in[dim*i], grad );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs ( out_grad[dim*i+d] - grad[d] ) < 1e-10 );
			}
		}
	}
#endif
}

template <int dim>
void
global_basis_function<dim>::pu_local_evaluate_and_grad ( valarray<double> &                 in,
                                                         global_basis_function_scratch &    scratch,
                                                         valarray<double> &                 vals,
                                                         valarray<double> &                 grad ) const
{
	pu_function_->local_evaluate_and_grad ( in, scratch.pu_fun_scratch, vals, grad );
}

template <int dim>
void
global_basis_function<dim>::local_approx_local_evaluate_and_grad ( valarray<double> &                 in,
                                                                   global_basis_function_scratch &    scratch,
                                                                   valarray<double> &                 vals,
                                                                   valarray<double> &                 grad ) const
{
	assert ( pu_function_ );
	const box<dim>& patch = pu_function_->access_box();
	
	assert ( local_approx_function_ );
	local_approx_function_->local_evaluate ( patch, in, vals );
	local_approx_function_->local_evaluate_grad( patch, in, grad );
}

//

template class global_basis_function<2>;
