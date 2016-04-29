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

#include "partition_of_unity_function.hh"

#include "box_utils.hh"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#if 0
#undef TRACE
#define TRACE cout << __FILE__ << ":" << __LINE__ << "\n";
#endif

using std::abs;
using std::cout;
using std::copy;
using std::exp;
using std::fabs;
using std::log;
using std::pow;
using std::slice;

template <int dim>
partition_of_unity_function<dim>::partition_of_unity_function ()
{
}

template <int dim>
partition_of_unity_function<dim>::~partition_of_unity_function ()
{
}

template <int dim>
const box<dim>&
partition_of_unity_function<dim>::access_box() const
{
	assert( patch_ );
	return *patch_;
}

template <int dim>
void
partition_of_unity_function<dim>::set( const typename box<dim>::ptr& ptch,
                                       const typename box<dim>::vec_ptr& np,
                                       const typename differentiable_function<dim>::ptr& pw )
{
	assert( ptch && pw );
	
	patch_             = ptch;
	neighbour_patches_ = np;
	patch_weight_      = pw;
	
}

#if 0
template <int dim>
void
partition_of_unity_function<dim>::set( const typename box<dim>::ptr& ptch,
                                       const typename box<dim>::vec_ptr& np,
                                       const typename differentiable_function<dim>::ptr& pw,
                                       const typename differentiable_function<dim>::vec_ptr& nw )
{
	assert( ptch && pw );

	patch_             = ptch;
	neighbour_patches_ = np;
	patch_weight_      = pw;
	neighbour_weights_ = nw;
}
#endif


template <int dim>
void
partition_of_unity_function<dim>::get_patch_decomposition_by_neighbours ( vector<box<dim> >& out ) const
{
#if 0
	// FIXME: delaying design decision.

	int num = neighbour_patches_.size();
	vector<box<dim> > tmp ( num );

	for ( int i=0; i<num; ++i ) {
		tmp[i] = *(neighbour_patches_[i]);
		assert( !tmp[i].empty() );
	}
#endif

	decompose_box( *patch_, neighbour_patches_, out );
}

template <int dim>
void
partition_of_unity_function<dim>::get_patch_decomposition_by_neighbours ( const box<dim>& bx, vector<box<dim> >& out ) const
{
#if 0
	// FIXME: delaying design decision.

	int num = neighbour_patches_.size();
	vector<box<dim> > tmp ( num );

	for ( int i=0; i<num; ++i ) {
		tmp[i] = *(neighbour_patches_[i]);
		assert( !tmp[i].empty() );
	}
#endif

	decompose_box( bx, neighbour_patches_, out );
}


template <int dim>
void
partition_of_unity_function<dim>::get_patch_decomposition_by_neighbours ( const partition_of_unity_function<dim>& first,
	                                                                  const partition_of_unity_function<dim>& second,
	                                                                  vector<box<dim> >& out )
{
	box<dim> restricted = *(first.patch_);
	restricted.clip_against( *(second.patch_) );

	const typename box<dim>::vec_ptr& first_vec  = first.neighbour_patches_;
	const typename box<dim>::vec_ptr& second_vec = second.neighbour_patches_;

	int first_size = first_vec.size();
	int total_num  = first_size + second_vec.size();
	typename box<dim>::vec_ptr combined( total_num );
	
	copy( first_vec.begin(), first_vec.end(), combined.begin() );
	copy( second_vec.begin(), second_vec.end(), combined.begin() + first_size );

	decompose_box( restricted, combined, out );
}

template <int dim>
double
partition_of_unity_function<dim>::global_evaluate ( const double* co ) const
{

	assert( patch_ );

	assert ( patch_->closed_intersect_point( co ) );
	
	if ( !patch_->closed_intersect_point( co ) ) {
#if 0 //ndef NDEBUG // Keep for future debug.
		std::clog << "Call to partition_of_unity_function global_evaluate with point outside patch.\n";
#endif
		return 0.0;
	}

	double tmp = 0.0;

	int num = neighbour_patches_.size();

	for ( int i=0; i<num; ++i ) {
		if ( (neighbour_patches_[i])->closed_intersect_point( co ) ) {
			tmp += patch_weight_->global_evaluate( *(neighbour_patches_[i]), co );
		}
	}

	assert( patch_ && patch_weight_ );

	double val = patch_weight_->global_evaluate( *patch_, co );
	
#if 0 //ndef NDEBUG // Keep for future debug.
	std::clog << "global evaluate called on patch_weight_ to give " << val << ", val/tmp is " << val/tmp << "\n";
#endif

	// Major change so that a patch is considered its own neighbour
	// in order to facilitate ease of assembly.
	// return  val /( tmp + val );

	assert ( val/tmp < 1 + 1e-10 );
	
	return val/tmp;
}

template <int dim>
double
partition_of_unity_function<dim>::local_evaluate ( const double* co ) const
{
	assert( patch_ );

	int num = neighbour_patches_.size();

	double sum_weights =0.0;

	double tmp_pt[dim];

	for ( int i=0; i<num; ++i ) {
		patch_->map_local_to_global( co, tmp_pt );

		// Minimize calls to closed_intersect_point. In particular,
		// avoid leaving assert calls involving closed_intersect_point.
		
		if ( (neighbour_patches_[i])->closed_intersect_point( tmp_pt ) ) {

			(neighbour_patches_[i])->map_global_to_local( tmp_pt, tmp_pt );

			sum_weights += patch_weight_->local_evaluate( *(neighbour_patches_[i]), tmp_pt );
		}
	}

	double this_weight = patch_weight_->local_evaluate( *patch_, co );

	// N.B. a patch is considered its own neighbour
	// in order to facilitate ease of assembly.
	
	return this_weight/sum_weights;
}

template <int dim>
void 
partition_of_unity_function<dim>::local_evaluate ( valarray<double>&        in,
                                                   pu_function_scratch &    scratch,
                                                   valarray<double>&        out ) const
{
	assert( patch_ );
	
	int num_neigh = neighbour_patches_.size();
	
	int in_size = in.size();
	
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	out.resize ( num_pts );

	valarray<double> &    global_pts = scratch.global;
	patch_->map_local_to_global ( in, global_pts );
	
	// Evaluate values even if not needed. Callgrind
	// shows attempt at reducing separate calls may be worthwhile.
	valarray<double> &    local_pts = scratch.local;
	local_pts.resize ( in_size );
	
	valarray<double> &    sum_weights = scratch.weight_sum;
	sum_weights.resize ( num_pts, 0.0 );
	
	valarray<bool> &    pred = scratch.pred;
	pred.resize ( num_pts, false );
	
	//	double local_pt[dim];
	
	for ( int n=0; n<num_neigh; ++n ) {
		(neighbour_patches_[n])->closed_intersect_point ( global_pts, pred );
		
		(neighbour_patches_[n])->map_global_to_local ( global_pts, pred, local_pts );
		
		assert ( patch_weight_ );
		
		patch_weight_->local_evaluate ( *(neighbour_patches_[n]), local_pts, pred, out );
		
		assert ( sum_weights.size() == out.size() );
		
#ifndef NDEBUG
		for ( size_t g=0; g<out.size(); ++g ) {
			assert ( out[g] > -1e-100 );
		}
#endif
		
		sum_weights += out;
	}
	
	patch_weight_->local_evaluate ( *patch_, in, out );
	
	assert ( sum_weights.min() > -1e-10 );
	
	out /= sum_weights;
	
	assert ( (-1e-10<out.min()) && ( out.max()<1+1e-10));
	
	// Internal consistency testing.
#if 0//ndef NDEBUG	
	for ( int i=0; i<num_pts; ++i ) {
		assert ( fabs( out[i] - local_evaluate ( &in[dim*i] ) ) < 1e-10 );
	}
#endif
}

template <int dim>
void
partition_of_unity_function<dim>::global_evaluate_grad( const double* co, double* out ) const
{
	assert( co && out );

	double weight_sum = 0.0;
	valarray<double> grad_sum(0.0, dim);
	valarray<double> grad_tmp(dim);

	int num = neighbour_patches_.size();

	for ( int i=0; i<num; ++i ) {
		if ( (neighbour_patches_[i])->closed_intersect_point( co ) ) {
			weight_sum += patch_weight_->global_evaluate( *(neighbour_patches_[i]), co );

			patch_weight_->global_evaluate_grad( *(neighbour_patches_[i]), co, &grad_tmp[0] );

			grad_sum += grad_tmp;
		}
	}

	double main_weight = patch_weight_->global_evaluate( *patch_, co );

	valarray<double> main_grad;

	patch_weight_->global_evaluate_grad( *patch_, co, &main_grad[0] );

	assert( main_grad.size() == dim );

	double weight_sum_squared = weight_sum*weight_sum;

	for ( int d=0; d<dim; ++d ) {
		out[d] = (main_grad[d]*weight_sum - main_weight*grad_sum[d])/weight_sum_squared;
	}
}

template <int dim>
void
partition_of_unity_function<dim>::local_evaluate_grad( const double* co, double* out ) const
{
	assert( co && out );
	
	double weight_sum = 0.0;
	valarray<double> grad_sum(0.0, dim);
	valarray<double> grad_tmp(dim);

	int num_neigh = neighbour_patches_.size();

	double global[dim];
	double local[dim];

	for ( int i=0; i<num_neigh; ++i ) {
		
		patch_->map_local_to_global( co, global );
		if ( (neighbour_patches_[i])->closed_intersect_point( global ) ) {

			assert ( neighbour_patches_[i] && !neighbour_patches_[i]->empty() );
			
			// Lots of heartache with this line missing. Only found when lots of repeated
			// values were seen within this loop.
			// 2008-09-02 Mike Li.
			neighbour_patches_[i]->map_global_to_local ( global, local );
			
			weight_sum += patch_weight_->local_evaluate( *(neighbour_patches_[i]), local );

			patch_weight_->local_evaluate_grad( *(neighbour_patches_[i]), local, &grad_tmp[0] );

			grad_sum += grad_tmp;
		}
	}

	double main_weight = patch_weight_->local_evaluate( *patch_, co );

	valarray<double> main_grad ( dim );

	assert ( patch_weight_ );
	patch_weight_->local_evaluate_grad( *patch_, co, &main_grad[0] );

	assert( main_grad.size() == dim );
	
	double weight_sum_squared = weight_sum*weight_sum;
	
	for ( int d=0; d<dim; ++d ) {
		out[d] = main_grad[d]/weight_sum - (main_weight*grad_sum[d])/weight_sum_squared;
	}
}

/**
 * Application of the gradient quotient rule to the partition of unity function.
 * phi is
 * 
 * 	W_i / ( sum_j W_j ) so gradient is
 * 
 * 	( (grad W_i)(sum_j W_j) - (W_i)(sum_j grad W_j) )/( sum_j W_j)^2
 */
template <int dim>
void
partition_of_unity_function<dim>::local_evaluate_grad ( valarray<double>&        in,
                                                        pu_function_scratch &    scratch,
                                                        valarray<double>&        out_grad ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;

	valarray<double> &    weight_sum = scratch.weight_sum;
	weight_sum.resize ( num_pts, 0.0 );
	
	valarray<double> &    grad_sum = scratch.grad_sum;
	grad_sum.resize ( in_size, 0.0 );
	
	int num_neigh = neighbour_patches_.size();
	
	valarray<double> &    global = scratch.global;
	global.resize ( in_size );
	patch_->map_local_to_global( in, global );
	
	valarray<bool> &    pred = scratch.pred;
	pred.resize ( num_pts, false );
	
	valarray<double> &    local = scratch.local;
	local.resize ( in_size );
	
	valarray<double> &    weight_tmp = scratch.weight_tmp;
	weight_tmp.resize ( num_pts );
	
	valarray<double> &    grad_tmp = scratch.grad_tmp;
	grad_tmp.resize ( in_size );
	
	valarray<double> &    main_weight = scratch.main_weight;
	main_weight.resize ( num_pts );
	
	// If this is changed, be sure to change the bit below too.
	valarray<double> & main_grad = out_grad;

	
	for ( int i=0; i<num_neigh; ++i ) {
		neighbour_patches_[i]->closed_intersect_point( global, pred );
		
		neighbour_patches_[i]->map_global_to_local ( global, pred, local );
		
		patch_weight_->local_evaluate( *(neighbour_patches_[i]), local, pred, weight_tmp );
		
		assert ( weight_sum.size() == weight_tmp.size() );
		weight_sum += weight_tmp;
		
		patch_weight_->local_evaluate_grad( *(neighbour_patches_[i]), local, pred, grad_tmp );
		
		assert ( grad_sum.size() == grad_tmp.size() );
		grad_sum   += grad_tmp;
	}
	
	patch_weight_->local_evaluate_grad( *patch_, in, main_grad );
	patch_weight_->local_evaluate( *patch_, in, main_weight );
	
	
	// 10% time in function allocating and deallocating unsigned long
	// when using slice.
	for ( int pt=0; pt<num_pts; ++pt ) {
		
		for ( int d=0; d<dim; ++d ) {
			main_grad[dim*pt + d] *= weight_sum[pt];
		}
		for ( int d=0; d<dim; ++d ) {
			grad_sum[dim*pt + d]  *= main_weight[pt];
		}
		
		// main_grad[ slice(dim*pt,dim,1) ] *= valarray<double>(weight_sum[pt], dim);
		// grad_sum[ slice(dim*pt,dim,1) ]  *= valarray<double>(main_weight[pt], dim);
	}

	// main_grad is another name for out_grad.
	assert ( out_grad.size() == grad_sum.size() );
	
	// Regression caused by operator += here in place of operator -=.
	// Debugging involved doubting the many variants of weight evaluation
	// functions. More in-built asserts should be placed in to make the
	// code self-checking when NDEBUG is not defined.
	// 2008-09-05 Mike Li.
	out_grad -= grad_sum;
	
	 weight_sum = pow ( weight_sum, 2.0 );
	/*
	for ( int i=0; i<num_pts; ++i ) {
		weight_sum[i] = weight_sum[i]*weight_sum[i];
	}
	*/
	for ( int i=0; i<num_pts; ++i ) {
		
		for ( int d=0; d<dim; ++d ) {
			out_grad[dim*i+d] /= weight_sum[i];
		}
		
		// out_grad[ slice(dim*i, dim, 1) ] /= valarray<double>(weight_sum[i],dim);
	}

	// Internal consistency testing.
#if 0//ndef NDEBUG
	{
		double grad[dim];
		for ( int i=0; i<num_pts; ++i ) {
			local_evaluate_grad( &in[dim*i], grad );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs ( out_grad[dim*i+d] - grad[d] ) < 1e-10 );
			}
		}
	}
#endif
}

template <int dim>
void
partition_of_unity_function<dim>::local_evaluate_and_grad ( valarray<double> &       in,
                                                            pu_function_scratch &    scratch,
                                                            valarray<double> &       out_vals,
                                                            valarray<double> &       out_grad ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;

	valarray<double> &    weight_sum = scratch.weight_sum;
	weight_sum.resize ( num_pts, 0.0 );
	
	valarray<double> &    grad_sum = scratch.grad_sum;
	grad_sum.resize ( in_size, 0.0 );
	
	int num_neigh = neighbour_patches_.size();
	
	valarray<double> &    global = scratch.global;
	global.resize ( in_size );
	patch_->map_local_to_global( in, global );
	
	valarray<bool> &    pred = scratch.pred;
	pred.resize ( num_pts, false );
	
	valarray<double> &    local = scratch.local;
	local.resize ( in_size );
	
	valarray<double> &    weight_tmp = scratch.weight_tmp;
	weight_tmp.resize ( num_pts );
	
	valarray<double> &    grad_tmp = scratch.grad_tmp;
	grad_tmp.resize ( in_size );
	
	valarray<double> &    main_weight = scratch.main_weight;
	main_weight.resize ( num_pts );
	
	// If this is changed, be sure to change the bit below too.
	valarray<double> & main_grad = out_grad;

	
	for ( int i=0; i<num_neigh; ++i ) {
		neighbour_patches_[i]->closed_intersect_point( global, pred );
		
		neighbour_patches_[i]->map_global_to_local ( global, pred, local );

#if 1
		patch_weight_->local_evaluate_and_grad ( *(neighbour_patches_[i]), local, pred, weight_tmp, grad_tmp );
#else
		patch_weight_->local_evaluate( *(neighbour_patches_[i]), local, pred, weight_tmp );
		patch_weight_->local_evaluate_grad( *(neighbour_patches_[i]), local, pred, grad_tmp );
#endif
	
		assert ( weight_sum.size() == weight_tmp.size() );
		weight_sum += weight_tmp;
		
		assert ( grad_sum.size() == grad_tmp.size() );
		grad_sum   += grad_tmp;
	}
	
// 2009-01-22 ML.
#if 1
	patch_weight_->local_evaluate_and_grad ( *patch_, in, main_weight, main_grad );
#else
	patch_weight_->local_evaluate_grad( *patch_, in, main_grad );
	patch_weight_->local_evaluate( *patch_, in, main_weight );
#endif
	
	out_vals.resize ( main_weight.size() );
	out_vals = main_weight;
	out_vals /= weight_sum;
	
	assert ( (-1e-10 < out_vals.min()) && (out_vals.max()<1+1e-10));
	
	// 10% time in function allocating and deallocating unsigned long
	// when using slice.
	for ( int pt=0; pt<num_pts; ++pt ) {
		
		for ( int d=0; d<dim; ++d ) {
			main_grad[dim*pt + d] *= weight_sum[pt];
		}
		for ( int d=0; d<dim; ++d ) {
			grad_sum[dim*pt + d]  *= main_weight[pt];
		}
		
		// main_grad[ slice(dim*pt,dim,1) ] *= valarray<double>(weight_sum[pt], dim);
		// grad_sum[ slice(dim*pt,dim,1) ]  *= valarray<double>(main_weight[pt], dim);
	}

	// main_grad is another name for out_grad.
	assert ( out_grad.size() == grad_sum.size() );
	
	// Regression caused by operator += here in place of operator -=.
	// Debugging involved doubting the many variants of weight evaluation
	// functions. More in-built asserts should be placed in to make the
	// code self-checking when NDEBUG is not defined.
	// 2008-09-05 Mike Li.
	out_grad -= grad_sum;
	
	 weight_sum = pow ( weight_sum, 2.0 );
	/*
	for ( int i=0; i<num_pts; ++i ) {
		weight_sum[i] = weight_sum[i]*weight_sum[i];
	}
	*/
	for ( int i=0; i<num_pts; ++i ) {
		
		for ( int d=0; d<dim; ++d ) {
			out_grad[dim*i+d] /= weight_sum[i];
		}
		
		// out_grad[ slice(dim*i, dim, 1) ] /= valarray<double>(weight_sum[i],dim);
	}

	// Internal consistency testing.
#if 0//ndef NDEBUG
	{
		double grad[dim];
		for ( int i=0; i<num_pts; ++i ) {
			local_evaluate_grad( &in[dim*i], grad );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs ( out_grad[dim*i+d] - grad[d] ) < 1e-10 );
			}
		}
	}
#endif
}

//


template class partition_of_unity_function<2>;
