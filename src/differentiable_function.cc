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

#include "differentiable_function.hh"

#include <algorithm>
#include <cmath>
#include <iostream>

using std::abs;
using std::copy;
using std::cout;
using std::fabs;
using std::slice;

template<int dim>
differentiable_function<dim>::differentiable_function()
{
}

template<int dim>
differentiable_function<dim>::~differentiable_function()
{
}

template<int dim>
void
differentiable_function<dim>::evaluate_and_grad ( const double *    co,
                                                  double *          value,
                                                  double *          grad ) const
{
	*value = this->evaluate ( co );
	this->evaluate_grad ( co, grad );
}

template<int dim>
void
differentiable_function<dim>::evaluate_grad ( valarray<double>&    in,
                                              valarray<double>&    out ) const
{
  int in_size = in.size();
  out.resize( in_size );

  assert( in_size%dim == 0 );
  int num_pts = in_size/dim;

  for( int i=0; i<num_pts; ++i ) {
	this->evaluate_grad( &in[dim*i], &out[dim*i] );
  }
}

template <int dim>
void
differentiable_function<dim>::evaluate_grad ( valarray<double>&    in,
                                              valarray<bool>&      pred,
                                              valarray<double>&    out ) const
{
	int in_size = in.size();
	out.resize( in_size );
	
	assert( in_size%dim == 0 );
	int num_pts = in_size/dim;
	
	for( int i=0; i<num_pts; ++i ) {
		if ( pred[i] ) {
			this->evaluate_grad( &in[dim*i], &out[dim*i] );
		} else {
			// Important when summing gradients.
			for ( int d=0; d<dim; ++d ) {
				out[dim*i+d] = 0.0;
			}
			// out[ slice(dim*i,dim,1) ] = zero_grad;
		}
	}
}

/**
	``grad_{x_i} f(x)'' is easy for global functions. For f defined in
	terms of a local function it's trickier.
	grad_{x_i} f(e(x)) = (grad_{e_i} f).(e_i partial x_i)
*/
template<int dim>
void
differentiable_function<dim>::global_evaluate_grad ( const box<dim>& bx,
			                             const double* co,
			                             double* grad ) const
{
  assert( co && grad );
  if ( this->is_global_function() ) {
    // Global function is straight forward.
    this->evaluate_grad( co, grad );
  } else {
    double* tmp_one = grad;
    double tmp_two[dim];
    
    bx.map_global_to_local( co, tmp_one );      // tmp_one has local coordinates now
    this->evaluate_grad( tmp_one, tmp_two );    // tmp_two has grad_{e_i} f
    
    bx.local_partial_global( grad );         // grad now has ( e_i partial x_i )

	  // N.B. tmp_one is simply another name for grad.
    for ( int d=0; d<dim; ++d ) {
	    grad[d] *= tmp_two[d];
    }
  }
}

// Just added so not as well tested as other code here. 2008-12-04 Michael LI.
template<int dim>
void
differentiable_function<dim>::global_evaluate_grad ( const box<dim>&      bx,
                                                     valarray<double>&    in,
                                                     valarray<double>&    out ) const
{
	int in_size = in.size();
	out.resize ( in_size );
	
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	if ( this->is_local_function() ) {
		valarray<double> local(in_size);
		bx.map_global_to_local( in, local );
		this->evaluate_grad( local, out );
		
		double factor[dim];
		bx.local_partial_global( &factor[0] );
		for ( int i=0; i<num_pts; ++i ) {
			for ( int d=0; d<dim; ++d ) {
				out[dim*i+d] *= factor[d];
			}
			// out[ slice(dim*i, dim, 1) ] *= factor;
		}
	} else {
		this->evaluate_grad ( in, out );
	}
	
	// Internal consistency testing.
#if 0 //ndef NDEBUG
	for ( int i=0; i<num_pts; ++i ) {
		double grad[dim];
		
		// GDB found copy and pasted section without adjusting
		// local_evaluate_grad to global_evaluate_grad which
		// tripped an assert in box.cc
		// 2008-12-04 Michael LI.
		global_evaluate_grad ( bx, &in[dim*i], grad );
		for ( int d=0; d<dim; ++d ) {
			assert ( fabs(out[dim*i + d] - grad[d] ) < 1e-10 );
		}
	}
#endif
}

template<int dim>
void
differentiable_function<dim>::local_evaluate_grad ( const box<dim>& bx,
                                                    const double*   co,
                                                    double*         grad ) const
{
  assert( co && grad );
  if ( this->is_local_function() ) {
    double* tmp_one = grad;
    double  tmp_two[dim];
    bx.local_partial_global( tmp_one );
    this->evaluate_grad( co, tmp_two );
	  
	  // N.B. tmp_one is simply another name for grad.
    for ( int d=0; d<dim; ++d ) {
	    grad[d] *= tmp_two[d];
    }
  } else {
    // Global function is straight forward.
    double temp[dim];
    bx.map_local_to_global( co, temp );
    this->evaluate_grad( temp, grad );
  }
}



template<int dim>
void
differentiable_function<dim>::local_evaluate_grad ( const box<dim>&      bx,
                                                    valarray<double>&    in,
                                                    valarray<double>&    out ) const
{
	int in_size = in.size();
	out.resize ( in_size );
	
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	if ( this->is_local_function() ) {
		double factor[dim];
		bx.local_partial_global( &factor[0] );
		this->evaluate_grad ( in, out );
		for ( int i=0; i<num_pts; ++i ) {
			for ( int d=0; d<dim; ++d ) {
				out[dim*i+d] *= factor[d];
			}
			// out[ slice(dim*i, dim, 1) ] *= factor;
		}		
	} else {
		valarray<double> global(in_size);
		bx.map_local_to_global( in, global );
		this->evaluate_grad( global, out );
	}
	
	// Internal consistency testing.
#if 0//ndef NDEBUG
	for ( int i=0; i<num_pts; ++i ) {
		double grad[dim];
		local_evaluate_grad ( bx, &in[dim*i], grad );
		for ( int d=0; d<dim; ++d ) {
			assert ( fabs(out[dim*i + d] - grad[d] ) < 1e-10 );
		}
	}
#endif
}

template <int dim>
void
differentiable_function<dim>::local_evaluate_and_grad ( const box<dim>&   bx,
                                 valarray<double>& in,
                                 valarray<double>& vals,
                                 valarray<double>& grad ) const
{
	int in_size = in.size();
	
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	vals.resize ( num_pts, 0.0 );
	grad.resize ( in_size, 0.0 );
	
	if ( this->is_local_function() ) {
		double factor[dim];
		bx.local_partial_global( &factor[0] );
		
		for ( int i=0; i<num_pts; ++i ) {
			this->evaluate_and_grad ( &in[dim*i], &vals[i], &grad[dim*i] );
			for ( int d=0; d<dim; ++d ) {
				grad[dim*i+d] *= factor[d];
			}
		}
	} else {
		valarray<double> global(in_size);
		bx.map_local_to_global( in, global );
		
		for ( int i=0; i<num_pts; ++i ) {
			this->evaluate_and_grad ( &global[dim*i], &vals[i], &grad[dim*i] );
		}
	}

#if 0//ndef NDEBUG
	{
		for ( int i=0; i<num_pts; ++i ) {
			assert ( fabs( vals[i] - local_evaluate(bx,&in[dim*i]) ) < 1e-10 );
		}
	}
	for ( int i=0; i<num_pts; ++i ) {
		double single_grad[dim];
		local_evaluate_grad ( bx, &in[dim*i], single_grad );
		for ( int d=0; d<dim; ++d ) {
			assert ( fabs(grad[dim*i + d] - single_grad[d] ) < 1e-10 );
		}
	}
#endif
}

template<int dim>
void
differentiable_function<dim>::local_evaluate_and_grad ( const box<dim>&   bx,
                                 valarray<double>& in,
                                 valarray<bool>&   pred,
                                 valarray<double>& vals,
                                 valarray<double>& grad ) const
{
	int in_size = in.size();
	
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	// Debugging found reversal of these arguments causing a failure to allocate
	// memory in other places! 2009-01-22 ML.
	vals.resize ( num_pts, 0.0 );
	grad.resize ( in_size, 0.0 );
	
	if ( this->is_local_function() ) {
		double factor[dim];
		bx.local_partial_global( &factor[0] );
		
		for ( int i=0; i<num_pts; ++i ) {
			if ( pred[i] ) {
				this->evaluate_and_grad ( &in[dim*i], &vals[i], &grad[dim*i] );
				for ( int d=0; d<dim; ++d ) {
					grad[dim*i+d] *= factor[d];
				}
			}
		}
	} else {
		valarray<double> global(in_size);
		bx.map_local_to_global( in, global );
		
		for ( int i=0; i<num_pts; ++i ) {
			if ( pred[i] ) {
				this->evaluate_and_grad ( &global[dim*i], &vals[i], &grad[dim*i] );
			}
		}
	}

#if 0//ndef NDEBUG
	{
		for ( int i=0; i<num_pts; ++i ) {
			if ( pred[i] ) {
				assert ( fabs( vals[i] - local_evaluate(bx,&in[dim*i]) ) < 1e-10 );
			}
		}
	}
	for ( int i=0; i<num_pts; ++i ) {
		double single_grad[dim];
		if ( pred[i] ) {
			local_evaluate_grad ( bx, &in[dim*i], single_grad );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs(grad[dim*i + d] - single_grad[d] ) < 1e-10 );
			}
		}
	}
#endif
}

template<int dim>
void
differentiable_function<dim>::local_evaluate_grad ( const box<dim>&      bx,
                                                    valarray<double>&    in,
                                                    valarray<bool>&      pred,
                                                    valarray<double>&    out ) const
{
	int in_size = in.size();
	out.resize ( in_size );
	
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	if ( this->is_local_function() ) {
		double factor[dim];
		bx.local_partial_global( &factor[0] );
		this->evaluate_grad ( in, pred, out );
		for ( int i=0; i<num_pts; ++i ) {
			if ( pred[i] ) {
				for ( int d=0; d<dim; ++d ) {
					out[dim*i+d] *= factor[d];
				}
				// out[ slice(dim*i, dim, 1) ] *= factor;
			}
		}
	} else {
		valarray<double> global(in_size);
		bx.map_local_to_global( in, pred, global );
		this->evaluate_grad( global, pred, out );
	}
	
	// Internal consistency testing.
#if 0//ndef NDEBUG
	double grad[dim];
	for ( int i=0; i<num_pts; ++i ) {
		if ( pred[i] ) {
			local_evaluate_grad ( bx, &in[dim*i], grad );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs(out[dim*i + d] - grad[d] ) < 1e-10 );
			}
		} else {
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs(out[dim*i + d] ) < 1e-10 );
			}
		}
	}
#endif
}

/**
	grad_{x_i} f(e(x(r))) = \sum_i ( f partial e_i ).( e_i partial x_i )
	so (local_eval_grad) product (local.local_partial_global)
*/
template<int dim>
void
differentiable_function<dim>::restricted_local_evaluate_grad ( const box<dim>& restricted_local,
                                 const box<dim>& local,
                                 const double*   co,
                                 double*         re ) const
{
  assert( co && re );
  if ( this->is_local_function() ) {
    double tmp_one[dim];
    double tmp_two[dim];
    local.map_restricted_local_to_local( restricted_local, co, &tmp_two[0] );
    this->evaluate_grad( &tmp_two[0], &tmp_one[0] );
    local.local_partial_global( &tmp_two[0] );
    restricted_local.global_partial_local( &tmp_two[0] );
    for ( int d=0; d<dim; ++d ) {
	    re[d] = tmp_one[d]*tmp_two[d];
    }
  } else {
    // Global function is straight forward.
    double temp[dim];
    restricted_local.map_local_to_global( co, &temp[0] );
    this->evaluate_grad( &temp[0], re );
  }
}

//

template class differentiable_function<1>;
template class differentiable_function<2>;
template class differentiable_function<3>;
