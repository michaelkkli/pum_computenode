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

#include "function.hh"
#include "box_utils.hh"

#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>

using std::abs;
using std::copy;
using std::cout;
using std::fabs;

template<int dim>
function<dim>::function() : ftype( local_function )
{
}

template<int dim>
function<dim>::~function(){}

template<int dim>
void function<dim>::set_global_function()
{
  ftype = global_function;
}

template<int dim>
void function<dim>::set_local_function()
{
  ftype = local_function;
}

template<int dim>
bool function<dim>::is_global_function() const
{
  return ftype == global_function;
}

template<int dim>
bool
function<dim>::is_local_function() const
{
  return ftype == local_function;
}

// Multipoint evaluation adapter.
template<int dim>
void
function<dim>::evaluate ( valarray<double>&    in,
                          valarray<double>&    out ) const
{
  int in_size = in.size();
  assert( in_size % dim == 0 );
  int num_pts = in_size/dim;
	
  out.resize( num_pts );

  for( int i=0; i<num_pts; ++i ) {
	out[i] = this->evaluate( &in[dim*i] );
  }
}

template<int dim>
void
function<dim>::evaluate ( valarray<double>&    in,
                          valarray<bool>&      pred,
                          valarray<double>&    out ) const
{
  int in_size = in.size();
  assert( in_size % dim == 0 );
  int num_pts = in_size/dim;

  out.resize( num_pts );

  for( int i=0; i<num_pts; ++i ) {
	if ( pred[i] ) {
		out[i] = this->evaluate( &in[dim*i] );
	} else {
		// This is important as we may want to sum the results.
		out[i] = 0.0;
	}
  }
}

// Single point evaluations.

template<int dim>
double
function<dim>::global_evaluate ( const box<dim>&    bx,
                                 const double*      co ) const
{
  assert( co );
	// Do not check for intersection as an uninitialized
	// box may be passed in for a global function defined
	// everywhere.
	// assert( bx.closed_intersect_point( co ) );

  if ( ftype == global_function ) {
    return this->evaluate( co );
  } else {
    double temp[dim];
    bx.map_global_to_local( co, &temp[0] );
    return this->evaluate( &temp[0] );
  }
}

template<int dim>
double function<dim>::local_evaluate ( const box<dim>& bx,
				       const double* co ) const
{
  assert( co );
  if ( ftype == local_function ) {
    return this->evaluate( co );
  } else {
    double temp[dim];
    // Map global points to local points.
    bx.map_local_to_global( co, &temp[0] );
    return this->evaluate( &temp[0] );
  }
}

template<int dim>
double
function<dim>::restricted_local_evaluate ( const box<dim>& rl,
					   const box<dim>& local,
					   const double* co ) const
{
  assert( co );
  if ( local_function == this->ftype ) {
    double temp[dim];
    local.map_restricted_local_to_local( rl, co, &temp[0] );
    return this->evaluate( &temp[0] );
  } else {
    double temp[dim];
    rl.map_local_to_global( co, &temp[0] );
    return this->evaluate( &temp[0] );
  }
}

// Multiple point evaluations.

template<int dim>
void
function<dim>::global_evaluate ( const box<dim>&      bx,
				 valarray<double>&    in,
				 valarray<double>&    out ) const
{
  if ( ftype == global_function ) {
    this->evaluate( in, out );
  } else {
    valarray<double> local;
    bx.map_global_to_local( in, local );
    this->evaluate( local, out );
  }
	// Internal consistency testing.
#if 0 //ndef NDEBUG
	{
		int in_size = in.size();
		assert ( in_size%dim ==0 );
		int num_pts = in_size/dim;
		
		for ( int i=0; i<num_pts; ++i ) {
			assert ( fabs( out[i] - global_evaluate(bx,&in[dim*i]) ) < 1e-10 );
		}
	}
#endif
}

template<int dim>
void
function<dim>::local_evaluate ( const box<dim>&      bx,
				valarray<double>&    in,
				valarray<double>&    out ) const
{
  if ( ftype == local_function ) {
    this->evaluate( in, out );
  } else {
    valarray<double> global;
    bx.map_local_to_global( in, global );
    this->evaluate( global, out );
  }
	
	// Internal consistency testing.
#if 0 //ndef NDEBUG
	{
		int in_size = in.size();
		assert ( in_size%dim ==0 );
		int num_pts = in_size/dim;
		
		for ( int i=0; i<num_pts; ++i ) {
			assert ( fabs( out[i] - local_evaluate(bx,&in[dim*i]) ) < 1e-10 );
		}
	}
#endif
}

template <int dim>
void
function<dim>::local_evaluate ( const box<dim>& bx,
	                      valarray<double>& in,
	                      valarray<bool>&   pred,
	                      valarray<double>& out ) const
{
	if ( ftype == local_function ) {
		this->evaluate ( in, pred, out );
	} else {
		valarray<double> mapped_points;
		bx.map_local_to_global( in, pred, mapped_points );
		this->evaluate( mapped_points, pred, out );
	}
	
	// Internal consistency testing.
#if 0 //ndef NDEBUG
	{
		int in_size = in.size();
		assert ( in_size%dim ==0 );
		int num_pts = in_size/dim;
		
		for ( int i=0; i<num_pts; ++i ) {
			if ( pred[i] ) {
				assert ( fabs( out[i] - local_evaluate(bx,&in[dim*i]) ) < 1e-10 );
			} else {
				assert ( fabs(out[i]) < 1e-10 );
			}
		}
	}
#endif
}

template<int dim>
void
function<dim>::restricted_local_evaluate ( const box<dim>& rl,
			    const box<dim>& local,
			    valarray<double>& points,
			    valarray<double>& results ) const
{
  if ( local_function == this->ftype ) {
    valarray<double> local_points;
    local.map_restricted_local_to_local( rl, points, local_points );
    this->evaluate( local_points, results );
  } else {
    valarray<double> global_points;
    rl.map_local_to_global( points, global_points );
    this->evaluate( global_points, results );
  }
}

// Instantiations needed for the benefit of subclasses.
// The compiler *will* complain if instantiations aren't present.
// On the other hand, the definitions above and the instantiations
// below may only appear once in a single translation unit whereas
// the definitions may be needed again when inheriting ...
// Reconsidering above comment as polynomial seems to do ok with just the header.

template class function<1>;
template class function<2>;
template class function<3>;
