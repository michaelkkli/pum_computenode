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

#include "box.hh"
#include "box_utils.hh"

#include <algorithm>
#include <valarray>
#include <vector>

#include <cassert>
#include <cstdlib>

#include <iostream>

using std::copy;
using std::cout;
using std::abs;
using std::fabs;
using std::copy;
using std::min;
using std::max;
using std::pow;
using std::valarray;
using std::vector;

template<int dim>
void
box<dim>::expand_by_factor ( const box<dim>& first,
			     double factor,
			     box<dim>& other )
{
  double a, b, c;
  const double* ext       = first.get();
  double* other_ext = other.get();
  assert( ext );
  assert( other_ext );
  for ( int d=0; d<dim; ++d ) {
    a = ext[d*dim    ];
    b = ext[d*dim + 1];
    c = 0.5 * (b-a);
    other_ext[d*dim    ] = factor * (a-c) + c;
    other_ext[d*dim + 1] = factor * (b-c) + c;
  }
}

template<int dim>
void
box<dim>::translate_centre ( const box<dim>& first,
			     const double* c_vec,
			     box<dim>& other )
{
  assert( c_vec );
  const double* ext       = first.get();
  double* other_ext = other.get();
  assert( ext );
  assert( other_ext );
  for ( int d=0; d<dim; ++d ) {
    other_ext[ d*dim     ] = ext[ d*dim     ] + c_vec[d];
    other_ext[ d*dim + 1 ] = ext[ d*dim + 1 ] + c_vec[d];
  }
}

template<int dim>
box<dim>::box(){}

template<int dim>
box<dim>::~box(){}

template<int dim>
void
box<dim>::set( const double* ext )
{
	assert( ext );
	copy( ext, ext+2*dim, extents );
#if 0
	for ( int i=0; i<2*dim; ++i ) {
		extents[i] = ext[i];
	}
#endif
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::translate( const double* disp )
{
	assert( disp );
	for ( int d=0; d<dim; ++d ) {
		extents[ 2*d     ] += disp[d];
		extents[ 2*d + 1 ] += disp[d];
	}
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::scale( const double* s )
{
	assert( s );
	double centre;
	for ( int d=0; d<dim; ++d ) {
		centre = 0.5*(extents[ 2*d ] + extents[ 2*d + 1 ]);
		extents[ 2*d     ] = (extents[ 2*d     ]-centre) * s[d] + centre;
		extents[ 2*d + 1 ] = (extents[ 2*d + 1 ]-centre) * s[d] + centre;
	}
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::scale( double s )
{
	double centre;
	for ( int d=0; d<dim; ++d ) {
		centre = 0.5*(extents[ 2*d ] + extents[ 2*d + 1 ]);
		extents[ 2*d     ] = (extents[ 2*d     ]-centre) * s + centre;
		extents[ 2*d + 1 ] = (extents[ 2*d + 1 ]-centre) * s + centre;
	}
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::recentre( const double* new_centre )
{
	assert( new_centre );
	double old_centre, add_disp;
	for ( int d=0; d<dim; ++d ) {
		old_centre = 0.5*(extents[ 2*d ] + extents[ 2*d + 1 ]);
		add_disp = new_centre[d] - old_centre;
		extents[ 2*d     ] += add_disp;
		extents[ 2*d + 1 ] += add_disp;
	}
}


// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
bool
box<dim>::empty() const
{
  for ( int i=0; i<dim; ++i ) {
    if ( extents[2*i] >= extents[2*i+1] ) {
      return true;
    }
  }
  return false;
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::get_centre_point ( double* co ) const
{
	assert( co );
	
	// Minus found instead of plus here. Messed up partitioning
	// of patches into interior and exterior.
	for ( int d=0; d<dim; ++d ) {
		co[d] = 0.5*(extents[ 2*d + 1 ] + extents[ 2*d ]);
	}
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
bool
box<dim>::closed_intersect_point ( const double* co ) const
{
	// Do not check for emptiness of box - far too expensive.
	
  for ( int d=0; d<dim; ++d ) {
    if ( co[d] < extents[2*d] || extents[2*d+1] < co[d] ) {
      return false;
    }
  }
  return true;
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template <int dim>
void
box<dim>::closed_intersect_point ( valarray<double>& in, valarray<bool>& out ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	out.resize( num_pts );
	
	bool continue_outer = false;
	
	for ( int i=0; i<num_pts; ++i ) {
		
		double* co = &in[dim*i];
		
		for ( int d=0; d<dim; ++d ) {
			if ( co[d] < extents[2*d] || extents[2*d+1] < co[d] ) {
				out[i] = false;
				
				assert ( out[i] == closed_intersect_point( co ) );
				
				continue_outer = true;
				break;
			}
		}
		if ( continue_outer ) {
			continue_outer = false;
			continue;
		}
		out[i] = true;
		assert ( out[i] == closed_intersect_point( co ) );
	}
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
bool
box<dim>::open_intersect_point ( const double* co ) const
{
  if ( this->empty() ) {
    return false;
  }
  for ( int d=0; d<dim; ++d ) {
    if ( co[d] <= extents[2*d] || extents[2*d+1] <= co[d] ) {
      return false;
    }
  }
  return true;
}

/**
	The strategy is to parametrize the line with parameter `t'. We project
	the box and line to each axis in turn and see which range of parameters
	give intersection of the projections.

	Special cases if line is very nearly vertical or very nearly horizontal.
*/
template<int dim>
void
box<dim>::intersect_line ( const double* line,
                           box_intersect_line_scratch & scratch,
                           bool* intersects_closed,
                           bool* intersects_open ) const
{
	assert( line && !empty() );
	
	const double* ext = this->get();
	
	// Easy exit loop to avoid creating array t.
	for ( int d=0; d<dim; ++d ) {
		double a = ext[ d*dim     ]; // Lower endpoint of interval d of box.
		double b = ext[ d*dim + 1 ]; // Upper endpoint of interval d of box.
		
		double r = line[         d ]; // Coordinate d of first point of line.
		double s = line[ dim*1 + d ]; // Coordinate d of second point of line.
		
		if ( (r<a && s<a) || (b<r && b<s) ) {
			// Intersection not possible.
			
			/**
				Very strange bug to find. The dereferencing operators
				were not present and false was being assigned to a pointer.
			*/
			if ( intersects_closed ) {
				*intersects_closed = false;
			}
			if ( intersects_open ) {
				*intersects_open = false;
			}
			return;
		}
	}
	
	// Parametrization variable `t'.
	vector<box<1> > &    t = scratch.box_vector;
	t.resize ( dim );
	
	for ( int d=0; d<dim; ++d ) {
		double a = ext[ d*dim     ]; // Lower endpoint of interval d of box.
		double b = ext[ d*dim + 1 ]; // Upper endpoint of interval d of box.
		
		double r = line[         d ]; // Coordinate d of first point of line.
		double s = line[ dim*1 + d ]; // Coordinate d of second point of line.
		
		if ( a<=r && a<=s && r<=b && s<= b ) {
			// Any parameter will do for this dimension.
			double other_ext[2] = { -1.0, 1.0 };
			t[d].set( &other_ext[0] );
			continue;
		}
		
		// We want parameter t such that
		// a <= r + t.(s-r) <= b
		// so t=(a-r)/(s-r) at lower end
		// and t=(b-r)/(s-r) at upper end.
		
		double smr = s-r; // s minus r.
		assert( fabs( static_cast<long double>(smr) ) > 1e-16 );
		
		double new_ext[2] = { (a-r)/smr, (b-r)/smr };
		if ( new_ext[0]>new_ext[1] ) {
			double tmp = new_ext[0];
			new_ext[0] = new_ext[1];
			new_ext[1] = tmp;
		}
		assert( new_ext[0] <= new_ext[1] );
		
		t[d].set( &new_ext[0] );
	}
	
	box<1> &    all_intersected = scratch.all_intersected;
	all_intersected = t[0];
	
	
	for ( int i=1; i<dim; ++i ) {
		all_intersected.clip_against( t[i] );
	}
	
	// If there is a range of parameter values, we retrieve
	// the values wanted and return from member function.
	if ( !all_intersected.empty() ) {
		if ( intersects_open ) {
			*intersects_open   = true;
		}
		if ( intersects_closed ) {
			*intersects_closed = true;
		}
		return;
	} else {
		if ( intersects_open ) {
			*intersects_open = false;
		}
	}
	
	// Continue if there isn't a range of parameter values.
	
	if ( !intersects_closed ) {
		// Finish if we don't care about a intersections of
		// the closed box and the given line.
		return;
	}
	
	// Now for the case that there isn't an interval where the
	// parameter t remains inside the box. There might be a single
	// value of t that is suitable.
	
	double candidate_t[dim];
	const box<1>& const_box_ref = all_intersected;
	const double* cand_ext = const_box_ref.get();
	for ( int d=0; d<dim; ++d ) {
		candidate_t[d] = cand_ext[2*d]; // Pick lower end point of parameter interval.
	}
	
	for ( int d=0; d<dim; ++d ) {
		if ( !t[d].closed_intersect_point( &candidate_t[0] ) ) {
			*intersects_closed = false;
			return;
		}
	}
	*intersects_closed = true;
}

template<int dim>
bool
box<dim>::closed_intersect_line ( const double* line, box_intersect_line_scratch & scratch ) const
{
	/**
		Very tricky bug to find. Valgrind memcheck tool
		found an if-statement depended on something uninitialized
		in boundary<dim>::closed_intersect_box.
		intersects_closed was uninitialized.
		Seems this was also responsible for the incorrect
		behaviour compiling with all debug off, profiling off,
		NDEBUG defined and -O3 on hilbert/archangel.
		2008-09-01 Mike Li.
	
		This was actually a symptom. The real problem was with intersect_line.
		Removing initialization here (initialization of intersects_closed suppressed
		the warning without solving the root cause.
		2008-10-03 Mike Li.
	*/
	bool intersects_closed;
	
	this->intersect_line( line, scratch, &intersects_closed, 0 );
	return intersects_closed;
}


/**
 * We drop a perpendicular from the centre of the box onto the line
 * to find the point of intersection with the line. We then test to see
 * if the line is inside the box.
 * 
 * If a and b are the end points, c is the centre and p is the foot of
 * the perpendicular.
 * 
 * 	(p-c).(b-a) = 0
 * 
 * 	p = a + k ( b - a )
 * 
 * 	k = [ (c-a).(b-a) ]/[(b-a).(b-a)]
 * 
 * 	
 */
template<int dim>
bool box<dim>::open_intersect_line ( const double* line, box_intersect_line_scratch & scratch ) const
{
#if 0 // Remove regression. Attempt at optimization introduces regression.
	assert ( line );
	
	valarray<double> a(dim), b(dim), c(dim), p(dim), tmp(dim);
	
	copy ( line,     line+dim, &a[0] );
	copy ( line+dim, line+2*dim, &b[0] );
	
	get_centre_point( &c[0] );
	
	p = c;
	p -= a;
	tmp = b;
	tmp -= a;
	p *= tmp;
	
	double k = p.sum();
	
	p = tmp;
	p *= tmp;
	
	k /= p.sum();
	
	if ( k < 0.0 || k > 1.0 ) {
		return false;
	}
	
	p = a;
	tmp *= k;
	p += tmp;
	
	cout << "Point is " << p[0] << ", " << p[1] << " ";
	if ( open_intersect_point( &p[0] ) ) {
		cout << "and intersects.\n";
	} else {
		cout << "and does not intersect.\n";
	}
	return open_intersect_point( &p[0] );
#endif

	/**
		Region intersection code for the petsc_solver-3 simulation was different
		for the cover set-up and the visualization cover set-up. The first used
		closed_intersect_line in which the uninitialized variable was already found.
		The second intersection test ends up using open_intersect_line and intersects_open
		was uninitialized.
		2008-10-03 Mike Li.
		This was actually a symptom. The real problem was with intersect_line.
		Removing initialization here.
		2008-10-03 Mike Li.
	*/
	bool intersects_open;
	this->intersect_line( line, scratch, 0, &intersects_open );
	return intersects_open;
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::global_partial_local ( double* re ) const
{
	assert( re );
	for ( int d=0; d<dim; ++d ) {
		re[d] = 0.5*(extents[ 2*d + 1 ] - extents[ 2*d ]);
	}
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::local_partial_global ( double* re ) const
{
	assert( re );
	for ( int d=0; d<dim; ++d ) {
		re[d] = 2.0/(extents[ 2*d + 1 ] - extents[ 2*d ]);
	}
}


template <int dim>
double
box<dim>::local_line_length ( double *    line ) const
{
	double total = 0.0;
	for ( int d=0; d<dim; ++d ) {
		double tmp = 2.0*(line[dim + d]-line[d])/(extents[dim*d+1]-extents[dim*d]);
		total += tmp*tmp;
	}
	return sqrt ( total );
}


// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::scale_and_translate ( double a0, double a1, double a2,
				double d0, double d1, double d2 )
{
  if ( this->empty() ) {
    return;
  }
  double scale_a[3] = { a0, a1, a2 };
  double displace[3] = {d0, d1, d2 };
  double tmp[2];
  double centre;
  for ( int d=0; d<dim; ++d ) {
      centre = 0.5*(extents[2*d+1] - extents[2*d]);
      tmp[0] = scale_a[d]*(extents[2*d] - centre) + centre + displace[d];
      tmp[1] = scale_a[d]*(extents[2*d+1] - centre) + centre + displace[d];
      extents[2*d] = tmp[0];
      extents[2*d+1] = tmp[1];
    }
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::map_local_to_global ( const double* co, double* out ) const
{
  assert( co && out );
  for ( int i=0; i<dim; ++i ) {
  	  if ( -1.-1e-10 > co[i] || co[i] > 1.+1e-10 ) {
	  	cout << "Problem!! i is " << i  << ", co[i] is " << co[i] << "\n";
	  	  abort();
	  }
	  assert( -1.-1e-10 < co[i] && co[i] < 1.+1e-10 );
	  out[i] = 0.5*(extents[2*i] + extents[2*i+1] + co[i]*(extents[2*i+1] - extents[2*i]) );
  }
	// Redundant and expensive paranoia.
	// assert ( closed_intersect_point( out ) );
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::map_global_to_local ( const double* co, double* out ) const
{
	assert ( co && out );
	
	// Point must be in box.
	assert( closed_intersect_point( co ) );
  for ( int i=0; i<dim; ++i ) {
	// Erroneous assert found by debugger. assert( -1. <= co[i] && co[i] <= 1. );
	  
	  // Debugging found this erroneous assert which prevents the output
	  // array from being the same as the input array. We wish to allow co == out
	  // for efficiency. 2008-08-19 Mike Li.
	  // assert( closed_intersect_point( co ) );
	
	  
	out[i] = ( 2.0*co[i] - extents[2*i] - extents[2*i+1] )/( extents[2*i+1] - extents[2*i] );
	assert( -1. <= out[i] && out[i] <= 1 );
  }
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void
box<dim>::map_restricted_local_to_local ( const box<dim>& rl,
					  const double* co,
					  double* out ) const
{
  //e(r) = ( (c+d-a-b) + r*(d-c) )/(b-a)
  assert( co && out );
  assert( closed_contains_box( rl ) );
  for ( int i=0; i<dim; ++i ) {
    assert ( -1.0 <= co[i] && co[i] <= 1.0 );

    double a=this->extents[2*i], b=this->extents[2*i+1];
    double c=rl.extents[2*i],    d=rl.extents[2*i+1];
    
    out[i] = ( (c+d-a-b) + co[i]*(d-c) )/(b-a);
    
#if 0 // Debug small number.
    out[i] = ( rl.extents[2*i]
	       + rl.extents[2*i+1]
	       - this->extents[2*i]
	       - this->extents[2*i+1]
	       + co[i]*( rl.extents[2*i+1] - rl.extents[2*i] )
	       )
      / ( this->extents[2*i+1]-this->extents[2*i] );
#endif
    
    assert( -1.0 <= out[i] && out[i] <= 1.0 );
  }
}

template<int dim>
void
box<dim>::map_local_to_global ( valarray<double>& points,
			   valarray<double>& result ) const
{
	
  int size = points.size();

  assert( size % dim == 0 );

  int num_pts = size / dim;

  // GDB helped find that the resize was incorrectly set to num_pts.
  // Corrected the related functions too.
  result.resize( size );

  for ( int i=0; i<num_pts; ++i ) {
	  // GDB helped find that the second argument was incorrectly
	  // set to &result[i]. Corrected the related functions too.
	  this->map_local_to_global( &points[dim*i], &result[dim*i] );
  }
}

template <int dim>
void
box<dim>::map_local_to_global ( valarray<double>& in, valarray<bool>& pred, valarray<double>& out ) const
{
  int in_size = in.size();

  assert( in_size % dim == 0 );

  int num_pts = in_size / dim;

  out.resize( in_size );

  for ( int i=0; i<num_pts; ++i ) {
	  if ( pred[i] ) {
		this->map_local_to_global( &in[dim*i], &out[dim*i] );
	  }
  }
}

template<int dim>
void
box<dim>::map_global_to_local ( valarray<double>& points,
				valarray<double>&       result ) const
{
  int size = points.size();
  assert( size % dim == 0 );

  int num_pts = size / dim;

  result.resize( size );
	
#if 0//ndef NDEBUG
	for ( int i=0; i<num_pts; ++i ) {
		assert ( this->closed_intersect_point( &points[dim*i] ) );
	}	
#endif

  for ( int i=0; i<num_pts; ++i ) {
	this->map_global_to_local( &points[dim*i], &result[dim*i] );
  }
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
// Callgrind shown heavily used.
template <int dim>
void
box<dim>::map_global_to_local ( valarray<double>& in,
                                valarray<bool>& pred,
                                valarray<double>& out ) const
{
  int in_size = in.size();

  assert( in_size % dim == 0 );

  int num_pts = in_size / dim;
	
	assert ( pred.size() == num_pts );

  out.resize( in_size );
  for ( int i=0; i<num_pts; ++i ) {
	  if ( pred[i] ) {
		double * global = &in[dim*i];
		double * local  = &out[dim*i];
		
		assert( closed_intersect_point( global ) );
		for ( int d=0; d<dim; ++d ) {
			const double & aa = extents[2*d];
			const double & bb = extents[2*d+1];
			
			local[d] = ( 2.0*global[d] - aa - bb )/( bb - aa );
			assert( -1. <= local[d] && local[d] <= 1 );
		}
		  // this->map_global_to_local( &in[dim*i], &out[dim*i] );
	  }
  }
	
#if 0//ndef NDEBUG
	for ( int i=0; i<num_pts; ++i ) {
		valarray<double> compare_in(dim), compare_out(dim);
		if ( pred[i] ) {
			copy ( &in[dim*i], &in[dim*i]+dim, &compare_in[0] );
			this->map_global_to_local( &compare_in[0], &compare_out[0] );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs( compare_out[d] - out[dim*i+d] ) < 1e-10 );
			}
		}
	}
#endif

	
}

template<int dim>
void
box<dim>::map_restricted_local_to_local ( const box<dim>& rl,
				valarray<double>& points,
				valarray<double>&       result ) const
{

#if 0 //ndef NDEBUG // Keep for future debug.
	if ( !closed_contains_box( rl ) ) {
		cout << "Local box\n";
		output_description( *this, cout );
		cout << "does not contain restricted local box\n";
		output_description( rl, cout );
	}
#endif
	
  assert( closed_contains_box( rl ) );
	
  int size = points.size();

  assert( size % dim == 0 );

  int num_pts = size / dim;

  result.resize( size );


  for ( int i=0; i<num_pts; ++i ) {
	  this->map_restricted_local_to_local( rl, &points[dim*i], &result[dim*i] );
  }
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
bool
box<dim>::open_intersect_box( const box<dim>& other ) const
{
  for ( int i=0; i<dim; ++i ) {
    if ( this->extents[2*i] > other.extents[2*i+1]
	 || this->extents[2*i+1] < other.extents[2*i] )
      {
	return false;
      }
  }
  return true;
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
bool
box<dim>::closed_intersect_box( const box<dim>& other ) const
{
  for ( int i=0; i<dim; ++i ) {
    if ( this->extents[2*i] >= other.extents[2*i+1]
	 || this->extents[2*i+1] <= other.extents[2*i] )
      {
	return false;
      }
  }
  return true;
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
bool
box<dim>::closed_contains_box( const box<dim>& other ) const
{
	assert( !empty() );
	assert( !other.empty() );
	
	// Debugging found variable clash with 'd' used twice.
	
	// Want p <= r and s <= q
	double p, q, r, s;
	for ( int d=0; d<dim; ++d ) {
		p = extents[2*d];
		q = extents[2*d+1];
		r = other.extents[2*d];
		s = other.extents[2*d+1];
		if ( p > r || s > q ) {
			
			return false;
		}
	}
	return true;
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
bool
box<dim>::open_contains_box( const box<dim>& other ) const
{
	assert( !empty() );
	assert( !other.empty() );

	// Want p < r and s < q
	double p, q, r, s;
	for ( int d=0; d<dim; ++d ) {
		p = extents[2*d];
		q = extents[2*d+1];
		r = other.extents[2*d];
		s = other.extents[2*d+1];
		if ( p >= r || s >= q ) {
			
			return false;
		}
	}
	return true;
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
void box<dim>::clip_against( const box<dim>& other )
{
	for ( int i=0; i<dim; ++i ) {
		extents[2*i]   = max( extents[2*i],   other.extents[2*i]   );
		extents[2*i+1] = min( extents[2*i+1], other.extents[2*i+1] );
	}
}

// TODO: Check changing 2 to dim.
// Not doing at the moment in order to avoid introducing regressions.
// 2009-04-16 ML.
template<int dim>
double
box<dim>::measure() const
{
	assert( dim >= 1 );
	assert( !this->empty() );

	double tmp = extents[1]-extents[0];

	for ( int d=1; d<dim; ++d ) {
		tmp *= extents[2*d+1] - extents[2*d];
	}
	return tmp;
}


// Explicit template initialization.

template class box<1>;
template class box<2>;
template class box<3>;
