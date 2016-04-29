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

#include <algorithm>
#include <fstream>
#include <list>
#include <set>
#include <stack>
#include <vector>
#include <cassert>
#include <cstdlib>


#include <iostream>
using std::cout;


using std::copy;
using std::div;
using std::list;
using std::set;
using std::vector;

template<int dim>
void
local_to_restricted_local ( valarray<double>& points,
			    const box<dim>& restricted_local,
			    const box<dim>& local_box,
			    valarray<double>& out )
{
  assert( points.size()%dim == 0 );
  assert( local_box.open_intersect_box( restricted_local ) );

  out.resize( points.size() );

  restricted_local.map_local_to_global( points, out );
  local_box.map_global_to_local( out, out );
}

template<int dim>
void
output_description ( const box<dim>& bx, ostream& out )
{
  const double* vals = bx.get();
  out << "box (";
  for ( int d=0; d<dim; ++d ) {
    out << vals[2*d] << ", " << vals[2*d+1] << ", ";
  }
  out << "\b\b)\n";
}

template<int dim>
void make_bounding_box ( const vector<box<dim> >& boxes,
			 box<dim>& bbox )
{
  int vec_size = boxes.size();
  if ( vec_size == 0 ) {
    return;
  }

  double extents[2*dim];

  const double* bget = 0;

  bget = boxes[0].get();

  // Set up first guess.
  for ( int i=0; i<2*dim; ++i ) {
    extents[i] = bget[i];
  }

  for ( int i=1; i<vec_size; ++i ) {
    bget = boxes[i].get();
    for ( int d=0; d<dim; ++d ) {
      extents[ 2*d ]   = std::min( extents[ 2*d ],   bget[ 2*d ] );
      extents[ 2*d+1 ] = std::max( extents[ 2*d+1 ], bget[ 2*d+1 ] );
    }
  }
  
  bbox.set( &extents[0] );
#if 0 // The above line should be equivalent.
  if ( dim == 1 ) {
    bbox.set( extents[0], extents[1] );
  } else if ( dim == 2 ) {
    bbox.set( extents[0], extents[1],
	      extents[2], extents[3] );
  } else if ( dim == 3 ) {
    bbox.set( extents[0], extents[1],
	      extents[2], extents[3],
	      extents[4], extents[5] );
  }
#endif
}

#if 0
template <int dim, typename container_t>
void split_box ( const box<dim>& bx, container_t& container )
{
  if ( bx.empty() ) {
    return;
  }
  double middle[dim];
  const double* ext = bx.get();
  for ( int d=0; d<dim; ++d ) {
      middle[d] = 0.5*(ext[2*d+1]-ext[2*d]);
    }
	  if ( 1 == dim ) {
	    container.push( box<dim>( ext[0],    middle[0] ) );
	    container.push( box<dim>( middle[0], ext[1] ) );
	  } else if ( 2 == dim ) {
	    // SW
	    container.push( box<dim>( ext[0], middle[0],
					   ext[2], middle[1]) );
	    // SE
	    container.push( box<dim>( middle[0], ext[1],
					   ext[2], middle[1]) );
	    // NW
	    container.push( box<dim>( ext[0], middle[0],
					   middle[1], ext[2+1]) );
	    // NE
	    container.push( box<dim>( middle[0], ext[1],
					   middle[1], ext[2+1]) );
	  } else if ( 3 == dim ) {
	    // Bottom (lower z)
	    // SW
	    container.push( box<dim>( ext[0], middle[0],
				      ext[2], middle[1],
				      ext[4], middle[2] ) );
	    // SE
	    container.push( box<dim>( middle[0], ext[1],
				      ext[2], middle[1],
				      ext[4], middle[2] ) );
	    // NW
	    container.push( box<dim>( ext[0], middle[0],
				      middle[1], ext[2+1],
				      ext[4], middle[2] ) );
	    // NE
	    container.push( box<dim>( middle[0], ext[1],
				      middle[1], ext[2+1],
				      ext[4], middle[2] ) );
	    // Top (lower z)
	    // SW
	    container.push( box<dim>( ext[0], middle[0],
				      ext[2], middle[1],
				      middle[2], ext[4+1] ) );
	    // SE
	    container.push( box<dim>( middle[0], ext[1],
				      ext[2], middle[1],
				      middle[2], ext[4+1] ) );
	    // NW
	    container.push( box<dim>( ext[0], middle[0],
				      middle[1], ext[2+1],
				      middle[2], ext[4+1]) );
	    // NE
	    container.push( box<dim>( middle[0], ext[1],
					   middle[1], ext[2+1],
				      middle[2], ext[4+1] ) );
	  }
}
#endif

template <int dim, typename coords>
int count_closed_intersect ( const box<dim>& bx,
			     coords& cds )
{
  const int size = cds.size();
  if ( size == 0 ) {
    return 0;
  }
  if ( size%dim != 0 ) {
    return 0;
  }
  int num_pts = size/dim;
  int count = 0;

  for ( int i=0; i<num_pts; ++i ) {
    if ( bx.closed_intersect_point( &cds[dim*i] ) ) {
      ++count;
    }
  }
  return count;
}

template <int dim, typename coords>
int count_open_intersect ( const box<dim>& bx,
			   coords& cds )
{
  const int size = cds.size();
  if ( size == 0 ) {
    return 0;
  }
  if ( size%dim != 0 ) {
    return 0;
  }
  int num_pts = size/dim;
  int count = 0;
  for ( int i=0; i<num_pts; ++i ) {
    if ( bx.open_intersect_point( &cds[dim*i] ) ) {
      ++count;
    }
  }
  return count;
}

#if 0 // Keep code but unused.
template <int dim>
void
closed_intersect_point ( const vector<box<dim> > & in, const double* co, valarray<bool>& out )
{
	int in_size = in.size();
	
	out.resize ( in_size );
	
#if 0
	bool continue_outer = false;
	
	// It turns out that this is an expensive loop to use and worse
	// better than to call closed_intersect_point on each box individually.
	// 2008-09-10 Mike Li.
	for ( int i=0; i<in_size; ++i ) {
		
		const double* ext = in[i].get();
		
		for ( int d=0; d<dim; ++d ) {
			if ( co[d] < ext[2*d] || ext[2*d+1] < co[d] ) {
				out[i] = false;
				
				assert ( out[i] == in[i].closed_intersect_point( co ) );
				
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
#endif
	

	for ( int i=0; i<in_size; ++i ) {
		out [i] = in[i].closed_intersect_point ( co );
	}

}
#endif

template<int dim>
void
decompose_box( const box<dim>& restricted,
	       const vector<box<dim> >& vec_in,
               vector<box<dim> >& vec_out )
{
	assert( !restricted.empty() );

	const double* ext = restricted.get();

	const double* tmp_ext = 0;
	
	set<double> points[dim];
	
	double lower, upper, tmp_d;
	for ( int d=0; d<dim; ++d ) {
		lower = ext[2*d];
		upper = ext[2*d+1];

		// TODO: remove
		// cout << "lower " << lower << ", upper " << upper << "\n";
		
		set<double>& pts = points[d];
		
		pts.insert( lower );
		pts.insert( upper );

		typename vector<box<dim> >::const_iterator it( vec_in.begin() ), end( vec_in.end() );
		for ( ; it!=end; ++it ) {
			assert( !it->empty() );
			
			tmp_ext = it->get();
			tmp_d = tmp_ext[2*d];
			if ( lower < tmp_d && tmp_d < upper ) {
				pts.insert( tmp_d );
				
				// TODO: remove
				//cout << "inserting " << tmp_d << "\n";
				
			}
			tmp_d = tmp_ext[2*d+1];
			if ( lower < tmp_d && tmp_d < upper ) {
				pts.insert( tmp_d );
				
				// TODO: remove
				//cout << "inserting " << tmp_d << "\n";
				
			}
		}
	}
	
	vector<vector<double> > vec_points( dim );
	
	for ( int d=0; d<dim; ++d ) {
		vec_points[d].resize( points[d].size() );
		copy( points[d].begin(), points[d].end(), vec_points[d].begin() );
	}
	
	make_boxes( vec_points, vec_out );
}

// The output is vector<box<dim> > rather than a vector of pointers
// as we're not likely to reuse the output after integration.
template <int dim>
void
decompose_box( const box<dim>&                   restricted,
               const typename box<dim>::vec_ptr& vec_in,
               vector<box<dim> >&                vec_out )
{
	assert( !restricted.empty() );

	const double* ext = restricted.get();

	const double* tmp_ext = 0;
	
	set<double> points[dim];
	
	double lower, upper, tmp_d;
	for ( int d=0; d<dim; ++d ) {
		lower = ext[2*d];
		upper = ext[2*d+1];

		// TODO: remove
		// cout << "lower " << lower << ", upper " << upper << "\n";
		
		set<double>& pts = points[d];
		
		pts.insert( lower );
		pts.insert( upper );

		typename box<dim>::vec_ptr::const_iterator it( vec_in.begin() ), end( vec_in.end() );
		for ( ; it!=end; ++it ) {
			assert( !(*it)->empty() );

			// Ensure we call the const version of get by calling.
			// through a const pointer held in a shared_ptr.
			shared_ptr<const box<dim> > tmp_sp_const = *it;
			tmp_ext = tmp_sp_const->get();

			tmp_d = tmp_ext[2*d];
			if ( lower < tmp_d && tmp_d < upper ) {
				pts.insert( tmp_d );
				
				// TODO: remove
				//cout << "inserting " << tmp_d << "\n";
				
			}
			tmp_d = tmp_ext[2*d+1];
			if ( lower < tmp_d && tmp_d < upper ) {
				pts.insert( tmp_d );
				
				// TODO: remove
				//cout << "inserting " << tmp_d << "\n";
				
			}
		}
	}
	
	vector<vector<double> > vec_points( dim );
	
	for ( int d=0; d<dim; ++d ) {
		vec_points[d].resize( points[d].size() );
		copy( points[d].begin(), points[d].end(), vec_points[d].begin() );
	}
	
	make_boxes( vec_points, vec_out );
}

#if 0
template <typename container_t, typename container2>
void
make_boxes ( int dim, const container_t* pts, container2& boxes )
{
	assert( 1 <= dim && dim <=3 );
	if ( dim == 1 ) {
		int num_boxes = pts[0].size() - 1;
		assert( num_boxes > 0 );
		boxes.resize( num_boxes );
		typename container2::iterator box_it( boxes.begin() );
		typename container_t::iterator it( pts[0].begin() ), end( pts[0].end() );
		double lower, upper = *it;
		++it;
		for ( ; it != end; ++it ) {
			lower = upper;
			upper = *it;
			double ext[2] = { lower, upper };
			box_it->set( &ext[0]  );
			++box_it;
		}
		return;
	} else if ( dim == 2 ) {
		int num_boxes = (pts[0].size()-1)*(pts[1].size()-1);
		assert( num_boxes > 0 );
		boxes.resize( num_boxes );
		typename container2::iterator box_it( boxes.begin() );
		typename container_t::iterator it( pts[0].begin() ), end( pts[0].end() );
		typename container_t::iterator outer_it( pts[1].begin() ), outer_end( pts[1].end() );
		double lower, upper = *it;
		++it;
		double outer_lower, outer_upper = *outer_it;
		++outer_it;
		for ( ; outer_it != outer_end; ++outer_it ) {
			outer_lower = outer_upper;
			outer_upper = *outer_it;
			for ( ; it != end; ++it ) {
				lower = upper;
				upper = *it;
				double ext[4] = { lower, upper, outer_lower, outer_upper };
				
				assert( box_it != boxes.end() ); // Don't run off the end.
				
				box_it->set( &ext[0] );
				++box_it;
			}
		}
		return;
	} else if ( dim == 3 ) {
		int num_boxes = (pts[0].size()-1)*(pts[1].size()-1)*(pts[3].size()-1);
		assert( num_boxes >0 );
		boxes.resize( num_boxes );
		typename container2::iterator box_it( boxes.begin() );
		typename container_t::iterator it( pts[0].begin() ), end( pts[0].end() );
		typename container_t::iterator middle_it( pts[1].begin() ), middle_end( pts[1].end() );
		typename container_t::iterator outer_it( pts[2].begin() ), outer_end( pts[2].end() );
		double lower, upper = *it;
		++it;
		double middle_lower, middle_upper = *middle_it;
		++middle_it;
		double outer_lower, outer_upper = *outer_it;
		++outer_it;
		for ( ; outer_it != outer_end; ++outer_it ) {
			outer_lower = outer_upper;
			outer_upper = *outer_it;
			for ( ; middle_it != middle_end; ++middle_it ) {
				middle_lower = middle_upper;
				middle_upper = * middle_it;
				for ( ; it != end; ++it ) {
					lower = upper;
					upper = *it;
					double ext[6] = { lower, upper, middle_lower, middle_upper, outer_lower, outer_upper };
					box_it->set( &ext[0] );
					++box_it;
				}
			}
		}
		return;
	}
}
#endif

template<int dim>
void
make_boxes ( const vector<vector<double> >& in, vector<box<dim> >& out )
{
	assert( in.size() == dim );
	
	int num_intervals[dim];
	int total_num = 1;
	for ( int d=0; d<dim; ++d ) {
		num_intervals[d] =  in[d].size() -1;
		total_num        *= in[d].size() -1;
		
		assert( in[d].size() >= 2 );
		
		// TODO: REMOVE
		//cout << "in[" << d << "].size() is " << in[d].size() << "\n";
		

	}
	out.resize( total_num );
	
	int running_product[dim];
	
	running_product[0] = num_intervals[0];
	// Loop deliberately starts at 1.
	for ( int d=1; d<dim; ++d ) {
		running_product[d] = running_product[d-1]*num_intervals[d];
		assert( running_product[d] > 0 );
	}
	
	double tmp_ext[2*dim];
	
	int indices[dim];
	
	int remainder;
	div_t divresult;
	
	// Go through all the boxes we want to do and work out the
	// indices into the vectors of double so we can make the boxes.
	for ( int i=0; i<total_num; ++i ) {
		
#if 0 // ndef NDEBUG // Keep for future debug.
		cout << "mb Total_num " << total_num << "\n";
#endif
		
		remainder = i;
		for ( int d=dim-1; d>0; --d ) {
			
			// GDB helped find that I had d in place of [d-1] here.
			
			divresult = div( remainder, running_product[d-1] );
			remainder  = divresult.rem;
			indices[d] = divresult.quot;
		}
		indices[0] = remainder;
		for ( int d=0; d<dim; ++d ) {
			
			// TODO: remove
			//cout << "make boxes indices["<<d<<"] is "<<indices[d]<<"\n";
			
			tmp_ext[2*d     ] = in[d][ indices[d]     ];
			tmp_ext[2*d + 1 ] = in[d][ indices[d] + 1 ];
			
			assert( tmp_ext[2*d] < tmp_ext[2*d + 1 ] );
		}
		out[i].set( &tmp_ext[0] );
	}
}

void make_box_with_radius ( double * centre, double radius, box<2>& bx )
{
	radius = fabs( radius );
	
	double    ext[2*2];
	
	for ( int d=0; d<2; ++d ) {
		ext[2*d] = -radius;  ext[2*d+1] = radius;
	}
	
	bx.set ( ext );
	bx.translate( centre );
}

template<int dim>
void box_project_axis ( const box<dim>& in, int axis, const box<dim-1>& out )
{
	assert( 0<=axis && axis<dim );
	const double* old_extents = in.get();
	double new_extents[2*(dim-1)];
	for ( int d=0, new_d=0; d<dim; ++d ) {
		if ( d == axis ) {
			continue;
		} else {
			new_extents[ new_d     ] = old_extents[ d     ];
			new_extents[ new_d + 1 ] = old_extents[ d + 1 ];
			++new_d;
		}
	}
	out.set( &new_extents[0] );
}

void gp_draw_single( const box<1>& box, std::ostream& str, double level )
{
  const double* ext = box.get();
  str << ext[0] << " " << level <<"\n";
  str << ext[1] << " " << level <<"\n";
}

void gp_draw_single( const box<2>& box, std::ostream& str, double level )
{
  const double* ext = box.get();
  // bottom left point
  str << ext[0] << " " << ext[2] << " " << level
      << "\n";
  // bottom right point
  str << ext[1] << " " << ext[2] << " " << level
      << "\n";
  // top right point
  str << ext[1] << " " << ext[3] << " " << level
      << "\n";
  // top left point
  str << ext[0] << " " << ext[3] << " " << level
      << "\n";
  // return to bottom left
  str << ext[0] << " " << ext[2] << " " << level
      << "\n";
  // leave line to end square
  str << "\n";
}


void gp_draw_single( const box<3>& box, std::ostream& str, double ignored )
{
  const double* ext = box.get();
  // first z-coord bottom left point
  str << ext[0] << " " << ext[2] << " " << ext[4]
      << "\n";
  // first z-coord bottom right point
  str << ext[1] << " " << ext[2] << " " << ext[4]
      << "\n";
  // first z-coord top right point
  str << ext[1] << " " << ext[3] << " " << ext[4]
      << "\n";
  // first z-coord top left point
  str << ext[0] << " " << ext[3] << " " << ext[4]
      << "\n";
  // first z-coord return to bottom left
  str << ext[0] << " " << ext[2] << " " << ext[4]
      << "\n";
  // leave line to end square
  str << "\n";

  // second z-coord bottom left point
  str << ext[0] << " " << ext[2] << " " << ext[5]
      << "\n";
  // second z-coord bottom right point
  str << ext[1] << " " << ext[2] << " " << ext[5]
      << "\n";
  // second z-coord top right point
  str << ext[1] << " " << ext[3] << " " << ext[5]
      << "\n";
  // second z-coord top left point
  str << ext[0] << " " << ext[3] << " " << ext[5]
      << "\n";
  // second z-coord return to bottom left
  str << ext[0] << " " << ext[2] << " " << ext[5]
      << "\n";
  // leave line to end square
  str << "\n";

  // leave line to end cube
  str << "\n";
}

template<int dim>
void
gp_draw( const std::vector<box<dim> >& vec,
	 std::ostream& str, double level )
{
  for ( size_t i=0; i<vec.size(); ++i ) {
    gp_draw_single( vec[i], str, level );
    str << "\n";
  }
}

// 

template void decompose_box<1>( const box<1>&, const vector<box<1> >&, vector<box<1> >& );
template void decompose_box<2>( const box<2>&, const vector<box<2> >&, vector<box<2> >& );
template void decompose_box<3>( const box<3>&, const vector<box<3> >&, vector<box<3> >& );

template void decompose_box( const box<1>&, const box<1>::vec_ptr&, vector<box<1> >& );
template void decompose_box( const box<2>&, const box<2>::vec_ptr&, vector<box<2> >& );
template void decompose_box( const box<3>&, const box<3>::vec_ptr&, vector<box<3> >& );

template void local_to_restricted_local<1> ( valarray<double>& points,
					     const box<1>& restricted_local,
					     const box<1>& local_box,
					     valarray<double>& out );

template void local_to_restricted_local<2> ( valarray<double>& points,
					     const box<2>& restricted_local,
					     const box<2>& local_box,
					     valarray<double>& out );

template void local_to_restricted_local<3> ( valarray<double>& points,
					     const box<3>& restricted_local,
					     const box<3>& local_box,
					     valarray<double>& out );

template void output_description<1> ( const box<1>&,
				      ostream& );

template void output_description<2> ( const box<2>&,
				      ostream& );

template void output_description<3> ( const box<3>&,
				      ostream& );

template void make_bounding_box<1> ( const vector<box<1> >&,
				     box<1>& );

template void make_bounding_box<2> ( const vector<box<2> >&,
				     box<2>& );

template void make_bounding_box<3> ( const vector<box<3> >&,
				     box<3>& );

#if 0
  template void split_box<1, std::stack<box<1> > > ( const box<1>&,
						     std::stack<box<1> >& );

  template void split_box<2, std::stack<box<2> > > ( const box<2>&,
						       std::stack<box<2> >& );

  template void split_box<3, std::stack<box<3> > > ( const box<3>&,
						       std::stack<box<3> >& );
#endif

template int count_closed_intersect( const box<1>&, valarray<double>& );
template int count_closed_intersect( const box<2>&, valarray<double>& );
template int count_closed_intersect( const box<3>&, valarray<double>& );

template int count_open_intersect( const box<1>&, valarray<double>& );
template int count_open_intersect( const box<2>&, valarray<double>& );
template int count_open_intersect( const box<3>&, valarray<double>& );

template void make_boxes ( const vector<vector<double> >& in, vector<box<1> >& out );
template void make_boxes ( const vector<vector<double> >& in, vector<box<2> >& out );
template void make_boxes ( const vector<vector<double> >& in, vector<box<3> >& out );

#if 0
template void make_boxes( int, const set<double>*, vector<box<1> >& );
template void make_boxes( int, const set<double>*, vector<box<2> >& );
template void make_boxes( int, const set<double>*, vector<box<3> >& );
#endif

template void gp_draw( const vector<box<1> >&, std::ostream&, double );
template void gp_draw( const vector<box<2> >&, std::ostream&, double );
template void gp_draw( const vector<box<3> >&, std::ostream&, double );

// template void closed_intersect_point<2> ( const vector<box<2> > & in, const double* co, valarray<bool>& out );
