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

#include "basic_dbinary_tree.hh"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <cassert>
#include <cmath>
#include <cstring>

using std::clog;
using std::cout;
using std::copy;
using std::find;
using std::inserter;
using std::pow;
using std::strcmp;
using std::strcpy;


template <int dim>
int
key_level ( const string& key )
{
	if ( key == "" ) {
		return -1;
	}
	return (key.size() - 1)/2;
}



#if 0
template <int dim>
string
project_key ( const string& key, int coord )
{
	assert( 0<=coord && coord<dim );
	
	string ret = key;
//    cout << "ret is "<<ret<<"\n";
	if ( coord == 0 ) {
		for ( size_t i=0; i<ret.length(); ++i ) {
			if ( ret[i] == '-' ) {
				continue;
			}
			if ( ret[i] == '0' || ret[i] == '2' || ret[i] == '4' || ret[i] == '6' ) {
				ret[i] = '0';
			} else {
				ret[i] = '1';
			}
		}
	} else if ( coord == 1 ) {
		for ( size_t i=0; i<ret.length(); ++i ) {
			if ( ret[i] == '-' ) {
				continue;
			}
			if ( ret[i] == '0' || ret[i] == '1' || ret[i] == '4' || ret[i] == '5' ) {
				ret[i] = '0';
			} else {
				ret[i] = '1';
			}
		}
	} else if ( coord == 2 ) {
		for ( size_t i=0; i<ret.length(); ++i ) {
			if ( ret[i] == '-' ) {
				continue;
			}
			if ( ret[i] == '0' || ret[i] == '1' || ret[i] == '2' || ret[i] == '3' ) {
				ret[i] = '0';
			} else {
				ret[i] = '1';
			}
		}
	}
	return ret;
}

template <>
string
project_key<1> ( const string& key, int coord )
{
	assert( coord==0 );
	return key;
}
#endif

template <>
void
project_key<1> ( const string& in, int coord, string& out )
{
	assert( coord == 0 );
	
	out = in;
}

template <>
void
project_key<2> ( const string& in, int coord, string& out )
{
	assert( 0<=coord && coord<2 );
	

	int size = in.size();
	
	// Valgrind tool memcheck showed invalid write of size 1.
	// Added space for terminating null character.
	// 2008-09-01 Mike Li.
	valarray<char> tmp(size + 1);
	{
		const char* c_str = in.c_str();
		strcpy ( &tmp[0], c_str );
	}
	
	// Regression caused by two things. Last entry of key was not projected.
	// Key character equal to '0' was unhandled. 2008-08-20 Mike Li.
	
	// Make sure we reach the last entry of key.
	int end = (size-1)/2 + 1;
	
	if ( coord == 0 ) {
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			
			// Regression caused by missing c == '0' check. 2008-08-20 Mike Li.
			
			if ( c == '0' || c == '2' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}
	} else {
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			
			// Regression caused by missing c == '0' check. 2008-08-20 Mike Li.
			
			if ( c == '0' || c == '1' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}	
	}
	
	out.assign ( &tmp[0], size );
	
#if 0 //ndef NDEBUG // Keep for future debug.
	std::clog << "\t\t\t" << in << " projected to coordinate " << coord << " is " << out << "\n";
#endif
}

template <>
void
project_key<3> ( const string& in, int coord, string& out )
{
	assert( 0<=coord && coord<3 );
	
	int size = in.size();
	
	char tmp[size+1];
	{
		const char* c_str = in.c_str();
		strcpy ( tmp, c_str );
	}
	
	// Make sure we reach last entry of key.
	int end = (size-1)/2 + 1;
	
	if ( coord == 0 ) {
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			if ( c == '0' || c == '2' || c == '4' || c== '6' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}
	} else if ( coord == 1 ) {
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			if ( c == '0' || c == '1' || c == '4' || c == '5' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}	
	} else {
		assert ( coord == 2 );
		
		for ( int i=0; i<end; ++i ) {
			char& c = tmp[2*i];
			if ( c == '0' || c == '1' || c == '2' || c == '3' ) {
				// Equality to '0' is left unchanged.
				c = '0';
			} else {
				c = '1';
			}
		}	
	}
	out.assign ( tmp, size );
}



template <>
void reconstruct_from_projected_keys<1> ( const vector<string>& in, string & out )
{
	assert ( in.size() == 1 );
	out = in[0];
}

template <>
void reconstruct_from_projected_keys<2> ( const vector<string>& in, string & out )
{
	int in_size = in.size();
	assert ( in_size == 2 );
	
	int str_size = in[0].size();
	int level    = (str_size - 1)/2;
	
	// Only safe thing to do is to initialize inside even though it may be done
	// outside already. Perhaps provide bool argument to allow choice.
	out = in[0];
	
	char tmp[] = "XX";
	
	//	clog << "Projections are " << in[0] << " and " << in[1] << ".\n";
	
	for ( int i=0; i<=level; ++i ) {	
		tmp[0] = in[0][2*i];
		tmp[1] = in[1][2*i];
		
		if ( strcmp(tmp,"00") == 0 ) {
			out[2*i] = '0';
		} else if ( strcmp(tmp,"10") == 0 ) {
			out[2*i] = '1';
		} else if ( strcmp(tmp,"01") == 0 ) {
			out[2*i] = '2';
		} else if ( strcmp(tmp,"11") == 0 ) {
			out[2*i] = '3';
		} else {
			abort();
		}
	}
	//	clog << "Reconstructed key is " << out << "\n";
	
}

template <>
void reconstruct_from_projected_keys<3> ( const vector<string>& in, string & out )
{
	int in_size = in.size();
	assert ( in_size == 3 );
	
	int level = key_level< 3 >(in[0]);
	
	out = in[0];
	
	string tmp ( "XXX" );
	
	for ( int i=0; i<=level; ++i ) {
		tmp[0] = in[0][2*i];
		tmp[1] = in[1][2*i];
		tmp[2] = in[2][2*i];
		
		if ( "000" == tmp ) {
			out[2*i] = '0';
		} else if ( "100" == tmp ) {
			out[2*i] = '1';
		} else if ( "010" == tmp ) {
			out[2*i] = '2';
		} else if ( "110" == tmp ) {
			out[2*i] = '3';
		} else if ( "001" == tmp ) {
			out[2*i] = '4';
		} else if ( "101" == tmp ) {
			out[2*i] = '5';
		} else if ( "011" == tmp ) {
			out[2*i] = '6';
		} else if ( "111" == tmp ) {
			out[2*i] = '7';
		} else {
			abort();
		}
	}
}




template <>
void reconstruct_from_projected_keys<1> ( const vector<string> (&in)[1], vector<string> & out )
{
	out = in[0];
}

template <>
void reconstruct_from_projected_keys<2> ( const vector<string> (&in)[2], vector<string> & out )
{
	
	const int num_x = in[0].size ();
	const int num_y = in[1].size ();
	
	int level = ( in[0][0].size() - 1 )/2;
	
	string base ( in[0][0].size(), '-' );
	
	for ( int i=0; i<level; ++i ) {
		base.at ( 2*i ) = 0;
	}
	
	out.resize ( num_y * num_x, base );
	
	// vector<string> to_reconstruct(2);
	
	char tmp[] = "XX";
	
	for ( int y=0; y<num_y; ++y ) {
		// to_reconstruct[1] = in[1][y];
		for ( int x=0; x<num_x; ++x ) {
			// to_reconstruct[0] = in[0][x];
			for ( int i=0; i<=level; ++i ) {	
				tmp[0] = in[0][x][2*i];
				tmp[1] = in[1][y][2*i];
				
				string& ref_out = out[y*num_x+x];
				
				if ( strcmp(tmp,"00") == 0 ) {
					ref_out[2*i] = '0';
				} else if ( strcmp(tmp,"10") == 0 ) {
					ref_out[2*i] = '1';
				} else if ( strcmp(tmp,"01") == 0 ) {
					ref_out[2*i] = '2';
				} else if ( strcmp(tmp,"11") == 0 ) {
					ref_out[2*i] = '3';
				} else {
					abort();
				}
			}
			
			// reconstruct_from_projected_keys<2>( to_reconstruct, out[ y*num_x + x ] );
		}
	}
}

template <>
void reconstruct_from_projected_keys<3> ( const vector<string> (&in)[3], vector<string> & out )
{
	const int num_x = in[0].size ();
	const int num_y = in[1].size ();
	const int num_z = in[2].size ();
	
	string base ( 2*key_level<3>(in[0][0])+1, 'X' );
	out.resize ( num_x * num_y * num_z, base );
	
	vector<string> to_reconstruct(3);
	
	for ( int z=0; z<num_z; ++z ) {
		for ( int y=0; y<num_y; ++y ) {
			to_reconstruct[1] = in[1][y];
			for ( int x=0; x<num_x; ++x ) {
				to_reconstruct[0] = in[0][x];
				reconstruct_from_projected_keys<2>( to_reconstruct,
				                                    out[ z*num_y*num_x + y*num_x + x ] );
			}
		}
	}
}


template <int dim>
basic_dbinary_tree<dim>::basic_dbinary_tree ()
: max_level_(0)
{
	double tmp[2*dim];
	for ( int d=0; d<dim; ++d ) {
		tmp[2*d]   = -1.0;
		tmp[2*d+1] =  1.0;
	}
	box<dim> tmp_box;
	tmp_box.set( tmp );
	
	initialize( tmp_box, 1 );
}

template <int dim>
void
basic_dbinary_tree<dim>::initialize ( const box<dim>& par_box, int max )
{
	assert( !par_box.empty() );
	parent_box_ = par_box;
	
	assert( 0< max && max <= 1000000 );
	max_level_ = max>0 ? max : 4 ;
	
	// Call of overload is ambiguous and (int,int) is not catered for -
	// must use some float variant if want to use pow.
	// num_1d_intervals_ = pow( 2, max );
	num_1d_intervals_ = ( 1 << max_level_ );
	
	num_1d_endpoints_ = num_1d_intervals_ + 1;
	
	succ_halving_num_1d_intervals_.resize( max_level_ + 1 );
	succ_halving_num_1d_intervals_[0] = num_1d_intervals_;
	for ( int i=1; i<max_level_+1; ++i ) {
		succ_halving_num_1d_intervals_[i] =
				succ_halving_num_1d_intervals_[i-1]/2;
	}
	assert ( succ_halving_num_1d_intervals_[max_level_] == 1 );
	
	for ( int i=0; i<dim; ++i ) {
		all_endpoints_[i].resize( num_1d_endpoints_ );
	}
	
	const box<dim>& tmp_const_bx_ref = parent_box_;
	const double* bx_pts = tmp_const_bx_ref.get();
	double A, B, step, length;
	for ( int d=0; d<dim; ++d ) {
		A = bx_pts[2*d];
		B = bx_pts[2*d+1];
		
		length = B-A;
		
		step = length/num_1d_intervals_;
		
		level_to_interval_length_[d].resize ( max_level_ + 1 );
		for ( int lev=0; lev<max_level_+1; ++lev ) {
			level_to_interval_length_[d][lev] = length/succ_halving_num_1d_intervals_[lev];
		}
		
		for ( int i=0; i<num_1d_endpoints_; ++i ) {
			// Outputting using vtp format found that [i]
			// was missing. 2008-08-22 Mike Li.
			all_endpoints_[d][i] = A + i*step;
		}
	}
	
	// Increasing powers of the number of end points.
	num_1d_endpoints_pows_.resize( dim+1 );
	num_1d_endpoints_pows_[0] = 1;
	for ( int i=1; i<=dim; ++i ) {
		// Debugging segfault found + in place of * here.
		// Modifications were being made separately on night and archangel/hilbert.
		// Night was still using pow with a static_cast<long double> to pacify
		// icpc and pgCC on francesca and bluebird. archangel/hilbert had an
		// incorrect implementation here hence non-working on archangel/hilbert.
		// 2008-08-27 Mike Li.
		num_1d_endpoints_pows_[i] = num_1d_endpoints_pows_[i-1] * num_1d_endpoints_;
	}
}

template <int dim>
int
basic_dbinary_tree<dim>::lower_1d_interval ( const string& k, int coord ) const
{
	assert ( 0<=coord && coord<dim );
	
	string tmp;
	project_key<dim>( k, coord, tmp );
	
	int lev = key_level<dim>(k);
	
	int count = 0;
	
	// Upper limit allow equal to do the level itself.
	for ( int i=0; i<=lev; ++i ) {
		if ( tmp[2*i] == '1' ) {
			count += succ_halving_num_1d_intervals_[i];
		}
	}
	return count;
}

template <int dim>
int
basic_dbinary_tree<dim>::num_intervals_covered ( const string& k ) const
{
	int lev = key_level<dim>( k );
	return succ_halving_num_1d_intervals_[lev];
}


/**
 * A simple optimization would be to copy the returned expression directly into
 * code. This function should definitely be left here for reference, regardless.
 */
template <int dim>
int
basic_dbinary_tree<dim>::num_intervals_covered ( int lev ) const
{
	return succ_halving_num_1d_intervals_[lev];
}

template <int dim>
int
basic_dbinary_tree<dim>::upper_1d_interval ( const string& k, int coord ) const
{
	assert ( 0<=coord && coord<dim );
	
	int lower = lower_1d_interval( k, coord );
	
	int covered = num_intervals_covered( k );
	
	return lower + (covered - 1);
}

template <int dim>
void
basic_dbinary_tree<dim>::get_key_1d_with_lower_end ( int level,
                                                     int desired_lower_end,
                                                     string & out ) const
{
	assert ( level >= 0 );
	
	assert ( desired_lower_end < num_1d_intervals_ );
	
	if ( level == 0 ) {
		assert ( desired_lower_end == 0 );
		out = "0";
	}
	
  out.assign( 2*level+1, '0');
  for ( int i=0; i<level; ++i ) {
    out.at( 1+2*i ) = '-';
  }
	
	int left_limit = 0;
	int count      = 1;
	while ( left_limit < desired_lower_end ) {
		
		assert ( count < succ_halving_num_1d_intervals_.size() );
		
		if ( desired_lower_end - left_limit >= succ_halving_num_1d_intervals_[count] ) {
			left_limit += succ_halving_num_1d_intervals_[count];
			assert ( count <= level );
			out[ 2*count] = '1';
		}
		++count;
	}
	
	// Ensure we were given a task that was possible.
	assert ( left_limit == desired_lower_end );
}

template <int dim>
void
basic_dbinary_tree<dim>::get_key_1d_with_upper_end ( int level, int upper_end, string & out ) const
{
	get_key_1d_with_lower_end ( level, upper_end - num_intervals_covered(level), out );
}


template <>
void
basic_dbinary_tree<1>::get_key_intersect_point ( int               level,
                                                 const double *    co,
                                                 string &          out ) const
{
	assert ( co );
	get_key_1d_intersect_point ( level, co[0], 0, out );
	
#ifndef NDEBUG
	box<1> containing;
	get_box ( out, containing );
	assert ( containing.closed_intersect_point( co ) );
#endif
}

template <>
void
basic_dbinary_tree<2>::get_key_intersect_point ( int               level,
                                                 const double *    co,
                                                 string &          out ) const
{
	vector<string> projections(2);
	
	get_key_1d_intersect_point ( level, co[0], 0, projections[0] );
	get_key_1d_intersect_point ( level, co[1], 1, projections[1] );
	
	reconstruct_from_projected_keys<2> ( projections, out );
	
#if 0//ndef NDEBUG
	box<2> containing;
	get_box ( out, containing );
	assert ( containing.closed_intersect_point( co ) );
#endif
}

template <>
void
basic_dbinary_tree<3>::get_key_intersect_point ( int               level,
                                                 const double *    co,
                                                 string &          out ) const
{
	vector<string> projections(3);
	
	get_key_1d_intersect_point ( level, co[0], 0, projections[0] );
	get_key_1d_intersect_point ( level, co[1], 1, projections[1] );
	get_key_1d_intersect_point ( level, co[2], 2, projections[2] );
	
	reconstruct_from_projected_keys<3> ( projections, out );
	
#if 0//ndef NDEBUG
	box<3> containing;
	get_box ( out, containing );
	assert ( containing.closed_intersect_point( co ) );
#endif
}

template <int dim>
void
basic_dbinary_tree<dim>::get_key_1d_intersect_point ( int         level,
                                                      double      co,
                                                      int         d,
                                                      string &    out ) const
{
	int num_interval_to_jump = succ_halving_num_1d_intervals_[level];
	
	assert ( 0<= d && d <=3 );
	const valarray<double> & endpts = all_endpoints_[d];
	
	assert ( !( co < endpts[0] ) );
	
	int upper = num_interval_to_jump;
	
	// Attempting to replace dbinary_tree with basic_dbinary_tree
	// causes the assert to be triggered when co is right at the edge of the
	// bounding box. Attempting to avoid by using < rather than <=.
	// 2009-02-13 ML.
	while ( endpts[ upper ] < co ) {
		assert ( num_interval_to_jump > 0 ); // Avoid infinite loop. 2009-06-27 ML.
		
		upper += num_interval_to_jump;
		
		assert ( upper < endpts.size() );
	}
	
	assert ( upper > 0 );
	assert ( upper <= num_1d_intervals_ );
	get_key_1d_with_upper_end ( level, upper, out );
	
#ifndef NDEBUG
	assert ( all_endpoints_[d][upper-num_interval_to_jump] <= co && co <= all_endpoints_[d][upper] );
#endif
}

template <int dim>
void
basic_dbinary_tree<dim>::get_key_1d_expanded_intersect_point ( int                 level,
                                                               double              co,
                                                               int                 d,
                                                               double              expansion,
                                                               vector<string> &    out ) const
{
	out.clear ();
	
	// Reserve has been shown to be useful (callgrind).
	// 2008-10-01 Mike Li.
	out.reserve ( 2 );
	out.resize ( 1 );
	
	assert ( level < succ_halving_num_1d_intervals_.size() );
	int num_interval_to_jump = succ_halving_num_1d_intervals_[level];
	
	assert ( 0<= d && d <=3 );
	const valarray<double> & endpts = all_endpoints_[d];
	
	int upper = num_interval_to_jump;
	while ( endpts[ upper ] < co ) {
		upper += num_interval_to_jump;
		assert ( upper < endpts.size() );
	}
	
	
	assert ( upper > 0 );
	assert ( upper <= num_1d_intervals_ );
	get_key_1d_with_upper_end ( level, upper, out[0] );
	
	// Weird things happen with too much overlap.
	assert ( expansion < 2.0 );
	
	assert ( endpts[upper] - co >  0.0 - 1e-10 ); // Allow zero.
	
	double percentage_from_interval_end = (endpts[upper] - co) / level_to_interval_length_[d][level];
	
	double percentage_overlap_of_expansion = expansion - 1.0;
	
	if ( 1.0 - percentage_from_interval_end < percentage_overlap_of_expansion ) {
		if ( upper == num_interval_to_jump ) {
			// Would overlap with interval to the left except
			// there is no interval to the left.
			return;
		}
		out.resize ( 2 );
		out[1] = out[0];
		get_key_1d_with_upper_end( level, upper-num_interval_to_jump, out[0] );
		return; // Assumes impossible to overlap with both left and right.
	} else if ( percentage_from_interval_end < percentage_overlap_of_expansion ) {
		if ( upper == num_1d_intervals_ ) {
			// Would overlap with interval to the right except
			// there is no interval to the right.
			return;
		}
		out.resize ( 2 );
		get_key_1d_with_lower_end( level, upper, out[1] );
	}
	
}


template <int dim>
void
basic_dbinary_tree<dim>::get_key_expanded_intersect_point ( int                 level,
                                                            const double *      co,
                                                            double              expansion,
                                                            vector<string> &    out ) const
{
	vector<string> intersects [dim];
	
	for ( int d=0; d<dim; ++d ) {
		get_key_1d_expanded_intersect_point( level,
		                                     co[d],
		                                     d,
		                                     expansion,
		                                     intersects[d] );
	}
	
	reconstruct_from_projected_keys<dim>( intersects, out );
}


template <int dim>
void
basic_dbinary_tree<dim>::get_1d_neighbours ( const string & key, vector<string> & out ) const
{
	int lower, upper;
	get_1d_endpoint_indices( key, 0, lower, upper );
	
	int level = key_level<dim>(key);
	
	if ( lower > 0 && upper < num_1d_intervals_ ) {
		out.resize(3);
		
		get_key_1d_with_upper_end( level, lower, out[0] );
		out[1] = key;
		get_key_1d_with_lower_end( level, upper, out[2] );
		return;
	}
	
	out.reserve(3);
	
	string tmp;
	
	if ( lower != 0 ) {
		get_key_1d_with_upper_end( level, lower, tmp );
		out.push_back(tmp);
	}

	out.push_back(key);
	
	if ( upper < num_1d_intervals_ ) {
		get_key_1d_with_lower_end( level, upper, tmp );
		out.push_back(tmp);
	}
}
template <int dim>
void
basic_dbinary_tree<dim>::get_neighbours ( const string & key, set<string> & out ) const
{
	vector<string> tmp;
	get_neighbours( key, tmp );
	
	out.clear();
	copy ( tmp.begin(), tmp.end(), inserter( out, out.end() ) );
}

template <>
void
basic_dbinary_tree<1>::get_neighbours ( const string & key, vector<string> & out ) const
{
	get_1d_neighbours( key, out );
}

template <>
void
basic_dbinary_tree<2>::get_neighbours ( const string & key, vector<string> & out ) const
{
	string projections[2];
	vector<string> neighbours_1d[2];
	
	for ( int i=0; i<2; ++i ) {
		project_key<2>( key, i, projections[i] );
		get_1d_neighbours( projections[i], neighbours_1d[i] );
	}
	
	int num_x = neighbours_1d[0].size();
	int num_y = neighbours_1d[1].size();
	
	out.resize ( num_x*num_y );
	
	vector<string> to_reconstruct(2);
	
	//	clog << "Expect " << num_x * num_y << " neighbours.\n";
	for ( int y=0; y<num_y; ++y ) {
		to_reconstruct[1] = neighbours_1d[1][y];
		for ( int x=0; x<num_x; ++x ) {
			to_reconstruct[0] = neighbours_1d[0][x];
			reconstruct_from_projected_keys<2>( to_reconstruct, out[ y*num_x + x ] );
			//			clog << "Neighbour found " << out[ y*num_x + x ] << "\n";
		}
	}
	
	assert ( find( out.begin(), out.end(), key) != out.end() );
}

template <>
void
basic_dbinary_tree<3>::get_neighbours ( const string & key, vector<string> & out ) const
{
	string projections[3];
	vector<string> neighbours_1d[3];
	
	for ( int i=0; i<3; ++i ) {
		project_key<3>( key, i, projections[i] );
		get_1d_neighbours( projections[i], neighbours_1d[i] );
	}
	
	const int num_x = neighbours_1d[0].size();
	const int num_y = neighbours_1d[1].size();
	const int num_z = neighbours_1d[2].size();
	
	out.resize ( num_x*num_y*num_z );
	
	vector<string> to_reconstruct(3);
	
	for ( int z=0; z<num_z; ++z ) {
		to_reconstruct[2] = neighbours_1d[2][z];
		for ( int y=0; y<num_y; ++y ) {
			to_reconstruct[1] = neighbours_1d[1][y];
			for ( int x=0; x<num_x; ++x ) {
				to_reconstruct[0] = neighbours_1d[0][x];
				reconstruct_from_projected_keys<3>( to_reconstruct,
				                                 out[ z*num_y*num_z + y*num_x + x ] );
			}
		}
	}
}
template <int dim>
void
basic_dbinary_tree<dim>::get_neighbours ( const set<string> & in, map<string, set<string> >& out ) const
{
	out.clear();
	
	typename set<string>::const_iterator set_it ( in.begin() );
	typename set<string>::const_iterator set_end ( in.end() );
	for ( ; set_it!=set_end; ++set_it ) {
		get_neighbours( *set_it, out[*set_it] );
	}
}

/**
 * Interval n has lower endpoint n and upper endpoint n+1.
 */
template <int dim>
void
basic_dbinary_tree<dim>::get_1d_endpoint_indices ( const string& k,
                                                   int           coord,
                                                   int&          lower,
                                                   int&          upper  ) const
{
	assert ( 0<=coord && coord<dim );
	
	string tmp;
	project_key<dim>( k, coord, tmp );
	
	int lev = key_level<dim>(k);
	assert ( lev > 0 );

	int count = 0;
	
	// Upper limit allow equal to do the level itself.
	for ( int i=0; i<=lev; ++i ) {
		if ( tmp[2*i] == '1' ) {
			assert ( i<succ_halving_num_1d_intervals_.size() );
			count += succ_halving_num_1d_intervals_[i];
		}
	}

	
	lower = count;

	assert ( lower == lower_1d_interval(k, coord) );
	
	int covered = succ_halving_num_1d_intervals_[lev];
	assert ( lev < succ_halving_num_1d_intervals_.size() );
	
	upper = lower + covered;
	
	assert ( upper == upper_1d_interval(k, coord) + 1 );

	// Testing get_key_1d_with_lower_end and get_key_1d_with_upper_end.
#ifndef NDEBUG // Keep for future debug.
	string should_match_tmp;
	get_key_1d_with_lower_end( lev, lower, should_match_tmp );
	assert ( should_match_tmp == tmp );
	get_key_1d_with_upper_end( lev, upper, should_match_tmp );
	assert ( should_match_tmp == tmp );
#endif
	
#if 0 // Vital for get_box so needs optimizing.
	lower = lower_1d_interval(k, coord);
	upper = upper_1d_interval(k, coord) + 1;
#endif
}

template <int dim>
void
basic_dbinary_tree<dim>::get_endpoint_indices ( const string& k, int* ind ) const
{
	assert ( ind );
	for ( int d=0; d<dim; ++d ) {
		get_1d_endpoint_indices( k, d, ind[2*d], ind[2*d+1] );
	}
}

template <int dim>
void
basic_dbinary_tree<dim>::get_box ( const string& k, box<dim>& bx ) const
{
	assert ( key_level<dim>(k) <= max_level_ );
	
	int ind[ 2*dim ];
	get_endpoint_indices( k, ind );
	
	double ext[2*dim];
	
	for ( int d=0; d<dim; ++d ) {
		ext[2*d  ] = (all_endpoints_[d])[ ind[2*d]   ];
		ext[2*d+1] = (all_endpoints_[d])[ ind[2*d+1] ];
	}
	
	bx.set(ext);
}

template <int dim>
void
basic_dbinary_tree<dim>::get_box ( const string& k, box<dim>& bx, double expansion ) const
{
	assert ( key_level<dim>(k) <= max_level_ );
	
	int ind[ 2*dim ];
	get_endpoint_indices( k, ind );
	
	double ext[2*dim];
	
	for ( int d=0; d<dim; ++d ) {
		ext[2*d  ] = (all_endpoints_[d])[ ind[2*d]   ];
		ext[2*d+1] = (all_endpoints_[d])[ ind[2*d+1] ];
	}
	
	bx.set(ext);
	bx.scale(expansion);
}


// Closed intersection.
template <int dim>
bool
basic_dbinary_tree<dim>::box_intersect_point ( const string& k, const double* co, double expansion ) const
{
	assert ( co );
	
	box<dim>    bx;
	
	get_box ( k, bx );
	bx.scale ( expansion );
	return bx.closed_intersect_point ( co );
	
#if 0 // Won't work due to lack of expansion.
	assert ( key_level<dim>(k) <= max_level_ );
	
	int ind[ 2*dim ];
	get_endpoint_indices( k, ind );
	
	for ( int d=0; d<dim; ++d ) {
		if ( co[d] < (all_endpoints_[d])[ ind[2*d]   ] ||
		     co[d] > (all_endpoints_[d])[ ind[2*d+1] ] )
		{
			return false;
		}
	}
	return true;
#endif

}

template <int dim>
bool
basic_dbinary_tree<dim>::are_neighbours ( const string& first, const string& second, double expansion ) const
{
#if 0 // Won't work due to lack of expansion.
	int firstind[2*dim];
	int secondind[2*dim];
	
	get_endpoint_indices( first, firstind );
	get_endpoint_indices( second, secondind );
	

	for ( int d=0; d<dim; ++d ) {
		if ( (all_endpoints_[d])[ firstind[2*d] ] > (all_endpoints_[d])[ secondind[2*d+1] ]
		     || (all_endpoints_[d])[ firstind[2*d+1] ] < (all_endpoints_[d])[ secondind[2*d] ] )
		{
			return false;
		}
	}
	return true;
#endif
	
	box<dim> firstbx, secondbx;
	
	get_box ( first, firstbx );
	firstbx.scale ( expansion );
	
	get_box ( second, secondbx );
	secondbx.scale ( expansion );
	
	return firstbx.open_intersect_box ( secondbx );

}

template <int dim>
void
basic_dbinary_tree<dim>::get_centre_point ( const string&    k,
                                            double *         co ) const
{
	// Based on get_box(const string&,box<dim>&);
	
	assert ( co );
	assert ( key_level<dim>(k) <= max_level_ );
	
	int ind[ 2*dim ];
	get_endpoint_indices( k, ind );
	
	for ( int d=0; d<dim; ++d ) {
		co[d] = 0.5*( (all_endpoints_[d])[ ind[2*d] ] + (all_endpoints_[d])[ ind[2*d+1] ] );
	}
}

template <int dim>
void
basic_dbinary_tree<dim>::get_box ( const vector<string> & in, vector<box<dim> > & out ) const
{
	int in_size = in.size();
	out.resize ( in_size );
	
	int ind[ 2*dim ];
	
	for ( int i=0; i<in_size; ++i ) {
		assert ( key_level<dim>(in[i]) <= max_level_ );
		get_endpoint_indices( in[i], ind );
		
		double ext[2*dim];
		
		for ( int d=0; d<dim; ++d ) {
			assert ( ind[2*d] < ind[2*d+1] );
			
			assert ( ind[2*d] < all_endpoints_[d].size() );
			ext[2*d  ] = (all_endpoints_[d])[ ind[2*d]   ];
			
			assert ( ind[2*d+1] < all_endpoints_[d].size() );
			ext[2*d+1] = (all_endpoints_[d])[ ind[2*d+1] ];
			
			assert ( ext[2*d  ] < ext[2*d+1] );
		}
		
		out[i].set(ext);
		assert ( !out[i].empty() );
	}
}

template <int dim>
void
basic_dbinary_tree<dim>::get_box ( const set<string>& in, vector<box<dim> >& res, double expansion ) const
{
	res.resize( in.size() );
	
	typename vector<box<dim> >::iterator res_it( res.begin() );
	
	typename set<string>::const_iterator it( in.begin() );
	typename set<string>::const_iterator end( in.end() );
	for ( ; it!=end; ++it, ++res_it ) {
		get_box( *it, *res_it );
		res_it->scale ( expansion );
	}
}

#if 0
template <int dim>
void
basic_dbinary_tree<dim>::get_box ( const vector<string>& in, vector<box<dim> >& res, double expansion ) const
{
	get_box ( in, res );
	
	int in_size = in.size();
	for ( int i=0; i<in_size; ++i ) {
		res[i].scale(expansion);
	}
}
#endif

template <>
void
basic_dbinary_tree<1>::get_grid_points ( valarray<double>& out ) const
{
	out.resize( num_1d_endpoints_ );
	
	out = all_endpoints_[0];
}

template <>
void
basic_dbinary_tree<1>::get_grid_points_3D ( valarray<double>& out, double init ) const
{
	// Compiler warning suggests the arguments were the wrong way around.
	//out.resize( init, 3*num_1d_endpoints_ );
	
	out.resize( 3*num_1d_endpoints_, init );
	assert ( out.size() == 3*num_1d_endpoints_ );
	
	for ( int i=0; i<num_1d_endpoints_; ++i ) {
		out[3*i   ] = all_endpoints_[0][i];
		// Other two entries already initialized to init.
	}
}

template <>
void
basic_dbinary_tree<2>::get_grid_points ( valarray<double>& out ) const
{
	out.resize( 2*num_1d_endpoints_pows_[2] );
	
	int ind;
	
	for ( int y=0; y<num_1d_endpoints_; ++y ) {
		for ( int x=0; x<num_1d_endpoints_; ++x ) {
			ind = 2*(y*num_1d_endpoints_ + x); 
			out[ ind     ] = (all_endpoints_[0])[x];
			out[ ind + 1 ] = (all_endpoints_[1])[y];
		}
	}
}

template <>
void
basic_dbinary_tree<2>::get_grid_points_3D ( valarray<double>& out, double init ) const
{
	out.resize( 3*num_1d_endpoints_pows_[2] );
	
	int ind;
	
	for ( int y=0; y<num_1d_endpoints_; ++y ) {
		for ( int x=0; x<num_1d_endpoints_; ++x ) {
			ind = 3*(y*num_1d_endpoints_ + x); 
			out[ ind     ] = (all_endpoints_[0])[x];
			out[ ind + 1 ] = (all_endpoints_[1])[y];
			out[ ind + 2 ] = init;
		}
	}
	
	assert ( out.size() % 3 == 0 );
}

template <>
void
basic_dbinary_tree<3>::get_grid_points ( valarray<double>& out ) const
{
	int N   = num_1d_endpoints_;
	int NN  = num_1d_endpoints_pows_[2];
	int NNN = num_1d_endpoints_pows_[3];
	
	out.resize( 3*NNN );
	
	int ind;
	
	for ( int z=0; z<N; ++z ) {
		for ( int y=0; y<N; ++y ) {
			for ( int x=0; x<N; ++x ) {
				ind = 3*(z*NN + y*N + x);
				out[ ind     ] = (all_endpoints_[0])[x];
				out[ ind + 1 ] = (all_endpoints_[1])[y];
				out[ ind + 2 ] = (all_endpoints_[2])[z];
			}
		}
	}
}

template<>
int
basic_dbinary_tree<1>::get_num_grid_points () const
{
	return num_1d_endpoints_;
}

template<>
int
basic_dbinary_tree<2>::get_num_grid_points () const
{
	return num_1d_endpoints_pows_[2];
}

template<>
int
basic_dbinary_tree<3>::get_num_grid_points () const
{
	return num_1d_endpoints_pows_[3];
}
template <>
void
basic_dbinary_tree<3>::get_grid_points_3D ( valarray<double>& out, double ignored ) const
{
	get_grid_points ( out );
}

template <>
void
basic_dbinary_tree<1>::get_box_grid_indices_stacked ( const string& k, int* gi ) const
{
	get_endpoint_indices ( k, gi );	
}

template <>
void
basic_dbinary_tree<1>::get_box_grid_indices_winding ( const string& k, int* gi ) const
{
	get_endpoint_indices ( k, gi );	
}

template <>
void
basic_dbinary_tree<2>::get_box_grid_indices_stacked ( const string& k, int* gi ) const
{
	int ind[4];
	get_endpoint_indices( k, ind );
	
	int N   = num_1d_endpoints_;
	
	// Clockwise from top left corner of square.
	// ind[0] is x-left, ind[1] is x-right
	// ind[2] is y-botton, ind[3] is y-top
	
	gi[2] = ind[3]*N + ind[0];    gi[3] = ind[3]*N + ind[1];
	gi[0] = ind[2]*N + ind[0];    gi[1] = ind[2]*N + ind[1];
}

template <>
void
basic_dbinary_tree<2>::get_box_grid_indices_winding ( const string& k, int* gi ) const
{
	int ind[4];
	get_endpoint_indices( k, ind );
	
	int N   = num_1d_endpoints_;
	
	// Clockwise from top left corner of square.
	// ind[0] is x-left, ind[1] is x-right
	// ind[2] is y-botton, ind[3] is y-top
	
	gi[3] = ind[3]*N + ind[0];    gi[2] = ind[3]*N + ind[1];
	gi[0] = ind[2]*N + ind[0];    gi[1] = ind[2]*N + ind[1];
}

template <>
void
basic_dbinary_tree<3>::get_box_grid_indices_stacked ( const string& k, int* gi ) const
{
	int ind[6];
	get_endpoint_indices( k, ind );
	
	int N   = num_1d_endpoints_;
	int NN  = num_1d_endpoints_pows_[2];
	
	// Clockwise from top left corner of square.
	// ind[0] is x-left, ind[1] is x-right
	// ind[2] is y-botton, ind[3] is y-top
	// ind[4] is z-down, ind[5] is z-up
	
	gi[6] = ind[5]*NN + ind[3]*N + ind[0];    gi[7] = ind[5]*NN + ind[3]*N + ind[1];
	gi[4] = ind[5]*NN + ind[2]*N + ind[0];    gi[5] = ind[5]*NN + ind[2]*N + ind[1];
	
	gi[2] = ind[4]*NN + ind[3]*N + ind[0];    gi[3] = ind[4]*NN + ind[3]*N + ind[1];
	gi[0] = ind[4]*NN + ind[2]*N + ind[0];    gi[1] = ind[4]*NN + ind[2]*N + ind[1];
}

template <>
void
basic_dbinary_tree<3>::get_box_grid_indices_winding ( const string& k, int* gi ) const
{
	int ind[6];
	get_endpoint_indices( k, ind );
	
	int N   = num_1d_endpoints_;
	int NN  = num_1d_endpoints_pows_[2];
	
	// Clockwise from top left corner of square.
	// ind[0] is x-left, ind[1] is x-right
	// ind[2] is y-botton, ind[3] is y-top
	// ind[4] is z-down, ind[5] is z-up
	
	gi[7] = ind[5]*NN + ind[3]*N + ind[0];    gi[6] = ind[5]*NN + ind[3]*N + ind[1];
	gi[4] = ind[5]*NN + ind[2]*N + ind[0];    gi[5] = ind[5]*NN + ind[2]*N + ind[1];
	
	gi[3] = ind[4]*NN + ind[3]*N + ind[0];    gi[2] = ind[4]*NN + ind[3]*N + ind[1];
	gi[0] = ind[4]*NN + ind[2]*N + ind[0];    gi[1] = ind[4]*NN + ind[2]*N + ind[1];
}

template <>
void
basic_dbinary_tree<1>::get_box_grid_connectivity_offsets ( const vector<string>&    keys,
                                                           valarray<int>&           con,
                                                           valarray<int>&           off,
                                                           int                      con_start,
                                                           int                      off_start )
{
	int num_keys = keys.size();
	assert ( num_keys > 0 );
	
	con.resize ( 2*num_keys );
	off.resize ( num_keys );
	
	for ( int i=0; i<num_keys; ++i ) {
		get_box_grid_indices_winding ( keys[i], &( static_cast<int&>(con[2*i]) ) );
	}
	
	if ( con_start != 0 ) {
		con += con_start;
	}
	
	for ( int i=0; i<num_keys; ++i ) {
		off[i] = 2*(i+1);
	}
	
	if ( off_start != 0 ) {
		off += off_start;
	}
}

template <>
void
basic_dbinary_tree<2>::get_box_grid_connectivity_offsets ( const vector<string>&    keys,
                                                           valarray<int>&           con,
                                                           valarray<int>&           off,
                                                           int                      con_start,
                                                           int                      off_start )
{
	int num_keys = keys.size();
	assert ( num_keys > 0 );
	
	con.resize ( 4*num_keys );
	off.resize ( num_keys );
	
	for ( int i=0; i<num_keys; ++i ) {
		get_box_grid_indices_winding ( keys[i], &( static_cast<int&>(con[4*i]) ) );
	}
	
	if ( con_start != 0 ) {
		con += con_start;
	}
	
	for ( int i=0; i<num_keys; ++i ) {
		off[i] = 4*(i+1);
	}
	
	if ( off_start != 0 ) {
		off += off_start;
	}
}

template <>
void
basic_dbinary_tree<3>::get_box_grid_connectivity_offsets ( const vector<string>&    keys,
                                                           valarray<int>&           con,
                                                           valarray<int>&           off,
                                                           int                      con_start,
                                                           int                      off_start )
{
	int num_keys = keys.size();
	assert ( num_keys > 0 );
	
	con.resize ( 8*num_keys );
	off.resize ( num_keys );
	
	for ( int i=0; i<num_keys; ++i ) {
		get_box_grid_indices_winding ( keys[i], &( static_cast<int&>(con[8*i]) ) );
	}
	
	if ( con_start != 0 ) {
		con += con_start;
	}
	
	for ( int i=0; i<num_keys; ++i ) {
		off[i] = 8*(i+1);
	}
	
	if ( off_start != 0 ) {
		off += off_start;
	}
}

template <int dim>
valarray<double>*
basic_dbinary_tree<dim>::access_endpoints_grid ()
{
	return all_endpoints_;
}

//

template int key_level<1> ( const string& key );
template int key_level<2> ( const string& key );
template int key_level<3> ( const string& key );

// template string project_key<1> ( const string& key, int coordinate );
// template string project_key<2> ( const string& key, int coordinate );
// template string project_key<3> ( const string& key, int coordinate );

template class basic_dbinary_tree<1>;
template class basic_dbinary_tree<2>;
template class basic_dbinary_tree<3>;
