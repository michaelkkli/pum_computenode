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

#include "extract_geometry_2d.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <valarray>
#include <vector>
using std::cerr;
using std::cout;
using std::equal;
using std::getline;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::stringstream;
using std::valarray;
using std::vector;

template<typename T>
void simple_off_output_geometry_2d( const valarray<T>& roi, ostream& out, T level=0 )
{
	assert( roi.size()%2 == 0 );
	
	int num_pts = roi.size()/2;
	
	out << "OFF\n" << num_pts << " 0 0\n";
	
	for ( int p=0; p<num_pts; ++p ) {
		out << roi[2*p] << " " << roi[2*p+1] << " " << level << "\n";
	}
	out << "\n";
}



template<typename T>
void take_points_subset_2d( const valarray<T>& in, int every, valarray<T>& out )
{
	assert( in.size() > 0 );
	assert( every >= 1 );
	assert( in.size()%2 == 0 );
	
	const int num_pts = in.size()/2;
	
	cout << "There are " << num_pts << " points.\n";
	
	const int remainder = num_pts%every;
	const int num_keep_pts = (num_pts - remainder)/every;
	
	cout << "Keep " << num_keep_pts << " points.\n";
	
	out.resize(2*num_keep_pts);
	
	for ( int i=0; i<num_keep_pts; ++i ) {
//		cout << "Kept point " << i << " is " << every*i << " value " << in[ 2*every*i   ] << "\n";
		out[2*i]   = in[ 2*every*i   ];
		out[2*i+1] = in[ 2*every*i+1 ];
	}
}

int main ( int argc, char* argv[] ) {

	string    in_filename;
	
	if ( argc == 1 ) {
		in_filename = "F2Rep 1 FRAP Series11_t00_ch01.tif";
	} else {
		in_filename.assign ( argv[1] );
	}
	
	size_t    dir_sep    = in_filename.find_last_of ( '/' );
	size_t    suffix_beg = in_filename.find_last_of ( '.' );
	
	if ( suffix_beg == string::npos ) {
		cout << "No suffix found. Exiting.\n";
		return 1;
	}
	
	string    base_name;
	
	if ( dir_sep != string::npos ) {
		if ( dir_sep != string::npos ) {
			base_name = in_filename.substr ( dir_sep + 1, suffix_beg );
		} else {
			base_name = in_filename.substr ( 0, suffix_beg );
		}
	}
	
	stringstream    ss;
	
	ss.str ( "extract_geometry_2d-1_" );
	ss << base_name << ".off";

	valarray<double> contour;
	extract_geometry_2d( in_filename.c_str(), contour, 0.2*255.0 );
	{
	  // Broken.
	  // ofstream file( (ss.str()).c_str() );
	  ofstream    file ( "extract_geometry_2d-1_contour.off" );
		simple_off_output_geometry_2d( contour, file );
	}
	
	ss.str ( "extract_geometry_2d-1_" );
	ss << base_name << "_thinned.off";

	{
		valarray<double> thinned;
		
		// Reduce number of points.
		take_points_subset_2d( contour, 10, thinned );
		// Broken.
		//		ofstream file( (ss.str()).c_str() );
	  ofstream    file ( "extract_geometry_2d-1_contour_thinned.off" );
		simple_off_output_geometry_2d( thinned, file );
	}

	return 0;
}
