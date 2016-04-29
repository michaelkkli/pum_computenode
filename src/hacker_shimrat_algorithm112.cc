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

#include "hacker_shimrat_algorithm112.hh"
#include "line_segments.hh"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <valarray>

using std::copy;
using std::cout;
using std::memcpy;
using std::slice;
using std::valarray;

bool point_in_polygon ( int n, double* in_x, double* in_y, double x0, double y0 )
{
	assert( in_x && in_y );
	
	valarray<double> x(n+1), y(n+1);
	
	// copy( in_x, in_x+n, x.begin() );
	// copy( in_y, in_y+n, y.begin() );
	memcpy ( &x[0], in_x, n*sizeof(double) );
	memcpy ( &y[0], in_y, n*sizeof(double) );
	
	x[n]=x[0];   y[n]=y[0];
	
	bool b = true;
	
	// TODO: check this and comment properly.
	for ( int i=0; i<n; ++i ) {
		// Hacker corrected typo y -> y0 and changed "strictly less"
		// to "less than or equal"
		// if ( y<y[i] == y>y[i+1] ) || ...
		// 	to
		// if ( ( y0<=y[i] == y0>y[i+1] ) && ...
		if ( ( (y0<=y[i]) == (y0>y[i+1]) ) &&
		     ( (x0-x[i]-(y0-y[i])*(x[i+1]-x[i])/(y[i+1]-y[i]))<0.0 ) )
		{
			b = !b;
		}
	}
	
	return !b;
}

bool point_in_polygon ( line_segments & ls, double x0, double y0 )
{
	int    n = ls.segments.size()/2 - 1;
	
	// Profiling has shown this is very inefficient with
	// lots of allocation and deallocation.
	// 2009-01-21 ML.
	// const valarray<double>    x = ls.segments[ slice (0, n+1, 2) ];
	// const valarray<double>    y = ls.segments[ slice (1, n+1, 2) ];
	
	double *    segco = &(ls.segments[0]);
	// x[i] == segco[2*i]
	// y[i] == segco[2*i +1]
	
	bool b = true;
	
	for ( int i=0; i<n; ++i ) {
		// Hacker corrected typo y -> y0 and changed "strictly less"
		// to "less than or equal"
		// if ( y<y[i] == y>y[i+1] ) || ...
		// 	to
		// if ( ( y0<=y[i] == y0>y[i+1] ) && ...
		if ( ( (y0<=segco[2*i+1]) == (y0>segco[2*(i+1)+1]) ) &&
		     ( (x0-segco[2*i]-(y0-segco[2*i+1])*(segco[2*(i+1)]-segco[2*i])/(segco[2*(i+1)+1]-segco[2*i+1]))<0.0 ) )
		{
			b = !b;
		}
	}
	
	return !b;
}

