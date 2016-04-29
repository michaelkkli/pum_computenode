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

#ifndef _SINGLEIMAGE_HH_
#define _SINGLEIMAGE_HH_

#include "differentiable_function.hh"
#include "function.hh"

#include <boost/shared_ptr.hpp>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>

using boost::shared_ptr;
using std::cerr;
using std::copy;
using std::cout;
using std::fstream;
using std::getenv;
using std::getline;
using std::ifstream;
using std::max;
using std::min;
using std::ofstream;
using std::string;
using std::valarray;
using std::vector;

class singleimage : public differentiable_function<2> {
public:
	typedef shared_ptr<singleimage>    ptr;
	singleimage();
	~singleimage();
	singleimage( const singleimage& other );
	singleimage& operator= ( const singleimage& other );
private:
	virtual double evaluate ( const double* ) const;
	virtual void evaluate_grad ( const double* co, double* grad ) const;
public:
	void deep_copy( const singleimage& );
#if 0
	static void draw_rectangle_on_image( int channel,
	                                     const string& filein,
	                                     int x0,
	                                     int x1,
	                                     int y0,
	                                     int y1,
	                                     const string& fileout );
	static void draw_polygon_on_image( int channel,
	                                   const string& filein,
	                                   const valarray<double>&,
	                                   const string& fileout );
#endif


public:

	void load( const string& filename, int channel );
	void rescale ( double scaling );

	// Positive values are taken as the maximum.
	// Returns the maximum that has been assumed.
	double rescale_percentage ( double select_max=-1.0 );
	void set_roi( double x0, double x1, double y0, double y1 );
	int get_x_size() const;
	int get_y_size() const;
public:
	double evaluate_global( double x, double y ) const;
	double evaluate_local( double x, double y ) const; // x, y in [-1,1]
	
	/**
		For use with 3D points in vtk visualization.
	*/
	void evaluate_global_3D_pts ( const double* in, int num_pts, double* out ) const;
	
	void save_local_gnuplot( const string& fileout,
	                         int points ) const;

private:
	bool               have_data;
	bool               roi_set;
	int                x_size;
	int                y_size;
	valarray<double>    data;
	double             roi_box[4]; // x0, x1, y0, y1
	double             centre[2];
	double             half_extent[2];
};



#endif // _SINGLEIMAGE_HH_
