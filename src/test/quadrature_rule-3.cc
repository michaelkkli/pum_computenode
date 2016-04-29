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
#include "function.hh"
#include "geometry_utils.hh"
#include "line_segments.hh"
#include "point_utils.hh"
#include "polynomial.hh"
#include "quadrature_rule.hh"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <valarray>

using std::clog;
using std::copy;
using std::cout;
using std::exp;
using std::fstream;
using std::memcpy;
using std::ofstream;
using std::scientific;
using std::slice;
using std::stringstream;

class integrand : public function<2> {
public:
	integrand () {
		this->set_global_function ();
	}
private:
	double evaluate ( const double* co ) const {
		return 1.0;
		//return co[0];
		//return 700*(co[0]*co[1]+56.7*co[1]*exp(co[0]))/(50.0+co[0]*sin(co[1]));
	}
};



int main ( int argc, char* argv[] ) {

	bool normal_orientation = true;
	
	if ( argc > 3 ) {
		int    num = atoi ( argv[3] );
		if ( 0==num ) {
			normal_orientation = false;
		} else {
			normal_orientation = true;
		}
	}
	
	box<2>    box1;
	double box_ext[] = { 0.0, 1.0, 0.0, 1.0 };
	box1.set ( box_ext );

	line_segments    segs;

	bool circular_boundary = false;
	
	if ( circular_boundary ) {
		valarray<double> vcoords;

		double displace_centre[] = {-1e-10, -1e-10};

		generate_circular_boundary( 100000, 1, displace_centre, vcoords );

		int tmp_size = vcoords.size();

		segs.segments.resize ( tmp_size +2 );
		copy ( &vcoords[0], &(vcoords[0])+tmp_size, &(segs.segments[0]) );

		segs.segments[ tmp_size ]   = vcoords[0];
		segs.segments[ tmp_size+1 ] = vcoords[1];
	} else {
		//	int num_segs = (segs.segments.size()/2) - 1;
		//	int num_pts_inside = 0; // Number of segment end point lying inside.

		vector<double>    seg_endpts;

		int segment_test = 5;
		
		if ( argc > 2 ) {
			stringstream ss;
			ss << argv[2];
			ss >> segment_test;
			if ( !ss ) {
				cout << "Illegal argument given!\n";
				abort();
			}
		}
		
		switch ( segment_test ) {
		case 0:
			// Down to up.
			// double seg_endpts[] = { 0.0, -1.0, 0.75, 0.25, 1.0, 2.0 };
			seg_endpts.push_back ( 0.0 ); seg_endpts.push_back ( -1.0 );
			seg_endpts.push_back ( 0.75 ); seg_endpts.push_back ( 0.25 );
			seg_endpts.push_back ( 1.0 ); seg_endpts.push_back ( 2.0 );
			break;
		case 1:
			// Up to down.
			// double seg_endpts[] = { 1.0, 2.0, 0.75, 0.25, 0.0, -1.0 };
			seg_endpts.push_back ( 1.0 ); seg_endpts.push_back ( 2.0 );
			seg_endpts.push_back ( 0.75 ); seg_endpts.push_back ( 0.25 );
			seg_endpts.push_back ( 0.0 ); seg_endpts.push_back ( -1.0 );
			break;
		case 2:
			// Polar lower right.
			// double seg_endpts[] = { 2.0, 2.0, 0.75, 0.25, 0.0, -1.0 };
			seg_endpts.push_back ( 2.0 ); seg_endpts.push_back ( 2.0 );
			seg_endpts.push_back ( 0.75 ); seg_endpts.push_back ( 0.25 );
			seg_endpts.push_back ( 0.0 ); seg_endpts.push_back ( -1.0 );
			break;
		case 3:
			// Polar upper left.
			// double seg_endpts[] = { 0.0, -1.0, 0.75, 0.25, 2.0, 2.0 };
			seg_endpts.push_back ( 0.0 ); seg_endpts.push_back ( -1.0 );
			seg_endpts.push_back ( 0.75 ); seg_endpts.push_back ( 0.25 );
			seg_endpts.push_back ( 2.0 ); seg_endpts.push_back ( 2.0 );
			break;
		case 4:
			seg_endpts.push_back ( 0.5 ); seg_endpts.push_back ( -0.5 );
			seg_endpts.push_back ( 0.5 ); seg_endpts.push_back ( 0.5 );
			seg_endpts.push_back ( 1.5 ); seg_endpts.push_back ( 0.5 );
			break;
		case 5:
			// Left to right.
			seg_endpts.push_back ( -1.0 ); seg_endpts.push_back ( 0.0 );
			seg_endpts.push_back (  2.0 ); seg_endpts.push_back ( 1.0 );
			seg_endpts.push_back ( -1.0 ); seg_endpts.push_back ( 0.0 );
			break;
		case 6:
			// Polar lower left.
			seg_endpts.push_back ( 1.5 ); seg_endpts.push_back ( 0.0 );
			seg_endpts.push_back ( 0.0 ); seg_endpts.push_back ( 1.5 );
			seg_endpts.push_back ( 1.5 ); seg_endpts.push_back ( 0.0 );
			break;
		default:
			clog << "Switch default: no segment test chosen.\n";
			abort();
		}

		
		
		segs.segments.resize ( 6 );
		memcpy ( &segs.segments[0], &seg_endpts[0], 6*sizeof(double) );
	}

	

	int num_quadpts = 4;

	if ( argc > 1 ) {
		stringstream ss;
		ss << argv[1];
		ss >> num_quadpts;
		if ( !ss ) {
			cout << "Illegal argument given!\n";
			abort();
		}
	}

	int entry_segment = 0;

	simple_quadrature_rule<2>    box_rule;
	simple_quadrature_rule<2>    inside_rule;
	simple_quadrature_rule<2>    outside_rule;

	generate_gauss_legendre_rule ( num_quadpts, box_rule );
	
	simple_quadrature_rule<1>    base_rule;
	generate_gauss_legendre_rule ( num_quadpts, base_rule );
	
	simple_quadrature_rule<2>    boundary_rule;
	
	generate_gauss_legendre_rule_local( base_rule, box1, segs, entry_segment, normal_orientation ? true : false, inside_rule, &boundary_rule );
	
	generate_gauss_legendre_rule_local( base_rule, box1, segs, entry_segment, normal_orientation ? false : true, outside_rule );
	


	{
		// Integrate on part of standard box.
		valarray<double>    temp_val ( inside_rule.weights.size() );

		integrand    my_fun;

		valarray<double>    eval_points;

		double    box_final_factor = box1.measure()/4;
		box1.map_local_to_global ( box_rule.points, eval_points );
		my_fun.global_evaluate ( box1, eval_points, temp_val );
		temp_val *= box_rule.weights;
		double    box_integral = temp_val.sum()*box_final_factor;
		
		box1.map_local_to_global ( inside_rule.points, eval_points );
		my_fun.global_evaluate ( box1, eval_points, temp_val );
		temp_val *= inside_rule.weights;
		double    inside_integral = temp_val.sum();

		box1.map_local_to_global ( outside_rule.points, eval_points );
		my_fun.global_evaluate ( box1, eval_points, temp_val );
		temp_val *= outside_rule.weights;
		double    outside_integral = temp_val.sum();

		double    sum_partition    = inside_integral+outside_integral;

		cout.precision ( 15 );
		cout << scientific;
		cout << "Partitioned integral : " << sum_partition << " (inside " << inside_integral << ", outside "
			<< outside_integral << ")\n";

		if ( circular_boundary ) {
			cout << "fabs from pi/4       : " << fabs( inside_integral-M_PI/4 )
				<< " (pi/4   " << M_PI/4.0
				<< ", inside  " << inside_integral << ")\n";
			cout << "Rel. err. from pi/4  : " << fabs( (4.0*inside_integral/M_PI)-1.0 ) << "\n";
		}
		cout << "Box integral         : " << box_integral << "\n";
		cout << "fabs of error        : " << scientific << fabs( sum_partition-box_integral ) << "\n";
		cout << "Rel. error           : " << fabs( sum_partition-box_integral ) / box_integral << "\n";
	}
	
	{
		valarray<double>    global_points;
		
		ofstream    file ( "quadrature_rule-3_draw.gnuplot" );
		file << "#!/usr/bin/gnuplot\n";
		file << "set size square\n";
		file << "plot [-1:2] [-1:2] '-' w l t \"box\"";
		file << ", '-' w l t \"boundary\", '-' w dots t \"inside\"";
		file << ", '-' w dots t \"outside\"";
		file << ", '-' w p pt 3 t \"boundary integral\"";
		file << "\n";
		gp_draw_single ( box1, file );
		file << "e\n";
		int stmpnum = segs.segments.size();
		gp_draw_points<2> ( segs.segments[slice(0,stmpnum,1)],
		                    valarray<double>(0.0,stmpnum/2),
		                    file );
		file << "e\n";
		file.precision ( 15 );
		box1.map_local_to_global ( inside_rule.points, global_points );
		gp_draw_points<2> ( global_points, inside_rule.weights, file, num_quadpts );
		file << "e\n";
		box1.map_local_to_global ( outside_rule.points, global_points );
		gp_draw_points<2> ( global_points, outside_rule.weights, file, num_quadpts );
		file << "e\n";
		box1.map_local_to_global ( boundary_rule.points, global_points );
		gp_draw_points<2> ( global_points, boundary_rule.weights, file, num_quadpts );
		file << "e\n";
		
	}
}

