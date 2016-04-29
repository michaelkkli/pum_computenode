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

#include "line_segments.hh"
#include "quadrature_rule.hh"

#include <cassert>
#include <cmath>
#include <iostream>

using std::clog;
using std::sqrt;

#ifndef eF77
#define eF77(s) s##_
#endif

extern "C" void eF77(dstevd) ( char* jobz,
                    int*  n,
                    double* d,
                    double* e,
                    double* z,
                    int*    ldz,
                    double* work,
                    int*    lwork,
                    int*    iwork,
                    int*    liwork,
                    int*    info );

void generate_gauss_legendre_rule ( int N, simple_quadrature_rule<1> & rule )
{
	assert ( N>0 );
	if ( N<=0 ) {
		rule.points.resize(0);
		rule.weights.resize(0);
		return;
	}

	char jobz = 'V';

	// Debugging helped find missing '&' here. 2009-04-15 ML.
	valarray<double> &    d = rule.points;
	d.resize ( N, 0.0 );

	valarray<double>    e(0.5,N-1);
	// Notice the 0.5 initialization so only division remains.
	for ( int i=0; i<N-1; ++i ) {
		e[i] /= sqrt( 1. - 1./(4*(i+1)*(i+1)) ); // Notice the +1 to account for indices.
	}

	int n = N; // Allow taking of non-const pointer.

	int ldz = N;
	valarray<double>    z( ldz*N );

	int lwork = 1 + 4*N + N*N;
	valarray<double>    work( lwork );

	int liwork = 3 + 5*N;
	valarray<int>    iwork( liwork );

	int info;
	eF77(dstevd) ( &jobz, &n, &d[0], &e[0], &z[0], &ldz, &work[0], &lwork, &iwork[0], &liwork, &info );

	// valarray<double>    quad_weights(N);
	valarray<double> &    quad_weights = rule.weights;
	quad_weights.resize ( N );
	
	for ( int i=0; i<N; ++i ) {
		quad_weights[i] = z[i*N];
	}
	quad_weights *= quad_weights;
	quad_weights *= 2.0;
}

void generate_gauss_legendre_rule ( int N, simple_quadrature_rule<2> & out )
{
	
	simple_quadrature_rule<1> quad1D;
	generate_gauss_legendre_rule( N, quad1D );
	const int numpts = quad1D.weights.size();
	
	out.points.resize( 2*numpts*numpts );
	out.weights.resize( numpts*numpts );
	
	for ( int y=0; y<numpts; ++y ) {
		for ( int x=0; x<numpts; ++x ) {
			out.points[ 2*(y*numpts + x) ] = quad1D.points[ x ];
			out.points[ 2*(y*numpts + x)+1 ] = quad1D.points[ y ];
		}
	}
	for ( int y=0; y<numpts; ++y ) {
		for ( int x=0; x<numpts; ++x ) {
			out.weights[ y*numpts + x ] =
				quad1D.weights[ x ]*quad1D.weights[ y ];
		}
	}
	
#if 0 // Check if broken.
	simple_quadrature_rule<1>    tmp1;
	generate_gauss_legendre_rule( N, tmp1 );

	out.points.resize  ( 2*N*N );
	out.weights.resize ( N*N );
	
	for ( int i=0; i<N; ++i ) {
		out.points[ slice(2*N*i,  N,2) ] = tmp1.points;    // Notice no index here.
		out.points[ slice(2*N*i+1,N,2) ] = tmp1.points[i]; // Notice index present.
	}

	for ( int i=0; i<N; ++i ) {
		out.weights[ slice(N*i,N,1) ] = tmp1.weights * tmp1.weights[i]; // Notice no index on first, index present on second.
	}
#endif
}

void generate_gauss_legendre_rule_local ( simple_quadrature_rule<1>&   base_rule,
                                          const box<2> &                 box1,
                                          line_segments &                segs,
                                          int                            entry_segment,
                                          bool                           do_inside,
                                          simple_quadrature_rule<2> &    out_rule,
                                          simple_quadrature_rule<2> *    bdry_rule )
{
	const double *    box_ext = box1.get ();
	
	int num_segs = (segs.segments.size()/2) - 1;
	int num_pts_inside = 0;
#if 0
	simple_quadrature_rule<1>    base_rule;
	generate_gauss_legendre_rule ( num_quadpts, base_rule );
#endif
	
	int num_quadpts = base_rule.weights.size();
	
	out_rule.points.resize  ( 2*num_quadpts*num_quadpts, 0.0 );
	out_rule.weights.resize (   num_quadpts*num_quadpts, 0.0 );
	
	if ( bdry_rule ) {
		bdry_rule->points.resize  ( 2*num_quadpts, 0.0 );
		bdry_rule->weights.resize (   num_quadpts, 0.0 );
	}

	int num_sections = num_quadpts;
	
	enum parameterization { cartesian_coords, polar_coords };
	int param_subcase = -1;

	vector<int>   inside_corners;
	vector<int>   outside_corners;

	// Both used for boundary integral. Values to allow fool-proof debugging.
	double entry_pt_exit_pt[4] = { 2345.0, 6789.0, 2468.0, 778899.0 };
	
	split_by_line_segments ( box1,
	                         segs,
	                         entry_segment, num_pts_inside,
	                         &entry_pt_exit_pt[0],
	                         &entry_pt_exit_pt[2], inside_corners, outside_corners );

	int     in_c_size = inside_corners.size();
	int     out_c_size = outside_corners.size();
	bool    sections_intersect_boundary = false;

	// cout << "in_c_size is " << in_c_size << ", out_c_size is " << out_c_size << "\n";

	parameterization    parameterization_to_use;

	// We assume the corners are given in an anticlockwise winding fashion.
	if ( (in_c_size==2) && (out_c_size==2) ) {
		parameterization_to_use = cartesian_coords;

		/*
			   2
			_______
			|     |
		3	|     |     1
			|_____|
			   0
		*/
		switch ( inside_corners[0] ) {
			case 0:
				param_subcase = do_inside ? 0 : 2;
				if ( bdry_rule ) {
					bdry_rule->points[slice(0,num_sections,2)] = base_rule.points;
				}
				break;
			case 1:
				param_subcase = do_inside ? 1 : 3;
				if ( bdry_rule ) {
					bdry_rule->points[slice(1,num_sections,2)] = base_rule.points;
				}
				break;
			case 2:
				param_subcase = do_inside ? 2 : 0;
				if ( bdry_rule ) {
					bdry_rule->points[slice(0,num_sections,2)] = base_rule.points;
				}
				break;
			case 3:
				param_subcase = do_inside ? 3 : 1;
				if ( bdry_rule ) {
					bdry_rule->points[slice(1,num_sections,2)] = base_rule.points;
				}
				break;
			default:
				clog << "Unsupported case.\n"; abort(); break;
		}
		sections_intersect_boundary = true;
	} else if ( (in_c_size==1) && (out_c_size==3) ) {

		/*
		      3 ______  2
		        |     |
		        |     |
		        |_____|
		      0         1
		*/
		
		parameterization_to_use     = polar_coords;
		param_subcase               = do_inside ? inside_corners[0] : outside_corners[1];
		sections_intersect_boundary = do_inside ? true : false;
	} else if ( (in_c_size==3) && (out_c_size==1) ) {
		parameterization_to_use     = polar_coords;
		param_subcase               = do_inside ? inside_corners[1] : outside_corners[0];
		sections_intersect_boundary = do_inside ? false : true;
	} else {
		clog << "Boundary integral case unsupported.\n"; abort();
	}
	
	if ( bdry_rule && (polar_coords == parameterization_to_use) && !sections_intersect_boundary ) {
		// The boundary integral scheme does not match up exactly
		// with the interior integral scheme so we must produce
		// custom sections in a more limited range of angles.
		
		double    pole_pt[2];
		
		/*
		      3 ______  2
		        |     |
		        |     |
		        |_____|
		      0         1
		*/
		
		switch ( param_subcase ) {
		case 0:
			pole_pt[0] = box_ext[0];    pole_pt[1] = box_ext[2];
			break;
		case 1:
			pole_pt[0] = box_ext[1];    pole_pt[1] = box_ext[2];
			break;
		case 2:
			pole_pt[0] = box_ext[1];    pole_pt[1] = box_ext[3];
			break;
		case 3:
			pole_pt[0] = box_ext[0];    pole_pt[1] = box_ext[3];
			break;
		default:
			abort();
			break;
		}
		
		double pole_r_zero, pole_theta_zero;
		
		polar_line_params ( entry_pt_exit_pt, &pole_r_zero, &pole_theta_zero );

		double rebase_entry_exit[4];
		rebase_entry_exit[0] = entry_pt_exit_pt[0] - pole_pt[0];
		rebase_entry_exit[1] = entry_pt_exit_pt[1] - pole_pt[1];
		rebase_entry_exit[2] = entry_pt_exit_pt[2] - pole_pt[0];
		rebase_entry_exit[3] = entry_pt_exit_pt[3] - pole_pt[1];
		
		double theta_one = atan2 ( rebase_entry_exit[0], rebase_entry_exit[1] );
		double theta_two = atan2 ( rebase_entry_exit[2], rebase_entry_exit[3] );
		
		double    po4       = M_PI/4.0;
		double    po2       = M_PI/2.0;
		double    p3o2      = 3.0*M_PI/2.0;
		double    twopi     = 2.0*M_PI;
		
		if ( theta_one < 0.0 ) {
			theta_one += twopi;
		}
		
		if ( theta_two < 0.0 ) {
			theta_two += twopi;
		}
		
		clog << "theta one and two " << theta_one << ", " << theta_two << "\n";

		double pi_d2 = M_PI/2.0;
		
		// Get correct range for use of tan.
		while ( theta_one > po2 ) {
			theta_one -= po2;
		}
		while ( theta_two > po2 ) {
			theta_two -= po2;
		}
		clog << "correct for tan theta one and two " << theta_one << ", " << theta_two << "\n";
		
		double theta_plus_d2  = 0.5*(theta_one + theta_two);
		double theta_minus_d2 = 0.5*(theta_two - theta_one);
		
		// Sweep from entry to exit without swapping. Make sure this matches up with weights.
		valarray<double>    local_line  (4);
		valarray<double>    global_line (4);
		
		for ( int i=0; i<num_sections; ++i ) {
			// double    angle     = 0.25*M_PI*(1.0+base_rule.points[i]);
			double angle = theta_plus_d2 + base_rule.points[i]*theta_minus_d2;
			
			double    tan_angle = tan(angle);
			// Set the origin of the ray/section in 0,1.
			// Set extremal point of ray/section in 2,3.
			switch ( param_subcase ) {
			case 0:
				local_line[0] = -1.0;
				local_line[1] = -1.0;
				if ( angle < po4 ) {
					local_line[2] =  1.0;
					local_line[3] = -1.0 + 2.0*tan_angle;
				} else {
					local_line[2] = -1.0 + 2.0/tan_angle;
					local_line[3] =  1.0;
				}
				break;
			case 1:
				local_line[0] =  1.0;
				local_line[1] = -1.0;
				if ( angle < po4 ) {
					local_line[2] =  1.0 - 2.0*tan_angle;
					local_line[3] =  1.0;
				} else {
					local_line[2] = -1.0;
					local_line[3] = -1.0 + 2.0/tan_angle;
				}
				break;
			case 2:
				local_line[0] =  1.0;
				local_line[1] =  1.0;
				if ( angle < po4 ) {
					local_line[2] = -1.0;
					local_line[3] =  1.0 - 2.0*tan_angle;
				} else {
					local_line[2] =  1.0 - 2.0/tan_angle;
					local_line[3] = -1.0;
				}
				break;
			case 3:
				local_line[0] = -1.0;
				local_line[1] =  1.0;
				if ( angle < po4 ) {
					local_line[2] =  1.0;
					local_line[3] =  1.0 - 2.0*tan_angle;
				} else {
					local_line[2] = -1.0 + 2.0/tan_angle;
					local_line[3] = -1.0;
				}
				break;
			default:
				clog << "Invalid case.\n"; abort();
				break;
			}
			
			box1.map_local_to_global ( local_line, global_line );
			
			double param1 = -3.0, param2 = -3.0;
			double pt[2];
			
			int num_segs_to_check = 1 + num_pts_inside;
			
			bool cross = false;
			
			int current_seg = entry_segment;
			
			while ( (!cross) && (num_segs_to_check>0) ) {
				cross = intersect_lines_2d ( &global_line[0], &segs.segments[2*current_seg], param1, param2, pt );
				
				++current_seg;
				
				if ( current_seg > num_segs-1 ) {
					current_seg -= num_segs;
				}
				
				--num_segs_to_check;
			}
			
			if ( (!cross) && (num_segs_to_check==0) ) {
				clog << "Problem: unexpected situation.\n"; abort();
			}
			
			global_line[2] = pt[0];
			global_line[3] = pt[1];
			box1.map_global_to_local ( &global_line[2], &local_line[2] );
			
			bdry_rule->points[ 2*i     ] = local_line[2];
			bdry_rule->points[ 2*i + 1 ] = local_line[3];
		}
	}

	
	double final_quad_factor;

	if ( cartesian_coords==parameterization_to_use ) {
		switch (param_subcase) {
			case 0:
				// Allow drop-through. No break here.
			case 2:
			final_quad_factor = (box_ext[1]-box_ext[0])/2; // Debugging found incorrect 3 and 2. 2009-06-26 ML.
				break;
			case 1:
				// Allow drop-through. No break here.
			case 3:
			final_quad_factor = (box_ext[3]-box_ext[2])/2; // Debugging found incorrect 1 and 0. 2009-06-26 ML.
				break;
			default:
				clog << "Invalid case.\n"; abort();
				break;
		}
	} else {
		// Assuming polar coordinates parameterization here.

		final_quad_factor = M_PI/4.0; // Half a corner angle of a square.
	}

#if 0 // Delete. 2009-06-26 ML.
	if ( cartesian_coords==parameterization_to_use ) {
		clog << "Debug: Param to use is cartesian_coords and final_quad_factor is "
			<<final_quad_factor << "\n";
	} else {
		clog << "Debug: Param to use is polar_coords and final_quad_factor is "
			<<final_quad_factor << "\n";
	}
#endif
	
	valarray<double>    local_line (4);
	for ( int i=0; i<num_sections; ++i ) {
		if ( cartesian_coords==parameterization_to_use ) {
			switch ( param_subcase ) {
				case 0:
					// Allow drop-through: no break here.
				case 2:
					// Vertical line.
					local_line[0] = base_rule.points[i];
					local_line[1] = -1.0;
					local_line[2] = base_rule.points[i];
					local_line[3] = 1.0;
					break;
				case 1:
					// Allow drop-through: no break here.
				case 3:
					// Horizontal line.
					local_line[0] = -1.0;
					local_line[1] = base_rule.points[i];
					local_line[2] = 1.0;
					local_line[3] = base_rule.points[i];
					break;
				default:
					clog << "Invalid case.\n"; abort();
			}
		} else {
			// Assuming polar coordinates parameterization here.

			// Debugging helped find we were taking values in (0,pi)
			// rather than (0,pi/2). 2009-04-18 ML.
			double    angle     = 0.25*M_PI*(1.0+base_rule.points[i]);
			double    po4       = M_PI/4;

			double    tan_angle = tan(angle);
			// Set the origin of the ray/section in 0,1.
			// Set extremal point of ray/section in 2,3.
			switch ( param_subcase ) {
				case 0:
					local_line[0] = -1.0;
					local_line[1] = -1.0;
					if ( angle < po4 ) {
						local_line[2] =  1.0;
						local_line[3] = -1.0 + 2.0*tan_angle;
					} else {
						local_line[2] = -1.0 + 2.0/tan_angle;
						local_line[3] =  1.0;
					}
					break;
				case 1:
					local_line[0] =  1.0;
					local_line[1] = -1.0;
					if ( angle < po4 ) {
						local_line[2] =  1.0 - 2.0*tan_angle;
						local_line[3] =  1.0;
					} else {
						local_line[2] = -1.0;
						local_line[3] = -1.0 + 2.0/tan_angle;
					}
					break;
				case 2:
					local_line[0] =  1.0;
					local_line[1] =  1.0;
					if ( angle < po4 ) {
						local_line[2] = -1.0;
						local_line[3] =  1.0 - 2.0*tan_angle;
					} else {
						local_line[2] =  1.0 - 2.0/tan_angle;
						local_line[3] = -1.0;
					}
					break;
				case 3:
					local_line[0] = -1.0;
					local_line[1] =  1.0;
					if ( angle < po4 ) {
						local_line[2] =  1.0;
						local_line[3] =  1.0 - 2.0*tan_angle;
					} else {
						local_line[2] = -1.0 + 2.0/tan_angle;
						local_line[3] = -1.0;
					}
					break;
				default:
					clog << "Invalid case.\n"; abort();
					break;
			}
		}
		valarray<double>    global_line (4);

		box1.map_local_to_global ( local_line, global_line );

		double param1 = -3.0, param2 = -3.0;
		double pt[2];

		int num_segs_to_check = 1 + num_pts_inside;

		bool cross = false;

		int current_seg = entry_segment;
		
		while ( (!cross) && (num_segs_to_check>0) ) {
			// Debugging helped find that 2*current_seg was incorrectly missing the 2 multiplier.
			// Strangely it worked in main but showed up as incorrect when function was made stand-alone.
			// 2009-04-19 ML.
			cross = intersect_lines_2d ( &global_line[0], &segs.segments[2*current_seg], param1, param2, pt );

			++current_seg;
			// Debugging helped find missing -1
			// so was unable to wrap around to zero.
			// 2009-06-25 ML.
			if ( current_seg > num_segs-1 ) {
				current_seg -= num_segs;
			}
			
			--num_segs_to_check;
		}

		if ( (!cross) && (num_segs_to_check==0) ) {
			//cout << "No crossing found and all relevant segments checked.\n";
			if ( sections_intersect_boundary==true ) {
				clog << "Problem: unexpected situation.\n"; abort();
			}
		}

		if ( cartesian_coords==parameterization_to_use ) {
		/*
			   2
			_______
			|     |
		3	|     |     1
			|_____|
			   0
		*/
			switch ( param_subcase ) {
				case 0:
					out_rule.weights[ slice(num_quadpts*i, num_quadpts, 1) ] =
						base_rule.weights*(base_rule.weights[i]*(pt[1] - box_ext[2])/2);
					local_line[3] = param1;
					if ( bdry_rule ) {
						bdry_rule->points[ 2*i + 1 ] = param1;
					}
					break;
				case 1:
					// Notice no indexing on second weights only.
					out_rule.weights[ slice(num_quadpts*i, num_quadpts, 1) ] =
						base_rule.weights*(base_rule.weights[i]*(box_ext[1] - pt[0])/2);
					local_line[0] = param1;
					if ( bdry_rule ) {
						bdry_rule->points[ 2*i     ] = param1;
					}
					break;
				case 2:
					out_rule.weights[ slice(num_quadpts*i, num_quadpts, 1) ] =
						base_rule.weights*(base_rule.weights[i]*(box_ext[3] - pt[1])/2);
					local_line[1] = param1;
					if ( bdry_rule ) {
						bdry_rule->points[ 2*i + 1 ] = param1;
					}
					break;
				case 3:
					// Notice no indexing on second weights only.
					out_rule.weights[ slice(num_quadpts*i, num_quadpts, 1) ] =
						base_rule.weights*(base_rule.weights[i]*(pt[0] - box_ext[0])/2);
					local_line[2] = param1;
					if ( bdry_rule ) {
						bdry_rule->points[ 2*i     ] = param1;
					}
					break;
				default:
					clog << "Invalid case.\n"; abort();
					break;
			}
		} else { // Polar coordinates.
			if ( cross ) {
				global_line[2] = pt[0];
				global_line[3] = pt[1];
				box1.map_global_to_local ( &global_line[2], &local_line[2] );
			} else {
				box1.map_local_to_global ( &local_line[2],  &global_line[2] );
			}
			
			if ( bdry_rule && sections_intersect_boundary ) {
				bdry_rule->points[ 2*i     ] = local_line[2];
				bdry_rule->points[ 2*i + 1 ] = local_line[3];
			}

			double calc[2];
			calc[0] = global_line[0]-global_line[2];
			calc[1] = global_line[1]-global_line[3];

			double    r_dist_squared = calc[0]*calc[0] + calc[1]*calc[1];

			slice    section_slice(num_quadpts*i, num_quadpts, 1);
			out_rule.weights[ section_slice ] =  base_rule.points + 1.0;
			out_rule.weights[ section_slice ] *= base_rule.weights * base_rule.weights[i]*(r_dist_squared/4.0);
		}
		valarray<double>    local_quad_points;
		map_line_parameters_to_points<2> ( base_rule.points, &local_line[0], local_quad_points );
		assert ( local_quad_points.size() == 2*num_quadpts );
		out_rule.points[ slice(2*num_quadpts*i, 2*num_quadpts, 1) ] = local_quad_points;
	}
	out_rule.weights *= final_quad_factor;
}

template<>
quadrature_rule<1>::quadrature_rule( quadrature_rule_type qt )
{
  // FIXME: Values from wikipedia so far.
  if ( gauss1 == qt ) {
	  quadrature_points.resize( 1 );
	  quadrature_weights.resize( 1 );
    quadrature_points[0]  = 0.0;

    quadrature_weights[0] = 2.0;
  } else if ( gauss2 == qt ) {
	  quadrature_points.resize( 2 );
	  quadrature_weights.resize( 2 );
	  
    double tmp = sqrt( 1.0/3.0 );
    quadrature_points[0]  = -tmp;
    quadrature_points[1]  = tmp;

    quadrature_weights[0] = 1.0;
    quadrature_weights[1] = 1.0;
  } else if ( gauss3 == qt ) {
	  quadrature_points.resize( 3 );
	  quadrature_weights.resize( 3 );
	  
    double tmp = sqrt( 3.0/5.0 );
    quadrature_points[0] = -tmp;
    quadrature_points[1] = 0.0;
    quadrature_points[2] = tmp;

    tmp = 5.0/9.0;
    quadrature_weights[0] = tmp;
    quadrature_weights[1] = 8.0/9.0;
    quadrature_weights[2] = tmp;
  } else if ( gauss4 == qt ) {
	  quadrature_points.resize( 4 );
	  quadrature_weights.resize( 4 );
	  
	  double common_pt = 2.0*sqrt(6.0/5.0);
	  double smaller_pt = sqrt( (3.0 - common_pt)/7.0 );
	  double larger_pt  = sqrt( (3.0 + common_pt)/7.0 );
	  double common_wg = sqrt( 30.0 );
	  double smaller_wg = (18.0 - common_wg)/36.0;
	  double larger_wg = (18.0 + common_wg)/36.0;
	  
	  quadrature_points[0] = -larger_pt;
	  quadrature_points[1] = -smaller_pt;
	  quadrature_points[2] =  smaller_pt;
	  quadrature_points[3] =  larger_pt;
	  
	  quadrature_weights[0] = smaller_wg;
	  quadrature_weights[1] = larger_wg;
	  quadrature_weights[2] = larger_wg;
	  quadrature_weights[3] = smaller_wg;
  } else if ( gauss5 == qt ) {
	  // Points.
	  {
		  quadrature_points.resize ( 5 );
		  double inner  = 2.0 * sqrt ( 10.0/7.0 );
		  double lesser = sqrt ( 5.0 - inner ) / 3.0;
		  double greater = sqrt ( 5.0 + inner ) / 3.0;
		  
		  quadrature_points[0] = -greater;
		  quadrature_points[1] = -lesser;
		  quadrature_points[2] =  0.0;
		  quadrature_points[3] =  lesser;
		  quadrature_points[4] =  greater;
	  }
	  // Weights.
	  {
		  quadrature_weights.resize ( 5 );
		  double tmp = 13.0 * sqrt ( 70.0 );
		  double smaller = ( 322.0 - tmp )/900.0;
		  double larger  = ( 322.0 + tmp )/900.0;

		  // Weights were found to be wrong. Smaller/larger mix-up. 2009-01-12 ML.
		  quadrature_weights[0] = smaller;
		  quadrature_weights[1] = larger;
		  quadrature_weights[2] = 128.0/225.0;
		  quadrature_weights[3] = larger;
		  quadrature_weights[4] = smaller;
	  }
  }
}

template<>
quadrature_rule<2>::quadrature_rule ( quadrature_rule_type qt )
{

  quadrature_rule<1> quad1D(qt);
  const int numpts = quad1D.num_points();
  
	quadrature_points.resize( 2*numpts*numpts );
	quadrature_weights.resize( numpts*numpts );
  
  for ( int y=0; y<numpts; ++y ) {
    for ( int x=0; x<numpts; ++x ) {
      quadrature_points[ 2*(y*numpts + x) ] = quad1D.quadrature_points[ x ];
      quadrature_points[ 2*(y*numpts + x)+1 ] = quad1D.quadrature_points[ y ];
    }
  }
  for ( int y=0; y<numpts; ++y ) {
    for ( int x=0; x<numpts; ++x ) {
      quadrature_weights[ y*numpts + x ] =
	quad1D.quadrature_weights[ x ]*quad1D.quadrature_weights[ y ];
    }
  }
  assert( quadrature_points.size()%2 == 0 );
}

template<>
quadrature_rule<3>::quadrature_rule ( quadrature_rule_type qt )
{
  quadrature_rule<1> quad1D(qt);
  const int numpts = quad1D.num_points();
	
	quadrature_points.resize( 3*numpts*numpts );
	quadrature_weights.resize( numpts*numpts );
	
  for ( int z=0; z<numpts; ++z ) {
    for ( int y=0; y<numpts; ++y ) {
      for ( int x=0; x<numpts; ++x ) {
        quadrature_points[ 3*(z*numpts*numpts + y*numpts + x)   ] = quad1D.quadrature_points[ x ];
        quadrature_points[ 3*(z*numpts*numpts + y*numpts + x)+1 ] = quad1D.quadrature_points[ y ];
	quadrature_points[ 3*(z*numpts*numpts + y*numpts + x)+2 ] = quad1D.quadrature_points[ z ];
      }
    }
    for ( int y=0; y<numpts; ++y ) {
      for ( int x=0; x<numpts; ++x ) {
        quadrature_weights[ z*numpts*numpts + y*numpts + x ] =
	   quad1D.quadrature_weights[ x ]*quad1D.quadrature_weights[ y ]*quad1D.quadrature_weights[ z ];
      }
    }
  }
}

template<int dim>
quadrature_rule<dim>::~quadrature_rule()
{
}

template<int dim>
int
quadrature_rule<dim>::num_points() const
{
  return quadrature_weights.size();
}
template<int dim>
valarray<double>&
quadrature_rule<dim>::access_quadrature_points()
{
	return quadrature_points;
}

template<int dim>
valarray<double>&
quadrature_rule<dim>::access_quadrature_weights()
{
	return quadrature_weights;
}

template<int dim>
double
quadrature_rule<dim>::integrate ( const box<dim>& bx,
				  const function<dim>& fn,
				  valarray<double>* in_va_ptr )
{
  valarray<double>* va_ptr = 0;
  valarray<double> scoped_va;
  if ( 0 == in_va_ptr ) {
    va_ptr = &scoped_va;
  } else {
    va_ptr = in_va_ptr;
  }

  fn.local_evaluate ( bx, quadrature_points, *va_ptr );

  (*va_ptr) *= quadrature_weights;

  return va_ptr->sum()*bx.measure()*std::pow(0.5,dim);

}

template<int dim>
double
quadrature_rule<dim>::integrate_product ( const box<dim>&      rl,
		    const box<dim>&      b0,
		    const function<dim>& f0,
		    const box<dim>&      b1,
		    const function<dim>& f1 )
{
  assert( b0.closed_intersect_box( rl ) && b1.closed_intersect_box( rl ) );
  valarray<double> tmp0;
  valarray<double> tmp1;
  f0.restricted_local_evaluate( rl, b0, this->quadrature_points, tmp0 );
  f1.restricted_local_evaluate( rl, b1, this->quadrature_points, tmp1 );

  assert( tmp0.size() == this->num_points() );
  assert( tmp1.size() == this->num_points() );

  tmp0 *= tmp1;
  tmp0 *= this->quadrature_weights;
  return tmp0.sum()*rl.measure()*std::pow(0.5,dim);
}

//

template class quadrature_rule<1>;
template class quadrature_rule<2>;
template class quadrature_rule<3>;

template struct simple_quadrature_rule<1>;
template struct simple_quadrature_rule<2>;
