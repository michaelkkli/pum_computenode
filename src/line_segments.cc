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
#include "box.hh"
#include "hacker_shimrat_algorithm112.hh"
#include "line_segments.hh"
#include "refinement_structure.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>

using std::atan2;
using std::back_inserter;
using std::cerr;
using std::clog;
using std::copy;
using std::cout;
using std::fabs;
using std::fill_n;
using std::reverse;
using std::rotate;
using std::sort;
using std::slice;
using std::swap;


void get_first_entry_num_inside ( line_segments& ls,
                                  box<2>& bx,
                                  int start_hint,
                                  int max_points_inside,
                                  int& actual_first_entry,
                                  int& num_found_inside )
{
	actual_first_entry = -1;
	num_found_inside  =  0;
	
	assert ( ls.segments.size()%2 == 0 );
	int num_segs = (ls.segments.size()/2)-1;
	
	if ( (start_hint<0) || (start_hint>num_segs-1) ) {
		start_hint = 0;
	}
	
	if ( (max_points_inside<0) || (max_points_inside>num_segs) ) {
		max_points_inside = num_segs;
	}
	
	// Make sure we start from a point and line segment that lies outside the box.
	// We don't expect to enter here if the start_hint is good.
	{
		int avoid_infinite_loop = 0;
		while ( bx.closed_intersect_point ( &ls.segments[2*start_hint] ) ) {
			if ( start_hint == 0 ) {
				start_hint = num_segs-1;
			} else {
				--start_hint;
			}
			
			++avoid_infinite_loop;
			if ( avoid_infinite_loop > num_segs ) {
				abort();
			}
		}
	}
	
	int checking = start_hint;
	
	{
		pair<int,int>    intersections; intersections.first = -1; intersections.second = -1;
	
		int num_tried = 0;
	
		// Keep advancing through line segments until we find one that intersects the box.
		while (1) {
			box_2d_intersect_line ( bx, &ls.segments[2*checking], 0, 0, intersections );
			
			if ( intersections.first != -1 ) {
				// Found first entry.
				actual_first_entry = checking;
				break;
			} else {
				++num_tried;
			
				// Advance to the next line segment.
				if ( checking == num_segs-1 ) {
					checking = 0;
				} else {
					++checking;
				}
			}
			
			if ( num_tried > max_points_inside ) {
				// No line segment intersects any edge of the box
				// within the specified 
				actual_first_entry = -1;
				num_found_inside  =  0;
				return;
			}
		}
	}
	
	// We have verified that the line segment with index `checking' does
	// intersect the box. We now find how many points lie inside the box.
	{
		int num_tried = 0;
		while (1) {
			if ( checking == num_segs-1 ) {
				checking = 0;
			} else {
				++checking;
			}
			
			if ( bx.closed_intersect_point ( &ls.segments[2*checking] ) ) {
				++num_found_inside;
			}
			
			++num_tried;
			
			if ( (num_tried>max_points_inside) || (num_found_inside > max_points_inside) ) {
				// Return if we find more points lying inside than there are segments.
				// This will never happen as we would have aborted further up when no point
				// is found to lie outside but this is a useful check to ensure we do not
				// have an infinite while-loop. 
				return;
			}
		}
	}
	
}

template <int dim>
void get_bounding_box ( line_segments& ls, box<dim>& bx )
{
	valarray<double>    segs = ls.segments;

	assert ( segs.size()%dim == 0 );
	int num_pts = (segs.size()-dim)/dim;
	assert ( num_pts > 0 );
	
	double ext[2*dim];
	
	for ( int d=0; d<dim; ++d ) {
		const valarray<double> & tmp = segs[slice(d,num_pts,dim)];
		ext[2*d]   = tmp.min();
		ext[2*d+1] = tmp.max();
	}
#if 0
	for ( int i=0; i<num_pts; ++i ) {
		if ( 0==i ) {
			for ( int d=0; d<dim; ++d ) {
				ext[2*d]   = segs[d];
				ext[2*d+1] = segs[d];
			}
		}
		
		for ( int d=0; d<dim; ++d ) {
			if ( segs[i*dim + d] < ext[2*d] ) {
				ext[2*d] = segs[i*dim+d];
			} else if ( ext[2*d+1] < segs[i*dim+d] ) {
				ext[2*d+1] = segs[i*dim+d];
			}
		}
	}
#endif
	bx.set ( ext );
}

void regularize_bounding_box ( box<2>& bx ) {
	const double* ext = static_cast<const box<2>&>(bx).get();
	
	double largest_length = ext[1]-ext[0];
	double tmp;
	
	for ( int d=1; d<2; ++d ) {
		tmp = ext[2*d+1] - ext[2*d];
		if ( tmp > largest_length ) {
			largest_length = tmp;
		}
	}
	
	double    centre[2*2];
	
	bx.get_centre_point( centre );
	
	make_box_with_radius ( centre, largest_length/2, bx );
}

/**
	Let A and B be the respective start and end points for the first line.
	Let C and D be the respective start and end points for the second line.

	We want to find parameters r, s in [-1, 1] such that
		A(1-r) + B(1+r) = C(1-s) + D(1+s)
	so
		A + B + r(B-A) = C + D + s(D-C)
		(B-A)r + (C-D)s = C + D - A - B

	[x_b - x_a	x_c - x_d] [r] = [ x_c + x_d - x_a - x_b ]
	[y_b - y_a	y_c - y_d] [s]   [ y_c + y_d - y_a - y_b ]

	[ M N ][r]=[U]
	[ O P ][s] [V]

	det = MP - ON

	[r] =  (1./det)[  P -N ][U]
	[s]            [ -O  M ][V]

	r = (PU-NV)/det
	s = (VM-OU)/det
*/
bool intersect_lines_2d ( const double* first,
                          const double* second,
                          double& param1,
                          double& param2,
                          double* inter_pt )
{
	assert ( first && second );

	const double* A = first;
	const double* B = first+2;
	const double* C = second;
	const double* D = second+2;

	double M = B[0] - A[0]; double N = C[0] - D[0];
	double O = B[1] - A[1]; double P = C[1] - D[1];

	double U = C[0] + D[0] - A[0] - B[0];
	double V = C[1] + D[1] - A[1] - B[1];

	double det = M*P - O*N;

	if ( fabs(det) < 1e-20 ) {
		param1 = param2 = -2.0;
		return false;
	}

	param1 = (P*U-N*V)/det;
	
	// Debugging found incorrect 1.0 in place of -1.0. 2008-12-01 Michael LI.
	if ( (param1 < -1.0-1e-10) || (param1 > 1.0+1e-10) ) {
		param1 = param2 = -2.0;
		return false;
	}
	
	param2 = (V*M-O*U)/det;
	
	// Debugging found incorrect 1.0 in place of -1.0. 2008-12-01 Michael LI.
	if ( (param2 < -1.0-1e-10) || (param2 > 1.0+1e-10) ) {
		param1 = param2 = -2.0;
		return false;
	}

	if ( inter_pt ) {
		inter_pt[0] = 0.5*( A[0] + B[0] + param1*(B[0]-A[0]) );
		inter_pt[1] = 0.5*( A[1] + B[1] + param1*(B[1]-A[1]) );
	}
	// Debugging found this return statement was missing. 2008-12-01 Michael Li.
	return true;
}

void box_2d_intersect_line ( const box<2>&    bx,
                            const double*    line,
                            double*          first_intrsct_pt,
                            double*          second_intrsct_pt,
                            pair<int,int>&   ordered_intersections )
{
	assert ( line );
	const double* ext = bx.get();
	double box_line[4];
	
	// Indicate no intersection.
	ordered_intersections.first  = -1;
	ordered_intersections.second = -1;
	
	// Values chosen to ensure they will be replaced.
	double entry_param =  2.0;
	double exit_param  = -2.0;
	
	double line_param, ignore_param;
	double tmp_intrsct[2];
	
	int num_intersects = 0;
	
	assert ( line );
	bool check;
	for ( int i=0; i<4; ++i ) {
		switch ( i ) {
			// x_start	x_end
			// y_start	y_end
		case 0:
			// Test bottom of box (y same - lowest).
			box_line[0] = ext[0];	box_line[2] = ext[1];
			box_line[1] = ext[2];   box_line[3] = ext[2];
			break;
		case 2:
			// Test top of box (y same - highest).
			box_line[0] = ext[1];	box_line[2] = ext[0];
			box_line[1] = ext[3];   box_line[3] = ext[3];
			break;
		case 3:
			// Test left edge (x same - lowest).
			box_line[0] = ext[0];	box_line[2] = ext[0];
			box_line[1] = ext[3];   box_line[3] = ext[2];
			break;
		case 1:
			// Test right edge ( x same - highest).
			box_line[0] = ext[1];	box_line[2] = ext[1];
			box_line[1] = ext[2];   box_line[3] = ext[3];
			break;
		default:
			abort();
			break;
		}
		
		// Debugging found incorrect use of the param arguments.
		// They were mistakenly thought to be entry exit params of
		// some sort. 2008-12-01 Michael Li.
		check = intersect_lines_2d ( line,
		                             &box_line[0],
		                             line_param,
		                             ignore_param,
		                             &tmp_intrsct[0] );
		if ( check ) {
			++num_intersects;
			if ( line_param < entry_param ) {
				entry_param = line_param;
				if ( first_intrsct_pt ) {
					first_intrsct_pt[0] = tmp_intrsct[0];
					first_intrsct_pt[1] = tmp_intrsct[1];
				}
				ordered_intersections.first = i;
			}
			if ( line_param > exit_param ) {
				// All this must be done because the
				// first seen intersection may also be the
				// exit intersection.
				// Equality of the pair denotes single intersection.
				exit_param = line_param;
				if ( second_intrsct_pt ) {
					second_intrsct_pt[0] = tmp_intrsct[0];
					second_intrsct_pt[1] = tmp_intrsct[1];
				}
				ordered_intersections.second = i;
			}
		}
	}
#if 0 // ndef NDEBUG // Not a valid test as line may be inside box. 2008-12-01 Michael Li.
	{
		box_intersect_line_scratch    bils;
		if ( bx.open_intersect_line ( line, bils ) ) {
			assert ( ordered_intersections.first != -1 );
		}
	}
#endif
}

template <int dim>
void gnuplot_output ( line_segments &    ls,
                      ostream &          out,
                      bool               with_labels,
                      int                max_num_labels )
{
	valarray<double> &    segments = ls.segments;
	assert ( segments.size()%dim == 0 );
	int num_segs = (segments.size()/dim)-1;
	assert(num_segs>0);
	
	out << "#!/usr/bin/gnuplot\n";
	
	int    num_to_skip = 1; 
	if ( num_segs > max_num_labels ) {
		num_to_skip = floor(0.5+(num_segs/max_num_labels));
	}
	// GDB helped to find case when num_to_skip was zero.
	// Previously coded without if-statement and the floor was
	// always carried out. 2009-02-20 ML.
	assert(num_to_skip>0);
	
	if ( with_labels ) {
		for ( int i=0; i<num_segs; i+=num_to_skip ) {
			if ( i>num_segs-1 ) {
				break;
			}
			out << "set label \" "
			    << i << "\" at ";
			for ( int d=0; d<dim; ++d ) {
				if ( d!=0 ) {
					out << ",";
				}
				out << segments[dim*i+d];
			}
			out << "\n";
		}
	}
	
	{
		box<dim> bx;
		get_bounding_box ( ls, bx );
		bx.scale ( 1.2 );
		const box<dim>& silly = bx;
		const double * ext = silly.get();
	
		out << "plot [" << ext[0] << ":" << ext[1] << "]"
		    << " ["<< ext[2*1] << ":" << ext[2*1+1] <<"] '-' w l\n";
	}

	gp_draw<dim> ( ls, out );
	

}

template <int dim>
void gp_draw ( line_segments &    ls,
               ostream &          out )
{
	valarray<double> &    segments = ls.segments;
	assert ( segments.size()%dim == 0 );
	int num_segs = (segments.size()/dim)-1;
	assert(num_segs>0);

	for ( int i=0; i<num_segs; ++i ) {
		for ( int d=0; d<dim; ++d ) {
			if ( d!=0 ) {
				out << "\t";
			}
			out << segments[dim*i+d];
		}
		out << "\n";
	}
	for ( int d=0; d<dim; ++d ) {
		if ( d!=0 ) {
			out << "\t";
		}
		out << segments[dim*(num_segs) + d];
	}
	out << "\n";
}

template <int dim>
void get_crossing_parameters ( valarray<double> *    grid,
                               int                   cd,
                               const double *        line,
                               valarray<double> &    params )
{
	assert ( grid && line );
	
	// Start and end points of line projected onto cd-axis.
	// Either one of a and b may be the larger.
	double    a = line[ cd ];
	double    b = line[ dim + cd ];
	
	if ( fabs( b-a ) < 1e-20 ) {
		// We need to divide by (b-a) later so will have to
		// ignore cases where the projection onto the coordinate cd
		// is very small and hopefully does not cross too many grid lines.
		params.resize(0);
		return;
	}
	
	bool    swapped = false;
	
	if ( a > b ) {
		swap ( a, b );
		swapped = true;
	}
	
	valarray<double> &    ref = grid[cd];
	
	int    refsize = ref.size();
	
	vector<double>    inside;
	double    tmp;
	for ( int i=0; i<refsize; ++i ) {
		tmp = ref[i];
		if ( (a<tmp) && (tmp<b) ) {
			inside.push_back( tmp );
		}
	}
	
	params.resize ( inside.size() );
	// Copy in order to make use of valarray functionality in place.
	copy ( inside.begin(), inside.end(), &params[0] );
	
	// 2*x = (a+b) + t(b-a)
	// so
	// t = (2*x - (a+b)) / (b-a)
	
	// Form parameters in place. params hold 'x' before transformation.
	params *= 2.0;
	params -= (a + b);
	params /= (b - a);
	
	if ( swapped ) {
		// The parameters are symmetric about zero
		// with respect to direction of the line.
		params *= -1.0;
		// Reverse to have parameters increasing in size.
		reverse ( &params[0], &params[0] + params.size() );
	}
}

void replace_parameters_by_midpoints ( vector<double>& params )
{
	int size = params.size();
	assert ( size >= 2 );
	
	for ( int i=0; i<size-1; ++i ) {
		params[i] = 0.5*(params[i]+params[i+1]);
	}
	params.pop_back();
}

template <int dim, class T>
void map_line_parameters_to_points ( T &    params,
                                     const double *      line,
                                     T &    points,
                                     double *            Jacobian_determinant )
{
	assert ( line );
	// For each dimension, x = 0.5*( (a+b) + t(b-a) )
	
	int num_pts = params.size();
	points.resize ( num_pts * dim );
	
	double a, b, apb, bma;
	for ( int d=0; d<dim; ++d ) {
		a = line[ d ];
		b = line[ dim + d ];
		apb = a + b;
		bma = b - a;
		for ( int i=0; i<num_pts; ++i ) {
			points[ i*dim + d ] = 0.5*(apb + params[i]*bma);
		}
	}
	
	if ( Jacobian_determinant ) {
		double &    jd = *Jacobian_determinant;
		jd = 1.0;
		
		for ( int d=0; d<dim; ++d ) {
			jd *= fabs ( 0.5 * bma );
		}
	}
}

template <int dim>
void get_grid_intersection_midpoints ( valarray<double>*     grid,
                                       line_segments &       segs,
                                       valarray<double> &    out_points,
                                       valarray<int> &       associated_segments,
                                       int                   hint,
                                       bool                  finish_after_first_exit )
{
	assert ( grid );
	
	// Get a bounding box for the grid.
	double grid_ext[2*dim];
	for ( int i=0; i<dim; ++i ) {
		assert ( grid[i].size() >= 2 );
		grid_ext[ i*dim     ] = grid[i][0];
		grid_ext[ i*dim + 1 ] = grid[i][grid[i].size()-1];
	}
	
	box<dim>    bx;
	bx.set( grid_ext );
	
	box_intersect_line_scratch    bils;
	
	valarray<double> &    seg_pts = segs.segments;
	
	int numsegs = (seg_pts.size()/dim)-1;
	
	// Don't allow a hint beginning at a non-existent segment.
	assert ( hint < numsegs );
	int    first_entry = hint;
	
	while ( !bx.open_intersect_line ( &seg_pts[dim*first_entry], bils ) ) {
		abort();
		++first_entry;
		// Make sure at least one segment intersects the box.
		assert ( first_entry < numsegs );
		if ( first_entry >= numsegs ) {
		  out_points.resize(0);
		  associated_segments.resize(0);
			return;
		}
	}
	
	vector<double>    tmp_out_pts;
	vector<int>       tmp_asegs;
	vector<double>    tmp_pts;
	
	valarray<double>    tmp_params;
	vector<double>      seg_params;
	
	bool    reached_first_exit = false;
	
	for ( int s=first_entry; s<numsegs; ++s ) {
		seg_params.clear();
		
		for ( int d=0; d<dim; ++d ) {
			get_crossing_parameters<dim> ( grid, d, &seg_pts[dim*s], tmp_params );
			copy ( &tmp_params[0],
			       &tmp_params[0]+tmp_params.size(),
			       back_inserter( seg_params ) );
		}
		
		if ( bx.open_intersect_point( &seg_pts[dim*s] ) ) {
			seg_params.push_back(-1.0);
		}
		
		sort ( seg_params.begin(), seg_params.end() );
		
		// Check if the forward end point is still in the grid box.
		if ( bx.open_intersect_point( &seg_pts[dim*(s+1)] ) ) {
			seg_params.push_back(1.0);
		} else {
			reached_first_exit = true;
		}
		
		if ( seg_params.size() < 2 ) {
			continue;
		}
		replace_parameters_by_midpoints( seg_params );
		
		// Debugging found missing multiply by dim here. 2008-12-01 Michael Li.
		map_line_parameters_to_points<dim> ( seg_params, &seg_pts[dim*s], tmp_pts );
		
		assert ( seg_params.size() == tmp_pts.size()/dim );
		
		copy ( tmp_pts.begin(),
		       tmp_pts.end(),
		       back_inserter( tmp_out_pts ) );
		
		// Debugging found missing division by dim here. 2008-12-01 Michael Li.
		fill_n ( back_inserter(tmp_asegs), tmp_pts.size()/dim, s );
		
		if ( (reached_first_exit==true)&&(finish_after_first_exit==true) ) {
			break;
		}
	}
	
	out_points.resize ( tmp_out_pts.size() );
	copy ( tmp_out_pts.begin(), tmp_out_pts.end(), &out_points[0] );
	
	associated_segments.resize ( tmp_asegs.size() );
	copy ( tmp_asegs.begin(), tmp_asegs.end(), &associated_segments[0] );

	assert ( tmp_out_pts.size()/2 == tmp_asegs.size() );
	if ( tmp_out_pts.size()/2 != tmp_asegs.size() ) {
		cerr << "Size mismatch.\n";
		abort();
	}

#if 0 // Mark for delete. 2008-12-10 ML.
	cout << "associated_segments are ";
	for ( int i=0; i<associated_segments.size(); ++i ) {
		if ( i!= 0 ) {
			cout << ", ";
		}
		cout << associated_segments[i];
	}
#endif
}

template <int dim>
void get_keys_intersect_point ( valarray<double> &             pts,
                                valarray<int> &                asso_segs,
                                line_segments &                linesg,
                                refinement_structure<dim> &    rs,
                                vector<string> &               ordered,
                                vector<int> &                  out_ent )
{
	assert ( rs.dtree );
	
	map<string,vector<int> >    incident_segments;
	
	basic_dbinary_tree<dim>& dtree = *(rs.dtree);
	
	int numpts = pts.size()/dim;
	assert ( asso_segs.size() == numpts );
	
	string    key;
	
	typename map<string,vector<int> >::iterator    mit,mend;
	
	ordered.clear();
	out_ent.clear();
	
	for ( int i=0; i<numpts; ++i ) {
		get_key_intersect_point<dim> ( rs, &pts[i*dim], key );

		mit = incident_segments.find(key);
		
		if ( mit == incident_segments.end() ) {
			// The points should go around the boundary and
			// so the patch keys will be in order.
			ordered.push_back(key);
		}
		incident_segments[key].push_back(asso_segs[i]);
#if 0 // Mark for delete. 2008-12-08 ML.
		cout << "point associated with "
		     << asso_segs[i]
		    << " lies in key " << key << "\n";
#endif
	}
	
	int num_segs = (linesg.segments.size()/2)-1;
	
	mit = incident_segments.begin();
	mend = incident_segments.end();
	for ( ; mit!=mend; ++mit ) {
		vector<int> &    inds = mit->second;
		
		// Debugging found repeats which we want to remove.
		// 2009-02-12 ML.
		vector<int>::iterator   uit = unique ( inds.begin(), inds.end() );
		inds.resize( uit - inds.begin() );
		
		if ( inds.size()==num_segs ) {
			// We don't handled all the segments in a single box yet.
			abort();
		}
		if ( (inds[0]==0) && (inds.back()==num_segs-1) ) {
#if 0	
			cout << "Entering wrap around of assoc_segs for " << mit->first <<"\n";
			cout << "Incident segments are ";
			for ( int i=0; i<inds.size(); ++i ) {
				if ( 0!=i) {
					cout << ", ";
				}
				cout << inds[i];
			}
			cout << "\n";
#endif
			
			// Due to wrapping around of the boundary segments,
			// it may happen that the first entry segments is not
			// the lowest index.
			int size = inds.size();
			int head = 0;
			
			// Debugging found assignment rather than equality check
			// here. Horrible. Horrible.
			// 2008-12-08 Michael LI.

			while ( (head < size) && (inds[head] == head) ) {
				++head;
			}
			
			rotate ( inds.begin(), inds.begin()+head, inds.end() );
			
#if 0
			cout << "After rotation, incident segments are ";
			for ( int i=0; i<inds.size(); ++i ) {
				if ( 0!=i) {
					cout << ", ";
				}
				cout << inds[i];
			}
			cout << "\n";
#endif
		}
		
	}
	
	
	int size = ordered.size();
	
	// Debugging found this line was incorrect and
	// was resizing to size giving the symptom that
	// out_ent was double the size it should have been.
	// 2008-12-07 Michael LI.
	out_ent.reserve(size);
	
	box<2>    bx;
	pair<int, int>    edge_ind;
	
#if 0 // Mark for delete 2008-12-08 Michael LI.	
	for ( int i=0; i<size; ++i ) {
		key = ordered[i];
		assert ( incident_segments.find(key) != incident_segments.end() );
		vector<int> &   v_incident = incident_segments[ key ];
		
		cout << "key " << key << " has incident segments ";
		for ( int j=0; j<v_incident.size(); ++j ) {
			if ( j!=0 ) {
				cout << ", ";
			}
			cout << v_incident[j];
		}
		cout << "\n";
	}
#endif

//	abort();

	// Using a parent box expansion factor of 1.0 causes a trip of the assert in here.
	// It appears that the job of the outer for loop is simply to assign the first entry
	// segment and after the use of unique above and rotation, this should simply be the first
	// number in v_incident. 2009-02-13 ML.
	for ( int i=0; i<size; ++i ) {
		key = ordered[i];
		assert ( incident_segments.find(key) != incident_segments.end() );
		vector<int> &   v_incident = incident_segments[ key ];
		
		assert (v_incident.size()>0 );
		out_ent.push_back ( v_incident[0] );
	}
#if 0 // Delete if the above is verified to be correct.
	for ( int i=0; i<size; ++i ) {
		key = ordered[i];
		assert ( incident_segments.find(key) != incident_segments.end() );
		vector<int> &   v_incident = incident_segments[ key ];
		
		dtree.get_box ( key, bx );
		int inner_size = v_incident.size();
		
		// abort();
		// For-loop below is incorrect. Will not properly find first entry.
		
		for ( int j=0; j<inner_size; ++j ) {
			edge_ind.first=-1; edge_ind.second=-1;
			box_2d_intersect_line ( bx,
			                        &(linesg.segments[2*v_incident[j]]),
			                        0,
			                        0,
			                        edge_ind );
			if ( edge_ind.first != -1 ) {
				// Found a line segment that intersects.
				out_ent.push_back( v_incident[j] );
				// Debugging found this break was missing.
				// 2008-12-07 Michael LI.
				break;
			} else if ( j == inner_size-1 ) {
				// No incident line segment intersects an edge of the box.
				abort();
			}
		}
	}
#endif
	
	assert ( ordered.size() == out_ent.size() );
}


void split_by_line_segments ( const box<2> &           bx,
                              line_segments &    ls,
                              int                ent_seg,
                              int &              num_pts_inside,
                              double *           out_entry_pt,
                              double *           out_exit_pt,
                              vector<int> &    inside_corners,
                              vector<int> &    outside_corners )
{

	double entry_pt[2] = {-10., -10.};
	double exit_pt[2]  = {-10., -10.};

	valarray<double> &    segs = ls.segments;
	int num_segs = (segs.size()/2) - 1;

	pair<int,int>    edge_ind;
	edge_ind.first  = edge_ind.second = -1;
#ifndef NDEBUG
	{
		box_intersect_line_scratch    bils;
		assert ( bx.open_intersect_line ( &segs[2*ent_seg] , bils) );
	}
#endif	
	box_2d_intersect_line ( bx, &segs[2*ent_seg], entry_pt, exit_pt, edge_ind );
	
	assert ( edge_ind.first >= 0 ); // Check there is an intersection.
	
	if ( out_entry_pt ) {
		out_entry_pt[0]=entry_pt[0]; out_entry_pt[1]=entry_pt[1];
	}

	if ( edge_ind.first == edge_ind.second ) {
	
	// A single edge is intersected so we know the edge does not go straight through.
	
	// Mark for delete. 2008-12-08 ML.
	//	cout << "There was one intersection point on segment " << ent_seg << " : looking for exit point.\n";
	
		// There was only a single intersection.
		// We need to find the exit point.
		// Special handling needed for exit
		// on the same edge.
		num_pts_inside = 0;
		int s=ent_seg;
		
		// While-loop allowing wrap-around.
		// Exit loop if all points are inside or if we find
		// a point that is outside.
		while ( (num_pts_inside<num_segs) ) {
			if ( bx.open_intersect_point( &segs[ ((s<num_segs-1) ? 2*(s+1) : 0) ] ) ) {
				++num_pts_inside;
			} else {
				break;
			}
			if ( s<num_segs-1 ) {
				++s;
			} else {
				s = 0;
			}
		}
		
		
		
		if ( num_pts_inside == num_segs ) {
			abort();
			// Currently not handled.
		}
		
		int exit_seg = ent_seg + num_pts_inside;
		
		// Should be harmless to allow this. 2009-07-06 ML.
		// assert ( ent_seg != exit_seg );
		
		// Debugging found there was a missing -1 here
		// so wrap around was not occurring properly.
		// 2008-12-08 Michael LI.
		if ( exit_seg > num_segs-1 ) {
			// Wrap around.
			exit_seg -= num_segs;
		}
		
		// Should be harmless to allow this. 2009-07-06 ML.
		// assert ( ent_seg != exit_seg );
		
		pair<int,int>    tmp_ind;
		tmp_ind.first = tmp_ind.second = -1;
		
		
		box_2d_intersect_line ( bx, &segs[2*exit_seg], exit_pt, 0, tmp_ind );
		

		assert ( (tmp_ind.first != -1) && (tmp_ind.second != -1) );
		
		if ( out_exit_pt ) {
			out_exit_pt[0]=exit_pt[0]; out_exit_pt[1]=exit_pt[1];
		}
		
		// There should be a single exit point intersection but
		// intersection with a corner could possibly produce two
		// in which case we take the second.
		edge_ind.second = tmp_ind.second;
		
		
		
		// Mark for delete. 2008-12-08 ML.
		// cout << "edge_ind " << edge_ind.first << ", " << edge_ind.second << "\n";
		
		if ( edge_ind.first == edge_ind.second ) {
			// Handle the case of entry and exit on the same edge.
			
// Test allow 2008-12-12 ML.
//			cout << "Should not enter here.\n";
//			abort();
			
			bool pos = false;
			
			inside_corners.resize(4);
			outside_corners.resize(0);
			
			vector<int> &    ic = inside_corners;
			
			switch ( edge_ind.first ) {
			case 0:
				ic[0]= 1; ic[1]=2; ic[2]=3; ic[3]=0;
				pos = (entry_pt[0] < exit_pt[0]);
				break;
			case 1:
				ic[0]= 2; ic[1]=3; ic[2]=0; ic[3]=1;
				pos = (entry_pt[1] < exit_pt[1]);
				break;
			case 2:
				ic[0]= 3; ic[1]=0; ic[2]=1; ic[3]=2;
				pos = (entry_pt[0] > exit_pt[0]);
				break;
			case 3:
				ic[0]= 0; ic[1]=1; ic[2]=2; ic[3]=3;
				pos = (entry_pt[1] > exit_pt[1]);
				break;
			default:
				abort();
			}
			
			if ( !pos ) {
				swap ( inside_corners, outside_corners );
			}
			return;
		}
		
		// If we enter and exit on distinct edges, continue as normal.
	} else {
		// The simplest case where the segment enters and
		// exits the box.
		if ( out_exit_pt ) {
			out_exit_pt[0]=exit_pt[0]; out_exit_pt[1]=exit_pt[1];
		}
		num_pts_inside = 0;
	}
	

	// pos denotes the positive sense as described below.
	vector<int> *    pos_inside;
	vector<int> *    pos_outside;
	
	int lower, upper;
	
	assert ( edge_ind.first  >= 0 && edge_ind.first  <= 3 );
	assert ( edge_ind.second >= 0 && edge_ind.second <= 3 );
	
	// The edges of the box are labelled
	// anticlockwise with the bottom edge being zero.
	// Our case-by-case handling assumes a lower index
	// is intersected before a higher index (positive sense). A simple
	// swapping takes care of the reverse case.
	// Actual single edge intersection is handled separately.
	if ( edge_ind.first < edge_ind.second ) {
		lower = edge_ind.first;
		upper = edge_ind.second;
		pos_inside  = &inside_corners;
		pos_outside = &outside_corners;
	} else {
		lower = edge_ind.second;
		upper = edge_ind.first;
		pos_inside  = &outside_corners;
		pos_outside = &inside_corners;
	}
	
	// Avoid dereferencing later.
	vector<int>& pos_in  = *pos_inside;
	vector<int>& pos_out = *pos_outside;
	
	switch (lower) {
	case 0:
	{
		switch ( upper ) {
		case 1:
		{
			pos_in.resize(3);
			pos_in[0]=2; pos_in[1]=3; pos_in[2]=0;
			pos_out.resize(1);
			pos_out[0]=1;
			break;
		}
		case 2:
		{
			pos_in.resize(2);
			pos_in[0]=3; pos_in[1]=0;
			pos_out.resize(2);
			pos_out[0]=1; pos_out[1]=2;
			break;
		}
		case 3:
		{
			pos_in.resize(1);
			pos_in[0]=0;
			pos_out.resize(3);
			pos_out[0]=1; pos_out[1]=2; pos_out[2]=3;
			break;
		}
		default:
			abort();
		}
		break;
	}
	case 1:
	{
		switch ( upper ) {
		case 2:
		{
			pos_in.resize(3);
			pos_in[0]=3; pos_in[1]=0; pos_in[2]=1;
			pos_out.resize(1);
			pos_out[0]=2;
			break;
		}
		case 3:
		{
			pos_in.resize(2);
			pos_in[0]=0; pos_in[1]=1;
			pos_out.resize(2);
			pos_out[0]=2; pos_out[1]=3;
			break;
		}
		default:
			abort();
		}
		// Break statement was missing. Misformed
		// vtk output. 2008-12-07 Michael LI.
		break;
		
	}
	case 2:
	{
		pos_in.resize(3);
		pos_in[0]=0; pos_in[1]=1; pos_in[2]=2;
		pos_out.resize(1);
		pos_out[0]=3;
		break;
	}
	default:
		cout << "Should not reach here. Claimed entry segment possibly does not intersect box.\n";
		abort();
	}
	
	// For debuggin 2008-12-07 Michael LI.
	assert ( inside_corners.size () != 0 );
	assert ( outside_corners.size () != 0 );
	
}

void get_intersect_keys_entry_indices ( refinement_structure<2> &    rs,
                                        line_segments &              ls,
                                        vector<string> &             keys,
                                        vector<int> &                entry_seg )
{
	valarray<double>    midpoints;
	valarray<int>       associated_segments;
	
	get_grid_intersection_midpoints<2> (
		(rs.dtree)->access_endpoints_grid(),
		ls,
		midpoints,
		associated_segments,
		0,
		false );
		
	
#if 0 // Mark for delete. 2008-12-10 ML.
	cout << "associated segments of "
	     << ((ls.segments.size()/2)-1) << " are :";
	for ( int i=0; i<associated_segments.size(); ++i ) {
		if ( i!= 0 ) {
			cout << ", ";
		}
		cout << associated_segments[i];
	}
	cout << "\n";
#endif
	
	assert ( midpoints.size()/2 == associated_segments.size() );
	
	get_keys_intersect_point<2> ( midpoints,
	                              associated_segments,
	                              ls,
	                              rs,
	                              keys,
	                              entry_seg );
	
#if 0 // Mark for delete. 2008-12-08 ML.
	cout << "Mike debug\n";
	int numpts = midpoints.size()/2;
	for ( int i=0; i<numpts; ++i ) {
		if ( i!=0 ) {
			cout << ", ";
		}
		cout << "(" << midpoints[2*i] << "," << midpoints[2*i+1] << ")[" << associated_segments[i] <<"]";
	}
	cout << "\n";
	for ( int i=0; i<keys.size(); ++i ) {
		if ( i!=0 ) {
			cout << ", ";
		}
		cout << keys[i];
	}
	cout << "\n";
	
#endif

	assert ( keys.size() == entry_seg.size() );
	
}

void get_inside_keys ( vector<string>& in, basic_dbinary_tree<2>& dtree, line_segments& ls, set<string>& out )
{
	box<2>    bx;
	double    centre[2];
	
	vector<string>::const_iterator    it ( in.begin() );
	vector<string>::const_iterator    end ( in.end() );
	
	out.clear();
	
	for ( ; it!=end; ++it ) {
		dtree.get_box ( *it, bx );
		bx.get_centre_point ( centre );
		if ( point_in_polygon ( ls, centre[0], centre[1] ) ) {
			// Assume input is ordered so
			// insert at the back.
			out.insert ( out.end(), *it );
		}
	}
}

void get_inside_points_3D ( const double * in, int num_pts, line_segments& ls, bool * out )
{
	assert ( in && out );
	
	for ( int i=0; i<num_pts; ++i ) {
		if ( point_in_polygon ( ls, in[3*i], in[3*i+1] ) ) {
			out[i] = true;
		} else {
			out[i] = false;
		}
	}
}

/*
	Line runs from (a1,a2) to (b1,b2) with equation of the form

	r = (a1,a2) + lambda (b1-a1,b2-a2)

	and polar coordinate form of line

	r = r0/cos(theta-theta0)
*/
void polar_line_params ( double* line, double* r0, double* theta0 )
{
	assert ( line );
	double a1 = line[0];    double a2 = line[1];
	double b1 = line[2];    double b2 = line[3];
	
	double bma1 = b1-a1; // (b minus a) 1
	double bma2 = b2-a2; // (b minus a) 2
	
	/*
		lambda = - a.(b-a)/(|b-a|^2)
	
		Location of closest point to origin/pole.
	*/
	assert ( fabs(bma1) > 0 );    assert ( fabs(bma2) > 0 );
	double lambda = - ( a1*bma1 + a2*bma2 ) / (bma1*bma1 + bma2*bma2);
	
	// Foot of perpendicular from pole to the line.
	double fp1 = a1 + lambda*bma1;
	double fp2 = a2 + lambda*bma2;
	
	double dist_line_from_pole = sqrt ( fp1*fp1 + fp2*fp2 );
	
	if ( r0 ) {
		*r0 = dist_line_from_pole;
		assert ( !(dist_line_from_pole <0.0) );
	}
	
	double foot_angle = atan2 ( fp1, fp2 );
	
	if ( foot_angle < 0 ) {
		foot_angle += 2.0*M_PI;
	}
	
	if ( theta0 ) {
		*theta0 = foot_angle;
	}
}

//

template void get_bounding_box<2> ( line_segments& ls, box<2>& bx );

template void get_crossing_parameters<1> ( valarray<double>*, int, const double*, valarray<double>& );
template void get_crossing_parameters<2> ( valarray<double>*, int, const double*, valarray<double>& );
template void get_crossing_parameters<3> ( valarray<double>*, int, const double*, valarray<double>& );

template void gnuplot_output<2> ( line_segments & , ostream&, bool, int );
template void gnuplot_output<3> ( line_segments & , ostream&, bool, int );

template void gp_draw<2> ( line_segments&, ostream & );

template void map_line_parameters_to_points<1,vector<double> > ( vector<double>&, const double*, vector<double>&, double* );
template void map_line_parameters_to_points<2,vector<double> > ( vector<double>&, const double*, vector<double>&, double* );
template void map_line_parameters_to_points<3,vector<double> > ( vector<double>&, const double*, vector<double>&, double* );

template void map_line_parameters_to_points<2,valarray<double> > ( valarray<double>&, const double*, valarray<double>&, double* );

template void get_grid_intersection_midpoints<1> ( valarray<double>*, line_segments&, valarray<double>&, valarray<int>&, int, bool );
template void get_grid_intersection_midpoints<2> ( valarray<double>*, line_segments&, valarray<double>&, valarray<int>&, int, bool );
template void get_grid_intersection_midpoints<3> ( valarray<double>*, line_segments&, valarray<double>&, valarray<int>&, int, bool );
