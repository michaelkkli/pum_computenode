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

#ifndef _LINE_SEGMENTS_HH_
#define _LINE_SEGMENTS_HH_

#include "box.hh"
#include "box_utils.hh"

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <valarray>
#include <vector>

using std::map;
using std::ostream;
using std::pair;
using std::set;
using std::string;
using std::valarray;
using std::vector;

template <int dim> class refinement_structure;
template <int dim> class basic_dbinary_tree;
class line_segments;

/**
	If there is an intersection, the function returns true.

	first_line is a pointer to the first of 2*dim doubles representing the start and end points of the first line

	second_line is a pointer to the first of 2*dim doubles representing the start and end points of the second line

	param1 is a double representing the intersection point in a [-1, 1] parameterization of the first line. Non-intersection is represented by a value less than -1.5.

	param2 is a double representing the intersection point in a [-1, 1] parameterization of the second line. Non-intersection is represented by a value less than -1.5.

	intersection_point is a pointer to the first of dim doubles to be used for returning the intersection point if an intersection occurs. A null pointer may be passed so that the intersection point is not calculated.
*/
bool intersect_lines_2d ( const double* first_line,
                          const double* second_line,
                          double& param1,
                          double& param2,
                          double* intersection_point );

/**
	Equality of the pair denotes single intersection.
	
	Return of intersection points is optional. Null pointers may be passed.
	
	     2
	    ____
	   |    |
	 3 |    | 1
	   |____|
	   
	     0
*/
void box_2d_intersect_line ( const box<2> &   bx,
                            const double *   line,
                            double *         first_intersection_pt,
                            double *         second_intrsctn_pt,
                            pair<int,int>&   ordered_intersections );

/**
	segments holds the coordinates x_0, .. , x_n points.

	If the segments represent a closed boundary, the last point is repeated.

	There are segments.size()/dim points.

	The first and last points are identical if the segments represent a closed boundary.

	Each adjacent pair of points represents one line segment and there are
	(segments.size()/dim)-1 line segments in total.
*/
struct line_segments {
	valarray<double>    segments;
};


/**
	Find the first entry segment and the number of line_segment points that
	lie inside the box. An invalid start_hint will result in starting from
	the first line segment.
	
	This function assumes that the boundary-box intersection is made up of
	no more than one piece.
	
	Max points inside is provided primarily for boundary patch finding when we know
	that a child patch can only intersect a subset of the line segments that
	its parent patch intersects.
	

	Use box_2d_intersect_line to find if an intersection exists.
	If there is no intersection with the edge of the box, the actual
	first entry is -1 and the actual_num_inside is zero.
*/
void get_first_entry_num_inside ( line_segments& ls,
                                  box<2>& bx,
                                  int start_hint,
                                  int max_points_inside,
                                  int& actual_first_entry,
                                  int& num_found_inside );

template <int dim>
void get_bounding_box ( line_segments& ls, box<dim>& bx );

void regularize_bounding_box ( box<2>& bx );

template <int dim>
void gnuplot_output ( line_segments &    ls,
                      ostream &          out,
                      bool               with_labels = true,
                      int                max_num_labels=50 );

template <int dim>
void gp_draw ( line_segments &    ls,
               ostream &          out );

/**
	Work out the parameters of grid line crossings.
*/
template <int dim>
void get_crossing_parameters ( valarray<double> *    grid,
                               int                   chosen_dimension,
                               const double *        line,
                               valarray<double> &    params );

void replace_parameters_by_midpoints ( vector<double>& params );

template <int dim, class T>
void map_line_parameters_to_points ( T &    params,
                                     const double *      line,
                                     T &    points,
                                     double *            Jacobian_determinant=0 );

/**
	To be used to identify boundary patches for a particular
	cover, and to be used in identifying boxes in a domain
	decomposition through which the boundary passes.
	
	grid points to the first element of an array of valarray<double> of dimension dim
	which is used to specify a grid.
*/
template <int dim>
void get_grid_intersection_midpoints ( valarray<double>*     grid,
                                       line_segments &       segs,
                                       valarray<double> &    out_points,
                                       valarray<int> &       associated_segments,
                                       int                   first_entry_hint=0,
                                       bool                  finish_after_first_exit=true );

/**
	Used to process the output from get_grid_intersection_midpoints.
	The vector is of intersecting keys without repetition and
	first entry gives the index of the first entry segment.
*/
template <int dim>
void get_keys_intersect_point ( valarray<double> &          points,
                                valarray<int> &             seg_indices,
                                line_segments &             ls,
                                refinement_structure<dim> & rs,
                                vector<string> &      intersecting,
                                vector<int> &         first_entry );

/**

*/
void split_by_line_segments ( const box<2> &           bx,
                              line_segments &    ls,
                              int                entry_segment,
                              int &              num_pts_inside,
                              double *           entry_pt,
                              double *           exit_pt,
                              vector<int> &    inside_corners,
                              vector<int> &    outside_corners );

/**
	Sufficient to get boundary keys for visualization. A little
	more work must be done to get boundary keys for simulation as
	the boxes are expanded which are then able to intersect the
	boundary even when the unexpanded boxes do not.
*/
void get_intersect_keys_entry_indices ( refinement_structure<2> &    rs,
                                              line_segments &              ls,
                                              vector<string> &             keys,
                                              vector<int> &                entry_seg );

/**
	Retrieve keys where the box centre lies within the region defined by line_segments.
	The retrieved keys may intersect the line_segments.
*/
void get_inside_keys ( vector<string>& in, basic_dbinary_tree<2>& dtree, line_segments& ls, set<string>& out );

/**
	Used to make a valarray<bool> mask for vtk_output points that are given in 3D
	coordinates.
*/
void get_inside_points_3D ( const double * in, int num_pts, line_segments& ls, bool * out );

void polar_line_params ( double* line, double* r0, double* theta0 );

#endif // _LINE_SEGMENTS_HH_
