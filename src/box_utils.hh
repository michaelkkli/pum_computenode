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

#ifndef _BOX_UTILS_HH_
#define _BOX_UTILS_HH_

#include "box.hh"

#include <iostream>
using std::ostream;

#include <valarray>
using std::valarray;

#include <vector>
using std::vector;

#include <algorithm>

template<int dim>
void
decompose_box( const box<dim>&,
	       const vector<box<dim> >&,
               vector<box<dim> >& );

template <int dim>
void
decompose_box( const box<dim>&,
               const typename box<dim>::vec_ptr&,
               vector<box<dim> >& );

#if 0
// A pointer to the start of a container_t[dim] must be passed.
// container_t is intended to be simliar to set<double>
// and container2 is intended to be similar to vector<box<dim> >.
// The implementation is the best and cleanest I can think of
// but doesn not count as nice.
template <typename container_t, typename container2>
void
make_boxes ( int, const container_t*, container2& );
#endif

template<int dim>
void
make_boxes ( const vector<vector<double> >&, vector<box<dim> >& );

void make_box_with_radius ( double * centre, double radius, box<2>& bx );

template<int dim>
void
local_to_restricted_local ( valarray<double>& points,
			    const box<dim>& restricted_local,
			    const box<dim>& local_box,
			    valarray<double>& out );

template<int dim>
void output_description ( const box<dim>&,
			  ostream& );

template<int dim>
void make_bounding_box ( const vector<box<dim> >&, box<dim>& );

#if 0
template <int dim, typename container_t>
void split_box ( const box<dim>&, container_t& container );
#endif

template <int dim, typename coords>
int count_closed_intersect ( const box<dim>&,
			     coords& );

template <int dim, typename coords>
int count_open_intersect ( const box<dim>&,
			   coords& );

//template <int dim>
//void closed_intersect_point ( const vector<box<dim> > &, const double*, valarray<bool>& );
			   
template<int dim>
void box_project_axis ( const box<dim>&, int, const box<dim-1>& );

// These are not templated because of the very different ways that
// the boxes are drawn.
void gp_draw_single( const box<1>&, std::ostream&, double level=0.0 );
void gp_draw_single( const box<2>&, std::ostream&, double level=0.0 );
void gp_draw_single( const box<3>&, std::ostream&, double ignored=0.0 );

template<int dim>
void gp_draw( const std::vector<box<dim> >&, std::ostream&, double level=0.0 );

#endif // _BOX_UTILS_HH_
