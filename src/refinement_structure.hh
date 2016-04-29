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

#ifndef _REFINEMENT_STRUCTURE_HH_
#define _REFINEMENT_STRUCTURE_HH_

#include "basic_dbinary_tree.hh"
#include "key_utils.hh"
#include "line_segments.hh"

#include <map>
#include <utility>
#include <set>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

template <int dim>
struct refinement_structure {
	typedef shared_ptr<refinement_structure<dim> >    ptr;

	map<string, int, levelwise_less>                  levels;
	typename basic_dbinary_tree<dim>::ptr             dtree;
	double                                            cover_factor;
};

template <int dim>
void generate_keys ( refinement_structure<dim> & rs, vector<string> & out );

template <int dim>
void generate_keys_and_neighbours ( refinement_structure<dim> &    rs,
                                   vector<string> &                keys,
                                   map<string, set<string> >&      nbrs,
                                   line_segments*                  boundary_segs=0,
                                   map<string, pair<int,int> >*    boundary_detail=0 );

template <int dim>
void get_key_intersect_point ( refinement_structure<dim> &    in,
                               const double*                  co,
                               string &                       out );

#endif // _REFINEMENT_STRUCTURE_HH_
