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

#ifndef _COVER_STRUCTURE_HH_
#define _COVER_STRUCTURE_HH_

#include "box.hh"
#include "refinement_structure.hh"

#include <boost/shared_ptr.hpp>

#include <map>
#include <set>
#include <string>
#include <vector>

using boost::shared_ptr;

using std::map;
using std::set;
using std::string;
using std::vector;

/**
 * Struct to minimize dependency of global_approximation_space on dbinary_tree.
 * A patch is conisidered its own neighbour and is explicitly present in the data
 * structures.
*/
template <int dim>
struct cover_structure {
	typedef shared_ptr<cover_structure<dim> >  ptr;
	
	string                                     support_name;
	
	typename refinement_structure<dim>::ptr    ref_structure;

	set<string>                                patch_keys;
	map<string, set<string> >                  patch_neighbour_keys;

	vector<box<dim> >                          patches;
	vector<vector<box<dim> > >                 neighbour_patches;
	
	vector<typename box<dim>::ptr>             ptr_patches;
	vector<vector<typename box<dim>::ptr> >    ptr_neighbour_patches;

	// Populated by generate_index_information.
	vector<string>                             index_to_patch_key;
	map<string, int>                           patch_key_to_index;
	vector<vector<int> >                       index_to_neighbour_indices;
};

/**
	Boundary information is now kept in cover_structure to allow
	restriction of integration. 2009-03-06 ML.
*/
template <>
struct cover_structure<2> {
	typedef shared_ptr<cover_structure<2> >  ptr;
	
	string                                     support_name;
	
	refinement_structure<2>::ptr    ref_structure;

	set<string>                                patch_keys;
	map<string, set<string> >                  patch_neighbour_keys;
	
	vector<box<2> >                          patches;
	vector<vector<box<2> > >                 neighbour_patches;
	
	vector<box<2>::ptr>             ptr_patches;
	vector<vector<box<2>::ptr> >    ptr_neighbour_patches;

	// Populated by generate_index_information.
	vector<string>                             index_to_patch_key;
	map<string, int>                           patch_key_to_index;
	vector<vector<int> >                       index_to_neighbour_indices;
	
	// To allow integration restriction.
	//	set<int>                                   boundary_indices;

	// These will allow restriction of domain integration.
	vector<bool>                               is_boundary_patch;
	//	vector<string>                             boundary_keys;
	//vector<int>                                entry_segment;
	
};


/**
	The index information can be generated purely from the patch set and
	neighbour information but we want to create fully valid cover_structures
	all in one go from dbinary_tree so all generation is done in dbinary_tree
	with correct decomposition of information separate from dbinary_tree.
*/
template <int dim>
void generate_index_information ( cover_structure<dim>& );

template <int dim>
void generate_patch_information ( cover_structure<dim> & );

template <int dim>
void generate_patch_boundary_information ( cover_structure<dim> &,
                                           set<string> & boundary_keys,
                                           line_segments& domain );

/**
 * Naive implementation making use only of cover_structure.
 * A more efficient implementation is possible with access to the
 * dbinary_tree but we want maximum separation for ease of coding
 * and testing. The better implementation should be made available
 * when the dbinary tree is supplied.
 */
template <int dim>
void get_patch_keys_incident_point ( cover_structure<dim>&, const double*, vector<string> & );

#endif // _COVER_STRUCTURE_HH_
