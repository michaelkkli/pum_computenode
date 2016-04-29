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

#ifndef _VTK_OUTPUT_HH_
#define _VTK_OUTPUT_HH_

#include "basic_dbinary_tree.hh"
#include "line_segments.hh"

#include <fstream>
#include <iostream>
#include <string>
#include <valarray>
#include <vector>

#include <boost/shared_ptr.hpp>

using boost::shared_ptr;

using std::string;
using std::ofstream;
using std::ostream;
using std::valarray;
using std::vector;

enum visualize_choice { inside_outside, inside, outside };

struct vtk_output {
	typedef shared_ptr<vtk_output> ptr;
	valarray<double>    points_3D;
	valarray<double>    scalars;
	valarray<int>       connectivity;
	valarray<int>       offsets;
	
	// 2009-09-05 New addition. Michael Li.
	valarray<double>    cell_scalars;
	
	vtk_output& operator= ( const vtk_output& other );
};

struct pvd_output {
	string head;
	vector<string> middle;
	string tail;
};

bool little_endian ();

/**
 * VTK binary format requires base64 encoding. Three chars
 * are encoded to four ASCII characters.
 */
void base64 ( unsigned char * in, ostream & out, int len );

// vtk_output& must not be const as we need access to valarray non-const operator []
void vtk_simple_output ( vtk_output&, const string &, bool binary=false );

void vtk_append_raw ( vtk_output&, const string &, bool do_cell_scalars=false );

void pvd_simple_output ( const pvd_output&, ofstream& );

void get_entry_connectivity_offsets ( vector<string> &         bdry_keys,
                                            vector<box<2> > &        boxes,
                                            vector<int> &            entry_seg,
                                            basic_dbinary_tree<2> &  dtree,
                                            line_segments &          ls,
                                            vector<double> &         entry,
                                            vector<int> &            connectivity,
                                            vector<int> &            offsets,
                                            int                      rebase_bdry,
                                            int                      rebase_first_ent,
                                            visualize_choice         viz_choice );

// Let's try to get just the boundary working first.
void make_boundary_vtk_output ( line_segments &              ls,
                                refinement_structure<2> &    rs,
                                vtk_output &                 vtkout,
                                visualize_choice             viz_choice=inside,
                                vector<string> *             optional_boundary_keys_out=0,
                                vector<int> *                opt_entry_set_out=0 );

#endif
