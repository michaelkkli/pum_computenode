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

#ifndef _BASIC_DBINARY_TREE_UTILS_HH_
#define _BASIC_DBINARY_TREE_UTILS_HH_

#include "basic_dbinary_tree.hh"

#include <iostream>
#include <string>
#include <vector>
using std::ostream;
using std::string;
using std::vector;

template <int dim>
void gnuplot_output ( vector<string> &             keys,
                      basic_dbinary_tree<dim> &    dtree,
                      ostream &                    out,
                      bool                         labels = true );


#endif // _BASIC_DBINARY_TREE_UTILS_HH_

