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

#ifndef _DOF_STRUCTURE_HH_
#define _DOF_STRUCTURE_HH_

#include <boost/shared_ptr.hpp>
#include <valarray>

using boost::shared_ptr;
using std::valarray;

struct dof_structure {
	typedef shared_ptr<dof_structure>   ptr;
	
	valarray<int>    patch_to_num_dof;             /// (1)
	valarray<int>    patch_to_start_index;         /// (2)
	valarray<int>    computenode_to_num_patches;   /// (3)
	valarray<int>    computenode_to_start_patch;   /// (4)
	valarray<int>    computenode_to_num_dof;       /// (5)
	valarray<int>    patch_indices_on_computenode; /// (6)
	valarray<int>    computenode_to_start_dof;     /// (7)
};

#endif // _DOF_STRUCTURE_HH_
