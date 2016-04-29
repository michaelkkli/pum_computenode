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

#ifndef _SOLUTION_EVALUATION_STRUCTURE_HH_
#define _SOLUTION_EVALUATION_STRUCTURE_HH_

#include <string>
#include <valarray>
#include <vector>

using std::string;
using std::valarray;
using std::vector;

struct solution_evaluation_structure {
	vector<vector<string> >    eval_keys;
	vector<vector<int> >       eval_indices;
};

void set_outside ( double val, solution_evaluation_structure &, valarray<double> & );

#endif // _SOLUTION_EVALUATION_STRUCTURE_HH_
