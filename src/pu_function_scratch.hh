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

#ifndef _PU_FUNCTION_SCRATCH_HH_
#define _PU_FUNCTION_SCRATCH_HH_

#include <valarray>

using std::valarray;

struct pu_function_scratch {
	// Function evaluation
	valarray<double>    vals_first;
	valarray<double>    vals_second;
	
	// Mapped points.
	valarray<double>    global;
	valarray<double>    local;
	
	valarray<bool>      pred;
	
	// Weight evaluation for global basis function
	// gradient evaluation.
	valarray<double>    weight_sum;
	valarray<double>    weight_tmp;
	valarray<double>    main_weight;
	
	// Gradient evaluation for global basis function.
	valarray<double>    grad_sum;
	valarray<double>    grad_tmp;
	valarray<double>    main_grad;
	
	//	valarray<double>    grad_dot_grad;
};

#endif // _PU_FUNCTION_SCRATCH_HH_
