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

#ifndef _INTEGRATION_SCRATCH_HH_
#define _INTEGRATION_SCRATCH_HH_

#include "global_basis_function_scratch.hh"
#include "quadrature_rule.hh"

#include <valarray>

using std::valarray;

struct integration_scratch {
	valarray<double>    restricted_local_quad_pts;
	valarray<double>    restricted_local_quad_weights;
	
	simple_quadrature_rule<2>    boundary_quad_rule;
	
	double              mapped_quadrature_final_factor;
	
	valarray<double>    global_fun_val;
	valarray<double>    global_fun_grad;
	
	// Used to keep a product that is used for function
	// integration.
	valarray<double>    global_fun_val_pu_val;
	valarray<double>    global_fun_val_squared_pu_val;
	valarray<double>    global_fun_grad_global_fun_grad;
	valarray<double>    global_fun_grad_squared_pu_val;
	
	valarray<double>    rhs_fun_val;
	valarray<double>    approx_val_rhs_fun_val;
	
	valarray<double>    approx_val_global_fun_val;
	valarray<double>    approx_grad_global_fun_grad;
	valarray<double>    approx_grad_dot_global_fun_grad;
	
	valarray<double>    first_approx_val_second_approx_val;
	valarray<double>    first_approx_grad_second_approx_grad;
	valarray<double>    first_approx_grad_dot_second_approx_grad;

	// Function evaluation
	valarray<double>    pu_vals_first;
	valarray<double>    pu_grad_first;
	valarray<double>    pu_vals_second;
	valarray<double>    pu_grad_second;
	valarray<double>    local_approx_vals_first;
	valarray<double>    local_approx_grad_first;
	valarray<double>    local_approx_vals_second;
	valarray<double>    local_approx_grad_second;
	valarray<double>    vals_first;
	valarray<double>    vals_second;
	valarray<double>    vals_third;

	valarray<double>    indicator_vals;
	valarray<double>    equilibrium_vals;
	valarray<double>    equilibrium_grad;
	
	valarray<double>    grad_first;
	valarray<double>    grad_second;

	// Mapped points.
	// valarray<double>    global;
	valarray<double>    local; // Use of this is deprecated.
	
	valarray<double>    local_pts_first;
	valarray<double>    local_pts_second;
	valarray<double>    global_pts_first;

	
	valarray<double>    grad_dot_grad;
	
	global_basis_function_scratch    global_basis_fun_scratch;
};

/*
template <int dim>
void resize ( integration_scratch &, int num_pts );
*/

#endif // _INTEGRATION_SCRATCH_HH_
