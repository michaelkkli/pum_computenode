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

#include "integration_scheme.hh"

//#include "interior_predicate.hh"
#include "point_utils.hh"
#include "box_utils.hh"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <set>

#include <vector>

#include <cmath>
			 
using std::copy;
using std::cout;
using std::find;
using std::map;
using std::ofstream;
using std::pow;
using std::set;
using std::slice;
using std::vector;


template<int dim>
integration_scheme<dim>::integration_scheme( quadrature_rule_type    qrt,
                                             quadrature_rule_type    bdry ) : 
	quad_rule( qrt ), boundary_quad_rule( bdry )
{
}
	
template<int dim>
integration_scheme<dim>::~integration_scheme()
{
}

template<int dim>
void
integration_scheme<dim>::get_quadrature_points ( valarray<double> &    out ) const
{
	// 2009-01-13 Originally missed out the dim* here. ML.
	out.resize ( dim*quad_rule.num_points() );
	out = quad_rule.access_quadrature_points();
}

// For integration!
template<int dim>
void
integration_scheme<dim>::get_mapped_quadrature_points_restricted_local_to_local (
	const box<dim>& restricted,
        const box<dim>& local,
	valarray<double>& out ) const
{
	// Debug
	assert( (quad_rule.access_quadrature_points()).size()%dim == 0 );
	local.map_restricted_local_to_local( restricted, quad_rule.access_quadrature_points(), out );
	
	// Debug
	assert( out.size() % dim == 0 );
}

template<int dim>
valarray<double>&
integration_scheme<dim>::access_quadrature_points() const
{
	return quad_rule.access_quadrature_points();
}

template<int dim>
valarray<double>&
integration_scheme<dim>::access_quadrature_weights() const
{
	return quad_rule.access_quadrature_weights();
}
	
template<int dim>
valarray<double>&
integration_scheme<dim>::access_boundary_quadrature_points() const
{
	return boundary_quad_rule.access_quadrature_points();
}

template<int dim>
valarray<double>&
integration_scheme<dim>::access_boundary_quadrature_weights() const
{
	return boundary_quad_rule.access_quadrature_weights();
}

template<int dim>
double
integration_scheme<dim>::get_mapped_quadrature_final_factor( const box<dim>& bx ) const
{
	return bx.measure()*pow(0.5, dim);
}

template <int dim>
double
integration_scheme<dim>::integrate_global_basis_function( const global_basis_function<dim>&    in,
                                                          integration_scratch &                scratch,
                                                          int                                  decomposition_level ) const
{
	integration_domain_decomposition<dim> tmp;

	while ( !in.get_integration_domain_decomposition( tmp, decomposition_level ) ) {
		--decomposition_level;
		assert( decomposition_level >= 0 );
	}

	int num = tmp.box_vector.size();

	assert( num >= 1 );

	valarray<double> &    pts  = scratch.local;
	valarray<double> &    vals = scratch.vals_first;

	double tmp_sum = 0.0;

	for ( int i=0; i<num; ++i ) {
		get_mapped_quadrature_points_restricted_local_to_local ( tmp.box_vector[i], in.access_box(), pts );

		in.local_evaluate( pts, scratch.global_basis_fun_scratch, vals );

		assert( vals.size() == (access_quadrature_weights()).size() );

		vals *= access_quadrature_weights();

		tmp_sum += get_mapped_quadrature_final_factor( tmp.box_vector[i] ) * vals.sum();
	}

	return tmp_sum;
}

template <int dim>
double
integration_scheme<dim>::integrate_product_global_basis_function( const global_basis_function<dim>& first,
                                                                  const global_basis_function<dim>& second,
                                                                  integration_scratch &             scratch,
                                                                  int                               decomposition_level ) const
{
	integration_domain_decomposition<dim> tmp;

	while ( !global_basis_function<dim>::get_integration_domain_decomposition( first, second, tmp, decomposition_level ) ) {
		--decomposition_level;
		assert( decomposition_level >= 0 );
	}

	int num = tmp.box_vector.size();
	
	assert( num >= 1 );

	valarray<double> &    pts         = scratch.local;
	valarray<double> &    vals_first  = scratch.vals_first;
	valarray<double> &    vals_second = scratch.vals_second;

	double tmp_sum = 0.0;

	for ( int i=0; i<num; ++i ) {
		
		get_mapped_quadrature_points_restricted_local_to_local ( tmp.box_vector[i], first.access_box(), pts );
		first.local_evaluate ( pts, scratch.global_basis_fun_scratch, vals_first  );
		
		get_mapped_quadrature_points_restricted_local_to_local ( tmp.box_vector[i], second.access_box(), pts );
		second.local_evaluate( pts, scratch.global_basis_fun_scratch, vals_second );
		
		assert( vals_first.size() == vals_second.size() );

		valarray<double> &    product = vals_first; // vals_first is recycled.
		product *= vals_second;

		assert( product.size() == (access_quadrature_weights()).size() );

		product *= access_quadrature_weights();

		double to_add = get_mapped_quadrature_final_factor( tmp.box_vector[i] ) * product.sum();
		
		assert ( to_add < 10.0 * tmp.box_vector[i].measure() );
		
		tmp_sum += to_add;
		
	}
	
	return tmp_sum;
}

template <int dim>
double
integration_scheme<dim>::integrate_product_global_basis_function_grad ( const global_basis_function<dim>& first,
                                                                        const global_basis_function<dim>& second,
                                                                        integration_scratch &             scratch,
                                                                        int                               decomposition_level ) const
{
	integration_domain_decomposition<dim> tmp;
	
	while ( !global_basis_function<dim>::get_integration_domain_decomposition( first, second, tmp, decomposition_level ) ) {
		--decomposition_level;
		assert( decomposition_level >= 0 );
	}
	
	int num = tmp.box_vector.size();	
	assert( num >= 1 );
	
	valarray<double> &    pts         = scratch.local;
	valarray<double> &    grad_first  = scratch.grad_first;
	valarray<double> &    grad_second = scratch.grad_second;
	
	double tmp_sum = 0.0;

	valarray<double> &    grad_dot_grad = scratch.grad_dot_grad;
	
	for ( int i=0; i<num; ++i ) {
		get_mapped_quadrature_points_restricted_local_to_local ( tmp.box_vector[i], first.access_box(), pts );
		assert ( pts.size() > 0 );
		
		int pts_size = pts.size();
		assert ( pts_size % dim == 0 );
		int num_pts = pts_size / dim;

		first.local_evaluate_grad ( pts, scratch.global_basis_fun_scratch, grad_first  );
		assert ( pts_size == grad_first.size() );
		
		get_mapped_quadrature_points_restricted_local_to_local ( tmp.box_vector[i], second.access_box(), pts );
		second.local_evaluate_grad ( pts, scratch.global_basis_fun_scratch, grad_second );
		assert ( pts_size == grad_second.size() );
		
		assert( grad_first.size() == grad_second.size() );
		
		valarray<double>& product = grad_second; // grad_second is recycled.
		product *= grad_first;

		grad_dot_grad.resize ( num_pts );
		grad_dot_grad = product[ slice(0, num_pts, dim) ];
		for ( int d=1; d<dim; ++d ) {
			grad_dot_grad += product[ slice( d, num_pts, dim ) ];
		}
		
		assert( grad_dot_grad.size() == (access_quadrature_weights()).size() );
		
		grad_dot_grad *= access_quadrature_weights();
		
		double to_add = get_mapped_quadrature_final_factor( tmp.box_vector[i] ) * grad_dot_grad.sum();
		
		//	assert ( to_add < 10.0 * tmp.box_vector[i].measure() );
		
		tmp_sum += to_add;
		
		assert ( tmp_sum < 1000.0 );
	}
	
	return tmp_sum;
}

template <int dim>
double
integration_scheme<dim>::integrate_product ( const function<dim>&                 global_fun,
                                             const global_basis_function<dim>&    gbasis_fun,
                                             integration_scratch &                scratch,
                                             int                                  decomposition_level ) const
{
	integration_domain_decomposition<dim> tmp;

	while ( !gbasis_fun.get_integration_domain_decomposition( tmp, decomposition_level ) ) {
		--decomposition_level;
		assert( decomposition_level >= 0 );
	}

	int num = tmp.box_vector.size();

	assert( num >= 1 );

	valarray<double> pts, vals_first, vals_second;

	double tmp_sum = 0.0;

	const function<dim>&              first  = global_fun;
	const global_basis_function<dim>& second = gbasis_fun;

	for ( int i=0; i<num; ++i ) {
		first.local_evaluate( tmp.box_vector[i], access_quadrature_points(), vals_first );
		
		get_mapped_quadrature_points_restricted_local_to_local ( tmp.box_vector[i], gbasis_fun.access_box(), pts );
		
		second.local_evaluate( pts, scratch.global_basis_fun_scratch, vals_second );

		assert( vals_first.size() == vals_second.size() );

		valarray<double>& product = vals_first; // vals_first is recycled.
		product *= vals_second;

		assert( product.size() == (access_quadrature_weights()).size() );

		product *= access_quadrature_weights();

		tmp_sum += get_mapped_quadrature_final_factor( tmp.box_vector[i] ) * product.sum();
	}

	return tmp_sum;
}

template<int dim>
void
integration_scheme<dim>::integrate_matrix_block ( vector<global_basis_function<dim> > &    test,
	                              vector<global_basis_function<dim> > &    trial,
                                      function<dim> &              indicator,
	                              valarray<double> &                     mass_matrix_vals,
                                      valarray<double> &                     bleach_matrix_vals,
                                      valarray<double> &                     stiffness_matrix_vals,
	                              integration_scratch &                  scratch,
	                              int                                    decomposition,
	                              differentiable_function<dim> *         equilibrium_fnp,
	                              valarray<double> *                     siggia_matrix_vals_ptr ) const
{
	integration_domain_decomposition<dim> tmp;

	while ( !global_basis_function<dim>::get_integration_domain_decomposition( test[0], trial[0], tmp, decomposition ) ) {
		--decomposition;
		assert( decomposition >= 0 );
	}

	int num_box = tmp.box_vector.size();
	
	assert( num_box >= 1 );
	
	valarray<double> &    restricted_local_quad_pts = scratch.restricted_local_quad_pts;

	valarray<double> &    local_pts_first  = scratch.local_pts_first;
	valarray<double> &    local_pts_second = scratch.local_pts_second;
	
	valarray<double> &    global_pts_first = scratch.global_pts_first;

	valarray<double> &    pu_vals_first = scratch.pu_vals_first;
	valarray<double> &    pu_grad_first = scratch.pu_grad_first;
	valarray<double> &    pu_vals_second = scratch.pu_vals_second;
	valarray<double> &    pu_grad_second = scratch.pu_grad_second;
	valarray<double> &    local_approx_vals_first = scratch.local_approx_vals_first;
	valarray<double> &    local_approx_grad_first = scratch.local_approx_grad_first;
	valarray<double> &    local_approx_vals_second = scratch.local_approx_vals_second;
	valarray<double> &    local_approx_grad_second = scratch.local_approx_grad_second;

	valarray<double> &    vals_first  = scratch.vals_first;
	valarray<double> &    vals_second = scratch.vals_second;
	valarray<double> &    vals_third  = scratch.vals_third;
	
	valarray<double> &    indicator_vals   = scratch.indicator_vals;
	valarray<double> &    equilibrium_vals = scratch.equilibrium_vals;
	valarray<double> &    equilibrium_grad = scratch.equilibrium_grad;
	
	valarray<double> &    grad_first  = scratch.grad_first;
	valarray<double> &    grad_second = scratch.grad_second;
	valarray<double> &    grad_dot_grad = scratch.grad_dot_grad;

	int test_size  = test.size();
	int trial_size = trial.size();
	
	assert( test_size  > 0 );
	assert( trial_size > 0 );

	int vals_size = test_size * trial_size;

	mass_matrix_vals.resize( vals_size, 0.0 );
	bleach_matrix_vals.resize( vals_size, 0.0 );
	stiffness_matrix_vals.resize( vals_size, 0.0 );
	if ( siggia_matrix_vals_ptr ) {
		siggia_matrix_vals_ptr->resize ( vals_size, 0.0 );
	}
	
	// Used to prevent too much copying.
	bool    last_was_interior = false;

	for ( int b=0; b<num_box; ++b ) {
// Mark for delete. 2009-01-14 ML.
//		(tmp.box_vector[b]).map_local_to_global ( access_quadrature_points(),
//                                                         global_pts_first );
		for ( int i=0; i<test_size; ++i ) {
			global_basis_function<dim> &    first = test[i];
			
			if ( i==0 ) {
				if ( !last_was_interior ) {
					get_quadrature_points ( restricted_local_quad_pts );
					last_was_interior = true;
				}
				
				(tmp.box_vector[b]).map_local_to_global ( restricted_local_quad_pts,
                                                                          global_pts_first );
			
				// This must be here, rather than the i-for-loop, due to the use of first.
				(first.access_box()).map_restricted_local_to_local( tmp.box_vector[b], restricted_local_quad_pts, local_pts_first );
				
// Mark for delete. 2009-01-14 ML.
				// 2009-01-13 Replacing with above in preparation for general handling of interior
				// and boundary integration. For some reason, the call below seems to execute quicker.
// 				get_mapped_quadrature_points_restricted_local_to_local ( tmp.box_vector[b],
// 				                                                         first.access_box(),
// 				                                                         local_pts_first );

				first.pu_local_evaluate_and_grad ( local_pts_first,
				                                   scratch.global_basis_fun_scratch,
				                                   pu_vals_first,
				                                   pu_grad_first );
			}
			
			first.local_approx_local_evaluate_and_grad ( local_pts_first,
			                                             scratch.global_basis_fun_scratch,
			                                             local_approx_vals_first,
			                                             local_approx_grad_first );
			
			int num_coord = local_pts_first.size();
			int num_pts   = num_coord/dim;
			
// Mark for delete. 2009-01-14 ML.
// 			first.local_evaluate_and_grad( local_pts_first,
// 			                               scratch.global_basis_fun_scratch,
// 			                               vals_first,
// 			                               grad_first );
			first.local_evaluate_and_grad( local_pts_first,
			                               pu_vals_first,
			                               pu_grad_first,
			                               local_approx_vals_first,
			                               local_approx_grad_first,
			                               scratch.global_basis_fun_scratch,
			                               vals_first,
			                               grad_first );
			
#if 0 // Global points are now calculated so we use those. 2008-12-04 Michael LI.
			// This was incorrectly set to global_evaluate
			// 2008-09-27 Mike Li.
			//indicator.local_evaluate( first.access_box(), pts, indicator_vals );
			indicator.local_evaluate ( tmp.box_vector[b], access_quadrature_points(), indicator_vals );
#endif
			indicator.global_evaluate ( tmp.box_vector[b], global_pts_first, indicator_vals );
			
			assert ( indicator_vals.min() > -1e-10 );
			assert ( indicator_vals.max() < 1 + 1e-10 );

			if ( equilibrium_fnp && siggia_matrix_vals_ptr )
			{
				equilibrium_fnp->global_evaluate ( tmp.box_vector[b], global_pts_first, equilibrium_vals );
				
				equilibrium_fnp->global_evaluate_grad ( tmp.box_vector[b], global_pts_first, equilibrium_grad );
			}

			for ( int j=0; j<trial_size; ++j ) {
				global_basis_function<dim> &    second = trial[j];
				
				if ( j==0 ) {
					(second.access_box()).map_restricted_local_to_local( tmp.box_vector[b],
					                                                restricted_local_quad_pts,
					                                                local_pts_second );
					                                                
// Mark for delete. 2009-01-14 ML.
					// 2009-01-13 Replacing with above in preparation for general handling of interior
					// and boundary integration.
// 					 get_mapped_quadrature_points_restricted_local_to_local ( tmp.box_vector[b],
// 					                                                          second.access_box(),
// 					                                                          local_pts_second );

					second.pu_local_evaluate_and_grad ( local_pts_second,
				                                   scratch.global_basis_fun_scratch,
				                                   pu_vals_second,
				                                   pu_grad_second );
				}

				second.local_approx_local_evaluate_and_grad ( local_pts_second,
			                                                     scratch.global_basis_fun_scratch,
			                                                     local_approx_vals_second,
			                                                     local_approx_grad_second );

// Mark for delete. 2009-01-14 ML.
// 				second.local_evaluate_and_grad( local_pts_second,
// 				                                scratch.global_basis_fun_scratch,
// 				                                vals_second,
// 				                                grad_second );
				second.local_evaluate_and_grad( local_pts_second,
				                                pu_vals_second,
				                                pu_grad_second,
				                                local_approx_vals_second,
				                                local_approx_grad_second,
				                                scratch.global_basis_fun_scratch,
				                                vals_second,
				                                grad_second );
				
				{
					
					
					assert( vals_first.size() == vals_second.size() );
			
					valarray<double> &    vals_product = vals_third; // Make sure to keep vals_second.
					
					// Debugging found that vals_third was not being initialized. 2008-12-04 Michael LI.
					vals_third.resize( vals_second.size() );
					vals_third = vals_second;
					
					vals_product *= vals_first;
			
					assert( vals_product.size() == (access_quadrature_weights()).size() );
			
					vals_product *= access_quadrature_weights();
			
					mass_matrix_vals[ i*trial_size + j ] += get_mapped_quadrature_final_factor( tmp.box_vector[b] ) * vals_product.sum();
					
					assert ( indicator_vals.size() == vals_product.size() );
					vals_product *= indicator_vals;
					bleach_matrix_vals[ i*trial_size + j ] += get_mapped_quadrature_final_factor( tmp.box_vector[b] ) * vals_product.sum();
				}

				{
					valarray<double>& grad_product = grad_second; // grad_second is recycled.
					grad_product *= grad_first;
					
					grad_dot_grad.resize ( num_pts );
					
					grad_dot_grad = grad_product[ slice(0, num_pts, dim) ]; // Initialization.
					for ( int d=1; d<dim; ++d ) {
						grad_dot_grad += grad_product[ slice( d, num_pts, dim ) ];
					}
					
					grad_dot_grad *= access_quadrature_weights();
			
					stiffness_matrix_vals[ i*trial_size + j ] += get_mapped_quadrature_final_factor( tmp.box_vector[b] ) * grad_dot_grad.sum();
				}
				
				if ( equilibrium_fnp && siggia_matrix_vals_ptr )
				{
					valarray<double> &    grad_product = grad_second;
					
					assert ( grad_first.size() == grad_product.size() );
					grad_product = grad_first;
					
					assert ( grad_product.size() == equilibrium_grad.size() );
					grad_product *= equilibrium_grad;
					
					assert ( grad_first.size() == equilibrium_grad.size() );
					
					grad_dot_grad = grad_product[ slice(0, num_pts, dim) ]; // Initialization.
					for ( int d=1; d<dim; ++d ) {
						grad_dot_grad += grad_product[ slice( d, num_pts, dim ) ];
					}
					
					valarray<double> &    equilibrium_integrand = grad_dot_grad;

					assert ( equilibrium_integrand.size() == vals_second.size() );
					equilibrium_integrand *= vals_second;

					assert ( equilibrium_integrand.size() == equilibrium_vals.size() );
					assert ( equilibrium_vals.min() > 0.0 );
					equilibrium_integrand /= equilibrium_vals;
					
					assert ( equilibrium_vals.min() > 0.0 );
					
					equilibrium_integrand *= access_quadrature_weights();
					
					(*siggia_matrix_vals_ptr)[ i*trial_size + j ] += get_mapped_quadrature_final_factor( tmp.box_vector[b] ) * equilibrium_integrand.sum();
				}
			}
		}
	}
	
#if 0//ndef NDEBUG
	for ( int i=0; i<test_size; ++i ) {
		for ( int j=0; j<trial_size; ++j ) {
			
			double a = stiffness_matrix_vals[ i*trial_size + j ];
			double b = integrate_product_global_basis_function_grad ( test[i],
				trial[j],
			                                                                        scratch,
				decomposition );
			
			// std::cout << " a is " << a << " and b is " << b << "\n";

			assert ( fabs(  a - b  ) < 1e-10 );
			
		}
	}
#endif
}

//

template class integration_scheme<2>;
