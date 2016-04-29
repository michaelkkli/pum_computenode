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

#ifndef _INTEGRATION_SCHEME_HH_
#define _INTEGRATION_SCHEME_HH_

// Caused a horrible cyclic dependency with unhelpful errors.
//#include "pum_discretization.hh"

#include "box.hh"
#include "dbinary_tree.hh"
#include "function.hh"
#include "global_basis_function.hh"
#include "indicator_function.hh"
#include "integration_scratch.hh"
#include "quadrature_rule.hh"

#include <boost/shared_ptr.hpp>
#include <functional>
#include <valarray>

using boost::shared_ptr;

using std::unary_function;
using std::valarray;



template<int dim>
class integration_scheme {
public:

	typedef shared_ptr<const integration_scheme<dim> > const_ptr;
	typedef shared_ptr<integration_scheme<dim> >       ptr;

public:
	integration_scheme( quadrature_rule_type = gauss4, quadrature_rule_type boundary_rule= gauss5 );
	~integration_scheme();
public:
//	void set_dbinary_tree ( const typename dbinary_tree<dim>::const_ptr & );

	/**
		Used to allow general handling of integration over boxes and
		polygons.
	*/
	void  get_quadrature_points ( valarray<double> & ) const;
	
	void  get_mapped_quadrature_points_restricted_local_to_local ( const box<dim>& restricted,
	                                                               const box<dim>& local,
	                                                               valarray<double>& ) const;
private:
	valarray<double>& access_quadrature_points() const;
	valarray<double>& access_quadrature_weights() const;
	
	valarray<double>& access_boundary_quadrature_points() const;
	valarray<double>& access_boundary_quadrature_weights() const;
public:
	double get_mapped_quadrature_final_factor( const box<dim>& ) const;
public:
	double integrate_global_basis_function( const global_basis_function<dim>&,
	                                        integration_scratch &,
	                                        int decomposition_level ) const;
	double integrate_product_global_basis_function( const global_basis_function<dim>&,
	                                                const global_basis_function<dim>&,
	                                                integration_scratch &,
	                                                int decomposition_level ) const;
	
	double integrate_product_global_basis_function_grad ( const global_basis_function<dim>&,
	                                                      const global_basis_function<dim>&,
	                                                      integration_scratch &,
	                                                      int decomposition_level ) const;
	
	double integrate_product ( const function<dim>&,
	                           const global_basis_function<dim>&,
	                           integration_scratch &,
	                           int decomposition_level ) const;
	
	void integrate_matrix_block ( vector<global_basis_function<dim> > &    patch_test,
	                              vector<global_basis_function<dim> > &    patch_trial,
	                              function<dim> &              indicator,
	                              valarray<double> &                     mass_matrix_block,
	                              valarray<double> &                     bleach_matrix_block,
	                              valarray<double> &                     stiffness_matrix_block,
	                              integration_scratch &                  scratch,
	                              int decomposition,
	                              differentiable_function<dim> *         equilibrium_fnp=0,
	                              valarray<double> *                     siggia_matrix_vals_ptr=0 ) const;
	
public:
	double integrate_product( const box<dim>&, const typename function<dim>::const_ptr&,
	                          const box<dim>&, const typename function<dim>::const_ptr&,
	                          const typename function<dim>::const_ptr& indicator_or_scaling );

private:
//	typename dbinary_tree<dim>::const_ptr dtree;
	mutable quadrature_rule<dim>    quad_rule;
	mutable quadrature_rule<1>      boundary_quad_rule;
private:
	integration_scheme<dim>( const integration_scheme<dim>& );            // Not implemented.
	integration_scheme<dim> operator= ( const integration_scheme<dim>& ); // Not implemented.
};

#endif // _INTEGRATION_SCHEME_HH_
