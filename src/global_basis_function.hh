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

#ifndef _GLOBAL_BASIS_FUNCTION_HH_
#define _GLOBAL_BASIS_FUNCTION_HH_


#include "box.hh"
#include "differentiable_function.hh"
#include "function.hh"
#include "global_basis_function_scratch.hh"
#include "integration_domain_decomposition.hh"
#include "partition_of_unity_function.hh"

#include <boost/shared_ptr.hpp>
#include <valarray>
#include <vector>

using boost::shared_ptr;
using std::valarray;
using std::vector;

//enum global_basis_function_type { gbf_interior, gbf_boundary };

template<int dim>
class global_basis_function {
public:
  typedef shared_ptr<global_basis_function<dim> >       ptr;
  typedef vector<ptr>                                   vec_ptr;
public:
	void set ( const typename box<dim>::vec_ptr&                 patches,
	           int                                               main_patch_index,
	           const typename differentiable_function<dim>::ptr& patch_weight,
	           const typename differentiable_function<dim>::ptr& local_approx );

	void set ( const typename partition_of_unity_function<dim>::ptr& pu_fn,
	           const typename differentiable_function<dim>::ptr&     local_approx );

	const box<dim>& access_box () const;
	
	/**
		Request integration domain decomposition of a certain level of refinement.
		If the level of refinement requested does not exist, false is returned and
		integration_domain_decomposition is not modified.
	*/
	bool get_integration_domain_decomposition( integration_domain_decomposition<dim>&,
	                                           int decomposition_level=0 ) const;

	/**
		Level zero is simple overlap. Level 1 is symmetric decomposition by neighbours.
		For some reason I can't make this a free function and declare it a friend.
	*/
	static bool get_integration_domain_decomposition( const global_basis_function& first,
	                                                  const global_basis_function& second,
	                                                  integration_domain_decomposition<dim>&, int level=0 );

	double global_evaluate ( const double* ) const;
	double local_evaluate ( const double* ) const;
	
	void global_evaluate_grad ( const double*, double* ) const;
	void local_evaluate_grad ( const double*, double* ) const;
	
	void local_evaluate_grad ( valarray<double>& points,
	                           global_basis_function_scratch &,
	                           valarray<double>& result ) const;
	
	void local_evaluate( valarray<double>& points,
	                     global_basis_function_scratch &,
	                     valarray<double>& result ) const;
	
	void local_evaluate_and_grad ( valarray<double> &                points,
	                               global_basis_function_scratch &   ,
	                               valarray<double> &                vals,
	                               valarray<double> &                grad ) const;
	
	/**
		We can keep the pu evaluations outside and iterate over
		local approx function evaluations. The older method would
		carry out redundant pu re-evaluations. The older method
		was for reduction of complexity at the basis function level
		but it's now clear we can get around this with not-too-major
		modifications. 2009-01-14 ML.
	*/
	void local_evaluate_and_grad ( valarray<double> &                in,
	                               valarray<double> &                pu_vals,
	                               valarray<double> &                pu_grad,
	                               valarray<double> &                local_approx_vals,
	                               valarray<double> &                local_approx_grad,
	                               global_basis_function_scratch &   scratch,
	                               valarray<double> &                vals,
	                               valarray<double> &                grad ) const;

	void pu_local_evaluate_and_grad ( valarray<double> &                 points,
	                                  global_basis_function_scratch &    scratch,
	                                  valarray<double> &                 vals,
	                                  valarray<double> &                 grad ) const;
	void local_approx_local_evaluate_and_grad ( valarray<double> &                 points,
	                                            global_basis_function_scratch &    scratch,
	                                            valarray<double> &                 vals,
	                                            valarray<double> &                 grad ) const;
	
private:
	typename partition_of_unity_function<dim>::ptr  pu_function_;
	typename differentiable_function<dim>::ptr      local_approx_function_;
};


#endif // _GLOBAL_BASIS_FUNCTION_HH_
