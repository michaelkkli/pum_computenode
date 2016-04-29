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

#ifndef _PARTITION_OF_UNITY_FUNCTION_HH_
#define _PARTITION_OF_UNITY_FUNCTION_HH_

#include "box.hh"
#include "differentiable_function.hh"
#include "pu_function_scratch.hh"

#include <valarray>
#include <vector>
using std::vector;

template <int dim>
class partition_of_unity_function {
public:
	typedef shared_ptr<partition_of_unity_function<dim> > ptr;
	typedef vector<ptr>                                   vec_ptr;
public:
	partition_of_unity_function ();
	~partition_of_unity_function ();
	const box<dim>& access_box() const;
	void set ( const typename box<dim>::ptr& ptch,
	           const typename box<dim>::vec_ptr& np,
	           const typename differentiable_function<dim>::ptr& pw );
#if 0
	// DEPRECATED.
	void set( const typename box<dim>::ptr& ptch,
	          const typename box<dim>::vec_ptr& np,
	          const typename differentiable_function<dim>::ptr& pw,
	          const typename differentiable_function<dim>::vec_ptr& nw );
#endif
	void get_patch_decomposition_by_neighbours ( vector<box<dim> >& ) const;
	void get_patch_decomposition_by_neighbours ( const box<dim>&, vector<box<dim> >& ) const;
public:
	static void get_patch_decomposition_by_neighbours ( const partition_of_unity_function<dim>&,
	                                                    const partition_of_unity_function<dim>&,
	                                                    vector<box<dim> >& );
public:
	double global_evaluate ( const double* ) const;
	double local_evaluate ( const double* ) const;
	void local_evaluate ( valarray<double> & in,
	                      pu_function_scratch &,
	                      valarray<double>& out ) const;
	
	void global_evaluate_grad( const double*, double* ) const;
	void local_evaluate_grad( const double*, double* ) const;
	void local_evaluate_grad ( valarray<double>& in,
	                           pu_function_scratch &,
	                           valarray<double>& out ) const;
	void local_evaluate_and_grad ( valarray<double> & in,
	                               pu_function_scratch &,
	                               valarray<double> & vals,
	                               valarray<double> & grad ) const;
	//double restricted_local_evaluate ( const box<dim>&, const double* ) const;
private:
	typename box<dim>::ptr                         patch_;
	typename box<dim>::vec_ptr                     neighbour_patches_;

	typename differentiable_function<dim>::ptr     patch_weight_;
	
	// typename differentiable_function<dim>::vec_ptr neighbour_weights_;
};

#endif // _PARTITION_OF_UNITY_FUNCTION_HH_
