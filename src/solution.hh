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

#ifndef _SOLUTION_HH_
#define _SOLUTION_HH_

#include "dof_structure.hh"
#include "function.hh"
#include "global_approximation_space.hh"
#include "petsc_solver.hh"

#include <boost/shared_ptr.hpp>
#include <valarray>
#include <vector>

using boost::shared_ptr;
using std::valarray;
using std::vector;

template <int dim>
class petsc_solver;

template <int dim>
class pum_convergence;

template<int dim>
class solution {
public:
  typedef shared_ptr<const solution<dim> > const_ptr;
  typedef shared_ptr<solution<dim> >       ptr;
	typedef vector<ptr>                vec_ptr;
	friend class petsc_solver<dim>;
	friend class pum_convergence<dim>;
public:
	solution ();
	
	// Size hint must be accurate.
	int set_global_approximation_space ( const typename global_approximation_space<dim>::ptr &,
	                                      int size_hint=0 );
	void set_dof_structure ( const typename dof_structure::ptr & );
	
	bool is_initialized () const;
	
	// Used for testing. Provided valarray must be of the right size.
	void set_coefficients ( valarray<double>& in );
	
	// DEPRECATED and removed to avoid having two ways of doing the same thing.
	// The version using a solution evaluation structure should always be used.
	// This will ensure better testing and avoid debugging confusion with two
	// evaluation methods moving out of step.
	// 2009-02-13 ML.
#if 0
	void global_evaluate ( valarray<double> &, valarray<double>& out );
#endif
	
	
	/**
		Evaluate solution without working out which patches cover the point.
	*/
	void global_evaluate ( valarray<double> &, solution_evaluation_structure &, valarray<double> & );

private:
  bool                                             initialized_;
  typename global_approximation_space<dim>::ptr    global_approx_space_;
  valarray<double>                                 global_coefficients_;
	typename dof_structure::ptr                dof_struct_;
};

#endif // _SOLUTION_HH_
