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

#ifndef _DIFFERENTIABLE_FUNCTION_HH_
#define _DIFFERENTIABLE_FUNCTION_HH_

#include "function.hh"
#include "box.hh"
#include <boost/shared_ptr.hpp>
#include <valarray>
#include <vector>
using boost::shared_ptr;
using std::valarray;
using std::vector;

template<int dim>
class differentiable_function : public function<dim> {
public:
	typedef shared_ptr<const differentiable_function<dim> > const_ptr;
	typedef shared_ptr<differentiable_function<dim> >       ptr;
	typedef vector<ptr>                                     vec_ptr;
public:
	differentiable_function ();
	~differentiable_function ();
private:
	differentiable_function ( const differentiable_function<dim>& );             // Not implemented.
	differentiable_function<dim>& operator= ( const differentiable_function<dim>& ); // Not implemented.
private:
	// The user provides this function.
	virtual void evaluate_grad ( const double* co, double* grad ) const=0;

	// The user may optionally override these member functions.
	virtual void evaluate_and_grad ( const double* co, double* value, double* grad ) const;
	virtual void evaluate_grad ( valarray<double>& in, valarray<double>& out ) const;
	virtual void evaluate_grad ( valarray<double>& in, valarray<bool>& pred, valarray<double>& out ) const;
public:
  void global_evaluate_grad ( const box<dim>&,
			        const double*, double*  ) const;
  void local_evaluate_grad ( const box<dim>&,
			       const double*, double* ) const;
  void restricted_local_evaluate_grad ( const box<dim>& restricted_local,
				          const box<dim>& local,
				          const double*, double* ) const;
public:
  void global_evaluate_grad ( const box<dim>&,
                              valarray<double>& points,
                              valarray<double>&       results ) const;
  void local_evaluate_grad (  const box<dim>&,
                              valarray<double>& points,
                              valarray<double>&       results ) const;

  void local_evaluate_and_grad ( const box<dim>&   bx,
                                 valarray<double>& points,
                                 valarray<double>& vals,
                                 valarray<double>& grad ) const;

  void local_evaluate_grad ( const box<dim>& bx,
                             valarray<double>& in,
                             valarray<bool>&   pred,
                             valarray<double>& out ) const;

  void local_evaluate_and_grad ( const box<dim>&   bx,
                                 valarray<double>& points,
                                 valarray<bool>&   pred,
                                 valarray<double>& vals,
                                 valarray<double>& grad ) const;


  void restricted_local_evaluate_grad ( const box<dim>&         restricted_local,
                                        const box<dim>&         local,
                                        valarray<double>& points,
                                        valarray<double>&       results ) const;
};

// Abstract class so need the definitions.
// Use both #include "function.hh" and #include "differentiable_function.cc"
// when inheriting from differentiable_function.

#endif // _DIFFERENTIABLE_FUNCTION_HH_
