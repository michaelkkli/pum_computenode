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

#ifndef _POLYNOMIAL_HH_
#define _POLYNOMIAL_HH_

#include "differentiable_function.hh"
#include "function.hh"

#include <boost/shared_ptr.hpp>
#include <algorithm>
#include <valarray>
#include <vector>
#include <cmath>

using boost::shared_ptr;
using std::valarray;
using std::vector;

template<int dim>
class polynomial : public differentiable_function<dim> {
public:
	typedef shared_ptr<const differentiable_function<dim> > const_ptr;
	typedef shared_ptr<differentiable_function<dim> >       ptr;
	typedef vector<ptr>                                     vec_ptr;
public:
  polynomial ();
  polynomial ( double factor, int xpow=0, int ypow=0, int zpow=0 );
  polynomial ( const polynomial& );
  polynomial& operator= ( const polynomial& );
  polynomial& operator+= ( const polynomial& );
  ~polynomial ();

  int get_num_terms() const;
  int get_order() const;

  void clear();

  //  void resize( int number_of_terms );
  void add_term ( double factor, int, int, int );
private:
  double evaluate ( const double* ) const;
  void evaluate_grad ( const double*, double* ) const;

  // Already implemented in function<dim>
  //void evaluate ( const valarray<double>&, valarray<double>& ) const;
private:
  void deep_copy ( const polynomial& );

private:
  int number_of_terms;
  int polynomial_order;

  double common_multiplier;
  vector<double> multipliers;
  // size is some multiple of dim.
  vector<int> powers;

};

#endif // _POLYNOMIAL_HH_
