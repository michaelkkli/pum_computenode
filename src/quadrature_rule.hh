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

#ifndef _QUADRATURE_RULE_HH_
#define _QUADRATURE_RULE_HH_

#include "function.hh"
#include "box.hh"
#include <boost/shared_ptr.hpp>
#include <valarray>
#include <vector>

using boost::shared_ptr;
using std::slice;
using std::valarray;
using std::vector;

template <int dim>
struct simple_quadrature_rule {
	valarray<double>    weights;
	valarray<double>    points;
};

class line_segments;

void generate_gauss_legendre_rule ( int N, simple_quadrature_rule<1> & );

void generate_gauss_legendre_rule ( int N, simple_quadrature_rule<2> & );

void generate_gauss_legendre_rule_local ( simple_quadrature_rule<1>&   base_rule,
                                          const box<2> &                 box1,
                                          line_segments &                segs,
                                          int                            entry_segment,
                                          bool                           do_inside,
                                          simple_quadrature_rule<2> &    out_rule,
                                          simple_quadrature_rule<2> *    boundary_rule=0);

enum quadrature_rule_type { gauss1, gauss2, gauss3, gauss4, gauss5 };

template<int dim>
class quadrature_rule {	
public:
	typedef shared_ptr<quadrature_rule<dim> >    ptr;
  quadrature_rule ( quadrature_rule_type = gauss3 );
  ~quadrature_rule();
public:
  int num_points() const;
  valarray<double>& access_quadrature_points();
  valarray<double>& access_quadrature_weights();
  double integrate ( const box<dim>&,
		     const function<dim>&,
		     valarray<double>* = 0 );

  double integrate_product ( const box<dim>&      restricted_local,
			     const box<dim>&      local1,
			     const function<dim>& f1,
			     const box<dim>&      local2,
			     const function<dim>& f2 );
public: // For bluebird
  valarray<double> quadrature_points;
  valarray<double> quadrature_weights;
public:
//  friend class quadrature_rule<2>;
//  friend class quadrature_rule<3>;
};

#endif // _QUADRATURE_RULE_HH_
