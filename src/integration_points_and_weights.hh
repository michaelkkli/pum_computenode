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

#ifndef _INTEGRATION_POINTS_AND_WEIGHTS_
#define _INTEGRATION_POINTS_AND_WEIGHTS_

#include "function.hh"
#include "box_utils.hh"
#include "box.hh"

#include <vector>
using std::vector;

#include <valarray>
using std::valarray;

#include <boost/shared_ptr.hpp>
using boost::shared_ptr;

template<int dim>
class integration_points_and_weights {
public:
  integration_points_and_weights( int n=1 );
  ~integration_points_and_weights();
  void evaluate_points( valarray<double>& points,
			const box<dim>& restricted_box,
			const box<dim>& local_box,
			const vector<shared_ptr<function<dim> > >& fn,
			vector<valarray<double> >& result,
			valarray<double>* temp=0 ) const;
private:
  integration_points_and_weights( const integration_points_and_weights& );
  integration_points_and_weights& operator= ( const integration_points_and_weights& );
public:
  valarray<double>& get_points();
  valarray<double>& get_weights();
protected:
  void set_number( int );

  valarray<double>& get_points();
  valarray<double>& get_weights();
private:
  int num_pts;
  valarray<double> points;
  valarray<double> weights;
};

#endif // _INTEGRATION_POINTS_AND_WEIGHTS_
