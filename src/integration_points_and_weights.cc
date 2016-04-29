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

#include "integration_points_and_weights.hh"

template<int dim>
integration_points_and_weights<dim>::integration_points_and_weights( int n ) :
  num_pts( n >= 1 ? n : 1 ),
  points( n >= 1 ? dim*n : 1 ),
  weights( n >= 1 ? n : 1 )
{
  if ( n == 1 ) {
    points[0] = 0.0;
    if ( dim > 1 ) {
      points[1] = 0.0;
    }
    if ( dim > 2 ) {
      points[2] = 0.0;
    }
    weights[0] = 2.0;
  }
}

template<int dim>
integration_points_and_weights<dim>::~integration_points_and_weights() { }

template<int dim>
void
integration_points_and_weights<dim>::evaluate_points( valarray<double>& points,
						      const box<dim>& restricted_box,
						      const box<dim>& local_box,
						      const vector<shared_ptr<function<dim> > >& fn,
						      vector<valarray<double> >& result,
						      valarray<double>* temp ) const
{
  valarray<double> va;
  valarray<double>* pva = &va;
  if ( temp != 0 ) {
    pva = temp;
  }
  valarray<double>& eval_pts = *pva;

  local_to_restricted_local<dim>( points,
				  restricted_box,
				  local_box,
				  eval_pts );

  int fsize = fn.size();
  if ( fsize == 0 ) {
    return;
  }
  result.resize( fsize );

  for ( int i=0; i<fsize; ++i ) {
    fn[i]->local_evaluate( local_box, eval_pts, result[i] );
  }
}

template<int dim>
const
valarray<double>& integration_points_and_weights<dim>::get_points() const
{
  return this->points;
}

template<int dim>
const
valarray<double>& integration_points_and_weights<dim>::get_weights() const
{
  return this->weights;
}

template<int dim>
valarray<double>& integration_points_and_weights<dim>::get_points()
{
  return this->points;
}

template<int dim>
valarray<double>& integration_points_and_weights<dim>::get_weights()
{
  return this->weights;
}

//

template class integration_points_and_weights<1>;
template class integration_points_and_weights<2>;
template class integration_points_and_weights<3>;
