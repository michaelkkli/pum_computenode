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

#ifndef _CUBIC_SPLINE_WEIGHT_HH_
#define _CUBIC_SPLINE_WEIGHT_HH_

#include "differentiable_function.hh"
#include "function.hh"

template<int dim>
class cubic_spline_weight : public differentiable_function<dim> {
public:
	cubic_spline_weight();
	~cubic_spline_weight();
private:
	cubic_spline_weight( const cubic_spline_weight<dim>& );
	cubic_spline_weight<dim>& operator= ( const cubic_spline_weight<dim>& );
private:
	double evaluate ( const double* ) const;

	void evaluate ( valarray<double>& in, valarray<double>& out ) const;
	void evaluate ( valarray<double>& in, valarray<bool>& pred, valarray<double>& out ) const;

	void evaluate_grad ( const double*, double* ) const;

	void evaluate_grad ( valarray<double>& in, valarray<double>& out ) const;
	void evaluate_grad ( valarray<double>& in, valarray<bool>& pred, valarray<double>& out ) const;

	void evaluate_and_grad ( const double* co, double* value, double* grad ) const;
};

#endif // _CUBIC_SPLINE_WEIGHT_HH_
