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

#ifndef _FUNCTION_UTILS_HH_
#define _FUNCTION_UTILS_HH_

#include "function.hh"
#include "box.hh"
#include <valarray>
using std::valarray;


template<int dim>
void restricted_local_function_product ( valarray<double>&                  in,
                                         const box<dim>&                          restricted_local,
                                         const box<dim>&                          first_bx,
                                         const typename function<dim>::const_ptr& first_fn,
                                         const box<dim>&                          second_bx,
                                         const typename function<dim>::const_ptr& second_fn,
                                         valarray<double>&                        out );


#endif // _FUNCTION_UTILS_HH_
