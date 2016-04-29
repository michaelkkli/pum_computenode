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

#include "function_utils.hh"

#include <cassert>

template<int dim>
void restricted_local_function_product ( valarray<double>&         in,
                                         const box<dim>&                 rl_bx,
                                         const box<dim>&                 first_bx,
                                         const function<dim>::const_ptr& first_fn,
                                         const box<dim>&                 second_bx,
                                         const function<dim>::const_ptr& second_fn,
                                         valarray<double>&               out )
{
	assert( in.size()>=dim );
	assert( in.size()%dim == 0 );
	
	assert( first_fn && second_fn );
	
	if ( !(first_fn && second_fn) ) {
		// We haven't been given two functions.
		out.clear();
		return;
	}
	
	first_fn->restricted_local_evaluate(in, rl_bx, first_bx, out);
	
	valarray<double> tmp;
	second_fn->restricted_local_evaluate(in, rl_bx, second_bx, out);
	
	assert( out.size() == tmp.size() );
	
	out *= tmp;
}

//

template void restricted_local_function_product<1>( valarray<double>&,
                                                    const box<1>&,
                                                    const box<1>&,
                                                    function<1>::const_ptr&,
                                                    const box<1>&,
                                                    function<1>::const_ptr&,
                                                    valarray<double>& );

template void restricted_local_function_product<2>( valarray<double>&,
                                                    const box<2>&,
                                                    const box<2>&,
                                                    function<2>::const_ptr&,
                                                    const box<2>&,
                                                    function<2>::const_ptr&,
                                                    valarray<double>& );

template void restricted_local_function_product<3>( valarray<double>&,
                                                    const box<3>&,
                                                    const box<3>&,
                                                    function<3>::const_ptr&,
                                                    const box<3>&,
                                                    function<3>::const_ptr&,
                                                    valarray<double>& );
