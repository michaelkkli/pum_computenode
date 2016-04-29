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

#include "indicator_function.hh"

#include <iostream>

using std::cout;

template<int dim>
indicator_function<dim>::indicator_function ( const string &                   r,
                                              const typename geometry<dim>::const_ptr & g ) : geom(g), region(r)
{
	assert( g );
	this->set_global_function();
}

template<int dim>
double
indicator_function<dim>::evaluate ( const double* co ) const
{
	// Keep: very useful when changing domains.
	//cout << "indicator_function asked to evaluate at point " << co[0] << ", " << co[1] << "\n";
	
	if ( geom->inside_region( region, co ) ) {
		return 1.0;
	} else {
		return 0.0;
	}
}

//

template class indicator_function<2>;
// template class indicator_function<3>;
