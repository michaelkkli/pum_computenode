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

#include "dbinary_tree_utils.hh"

#include <iostream>

using std::cout;

template<int dim>
bool
less_level<dim>::operator() ( const string& one, const string& two ) const
{
	assert( one.size()>0 && two.size()>0 && one[0]==two[0] );
	return key_level<dim>(one) < key_level<dim>(two);
}

template<int dim>
bool
greater_level<dim>::operator() ( const string& one, const string& two ) const
{
	assert( one.size()>0 && two.size()>0 && one[0]==two[0] );
	return key_level<dim>(one) > key_level<dim>(two);
}

template<int dim>
bool
is_sub_key<dim>::operator() ( const string& one, const string& two ) const
{
	return two.find( one ) == 0;
}

//



template class less_level<1>;
template class less_level<2>;
template class less_level<3>;

template class greater_level<1>;
template class greater_level<2>;
template class greater_level<3>;

template class is_sub_key<1>;
template class is_sub_key<2>;
template class is_sub_key<3>;

