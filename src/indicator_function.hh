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

#ifndef _INDICATOR_FUNCTION_HH_
#define _INDICATOR_FUNCTION_HH_

#include "function.hh"
#include "geometry.hh"
#include <boost/shared_ptr.hpp>
#include <string>
using boost::shared_ptr;
using std::string;

template <int dim>
class indicator_function : public function<dim> {
public:
	typedef shared_ptr<indicator_function<dim> >       ptr;
public:
	indicator_function( const string & region, const typename geometry<dim>::const_ptr & );
private:
	double evaluate ( const double* ) const;
private:
	string                               region;
	typename geometry<dim>::const_ptr    geom;
};

#endif // _INDICATOR_FUNCTION_HH_
