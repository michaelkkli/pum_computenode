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

#ifndef _SINGLE_STEP_HH_
#define _SINGLE_STEP_HH_

#include "pum_discretization.hh"
#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <valarray>
using boost::shared_ptr;
using std::map;
using std::string;
using std::valarray;

template<int dim>
class single_step {
public:
	single_step();
	~single_step();
private:
	single_step( const single_step<dim>& );                  // Not implemented.
	single_step<dim>& operator= ( const single_step<dim>& ); // Not implemented.
private:
	bool initialized;
	map<string,valarray<double> >        species_coefficients;
	shared_ptr<pum_discretization<dim> > pum_disc;
};

#endif // _SINGLE_STEP_HH_
