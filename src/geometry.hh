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

#ifndef _GEOMETRY_HH_
#define _GEOMETRY_HH_

#include "boundary.hh"
#include "box.hh"
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <valarray>
#include <vector>
using boost::shared_ptr;
using std::map;
using std::ostream;
using std::set;
using std::string;
using std::valarray;
using std::vector;

template<int dim>
class geometry {
public:
  typedef shared_ptr<const geometry<dim> > const_ptr;
  typedef shared_ptr<geometry<dim> >       ptr;
public:
	geometry();
	~geometry();
public:
	void create_boundary ( const string& name, valarray<double>& );
	void create_region ( const string& name, const string& oriented_boundary = "" );
	const map<string, boundary<dim> >& access_boundaries () const;
	const map<string, vector<string> >& access_regions () const;
	void get_bounding_box ( box<dim>& );
	void update(); // geometry doesn't itself need this but boundary does.
	bool inside_region ( const string& region, const double* co ) const;
	void intersect_region ( const string & region, const vector<box<dim> > &, valarray<bool>& ) const;
	bool is_boundary ( const string& name ) const;
	bool is_region ( const string& name ) const;
	bool is_support ( const string& name ) const;
	bool gp_draw ( const string& region_or_boundary, ostream& out, double level=0.0 ) const;
private:
	set<string>                  created;
	map<string, boundary<dim> >  boundaries;
	map<string, vector<string> > regions;
private:
  geometry ( const geometry<dim>& );                 // Not implemented.
  geometry<dim>& operator= ( const geometry<dim>& ); // Not implemented.
};

#endif // _GEOMETRY_HH_
