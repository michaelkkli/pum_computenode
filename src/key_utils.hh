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

#ifndef _KEY_UTILS_HH_
#define _KEY_UTILS_HH_

#include <string>
#include <vector>
using std::string;
using std::vector;

template <int dim>
struct key_traits {
  key_traits();
  ~key_traits();
  int          number_of_subdivisions;
  vector<char> subdivision_labels;
};

template <int dim>
void recursive_modify_keys ( int mod_level,
			     vector<string>::iterator it,
			     int number_to_operate_on );


/**
	Generate keys representing a complete refinement of the requested level.
*/
template <int dim>
void generate_keys ( int                 level,
                     vector<string> &    out );

template <int dim>
void generate_keys ( const string &      start,
                     int                 level,
                     vector<string> &    out );

/**
	Non-templated minimum-work function object class to
	order keys first by level and then naturally by string order.
*/
class levelwise_less {
	public:
		bool operator() ( const string &   a,
		                  const string &   b );
};

#endif // _KEY_UTILS_HH_
