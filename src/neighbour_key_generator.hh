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

#ifndef _NEIGHBOUR_KEY_GENERATOR_HH_
#define _NEIGHBOUR_KEY_GENERATOR_HH_

#include <string>
#include <vector>

using std::string;
using std::vector;

/**
 * For level 7, update_neighbours takes 94% of run time.
 * 91% is for creating boxes with 57% projecting keys.
 */
template <int dim>
class neighbour_key_generator {
public:
	neighbour_key_generator ();
	~neighbour_key_generator ();
	
	void generate_neighbours ( const string & key, vector<string> & out ) const;
	
private:
	void generate_neighbours ( const string & key,
	                           const vector<string> & tails,
	                           vector<string> & out ) const;

};

#endif // _NEIGHBOUR_KEY_GENERATOR_HH_