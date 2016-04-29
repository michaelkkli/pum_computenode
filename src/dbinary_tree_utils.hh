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

#ifndef _DBINARY_TREE_UTILS_HH_
#define _DBINARY_TREE_UTILS_HH_

#include "key_utils.hh"
#include "dbinary_tree.hh"
#include <cassert>
#include <functional>
#include <string>
using std::binary_function;
using std::string;

template<int dim>
class less_level : public binary_function<const string&, const string&, bool> {
public:
	bool operator() ( const string&, const string& ) const;
};

template<int dim>
class greater_level : public binary_function<const string&, const string&, bool> {
public:
	bool operator() ( const string&, const string& ) const;
};

template<int dim>
class is_sub_key : public binary_function<const string&, const string&, bool> {
public:
	bool operator() ( const string&, const string& ) const;
};

#endif // _DBINARY_TREE_UTILS_HH_
