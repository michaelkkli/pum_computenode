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

#ifndef _MY_CLASS_HH_

template <int dim>
class my_class {
public:
	my_class ();
	~my_class ();
private:
	my_class& my_class ( const my_class& );  // Not implemented.
	my_class& operator= ( const my_class& ); // Not implemented.
};

#endif // _MY_CLASS_HH_
