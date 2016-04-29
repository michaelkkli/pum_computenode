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

#ifndef _STROUSTRUP_HH_
#define _STROUSTRUP_HH_

template <class In, class Out, class Pred>
Out copy_if (In first, In last, Out res, Pred p)
{
	while (first != last) {
		if (p(*first)) *res++ = *first;
		++first;
	}
	return res;
}

#endif // _STROUSTRUP_HH_
