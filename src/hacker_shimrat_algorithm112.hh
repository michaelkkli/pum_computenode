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

#ifndef _HACKER_SHIMRAT_ALGORITHM112_HH_
#define _HACKER_SHIMRAT_ALGORITHM112_HH_

class line_segments;

/**
 * Point in polygon appears insensitive to boundary orientation.
 * 2009-07-25 ML.
 */

/**
	Algorithms H.J. Wegstein, Editor

	Algorithm 112 Position of point relative to polygon.
	M. Shimrat, University of Alberta, Calgary, Alberta, Canada.
	Comm. ACM, Aug. 1962.

	Certification of Algorithm 112 Position of point relative to polygon.
	Richard Hacker, The Boeing Co., Seattle, Wash.
	Comm. ACM, 5:606, 1962.
*/
bool point_in_polygon ( int n, double* x, double* y, double x0, double y0 );


/**
	Efficient version for line_segments data structure.
*/
bool point_in_polygon ( line_segments&, double x0, double y0 );


#endif // _HACKER_SHIMRAT_ALGORITHM112_HH_
