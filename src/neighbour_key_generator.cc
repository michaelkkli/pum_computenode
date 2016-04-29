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

#include "neighbour_key_generator.hh"

#include <algorithm>

using std::sort;

template <int dim>
neighbour_key_generator<dim>::neighbour_key_generator ()
{
}

template <int dim>
neighbour_key_generator<dim>::~neighbour_key_generator ()
{
}

template <>
void
neighbour_key_generator<1>::generate_neighbours ( const string & key, vector<string> & out ) const
{
	int key_size = key.size();
	
	// Key represents box of less than level 2.
	if ( key_size < 5 ) {
		if ( key == "" || key == "-1" ) {
			abort();
			return;
		}
		if ( key == "0" ) {
			out.resize(1);
			out[0] = "0";
			return;
		}
		// Assume key of level 1.
		if ( key_size == 3 ) {
			out.resize(2);
			out[0] = "0-0";
			out[1] = "0-1";
			return;
		}
	}
	
	vector<string> tails;
	tails.reserve ( 3 );
	
	char penultimate = key[ key_size - 3 ];
	char ultimate    = key[ key_size - 1 ];
	
	switch ( penultimate ) {
	case '0':
		tails.push_back( "0-0" );
		tails.push_back( "0-1" );
		switch ( ultimate ) {
		case '0':
			break;
		case '1':
			tails.push_back( "1-0" );
			break;
		default:
			abort();
			break;
		}
		break; // penultimate
	case '1':
		switch ( ultimate ) {
		case '0':
			tails.push_back( "0-1" );
			break;
		case '1':
			break;
		default:
			abort();
			break;
		}
		tails.push_back( "1-0" );
		tails.push_back( "1-1" );
		break;
	default:
		abort();
		break;
	}
	
	generate_neighbours( key, tails, out );
	
}

template <>
void
neighbour_key_generator<2>::generate_neighbours ( const string & key, vector<string> & out ) const
{
	int key_size = key.size();
	
	// Key represents box of less than level 2.
	if ( key_size < 5 ) {
		if ( key == "" || key == "-1" ) {
			abort();
			return;
		}
		if ( key == "0" ) {
			out.resize(1);
			out[0] = "0";
			return;
		}
		// Assume key of level 1.
		if ( key_size == 3 ) {
			out.resize(4);
			out[0] = "0-0";
			out[1] = "0-1";
			out[2] = "0-2";
			out[3] = "0-3";
			return;
		}
	}
	
	vector<string> tails;
	tails.reserve ( 9 );
	
	char penultimate = key[ key_size - 3 ];
	char ultimate    = key[ key_size - 1 ];
	
	switch ( penultimate ) {
	case '0':
		tails.push_back( "0-0" );
		tails.push_back( "0-1" );
		tails.push_back( "0-2" );
		tails.push_back( "0-3" );
		switch ( ultimate ) {
		case '0':
			break;
		case '1':
			tails.push_back( "1-0" );
			tails.push_back( "1-2" );
			break;
		case '2':
			tails.push_back( "2-0" );
			tails.push_back( "2-1" );
			break;
		case '3':
			tails.push_back( "1-0" );
			tails.push_back( "1-2" );
			tails.push_back( "2-0" );
			tails.push_back( "2-1" );
			tails.push_back( "3-0" );
		default:
			abort();
			break;
		}
		break;
	case '1':
		tails.push_back( "1-0" );
		tails.push_back( "1-1" );
		tails.push_back( "1-2" );
		tails.push_back( "1-3" );
		switch ( ultimate ) {
		case '0':
			tails.push_back( "0-1" );
			tails.push_back( "0-3" );
			break;
		case '1':
			break;
		case '2':
			tails.push_back( "0-1" );
			tails.push_back( "0-3" );
			tails.push_back( "2-1" );
			tails.push_back( "3-0" );
			tails.push_back( "3-1" );
			break;
		case '3':
			tails.push_back( "3-0" );
			tails.push_back( "3-1" );
		default:
			abort();
			break;
		}
		break;
	case '2':
		tails.push_back( "2-0" );
		tails.push_back( "2-1" );
		tails.push_back( "2-2" );
		tails.push_back( "2-3" );
		switch ( ultimate ) {
		case '0':
			tails.push_back( "0-2" );
			tails.push_back( "0-3" );
			break;
		case '1':
			tails.push_back( "0-2" );
			tails.push_back( "0-3" );
			tails.push_back( "1-2" );
			tails.push_back( "3-0" );
			tails.push_back( "3-2" );
			break;
		case '2':
			break;
		case '3':
			tails.push_back( "3-0" );
			tails.push_back( "3-2" );
			break;
		default:
			abort();
			break;
		}
		break;
	case '3':
		tails.push_back( "3-0" );
		tails.push_back( "3-1" );
		tails.push_back( "3-2" );
		tails.push_back( "3-3" );
		switch ( ultimate ) {
		case '0':
			tails.push_back( "0-3" );
			tails.push_back( "1-2" );
			tails.push_back( "1-3" );
			tails.push_back( "2-1" );
			tails.push_back( "2-3" );
			break;
		case '1':
			tails.push_back( "1-2" );
			tails.push_back( "1-3" );
			break;
		case '2':
			tails.push_back( "2-1" );
			tails.push_back( "2-3" );
			break;
		case '3':
			break;
		default:
			abort();
			break;
		}
		break;
	default:
		abort();
		break;
	}
	
	sort ( tails.begin(), tails.end() );
	
	generate_neighbours( key, tails, out );
}

template <int dim>
void
neighbour_key_generator<dim>::generate_neighbours ( const string &         key,
                                                    const vector<string> & tails,
                                                    vector<string> &       out ) const
{
	int tails_size = tails.size();
	int pos        = tails_size - 3;
	out.assign ( tails_size, key );
	
	for ( int i=0; i<tails_size; ++i ) {
		out[i].replace(pos, 3, tails[i]);
	}
}

//

template class neighbour_key_generator<1>;
template class neighbour_key_generator<2>;
