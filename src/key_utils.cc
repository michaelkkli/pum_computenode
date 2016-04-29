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

#include "key_utils.hh"
#include <cassert>

template <int dim>
key_traits<dim>::key_traits() {
    number_of_subdivisions = 1;
    for ( int i=0; i<dim; ++i ) {
      number_of_subdivisions *= 2;
    }
    subdivision_labels.resize( number_of_subdivisions );
    subdivision_labels[0] = '0';
    for ( int i=1; i<number_of_subdivisions; ++i ) {
      subdivision_labels[i] = subdivision_labels[i-1] + 1;
    }
}

template <int dim>
key_traits<dim>::~key_traits()
{
}

template <int dim>
void recursive_modify_keys ( int mod_level,
			     vector<string>::iterator it,
			     int num )
{
  assert( mod_level > 0 );

  key_traits<dim> kt;

  int number_to_modify = num/kt.number_of_subdivisions;

	vector<vector<string>::iterator> start_points( kt.number_of_subdivisions );

  for ( int i=0; i<kt.number_of_subdivisions; ++i ) {
    start_points[i] = it + i*number_to_modify;
  }

  vector<string>::iterator mod_it;
  for ( int i=0; i<kt.number_of_subdivisions; ++i ) {
    mod_it = start_points[i];
    for ( int j=0; j<number_to_modify; ++j ) {
      mod_it->at( 2*mod_level ) = kt.subdivision_labels[i];
      ++mod_it;
    }
  }

  if ( number_to_modify > 1 ) {
    for ( int i=0; i<kt.number_of_subdivisions; ++i ) {
      recursive_modify_keys<dim>( mod_level+1, start_points[i], number_to_modify );
    }
  }
}

template <int dim>
void generate_keys ( int level, vector<string>& out )
{
  assert( level > 0 );

  if ( level == 0 ) {
    out.resize(1);
    out[0] = '0';
    return;
  }

  key_traits<dim> kt;

  string base_key( 2*level+1, '0');
  for ( int i=0; i<level; ++i ) {
    base_key.at( 1+2*i ) = '-';
  }
//  cout << "Base key is " << base_key << "\n";

  int number_of_keys = 1;

  for ( int i=0; i<level; ++i ) {
    number_of_keys *= kt.number_of_subdivisions;
  }

  out.assign( number_of_keys, base_key );

  recursive_modify_keys<dim>( 1, out.begin(), out.size() );

#if 0 //ndef NDEBUG // Keep for future debug.
  cout << "There are " << number_of_keys << " keys :\n";
  for ( int i=0; i<number_of_keys; ++i ) {
    cout << out[i] << "\n";
  }
#endif

}

template <int dim>
void generate_keys ( const string &      start,
                     int                 lev,
                     vector<string> &    out )
{
	assert ( start.size()>0 && (start.size()%2 == 1) );
	int inlevel = (start.size()-1)/2;

  if ( lev == 0 ) {
    out.resize(1);
    out[0] = start;
    return;
  }

  key_traits<dim> kt;

  int level = inlevel + lev;

  string base_key( 2*level+1, '0');
  for ( int i=0; i<level; ++i ) {
    base_key.at( 1+2*i ) = '-';
  }
  base_key.replace ( 0, start.size(), start );
//  cout << "Base key is " << base_key << "\n";

  int number_of_keys = 1;

  for ( int i=0; i<level; ++i ) {
    number_of_keys *= kt.number_of_subdivisions;
  }

  out.assign( number_of_keys, base_key );

  recursive_modify_keys<dim>( inlevel+1, out.begin(), out.size() );
}

bool
levelwise_less::operator() ( const string &   a,
                             const string &   b )
{
	size_t    asize = a.size();
	size_t    bsize = b.size();
	
	if ( asize < bsize ) {
		return true;
	} else if ( asize > bsize ) {
		return false;
	} else {
		return a < b;
	}
}

//

template struct key_traits<1>;
template struct key_traits<2>;
template struct key_traits<3>;

template void recursive_modify_keys<1>( int, vector<string>::iterator, int );
template void recursive_modify_keys<2>( int, vector<string>::iterator, int );
template void recursive_modify_keys<3>( int, vector<string>::iterator, int );

template void generate_keys<1>( int, vector<string>& );
template void generate_keys<2>( int, vector<string>& );
template void generate_keys<3>( int, vector<string>& );

template void generate_keys<1>( const string&, int, vector<string>& );
template void generate_keys<2>( const string&, int, vector<string>& );
template void generate_keys<3>( const string&, int, vector<string>& );
