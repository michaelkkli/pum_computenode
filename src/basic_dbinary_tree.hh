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

#ifndef _BASIC_DBINARY_TREE_HH_
#define _BASIC_DBINARY_TREE_HH_

#include "box.hh"

#include <boost/shared_ptr.hpp>

#include <map>
#include <set>
#include <string>
#include <vector>

using boost::shared_ptr;

using std::map;
using std::set;
using std::string;
using std::vector;

template <typename T>
void copy_vector_to_vec_ptr ( const vector<T>& in, typename T::vec_ptr& out )
{
	int size = in.size();
	out.resize( size );
	for ( int i=0; i<size; ++i ) {
		out[i].reset( new T( in[i] ) );
	}
}

template <int dim>
int
key_level ( const string& key );

#if 0
template <int dim>
string
project_key ( const string& key, int coordinate );
#endif

/**
Efficient projection of keys.
*/
template <int dim>
void project_key( const string& s, int coordinate, string& );

template <int dim>
void reconstruct_from_projected_keys ( const vector<string>& in, string & out );

template <int dim>
void reconstruct_from_projected_keys ( const vector<string> (&in)[dim], vector<string> & out );

/**
 * A grid of points is created from which it is possible to create
 * the boxes both for computation and for visualization.
 */

template <int dim>
class basic_dbinary_tree {
public:
	typedef shared_ptr<basic_dbinary_tree<dim> > ptr;
public:
	basic_dbinary_tree ();
	void initialize ( const box<dim>&, int maximum_level );
	
	/**
	 * A key of max_level_ will give the exact index of
	 * 1D interval. A key of lower level will give the
	 * lowest index interval that is covered by the long
	 * interval represented by the key.
	 */
	int lower_1d_interval ( const string&, int coord ) const;
	
	int num_intervals_covered ( const string& ) const;
	int num_intervals_covered ( int lev ) const;
	int upper_1d_interval ( const string&, int coord ) const;
	
	void get_key_1d_with_lower_end ( int level, int lower_end, string & out ) const;
	void get_key_1d_with_upper_end ( int level, int upper_end, string & out ) const;
	
	void get_key_1d_intersect_point ( int         level,
	                                  double      co,
	                                  int         dimension,
	                                  string &    out ) const;
	
	/**
	 * This won't be used nearly as much as the one considering expanded
	 * boxes.
	 */
	void get_key_intersect_point ( int               level,
	                               const double *    co,
	                               string &          out ) const;
	
	/**
	 * Open intersection as kcachegrind/valgrind/callgrind show strict
	 * inequalities are cheaper than allowing equality.
	 */
	void get_key_1d_expanded_intersect_point ( int                 level,
	                                           double              co,
	                                           int                 dimension,
	                                           double              expansion,
	                                           vector<string> &    out ) const;
	
	void get_key_expanded_intersect_point ( int                 level,
	                                        const double *      co,
	                                        double              expansion,
	                                        vector<string> &    out ) const;
	
	/**
	 * Assumes a 1d key passed in. Project before call.
	 */
	void get_1d_neighbours ( const string & key, vector<string> & out ) const;
	
	void get_neighbours ( const string & key, vector<string> & out ) const;
	void get_neighbours ( const string & key, set<string> & out ) const;
	void get_neighbours ( const set<string> & in, map<string, set<string> >& out ) const;
	
	/**
	 * This method is used to create boxes rapidly.
	 */
	void get_1d_endpoint_indices ( const string&, int coord, int& lo, int& up ) const;
	
	void get_endpoint_indices ( const string&, int* ) const;
	
	void get_box ( const string&, box<dim>& ) const;
	
	void get_box ( const string&, box<dim>&, double expansion ) const;
	
	void get_box ( const vector<string> &, vector<box<dim> > & out ) const;

	void get_centre_point ( const string&, double * ) const;

	void get_box ( const set<string>& in, vector<box<dim> >& res, double expansion ) const;

	
	// Unused.
	// void get_box ( const vector<string> &, vector<box<dim> > & out, double expansion ) const;
	
	bool box_intersect_point ( const string& k, const double* co, double expansion ) const;
	bool are_neighbours ( const string& first, const string& second, double expansion ) const;
	
	void get_grid_points ( valarray<double> & ) const;
	void get_grid_points_3D ( valarray<double> &, double init=0.0 ) const;
	int  get_num_grid_points () const;
	
	
	/**
	 * For VTK_PIXEL and VTK_VOXEL.
	 */
	void get_box_grid_indices_stacked ( const string&, int* ) const;
	
	/**
	 * For VTK_QUAD and HEXAHEDRON.
	 */
	void get_box_grid_indices_winding ( const string&, int* ) const;
	
	void get_box_grid_connectivity_offsets ( const vector<string>&,
	                                         valarray<int>& con,
	                                         valarray<int>& off,
	                                         int con_start=0,
	                                         int off_start=0 );

	valarray<double>* access_endpoints_grid ();
private:
	box<dim>         parent_box_;
	int              max_level_;
	int              num_1d_intervals_;
	int              num_1d_endpoints_;
	valarray<int>    succ_halving_num_1d_intervals_;
	valarray<double> all_endpoints_[dim];
	valarray<int>    num_1d_endpoints_pows_;
	valarray<double> level_to_interval_length_[dim];
};

#endif
