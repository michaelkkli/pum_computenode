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

#ifndef _DBINARY_TREE_HH_
#define _DBINARY_TREE_HH_

#include "basic_dbinary_tree.hh"
#include "box.hh"
#include "cover_structure.hh"
#include "geometry.hh"

#include <boost/shared_ptr.hpp>

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <vector>
#include <cassert>

using boost::shared_ptr;

using std::map;
using std::ostream;
using std::set;
using std::stack;
using std::string;
using std::vector;

#if 0
template<int dim>
int key_level( const string& );
#endif
			 
#if 0 // Moved to basic_dbinary_tree. Delete. 2009-02-11 ML.
template <typename T>
void copy_vector_to_vec_ptr ( const vector<T>& in, typename T::vec_ptr& out )
{
	int size = in.size();
	out.resize( size );
	for ( int i=0; i<size; ++i ) {
		out[i].reset( new T( in[i] ) );
	}
}
#endif

/**
   The keys of a dbinary_tree are the most important. In 1D they
   are of the form "0-0-0-1" where 0 indicates left and 1 indicates
   a right interval.
 */
template<int dim>
class dbinary_tree {
public:
	typedef shared_ptr<const dbinary_tree<dim> > const_ptr;
	typedef shared_ptr<dbinary_tree<dim> >       ptr;
public:
  dbinary_tree ();
  ~dbinary_tree ();
public:
  void set_geometry ( shared_ptr<geometry<dim> >&, double expand_bounding_box=1.0 );
  void set_parent_box ( const box<dim>& par_box );
  void set_max_level ( int new_max_level );
  void update_neighbours();
  void update();
  bool is_updated() const;
  void set_initial_refinement( int level );
  void set_cover_factor( double cf );
  void dump_cout( int n = 20 ) const;
  void dump_gp_template_boxes( ostream& out );
  void dump_gp( ostream& out ) const;
  void gp_draw_parent_box( ostream& out, double level=0.0 ) const;
  void gp_draw_boundary_cover ( const string& boundary, ostream& out, double level ) const;
  void gp_draw_region_interior_cover ( const string& region, ostream& out, double level ) const;
  void gp_draw_region_boundary_cover ( const string& region, ostream& out, double level ) const;
  // Do not provide is interior key in general as the list of things
  // to check will often be much larger for a large number of small patches.
  // The exception will be for highly convoluted boundaries.
  bool region_key_is_boundary_key ( const string& region, const string& key ) const;
public: // Const access member functions.
  const geometry<dim>& access_geometry () const;
  const set<string>& access_region_boundary_keys( const string& region ) const;
  const set<string>& access_region_interior_keys( const string& region ) const;
  const set<string>& access_region_all_keys( const string& region ) const;
  const set<string>& access_boundary_keys ( const string& boundary ) const;
  void get_patch_keys( const string& region_or_boundary, set<string>& res ) const;
	void get_neighbour_keys( const string& region_or_boundary, map<string, set<string> >& ) const;
  const set<string>& access_patch_keys( const string& region_or_boundary ) const;
  int get_num_patches( const string& region_or_boundary ) const;
  bool support_has_key( const string& region_or_boundary, const string& key ) const;
  void get_patches( const set<string>& in, vector<box<dim> >& res ) const;
  void get_region_patches( const string& region, vector<box<dim> >& res ) const;
  void get_region_boundary_patches ( const string& region, vector<box<dim> >& out ) const;
  void get_region_interior_patches ( const string& region, vector<box<dim> >& out ) const;
  void get_boundary_patches( const string& boundary, vector<box<dim> >& res ) const;
  const box<dim>& access_parent_box() const;
	
  string find_patch_containing_point( const string& region, const double* co ) const;
	
	// Assumes all patches still of initial size.
	void    get_patch_containing_point ( const double* co, vector<string> & ) const;
	
  void get_boundary_neighbours( const string& boundary, const string& key, set<string>& out ) const;
  void get_region_neighbours( const string& region, const string& key, set<string>& res ) const;
  void get_region_neighbours_incident_point ( const string& region, const string& key, const double* co, set<string>& res ) const;
	
 	void get_global_box( const string&, box<dim>& ) const;
	void get_global_box( const vector<string> & in, vector<box<dim> >& out ) const;
#if 0 // Should not be necessary and horribly inefficient if used.
	void get_patch_keys_incident_point ( const string&      region,
	                                     const double*      co,
	                                     vector<string>&    out ) const;
#endif
	
	/// Used to get nice region for visualization.
	void get_intersect_region ( const string&, vector<string> & in, vector<string> & out );


	// Caching is inefficient.
	// const box<dim>& access_global_box( const string& ) const;
	
	void get_cover_structure ( const string& support_name, typename cover_structure<dim>::ptr& );
private:
#if 0
  void project_key( const string& s, int coordinate, string& ) const;
  void get_endpoints_for_projected_key ( string s, double& e0, double& e1 ) const;
	  void create_global_box ( const string& key, box<dim>& out ) const;
#endif
  double centre_for_projected_key ( const string& ) const;

//  void get_box ( string s, box<dim>& ) const;
  void split_dbinary_node_and_push ( string, stack<string>& );
	// int level ( string& s ) const;

// dbinary_tree is DEPRECATED so we make certain members public to facili
public:
	int              initial_level;
	int                    max_level;
private:
  bool                   updated_neighbours;
  bool                   updated;
  bool                   updated_global_boxes;
  bool                   have_parent_box;
  box<dim>               parent_box;
  double                 cover_factor;
  vector<box<dim> >      template_boxes;
	
	
	// Removing structure as profiling has shown its use to
	// be inefficient for large covers.
	//  map<string, box<dim> > global_boxes;
	
	
  valarray<double>       powers_of_one_half;
  valarray<int>          powers_of_two;
  shared_ptr<geometry<dim> > geom;
	typename basic_dbinary_tree<dim>::ptr basic_dtree;
	typename cover_structure<dim>::ptr    cover_struct_;
private:
  set<string>               present_tree;
  map<string, set<string> > present_neighbours;
  map<string, set<string> > present_boundary_keys; // Accessed by boundary name.
  map<string, set<string> > present_interior_keys; // Accessed by oriented boundary name.

  // We keep boundary and interior keys separate as
  // they are handled quite differently and we don't
  // want to mix them together only to separate them
  // when we need to use them.
  map<string, set<string> > present_region_boundary_keys;
  map<string, set<string> > present_region_interior_keys;
  map<string, set<string> > present_region_all_keys; // Accessed by region.
private:
  dbinary_tree ( const dbinary_tree<dim>& );                 // Not implemented.
  dbinary_tree<dim>& operator= ( const dbinary_tree<dim>& ); // Not implemented.
};

#endif // _DBINARY_TREE_HH_
