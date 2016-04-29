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

#ifndef _PUM_APPROX_HH_
#define _PUM_APPROX_HH_

#include "function.hh"
#include "box.hh"



#include <vector>
using std::vector;

#include <valarray>
using std::valarray;

#include <string>
using std::string;
	 
#include <map>
using std::map;

#include <set>
using std::set;

#include <stack>
using std::stack;

template <int dim, typename key_t=string>
class pum_approx {
public:
  typedef key_t                             key_type;
  typedef typename function<dim>::const_pfunction const_pfunction;
public:
  pum_approx();
  ~pum_approx();
private:
  pum_approx( const pum_approx<dim>& );            // Not implemented.
  pum_approx<dim>& operator=( const pum_approx<dim>& ); // Not implemented.
public:
  pum_approx<dim, key_t>&                 add_patch( const key_type&, const box<dim>& );
  pum_approx<dim, key_t>&          set_patch_weight( const key_type&, const_pfunction& );
  pum_approx<dim, key_t>& set_local_approx_function( const key_type&, const vector<const_pfunction>& );
  pum_approx<dim, key_t>&              remove_patch( const key_type& );
  int num_local_coefficients_needed () const;
  void form_global_mass_matrix(  ) const;
  void give_local_coefficients ( const valarray<double>& );
  double seq_global_evaluate ( const double* ) const;
  void output_gnuplot ( const char*,
			int points_per_dimension=100,
			bool on_boundary=true ) const;
  double integrate_function ( box<dim>&, const_pfunction& ) const;
public: // const so these can be re-entrant.
  void get_neighbours( const key_type&, set<key_type>& ) const;
  void get_patches_and_weights( const vector<key_type>&,
                                vector<box<dim> >&,
                                vector<const_pfunction>& ) const;
  void get_local_approx_functions( const key_type&, vector<const_pfunction>& ) const;
  double get_coefficient( const key_type&, int ) const;
private:
  void insert_petsc_block( const key_type&, const key_type&, int, int, double* ) const;
  key_type find_patch_containing_point ( const double*,
				    bool lowest_index=true ) const;

  // Evaluate local approximation within a particular patch.
  double evaluate_local_approx ( const key_type&, const double* ) const;

  void register_change();

  // It may be possible to keep all knowledge of computer cluster
  // architecture within update().
  void update();
  void update_neighbours();


private:
  box<dim> bounding_box;
  bool     updated;
  bool     coefficients_provided;
  
private:
  set<key_type>                     global_keys;
  set<key_type>                     local_keys;
  map<key_type, box<dim> >          patches;
  map<key_type, const_pfunction>          patch_weights;
  map<key_type, vector<const_pfunction> > local_approx_functions;
  map<key_type, set<key_type> >     neighbours;
  map<key_type, int>                global_index_start;
  map<key_type, int>                local_index_start;
private:
  int              num_local_patches;
  int              local_matrix_size;
  int              global_matrix_size;
  valarray<double> local_coefficients;
};

#endif // _PUM_APPROX_HH_
