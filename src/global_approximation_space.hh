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

#ifndef _GLOBAL_APPROXIMATION_SPACE_HH_
#define _GLOBAL_APPROXIMATION_SPACE_HH_

#include "cover_structure.hh"
#include "dbinary_tree.hh"
#include "function.hh"
#include "global_basis_function.hh"
#include "solution_evaluation_structure.hh"

#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <vector>

using boost::shared_ptr;
using std::map;
using std::string;
using std::vector;

enum local_basis_type { monomial1, monomial2, monomial3 };

/**
	Class representing the global approximation space.
	Notice that this class has no concept of integration.
*/
template<int dim>
class global_approximation_space {
public:
	typedef shared_ptr<const global_approximation_space<dim> > const_ptr;
	typedef shared_ptr<global_approximation_space<dim> >       ptr;
public:
	global_approximation_space( local_basis_type = monomial1 );
	~global_approximation_space();
public:
#if 0
  bool is_updated () const;
#endif

	void give_dbinary_tree ( typename dbinary_tree<dim>::ptr & );
	
	void set_support ( const string & );
	
	/**
		Build the cached_pu_functions_ in here.
	*/
	void give_cover_structure ( typename cover_structure<dim>::ptr& );

	/**
		Get number of degrees of freedom per patch.
	*/
	void get_patch_to_num_dof ( valarray<int>& ) const;
	int est_patch_num_dof () const;

	double get_patch_start_index ( int patch_index ) const;

	//	void get_global_basis_by_index ( int, global_basis_function<dim>& ) const;
	void get_global_basis_by_patch ( int, vector<global_basis_function<dim> >& ) const;
	
	double global_evaluate ( const double* co, const string & key ) const;

	/**
		Would require a quick way of finding a patch containing this point.
		Implies use of dbinary_tree.
	*/
	void global_evaluate ( const double* co ) const;
	
	/**
	 * Both keys and indices are retrieved to avoid a double lookup.
	 * global_approximation function knows about the cover_structure
	 * and solution knows about the dof_structure but both pieces are needed
	 * for successful evaluation. Information is passed back and
	 * forth in preference to a higher degree of coupling.
	 */
	void get_incident_patch_keys_and_indices ( const double*, vector<string>&, vector<int>& ) const;
	
	void get_solution_evaluation_structure ( valarray<double> &, solution_evaluation_structure & ) const;
	
	/**
	 * Allow evaluation given the patch keys for the incident patches.
	 * Intended to be called by solution as that is where the coefficients
	 * are known.
	 */
	double global_evaluate ( const double*,
	                         const vector< string > &            incident_patch_keys,
	                         const vector< const double* > &    coeff_starts ) const;
	
  /**
	Number of functions in global basis.
  */
	//  int size() const;
  
#if 0
  /**
	Check if two global basis functions have overlapping support.
  */
  bool overlapping_support( int, int ) const;
public:
  void set_dbinary_tree ( typename dbinary_tree<dim>::ptr & );
  void set_support ( const string& );
  void update ();
#endif
private:
	string                                                      support_;
	typename differentiable_function<dim>::ptr                  default_weight_;
	typename differentiable_function<dim>::vec_ptr              default_local_basis_;

	typename dbinary_tree<dim>::ptr                             dtree_;
	
	typename cover_structure<dim>::ptr                          cover_struct_;

	map<string, typename partition_of_unity_function<dim>::ptr> cached_pu_functions_;
	map<string, typename differentiable_function<dim>::vec_ptr> local_approx_functions_;
};

#endif // _GLOBAL_APPROXIMATION_SPACE_HH_
