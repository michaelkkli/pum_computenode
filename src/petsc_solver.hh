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

#ifndef _PETSC_SOLVER_HH_
#define _PETSC_SOLVER_HH_

#include "dof_structure.hh"
#include "function.hh"
#include "pum_discretization.hh"
#include "solution.hh"

#include <petscts.h>

#include <boost/shared_ptr.hpp>

#include <map>
#include <string>
#include <valarray>
#include <vector>

using boost::shared_ptr;

using std::map;
using std::slice;
using std::string;
using std::valarray;
using std::vector;

template <int dim>
class solution;

// Only pum_discretization knows how to draw things.
// We'll keep petsc_solver as dumb as possible.

template<int dim>
class petsc_solver {
public:
	typedef shared_ptr<petsc_solver<dim> > ptr;
public:
	petsc_solver();
	~petsc_solver();
public:
	/**
		Index distribution gives the number of degrees of freedom per patch.
		
		The patch_indices_on_computenode are returned by reference so we don't
		need to store information on the index_distribution and avoid need to
		synchronize with a member function that returns the indices.

		A sudden jump in execution time from level 7 to 8 for basis order 1 is caused
		by reallocation of memory by PETSc. The same thing happens going from level 5 to 6
		with basis order 2. The row_nonzero must be adjusted accordingly.
	*/
	PetscErrorCode create_petsc_objects ( const valarray<int>& patch_to_num_dof, int row_nonzero=27 );
	PetscErrorCode destroy_petsc_objects ();
	
	void get_patch_indices_on_computenode ( valarray<int>& ) const;

	PetscErrorCode set_matrix_entries ( int first_patch_index,
	                                    int second_patch_index,
	                                    valarray<double>& vals );
	
	PetscErrorCode set_matrix_entries ( int first_patch_index,
	                                    int second_patch_index,
	                                    valarray<double> & mass_matrix_vals,
	                                    valarray<double> & bleach_matrix_vals,
	                                    valarray<double> & stiffness_matrix_vals,
	                                    valarray<double> * siggia_matrix_vals=0 );
	
	PetscErrorCode matrix_assembly_begin ();
	PetscErrorCode matrix_assembly_end   ();
	
#if 0 // Appears unused - will not update for siggia_matrix. 2008-12-04 Michael LI.
	PetscErrorCode quick_setup_matrices_after_first_assembly ();
#endif
	
	PetscErrorCode quick_assemble_matrices ();
	
	PetscErrorCode load_mass_matrix ( double=1.0 );
	PetscErrorCode add_scaled_stiffness_matrix ( double );
	PetscErrorCode add_scaled_bleach_matrix ( double );
	
	PetscErrorCode add_scaled_siggia_matrix ( double );
	
	PetscErrorCode set_rhs_vector_entries ( int patch_index, valarray<PetscScalar>& vals );
	PetscErrorCode rhs_vector_assembly_begin ();
	PetscErrorCode rhs_vector_assembly_end   ();
	
	PetscErrorCode save_as_mass_matrix ();
	PetscErrorCode save_as_stiffness_matrix ();
	
	PetscErrorCode save_solution1 ();
	PetscErrorCode save_domain_integration_vec ();
	PetscErrorCode save_bleach_integration_vec ();
	
	PetscErrorCode view_petsc_matrix ( const char* = 0 );
	PetscErrorCode view_mass_matrix ( const char* = 0 );
	PetscErrorCode view_stiffness_matrix ( const char* = 0 );
	PetscErrorCode view_bleach_matrix ( const char* = 0 );
	
	PetscErrorCode view_solution ( const char* = 0 );
	
	double integrate_solution_over_domain ();
	double integrate_solution_over_bleach_region ();
	
	PetscErrorCode scale_current_and_add_mass_matrix ( double );
	
	PetscErrorCode zero_rhs_vec ();
	PetscErrorCode zero_petsc_solution ();
	
	PetscErrorCode rhs_vec_add_mass_matrix_solution1 ();
	PetscErrorCode rhs_vec_mass_mult_previous_solution ();
	
	void get_dof_structure ( dof_structure::ptr& ) const;
	
	PetscErrorCode solve ( solution<dim>&, bool output_iter=true );

private:
	PetscMPIInt      comm_rank_;
	PetscMPIInt      comm_size_;
	bool             petsc_objects_created_;
	
	bool             mass_matrix_created_;
	bool             bleach_matrix_created_;
	bool             stiffness_matrix_created_;
	bool             siggia_matrix_created_;
	
	bool             solution1_created_;
	bool             domain_integration_vec_created_;
	bool             bleach_integration_vec_created_;

	bool             ksp_operators_set_;
	
	valarray<int>    patch_to_num_dof_;             /// (1)
	valarray<int>    patch_to_start_index_;         /// (2)
	valarray<int>    computenode_to_num_patches_;   /// (3)
	valarray<int>    computenode_to_start_patch_;   /// (4)
	valarray<int>    computenode_to_num_dof_;       /// (5)
	valarray<int>    patch_indices_on_computenode_; /// (6)
	valarray<int>    computenode_to_start_dof_;     /// (7)

private:
	Vec              petsc_solution_;
	Mat              petsc_matrix_;
	
	Mat              mass_matrix_;
	Mat              stiffness_matrix_;
	Mat              bleach_matrix_;
	Mat              siggia_matrix_;
	
	Vec              solution1_;
	
	Vec              domain_integration_vec_;
	Vec              bleach_integration_vec_;
	
	Vec              rhs_vec_;
	KSP              ksp_;
private:
	petsc_solver<dim> ( const petsc_solver<dim>& );             // Not implemented.
	petsc_solver<dim>& operator= ( const petsc_solver<dim>& ); // Not implemented.
};


/*
	Removed from class.
*/
#if 0
public:
	void set_pum_discretization( shared_ptr<pum_discretization<dim> >& );
	void update();
	void assemble_for_initial_projection();
	void solve_initial_projection();

	void species_set_function( const string& species, const pfunction& species_func );
private:
	bool             have_discretization;
	bool             updated;
	vector<string>   patch_keys_on_computenode;
	map<string, const_pfunction> species_functions;
	shared_ptr<pum_discretization<dim> > dis;
#endif

#endif // _PETSC_SOLVER_HH_
