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

#ifndef _PUM_CONVERGENCE_HH_
#define _PUM_CONVERGENCE_HH_

#include "cover_structure.hh"
#include "differentiable_function.hh"
#include "line_segments.hh"
#include "refinement_structure.hh"
#include "solution.hh"
#include "vtk_output.hh"

#include <petscts.h>

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <string>
#include <valarray>

using boost::shared_ptr;
using std::cout;
using std::string;
using std::valarray;

struct species_system_data {
	double                          time_step_dt;
	valarray<PetscReal>             diffusion_coefficient;
	valarray<PetscReal>             spot_bleach_coefficient;
	valarray<PetscReal>             global_bleach_coefficient;
	vector<valarray<PetscReal> >    reaction_coefficients;
	vector<valarray<PetscReal> >    volume_adjust;
};

class circular_indicator_function : public function<2> {
public:
	circular_indicator_function ( double centre_x, double centre_y, double radius )
		: c_x_(centre_x), c_y_(centre_y), r_s_(radius*radius)
	{
		this->set_global_function();
	}
private:
	double evaluate ( const double* p ) const {
		if ( (c_x_-p[0])*(c_x_-p[0]) + (c_y_-p[1])*(c_y_-p[1]) < r_s_ ) {
			return 1.0;
		} else {
			return 0.0;
		}
	}
	double   c_x_;
	double   c_y_;
	double   r_s_;
};

// We allow restriction of the integration region at boundaries.
namespace {
	enum integral_type { integ_inside, integ_outside, integ_both };
};

template <int dim>
struct pum_details {
};

struct diagnostic_quad_details {
	vector<double>      all_interior_points;
	valarray<double>    interior_points;
};

struct pum_return_information {
	int   num_dof;
};

template <>
struct pum_details<2> {
	typedef shared_ptr<pum_details<2> >    ptr;
	line_segments                          domain;
	bool                                   anticlockwise_dom;
	
	basic_dbinary_tree<2>::ptr             dtree;
	refinement_structure<2>::ptr           ref_struct;
	
	// Optional - will use ref_struct if not set.
	refinement_structure<2>::ptr           viz_ref_struct;
	
	local_basis_type                       local_approx_space;
	
	differentiable_function<2>::ptr        fun_ptr;
	
	function<2>::ptr                       rhs_ptr;
	
	function<2>::ptr                       bleach_indicator_fun;
	differentiable_function<2>::ptr        equilibrium_fun;
	
	int                                    integration_decomposition;
	
	
	simple_quadrature_rule<2>                      interior_box_cubature;
	
	simple_quadrature_rule<1>                      boundary_base_rule;
#if 0
	// TODO: change to use quadrature_rule.
	valarray<double>                               box_cubature_points;
	valarray<double>                               box_cubature_weights;
#endif

};


struct convergence_results {
	valarray<int>       patch_to_num_dof;

	double    function_val_integral;
	double    function_val_squared_integral;
	double    function_grad_squared_integral;
	
	double    approx_val_integral;
	double    approx_val_function_val_integral;
	double    approx_grad_function_grad_integral;
		
	double    approx_val_squared_integral;
	double    approx_grad_squared_integral;
};

struct matrix_vector_particulars {
	double              pu_val_function_val_vec_vals;
	double              pu_val_function_val_squared_vec_vals;
	double              pu_val_function_grad_squared_vec_vals;
	
	valarray<double>    approx_val_rhs_fun_vec_vals;
	
	valarray<double>    approx_val_vec_vals;
	valarray<double>    approx_val_function_val_vec_vals;
	valarray<double>    approx_grad_function_grad_vec_vals;
	
	valarray<double>    mass_matrix_vals;
	valarray<double>    stiffness_matrix_vals;
	
	valarray<double>    bleach_matrix_vals;
	
	valarray<double>    siggia_matrix_vals;
};

template <int dim>
class pum_convergence {
public:
	typedef shared_ptr<pum_convergence<dim> >    ptr;
	pum_convergence ();
	~pum_convergence ();
public:
	// A scalar system has num_species = 0.
	// Species is only for a system of equations.
	// 2009-11-07 ML.
	void initialize ( typename pum_details<dim>::ptr &    details,
	                  pum_return_information *            pum_ret_info=0,
	                  diagnostic_quad_details * = 0,
	                  bool assemble_bleach_matrix = false,
	                  bool assemble_siggia_matrix = false,
	                  int  num_species = 0 );
	
	void draw_cover ( const string &, bool labels = true );
	
	void set_problem_L2_projection ();
	
	/**
		-\kappa \Delta u + u = f
		The stiffness matrix will be added with coefficient  -(-\kappa) = \kappa
	*/
	void set_problem_Helmholtz ( double kappa, double a );
	void set_problem_Poisson ( double kappa );
	PetscErrorCode solve_and_get_results ( convergence_results & );
	void output_vtk ( const string& solution, const string& error="", bool do_cell_scalars=false );
	void set_sum_all_pu_test ( int step=3 );
	
	PetscErrorCode load_mass_matrix ( double=1.0 );
	PetscErrorCode load_stiffness_matrix ( double=1.0 );
	
	PetscErrorCode add_scaled_stiffness_matrix ( double );
	PetscErrorCode add_scaled_bleach_matrix ( double );
	PetscErrorCode add_scaled_siggia_matrix ( double );
	
	PetscErrorCode solve ( bool output_iter );
	
	PetscErrorCode rhs_vec_mass_mult_previous_solution ();
	
	double integrate_solution_over_domain ();
	
	PetscErrorCode species_set_initial_proportions ( double* );
	
	PetscErrorCode species_set_matrices ( species_system_data& species_sys_dat,
	                                      bool                 spot_bleach_on,
	                                      valarray<bool>*      species_mask=0 );
	
	PetscErrorCode species_solve ( species_system_data& species_sys_dat );
	
	PetscErrorCode species_integrate_over_domain ( valarray<double> & );
	
	PetscErrorCode species_sum_output_vtk ( const string& solution, const string& error="", int sum_first_n_species=2 );
	
private:
	PetscErrorCode create_petsc_objects ( bool create_bleach_matrix=false,
	                                      bool create_siggia_matrix=false );
	PetscErrorCode destroy_petsc_objects ();
	void assemble_matrix_vector ( diagnostic_quad_details * = 0 );
	void assemble_matrix_vector_block ( int,
	                                    int,
	                                    integration_scratch& is,
	                                    matrix_vector_particulars& mvs,
	                                    diagnostic_quad_details * = 0 );
	
	void calculate_entries ( const box<dim>&                      restricted_local_box,
	                         vector<global_basis_function<dim> >& test,
	                         vector<global_basis_function<dim> >& trial,
	                         integration_scratch&                 is,
	                         matrix_vector_particulars&           mvs,
	                         bool                                 patch_indices_match/*,
	                         bool                                 this_is_box_interior*/ );
	
	PetscErrorCode set_matrix_entries ( int,
	                                    int,
	                                    matrix_vector_particulars& mvs );
	PetscErrorCode petsc_objects_all_assemble ();

	
	void get_dof_structure ( dof_structure::ptr& ) const;
	PetscErrorCode solve ( solution<dim>&, bool output_iter );
	
	
private:
	bool                               initialized_;
	
	typename pum_details<dim>::ptr     pum_d_;
	
	typename cover_structure<dim>::ptr   cov_struct_;
	
	// Currently must be a shared_ptr as the local spaces is set
	// on construction.
	typename global_approximation_space<dim>::ptr    global_approx_space_;

	valarray<int>                      patch_to_num_dof_;             /// (1)
	valarray<int>                      patch_to_start_index_;         /// (2)
	valarray<int>                      computenode_to_num_patches_;   /// (3)
	valarray<int>                      computenode_to_start_patch_;   /// (4)
	valarray<int>                      computenode_to_num_dof_;       /// (5)
	valarray<int>                      patch_indices_on_computenode_; /// (6)
	valarray<int>                      computenode_to_start_dof_;     /// (7)

	PetscMPIInt                        comm_rank_;
	PetscMPIInt                        comm_size_;
	
	bool                               petsc_objects_created_;
	bool                               ksp_matrix_populated_;
	bool                               ksp_operators_set_;
	
	//int                                integ_domain_decomposition_;
	
	// Length is number of patches.
	Vec                                pu_val_function_val_vec_;
	Vec                                pu_val_function_val_squared_vec_;
	Vec                                pu_val_function_grad_squared_vec_;
	
	// Length is size of PUM basis.
	Vec                                approx_val_vec_;
	Vec                                approx_val_function_val_vec_;
	Vec                                approx_grad_function_grad_vec_;
	
	Mat                                mass_matrix_;
	Mat                                stiffness_matrix_;

	KSP                                ksp_;

	Mat                                ksp_matrix_;
	Vec                                ksp_rhs_vec_;
	Vec                                ksp_solution_;
	
	Vec                                calculation_vec_;
	
	// Adding support for different species.
	// Please note that Mat and Vec are typedefs for
	// pointers to PETSc internal data structures.
	// 2009-11-01 ML.
	int                                num_species_;
	bool                               species_matrices_vectors_created_;
	vector<Mat>                        species_matrices_;
	vector<bool>                       species_matrices_populated_;
	vector<Vec>                        species_previous_vec_;
	vector<Vec>                        species_current_vec_;
	
	// Special optional matrices.
	Mat     bleach_matrix_;
	bool    bleach_matrix_created_;
	
	Mat     siggia_matrix_;
	bool    siggia_matrix_created_;
	
	vtk_output                         vtk_out_;
	vtk_output                         vtk_error_out_;
	
	typename solution<dim>::ptr        sol_;
	solution_evaluation_structure      sol_eval_struct_;
private:
	pum_convergence ( const pum_convergence & );           // Not implemented.
	pum_convergence& operator= ( const pum_convergence& ); // Not implemented.
};

#endif // _PUM_CONVERGENCE_HH_
