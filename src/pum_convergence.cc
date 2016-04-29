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

#include "basic_dbinary_tree_utils.hh"
#include "box_utils.hh"
#include "pum_convergence.hh"
#include "solution_evaluation_structure.hh"

//#include <valgrind/callgrind.h>

#include <petscksp.h>

#include <algorithm>
#include <cmath>
#include <iostream>

using std::clog;
using std::copy;
using std::fabs;

template <int dim>
pum_convergence<dim>::pum_convergence() :
initialized_(false),
petsc_objects_created_(false),
ksp_matrix_populated_(false),
ksp_operators_set_(false),
bleach_matrix_created_(false),
siggia_matrix_created_(false),
num_species_(0),
species_matrices_vectors_created_(false)
{
	MPI_Comm_rank( PETSC_COMM_WORLD, &comm_rank_ );
	MPI_Comm_size( PETSC_COMM_WORLD, &comm_size_ );
}

template <int dim>
pum_convergence<dim>::~pum_convergence() {
	destroy_petsc_objects();
}

template <>
void
pum_convergence<2>::initialize ( pum_details<2>::ptr &        details,
                                 pum_return_information *     pum_return_info,
                                 diagnostic_quad_details *    diag_quad_d,
                                 bool                         assemble_bleach_matrix,
                                 bool                         assemble_siggia_matrix,
                                 int                          num_species )
{
	num_species_ = num_species;
	
	if ( assemble_bleach_matrix ) {
		assert ( details->bleach_indicator_fun );
	}
	if ( assemble_siggia_matrix ) {
		assert ( details->equilibrium_fun );
	}
	
	initialized_ = false;
	
//	if ( petsc_objects_created_ ) {
//		destroy_petsc_objects();
//	}

	PetscLogDouble log_t0, log_t1;
	
	assert ( details );
	pum_d_ = details;
	
	bool viz_ref_same_as_sim_ref = true;
	if ( pum_d_->viz_ref_struct ) {
		viz_ref_same_as_sim_ref = false;
	} else {
		pum_d_->viz_ref_struct = pum_d_->ref_struct;
	}
	

	cov_struct_.reset ( new cover_structure<2> );
	cov_struct_->ref_structure = pum_d_->ref_struct;
	
	map<string, set<string> > &    nbrs = cov_struct_->patch_neighbour_keys;

	vector<string> tree_all_keys;
	
	generate_keys_and_neighbours< 2 >( *(cov_struct_->ref_structure), tree_all_keys, nbrs );

	//	PetscPrintf ( PETSC_COMM_WORLD, "Early abort for massif level 10 test.\n" );
	//	exit(0);

	
	set<string>    inside_poss_bdry_keys;
	
	// Centres lie inside. 2009-07-06 ML.
	// Point in polygon appears insensitive to boundary orientation
	// hence no reversal if clockwise winding. 2009-07-25 ML.
	get_inside_keys ( tree_all_keys,
	                  *(cov_struct_->ref_structure->dtree),
	                  pum_d_->domain,
	                  inside_poss_bdry_keys );

	tree_all_keys.clear();
	
	set<string>        bdry_keys;
	vector<string>    vec_bdry_keys;
	vector<int>       entry_seg;
	
	vector<string>    viz_vec_bdry_keys;
	vector<int>       viz_entry_seg;
	
	if ( viz_ref_same_as_sim_ref ) {
		// Take boundary keys made for visualization.
		make_boundary_vtk_output ( pum_d_->domain,
		                           *(pum_d_->ref_struct),
		                           vtk_out_,
		                           pum_d_->anticlockwise_dom ? inside : outside, &vec_bdry_keys, &entry_seg );
	} else {
		assert ( pum_d_->viz_ref_struct );
		make_boundary_vtk_output ( pum_d_->domain,
		                           *(pum_d_->viz_ref_struct),
		                           vtk_out_,
		                           pum_d_->anticlockwise_dom ? inside : outside, &viz_vec_bdry_keys, &viz_entry_seg );
	

		
		// We need to form correct boundary keys for the simulation - the visualization
		// boundary keys are inappropriate.
		get_intersect_keys_entry_indices ( *(pum_d_->ref_struct), pum_d_->domain, vec_bdry_keys, entry_seg );
	}
	
	copy ( vec_bdry_keys.begin(), vec_bdry_keys.end(),
	       inserter( bdry_keys, bdry_keys.begin() ) );

	
	{
		set<string>    strict_inside_keys;
		vector<string>    strict_inside_keys_vec;
		
		if ( viz_ref_same_as_sim_ref ) {
		set_difference ( inside_poss_bdry_keys.begin(),
				inside_poss_bdry_keys.end(),
				bdry_keys.begin(),
				bdry_keys.end(),
				inserter( strict_inside_keys,strict_inside_keys.begin() ) );
		
			strict_inside_keys_vec.resize ( strict_inside_keys.size() );
			
		copy ( strict_inside_keys.begin(),
                       strict_inside_keys.end(),
                       strict_inside_keys_vec.begin() );
		} else {
			vector<string> viz_tree_all_keys;

			generate_keys< 2 >( *(pum_d_->viz_ref_struct), viz_tree_all_keys );
			
			set<string>    viz_inside_poss_bdry_keys;
			
			get_inside_keys ( viz_tree_all_keys,
			                  *(pum_d_->viz_ref_struct->dtree),
			                  pum_d_->domain,
			                  viz_inside_poss_bdry_keys );
			
			set<string>    viz_bdry_keys;
			
			copy ( viz_vec_bdry_keys.begin(), viz_vec_bdry_keys.end(),
			       inserter ( viz_bdry_keys, viz_bdry_keys.begin() ) );
			
			set_difference ( viz_inside_poss_bdry_keys.begin(),
			                 viz_inside_poss_bdry_keys.end(),
			                 viz_bdry_keys.begin(),
			                 viz_bdry_keys.end(),
			                 inserter ( strict_inside_keys, strict_inside_keys.begin() ) );
			
			strict_inside_keys_vec.resize ( strict_inside_keys.size() );
			
			copy ( strict_inside_keys.begin(),
			       strict_inside_keys.end(),
			       strict_inside_keys_vec.begin() );
		}

		valarray<int>    boundary_connectivity ( vtk_out_.connectivity );
		valarray<int>    boundary_offsets ( vtk_out_.offsets );
		
		valarray<int>    strict_inside_connectivity;
		valarray<int>    strict_inside_offsets;
			
		if ( strict_inside_keys_vec.size() > 0 ) {
			pum_d_->ref_struct->dtree->get_box_grid_connectivity_offsets ( strict_inside_keys_vec,
			                                                               strict_inside_connectivity,
			                                                               strict_inside_offsets );

		
			// No modifications need to be made to connectivity as numbering is consistent.
			vtk_out_.connectivity.resize (
			strict_inside_connectivity.size() + boundary_connectivity.size() );
			vtk_out_.connectivity[ slice(0,strict_inside_connectivity.size(),1) ] =
			strict_inside_connectivity;
			vtk_out_.connectivity[ slice(strict_inside_connectivity.size(), boundary_connectivity.size(),1) ] =
			boundary_connectivity;
	
			// An adjustment will need to be made to offsets.
			size_t    strict_inside_connect_new_offset = strict_inside_offsets[strict_inside_offsets.size()-1];
		
			vtk_out_.offsets.resize (
			strict_inside_offsets.size() + boundary_offsets.size() );
			vtk_out_.offsets[ slice(0,strict_inside_offsets.size(),1) ] =
			strict_inside_offsets;
			vtk_out_.offsets[ slice(strict_inside_offsets.size(),boundary_offsets.size(),1) ] =
			boundary_offsets + valarray<int> ( strict_inside_connect_new_offset, boundary_offsets.size() );
		}
		
		
		
		PetscPrintf ( PETSC_COMM_WORLD,
		              "There are %i boundary patches, %i patches with their centre inside the domain and  %i patches lying strictly inside the domain.\n",
		              bdry_keys.size(),
		              inside_poss_bdry_keys.size(),
		              strict_inside_keys.size() );
	}
	
	set_union ( bdry_keys.begin(), bdry_keys.end(),
			            inside_poss_bdry_keys.begin(), inside_poss_bdry_keys.end(),
			            inserter( cov_struct_->patch_keys, cov_struct_->patch_keys.begin() ) );

	
//	CALLGRIND_START_INSTRUMENTATION;
	
	PetscPrintf ( PETSC_COMM_WORLD, "generate_index_information" ); PetscGetTime ( &log_t0 );
	generate_index_information ( *cov_struct_ );
	PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, " took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscPrintf ( PETSC_COMM_WORLD, "generate_patch_information" ); PetscGetTime ( &log_t0 );
	generate_patch_information ( *cov_struct_ );
	PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, "took %f minutes.\n", (log_t1-log_t0)/60.0 );

	// These boundary keys are insufficient as strictly interior boxes
	// may expand to patches that intersect the boundary. This must be
	// accounted for in order to guarantee integration within the domain.
	// 2009-07-06 ML.
	PetscPrintf ( PETSC_COMM_WORLD, "generate_patch_boundary_information" ); PetscGetTime ( &log_t0 );
	generate_patch_boundary_information ( *cov_struct_, bdry_keys, pum_d_->domain );
	PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, " took %f minutes.\n", (log_t1-log_t0)/60.0 );

#if 0 // Making sure we got all the boxes that expand to be boundary patches. 2009-07-06 ML.
	{
		clog << "MMM Debug bdry_keys are :\n";
		set<string>::iterator iit(bdry_keys.begin()), eend(bdry_keys.end());
		for ( ; iit!=eend; ++iit ) {
			clog << *iit << ", ";
		}
	}
#endif
	
	bdry_keys.clear();
	inside_poss_bdry_keys.clear();
	
	
	// PetscPrintf ( PETSC_COMM_WORLD, "Clearing cover_structure->patch_neighbour_keys." );
	// CANNOT clear either of these. 2009-06-13 ML.
	// cov_struct_->patch_keys.clear();
	// cov_struct_->patch_neighbour_keys.clear();
	// In particular, these are used in get_solution_evaluation_structure. 2009-06-24 ML.
	
	global_approx_space_.reset ( new global_approximation_space<2> ( details->local_approx_space ) );
	global_approx_space_->give_cover_structure( cov_struct_ );
	
	// Debugging found that a hint was being provided which was incorrect.
	// The global approximation space should be used to determine the
	// space required to store the coefficients. patch_to_num_dof_.sum()
	// was being provided which was probably zero on first use and then
	// nonzero when pum_convergence was being reused. valgrind/memcheck
	// found this one.
	//
	// Nasty little bug that threw up double deletion errors claiming to
	// be VecDestroy problems.
	// 2009-02-20 ML.
	sol_.reset ( new solution<2> );
	int num_dof = sol_->set_global_approximation_space ( global_approx_space_ );
	
	if ( pum_return_info ) {
		pum_return_info->num_dof = num_dof;
	}

	
	PetscPrintf ( PETSC_COMM_WORLD, "create_petsc_objects" ); PetscGetTime ( &log_t0 );
	create_petsc_objects( assemble_bleach_matrix, assemble_siggia_matrix );
	PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, " took %f minutes.\n", (log_t1-log_t0)/60.0 );
	assert ( petsc_objects_created_ == true );
	
	// Valgrind/massif shows that around 60% of 1.284 GB is taken up by
	// the three matrices (mass, stiffness, ksp matrix). This doesn't leave
	// sufficient memory to attempt a solve on night which only has 2 GB.
	// 2009-05-28 ML.
	//PetscPrintf ( PETSC_COMM_WORLD, "Early exit for massif." );
	//exit(0);
	
	PetscPrintf ( PETSC_COMM_WORLD, "assemble_matrix_vector" ); PetscGetTime ( &log_t0 );
	assemble_matrix_vector( diag_quad_d );
	PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, " took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscPrintf ( PETSC_COMM_WORLD, "petsc_objects_all_assemble" ); PetscGetTime ( &log_t0 );
	petsc_objects_all_assemble ();
	PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, " took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	initialized_ = true;
//	CALLGRIND_STOP_INSTRUMENTATION;
}

template <>
void pum_convergence<2>::draw_cover ( const string & filename, bool labels )
{
	ofstream out ( filename.c_str() );
	
	out << "";
	
	vector<string> &             keys = cov_struct_->index_to_patch_key;
	basic_dbinary_tree<2> &    dtree = *(pum_d_->dtree);
	
	{
		int num_keys = keys.size();
		
		vector<box<2> >    boxes;
		dtree.get_box ( keys, boxes );
		
		out << "#!/usr/bin/gnuplot -persist\n";
		out << "set size square\n";
		
		double centre[2];
		if ( labels ) {
			for ( int i=0; i<num_keys; ++i ) {
				out << "set label \""
					<< keys[i] << "\" at ";
				boxes[i].get_centre_point ( centre );
				for ( int d=0; d<2; ++d ) {
					if ( d!=0 ) {
						out << ",";
					}
					out << centre[d];
				}
				out << " font \"Times,5\"\n";
			}
		}
		
		out << "plot '-' w l, '-' w l\n";
		
		gp_draw<2> ( pum_d_->domain, out );
		out << "e\n";
		gp_draw ( boxes, out );
	}
}


template <int dim>
void
pum_convergence<dim>::set_problem_L2_projection ()
{
	assert ( initialized_ );
	
	if ( comm_rank_ == comm_size_-1 ) {
		int num_pts = (vtk_out_.points_3D).size() / 3;
		
		valarray<double> points_two_D ( 2*num_pts );
#if 0
	const valarray<double>& points_three_D = vtk_out_.points_3D;
	
	points_two_D[ slice(0, num_pts, 2) ] = points_three_D[ slice(0, num_pts, 3) ];
	points_two_D[ slice(1, num_pts, 2) ] = points_three_D[ slice(1, num_pts, 3) ];
#else
		points_two_D[ slice(0, num_pts, 2) ] = (vtk_out_.points_3D)[ slice(0, num_pts, 3) ];
		points_two_D[ slice(1, num_pts, 2) ] = (vtk_out_.points_3D)[ slice(1, num_pts, 3) ];
#endif
		vtk_out_.scalars.resize(points_two_D.size(), -1.0);
		assert ( pum_d_->rhs_ptr);
		(pum_d_->rhs_ptr)->global_evaluate ( box<2>(), points_two_D, vtk_out_.scalars );
		
		PetscPrintf( PETSC_COMM_SELF, "There are %i points, %i scalars and max scalar is %f.\n",
		             vtk_out_.points_3D.size()/3, vtk_out_.scalars.size(), vtk_out_.scalars.max() );
		
		vtk_out_.points_3D[ slice(2, vtk_out_.scalars.size(), 3) ] = vtk_out_.scalars;

		
		vtk_append_raw ( vtk_out_, "pum_convergence-1_function_to_project.vtp" );
	}

	load_mass_matrix ();
}

template <int dim>
void
pum_convergence<dim>::set_problem_Helmholtz ( double kappa, double a )
{
	assert ( initialized_ );
	
	load_mass_matrix ( a );
	add_scaled_stiffness_matrix ( -kappa );
}

template <int dim>
void
pum_convergence<dim>::set_problem_Poisson ( double kappa )
{
	assert ( initialized_ );
	
	load_stiffness_matrix ( kappa );
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::solve ( bool output_iter )
{
	solve ( *sol_, output_iter );
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::solve_and_get_results ( convergence_results &    cr )
{
	assert ( initialized_ );
	
	solve ( *sol_ , true );
	
	cr.patch_to_num_dof.resize ( patch_to_num_dof_.size() );
	cr.patch_to_num_dof = patch_to_num_dof_;
	
	PetscErrorCode ierr;
	
// 	PetscPrintf ( PETSC_COMM_WORLD, "Suspect: pu_val_function_val_vec_ last entry.\n" );
// 	VecView( pu_val_function_val_vec_, PETSC_VIEWER_STDOUT_WORLD );
	
	ierr = VecSum ( pu_val_function_val_vec_,          &cr.function_val_integral ); CHKERRQ(ierr);
	ierr = VecSum ( pu_val_function_val_squared_vec_,  &cr.function_val_squared_integral ); CHKERRQ(ierr);
	ierr = VecSum ( pu_val_function_grad_squared_vec_, &cr.function_grad_squared_integral ); CHKERRQ(ierr);
	
	ierr = VecDot ( ksp_solution_, approx_val_vec_                , &cr.approx_val_integral ); CHKERRQ(ierr);
	ierr = VecDot ( ksp_solution_, approx_val_function_val_vec_   , &cr.approx_val_function_val_integral ); CHKERRQ(ierr);
	ierr = VecDot ( ksp_solution_, approx_grad_function_grad_vec_ , &cr.approx_grad_function_grad_integral ); CHKERRQ(ierr);
	
	ierr = MatMult ( mass_matrix_,     ksp_solution_, calculation_vec_ ); CHKERRQ(ierr);
	ierr = VecDot  ( calculation_vec_, ksp_solution_, &cr.approx_val_squared_integral ); CHKERRQ(ierr);
	
	ierr = MatMult ( stiffness_matrix_, ksp_solution_, calculation_vec_ ); CHKERRQ(ierr);
	ierr = VecDot  ( calculation_vec_,  ksp_solution_, &cr.approx_grad_squared_integral ); CHKERRQ(ierr);
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::rhs_vec_mass_mult_previous_solution ()
{	
	PetscErrorCode ierr;
	ierr = MatMult ( mass_matrix_, ksp_solution_, ksp_rhs_vec_ ); CHKERRQ(ierr);
	return ierr;
}

template <int dim>
double
pum_convergence<dim>::integrate_solution_over_domain ()
{
	double tmp;
	PetscErrorCode ierr = VecDot ( ksp_solution_, approx_val_vec_, &tmp ); CHKERRQ(ierr);
	return tmp;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::species_integrate_over_domain ( valarray<double> & out )
{
	assert ( num_species_ > 0 );
	PetscErrorCode ierr = 0;
	for ( int i=0; i<num_species_; ++i ) {
		ierr = VecDot ( species_previous_vec_[i], approx_val_vec_, &out[i] );CHKERRQ(ierr);
	}
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::species_sum_output_vtk ( const string& solution, const string& error, int num )
{
	VecZeroEntries( ksp_solution_ );
	for ( int i=0; i<num; ++i ) {
		VecAXPY(ksp_solution_,1.0,species_previous_vec_[i]);
	}
	
	int sendcount = computenode_to_num_dof_[comm_rank_];
	
	// PetscScalar must be double.
	double* sendbuf = 0;
	
	PetscErrorCode ierr = VecGetArray ( ksp_solution_, &sendbuf ); CHKERRQ(ierr);
	
	if ( comm_size_ == 1 ) {
		assert ( comm_rank_ == 0 );
		assert ( sen_dbuf );
		assert ( sol_->global_coefficients_.size() == computenode_to_num_dof_[0] );
		
		// solve would otherwise pass sol in as *sol_.
		copy ( sendbuf, sendbuf + computenode_to_num_dof_[0], &sol_->global_coefficients_[0] );
	} else {
		// void* sendbuf, int sendcount, MPI_Datatype sendtype,
		// void* recvbuf
		// int* recvcounts
		// int* displs
		
		MPI_Allgatherv ( sendbuf, sendcount, MPI_DOUBLE,
		                 &sol_->global_coefficients_[0],
		                 &computenode_to_num_dof_[0],
		                 &computenode_to_start_dof_[0],
		                 MPI_DOUBLE,
		                 PETSC_COMM_WORLD );
	}
	
	ierr = VecRestoreArray ( ksp_solution_, &sendbuf ); CHKERRQ(ierr);
	
	{
		dof_structure::ptr tmp_dof_st;
		get_dof_structure( tmp_dof_st );
		assert ( tmp_dof_st );
		sol_->set_dof_structure( tmp_dof_st );
	}
	
	this->output_vtk(solution,error);
	
	return ierr;
}

template <int dim>
void
pum_convergence<dim>::output_vtk ( const string& name, const string& error_name, bool do_cell_scalars )
{
	const valarray<double>& points_three_D = vtk_out_.points_3D;
	
	int num_pts = points_three_D.size() / 3;
	
	valarray<double> points_two_D ( 2*num_pts );
	
	points_two_D[ slice(0, num_pts, 2) ] = points_three_D[ slice(0, num_pts, 3) ];
	points_two_D[ slice(1, num_pts, 2) ] = points_three_D[ slice(1, num_pts, 3) ];

	solution_evaluation_structure    sol_eval_struct;
	global_approx_space_->get_solution_evaluation_structure ( points_two_D, sol_eval_struct );
	
	sol_->global_evaluate ( points_two_D, sol_eval_struct, vtk_out_.scalars );
	vtk_out_.points_3D[ slice(2, vtk_out_.scalars.size(), 3) ] = vtk_out_.scalars;
	
	// Currently cell data is used to show which MPI processes own which patches.
	// Only to be used for small problem sizes. 2009-09-05 ML.
	if ( do_cell_scalars ) {
		int num_patches = computenode_to_num_patches_.sum();
		
		// The polygons are both the cell interiors and the boundary cells.
		vtk_out_.cell_scalars.resize( vtk_out_.offsets.size(), -10.0 );
		
		// Loop-carried dependency that can be eliminated
		// for large problems but it makes little sense to
		// optimize at this time. 2009-09-05 ML.
		int num_done = 0;
		
		for ( int i=0; i<comm_size_; ++i ) {
			int upper = computenode_to_num_patches_[i];
			
			for ( int pat=0; pat<upper; ++pat ) {
				vtk_out_.cell_scalars[num_done + pat] = i;
			}
			num_done += upper;
		}
	}
	
	vtk_append_raw ( vtk_out_, name, do_cell_scalars );
	
	if ( error_name != "" ) {
		vtk_error_out_ = vtk_out_;
		
		assert ( pum_d_->fun_ptr);
		// This is likely to cause evaluation outside the domain which makes the autoscaling
		// option of ParaView less useful. On option is to avoid global evaluations for the
		// known-solution global function but this may be brittle with regard to future changes.
		// The second option is to zero the points outside the domain at the last step.
		// 2009-06-16 ML.
		(pum_d_->fun_ptr)->global_evaluate ( box<2>(), points_two_D, vtk_error_out_.scalars );
		
		// L2 error integrand: take difference and square.
		vtk_error_out_.scalars -= vtk_out_.scalars;
		vtk_error_out_.scalars *= vtk_error_out_.scalars;
		
		set_outside ( 0.0, sol_eval_struct, vtk_error_out_.scalars );
		
		vtk_error_out_.points_3D[ slice(2, vtk_error_out_.scalars.size(), 3) ] = vtk_error_out_.scalars;
		
		vtk_append_raw ( vtk_error_out_, error_name, do_cell_scalars );
	}
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::create_petsc_objects ( bool create_bleach_matrix,
                                             bool create_siggia_matrix )
{
	if ( petsc_objects_created_ ) {
		destroy_petsc_objects();
	}
	
	/**
	 * Keep a copy of patch_to_num_dof. (1)
	 */
	assert (global_approx_space_ );
	global_approx_space_->get_patch_to_num_dof ( patch_to_num_dof_ );
	
	int row_nonzero = global_approx_space_->est_patch_num_dof() * static_cast<int>( pow(3,dim) );
	PetscPrintf ( PETSC_COMM_WORLD, "Estimated row nonzero is %i.\n", row_nonzero );
	
	/**
	 * Determine patch start indices (global dof numbering). (2)
	 */
	{
		int num_patches = patch_to_num_dof_.size();
		patch_to_start_index_.resize( num_patches );
		patch_to_start_index_[0] = 0;
		for ( int i=1; i<num_patches; ++i ) {
			patch_to_start_index_[i] = patch_to_start_index_[i-1] + patch_to_num_dof_[i-1];
		}
	}


	{
	/**
	 * Distribute the patches between the computenodes. (3)
	 */
		computenode_to_num_patches_.resize( comm_size_ );

		div_t res = div( patch_to_num_dof_.size(), comm_size_ );
		int guaranteed = res.quot;
		int cutoff = res.rem;

		for ( int i=0; i<comm_size_; ++i ) {
			if ( i<cutoff ) {
				computenode_to_num_patches_[i] = guaranteed + 1;
			} else {
				computenode_to_num_patches_[i] = guaranteed;
			}
		}

	/**
	 * Based on the distribution of patches, distribute the dofs. (4)
	 */
		computenode_to_start_patch_.resize( comm_size_ );

		for ( int i=0; i<comm_size_; ++i ) {
			if ( i<=cutoff ) {
				computenode_to_start_patch_[i] = i*(guaranteed+1);
			} else {
				computenode_to_start_patch_[i] = cutoff*(guaranteed+1)+(i-cutoff)*guaranteed;
			}
		}
	}
	
	/**
	 * Distribute the dof between the computenodes. (5)
	 */
	{
		computenode_to_num_dof_.resize( comm_size_ );
		for ( int i=0; i<comm_size_; ++i ) {
			slice s ( computenode_to_start_patch_[i], computenode_to_num_patches_[i], 1 );
			const valarray<int>& tmp = patch_to_num_dof_[s];

			computenode_to_num_dof_[i] = tmp.sum();
		}
	}
	
	/**
	 * Work out patches this computenode is responsible for. (6)
	 */
	{
		int num   = computenode_to_num_patches_[comm_rank_];
		int start = computenode_to_start_patch_[comm_rank_];
		patch_indices_on_computenode_.resize( num );
		for ( int i=0; i<num; ++i ) {
			patch_indices_on_computenode_[i] = start + i;
		}
	}
	
	/**
	 * Work out how the dof are split between the computenodes.
	 */
	{
		computenode_to_start_dof_.resize ( comm_size_ );
		computenode_to_start_dof_[0] = 0;
		for ( int i=1; i<comm_size_; ++i ) {
			computenode_to_start_dof_[i] = computenode_to_start_dof_[i-1] + computenode_to_num_dof_[i-1];
		}
	}
	
	/**
	 * Created petsc objects.
	 */
	{
		int num_patches = computenode_to_num_patches_.sum();
		int global_size = computenode_to_num_dof_.sum();
		
		PetscErrorCode ierr;
		
		// Length is number of patches.
		ierr = VecCreate(PETSC_COMM_WORLD, &pu_val_function_val_vec_); CHKERRQ(ierr);
		ierr = VecSetSizes(pu_val_function_val_vec_, computenode_to_num_patches_[comm_rank_], num_patches); CHKERRQ(ierr);
		ierr = VecSetFromOptions(pu_val_function_val_vec_); CHKERRQ(ierr);
		ierr = VecZeroEntries(pu_val_function_val_vec_); CHKERRQ(ierr);
		
		ierr = VecCreate(PETSC_COMM_WORLD, &pu_val_function_val_squared_vec_); CHKERRQ(ierr);
		ierr = VecSetSizes(pu_val_function_val_squared_vec_, computenode_to_num_patches_[comm_rank_], num_patches); CHKERRQ(ierr);
		ierr = VecSetFromOptions(pu_val_function_val_squared_vec_); CHKERRQ(ierr);
		ierr = VecZeroEntries(pu_val_function_val_squared_vec_); CHKERRQ(ierr);
		
		ierr = VecCreate(PETSC_COMM_WORLD, &pu_val_function_grad_squared_vec_); CHKERRQ(ierr);
		ierr = VecSetSizes(pu_val_function_grad_squared_vec_, computenode_to_num_patches_[comm_rank_], num_patches); CHKERRQ(ierr);
		ierr = VecSetFromOptions(pu_val_function_grad_squared_vec_); CHKERRQ(ierr);
		ierr = VecZeroEntries(pu_val_function_grad_squared_vec_); CHKERRQ(ierr);
		
		// The next three are size of PUM basis.
		ierr = VecCreate(PETSC_COMM_WORLD, &approx_val_vec_); CHKERRQ(ierr);
		ierr = VecSetSizes(approx_val_vec_, computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
		ierr = VecSetFromOptions(approx_val_vec_); CHKERRQ(ierr);
		ierr = VecZeroEntries(approx_val_vec_); CHKERRQ(ierr);
		
		ierr = VecCreate(PETSC_COMM_WORLD, &approx_val_function_val_vec_); CHKERRQ(ierr);
		ierr = VecSetSizes(approx_val_function_val_vec_, computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
		ierr = VecSetFromOptions(approx_val_function_val_vec_); CHKERRQ(ierr);
		ierr = VecZeroEntries(approx_val_function_val_vec_); CHKERRQ(ierr);
		
		ierr = VecCreate(PETSC_COMM_WORLD, &approx_grad_function_grad_vec_); CHKERRQ(ierr);
		ierr = VecSetSizes(approx_grad_function_grad_vec_, computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
		ierr = VecSetFromOptions(approx_grad_function_grad_vec_); CHKERRQ(ierr);
		ierr = VecZeroEntries(approx_grad_function_grad_vec_); CHKERRQ(ierr);
		
		ierr = MatCreateMPIAIJ ( PETSC_COMM_WORLD,
		                         computenode_to_num_dof_[comm_rank_],
		                         computenode_to_num_dof_[comm_rank_],
		                         global_size,
		                         global_size,
		                         row_nonzero,
		                         PETSC_NULL,
		                         0,
		                         PETSC_NULL,
		                         &mass_matrix_ );
		ierr = MatSetFromOptions(mass_matrix_); CHKERRQ(ierr);
		ierr = MatSeqSBAIJSetPreallocation ( mass_matrix_,
		                                     1,
		                                     row_nonzero,
		                                     PETSC_NULL ); CHKERRQ(ierr);
		
		ierr = MatCreateMPIAIJ ( PETSC_COMM_WORLD,
		                         computenode_to_num_dof_[comm_rank_],
		                         computenode_to_num_dof_[comm_rank_],
		                         global_size,
		                         global_size,
		                         row_nonzero,
		                         PETSC_NULL,
		                         0,
		                         PETSC_NULL,
		                         &stiffness_matrix_ );
		ierr = MatSetFromOptions(stiffness_matrix_); CHKERRQ(ierr);
		ierr = MatSeqSBAIJSetPreallocation ( stiffness_matrix_,
		                                     1,
		                                     row_nonzero,
		                                     PETSC_NULL ); CHKERRQ(ierr);
		
		
		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp_); CHKERRQ(ierr);
		
		// PetscPrintf ( PETSC_COMM_WORLD, "Guess nonzero in attempt to reduce its.\n" );
		// Nonzero guess doesn't help with Helmholtz Eqn. Perhaps the conditioning of the
		// matrix needs more consideration. Allowing initial nonzero guess definitely does
		// help time dependent problems when the new solution is expected to be close to
		// the previous solution.
		// 2009-02-25 ML.
		// Commenting out until use of instationary. 2009-03-13 ML.
		//ierr = KSPSetInitialGuessNonzero ( ksp_, PETSC_TRUE ); CHKERRQ(ierr);
		
		ierr = KSPSetFromOptions(ksp_); CHKERRQ(ierr);

		ierr = MatCreateMPIAIJ ( PETSC_COMM_WORLD,
		                         computenode_to_num_dof_[comm_rank_],
		                         computenode_to_num_dof_[comm_rank_],
		                         global_size,
		                         global_size,
		                         row_nonzero,
		                         PETSC_NULL,
		                         0,
		                         PETSC_NULL,
		                         &ksp_matrix_ );
		ierr = MatSetFromOptions(ksp_matrix_); CHKERRQ(ierr);
		ierr = MatSeqSBAIJSetPreallocation ( ksp_matrix_,
		                                     1,
		                                     row_nonzero,
		                                     PETSC_NULL ); CHKERRQ(ierr);
		
		if ( create_bleach_matrix ) {
			ierr = MatCreateMPIAIJ ( PETSC_COMM_WORLD,
			                         computenode_to_num_dof_[comm_rank_],
			                         computenode_to_num_dof_[comm_rank_],
			                         global_size,
			                         global_size,
			                         row_nonzero,
			                         PETSC_NULL,
			                         0,
			                         PETSC_NULL,
			                         &bleach_matrix_ );
			ierr = MatSetFromOptions(bleach_matrix_); CHKERRQ(ierr);
			bleach_matrix_created_ = true;
			// Notice no preallocation as there are expected to be very few entries.
		}
		
		if ( create_siggia_matrix ) {
			ierr = MatCreateMPIAIJ ( PETSC_COMM_WORLD,
			                         computenode_to_num_dof_[comm_rank_],
			                         computenode_to_num_dof_[comm_rank_],
			                         global_size,
			                         global_size,
			                         row_nonzero,
			                         PETSC_NULL,
			                         0,
			                         PETSC_NULL,
			                         &siggia_matrix_ );
			ierr = MatSetFromOptions(siggia_matrix_); CHKERRQ(ierr);
			siggia_matrix_created_ = true;
		}
		
		ierr = VecCreate(PETSC_COMM_WORLD, &ksp_rhs_vec_); CHKERRQ(ierr);
		ierr = VecSetSizes(ksp_rhs_vec_, computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
		ierr = VecSetFromOptions(ksp_rhs_vec_); CHKERRQ(ierr);
		ierr = VecZeroEntries(ksp_rhs_vec_); CHKERRQ(ierr);
		
		ierr = VecCreate(PETSC_COMM_WORLD, &ksp_solution_); CHKERRQ(ierr);
		ierr = VecSetSizes(ksp_solution_, computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
		ierr = VecSetFromOptions(ksp_solution_); CHKERRQ(ierr);
		
		ierr = VecCreate(PETSC_COMM_WORLD, &calculation_vec_); CHKERRQ(ierr);
		ierr = VecSetSizes(calculation_vec_, computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
		ierr = VecSetFromOptions(calculation_vec_); CHKERRQ(ierr);
		
		if ( num_species_ > 0 ) {
			species_matrices_.resize ( num_species_, PETSC_NULL );
			species_matrices_populated_.resize ( num_species_, false );
			species_previous_vec_.resize ( num_species_, PETSC_NULL );
			species_current_vec_.resize ( num_species_, PETSC_NULL );
			
			for ( int i=0; i<num_species_; ++i ) {
				ierr = MatCreateMPIAIJ ( PETSC_COMM_WORLD,
				                         computenode_to_num_dof_[comm_rank_],
				                         computenode_to_num_dof_[comm_rank_],
				                         global_size,
				                         global_size,
				                         row_nonzero,
				                         PETSC_NULL,
				                         0,
				                         PETSC_NULL,
				                         &species_matrices_[i] );
				ierr = MatSetFromOptions(species_matrices_[i]); CHKERRQ(ierr);
				ierr = MatSeqSBAIJSetPreallocation ( species_matrices_[i],
				                                     1,
				                                     row_nonzero,
				                                     PETSC_NULL ); CHKERRQ(ierr);
				
				ierr = VecCreate(PETSC_COMM_WORLD, &species_previous_vec_[i]); CHKERRQ(ierr);
				ierr = VecSetSizes(species_previous_vec_[i], computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
				ierr = VecSetFromOptions(species_previous_vec_[i]); CHKERRQ(ierr);
				
				ierr = VecCreate(PETSC_COMM_WORLD, &species_current_vec_[i]); CHKERRQ(ierr);
				ierr = VecSetSizes(species_current_vec_[i], computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
				ierr = VecSetFromOptions(species_current_vec_[i]); CHKERRQ(ierr);
			}
			species_matrices_vectors_created_ = true;
		}
#if 0
		{
			int first_row, last_row;
			MatGetOwnershipRange ( ksp_matrix_, &first_row, &last_row );
			PetscPrintf( PETSC_COMM_SELF,
			             "Computenode %i owns %i of a total of %i degrees of freedom : [%i, %i).\n",
			             comm_rank_,
			             computenode_to_num_dof_[comm_rank_],
			             global_size,
			             first_row,
			             last_row
			              );
		}
#endif

		petsc_objects_created_ = true;
		return ierr;
	}
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::destroy_petsc_objects ()
{
	PetscErrorCode ierr;
	if ( petsc_objects_created_ ) {
		ierr = VecDestroy ( pu_val_function_val_vec_ ); CHKERRQ(ierr);
		ierr = VecDestroy ( pu_val_function_val_squared_vec_ ); CHKERRQ(ierr);
		ierr = VecDestroy ( pu_val_function_grad_squared_vec_ ); CHKERRQ(ierr);
		
		ierr = VecDestroy ( approx_val_vec_ ); CHKERRQ(ierr);
		ierr = VecDestroy ( approx_val_function_val_vec_ ); CHKERRQ(ierr);
		ierr = VecDestroy ( approx_grad_function_grad_vec_ ); CHKERRQ(ierr);
		
		ierr = MatDestroy ( mass_matrix_ ); CHKERRQ(ierr);
		ierr = MatDestroy ( stiffness_matrix_ ); CHKERRQ(ierr);

		ierr = KSPDestroy ( ksp_ ); CHKERRQ(ierr);

		ierr = MatDestroy ( ksp_matrix_ ); CHKERRQ(ierr);
		ierr = VecDestroy ( ksp_rhs_vec_ ); CHKERRQ(ierr);
		ierr = VecDestroy ( ksp_solution_ ); CHKERRQ(ierr);
		
		ierr = VecDestroy ( calculation_vec_ ); CHKERRQ(ierr);
		
		ksp_operators_set_     = false;
		ksp_matrix_populated_  = false;
		petsc_objects_created_ = false;
		initialized_           = false;
	}
	if ( bleach_matrix_created_ ) {
		ierr = MatDestroy ( bleach_matrix_ ); CHKERRQ(ierr);
	}
	if ( siggia_matrix_created_ ) {
		ierr = MatDestroy ( siggia_matrix_ ); CHKERRQ(ierr);
	}
	if ( species_matrices_vectors_created_ ) {
		assert ( num_species_ > 0 );
		for ( int i=0; i<num_species_; ++i ) {
			ierr = MatDestroy ( species_matrices_[i] ); CHKERRQ(ierr);
			ierr = VecDestroy ( species_previous_vec_[i] ); CHKERRQ(ierr);
			ierr = VecDestroy ( species_current_vec_[i] ); CHKERRQ(ierr);
		}
		species_matrices_vectors_created_ = false;
	}
	return ierr;
}

template <int dim>
void
pum_convergence<dim>::assemble_matrix_vector ( diagnostic_quad_details * diag_quad_d )
{
	assert ( cov_struct_ );
	const vector<vector<int> >& index_to_neighbour_indices = cov_struct_->index_to_neighbour_indices;
	
	int first_patch_ind;
	int num_neigh;
	int second_patch_ind;
	
	int num = patch_indices_on_computenode_.size();
	
	integration_scratch          is;
	matrix_vector_particulars    mvs;
	
#pragma omp parallel for private( first_patch_ind,num_neigh,second_patch_ind, is, mvs)
	for ( int i=0; i<num; ++i ) {
		first_patch_ind = patch_indices_on_computenode_[i];
		
		const vector<int>& neighbours = index_to_neighbour_indices[first_patch_ind];
		num_neigh = neighbours.size();
		
		for ( int j=0; j<num_neigh; ++j ) {
			second_patch_ind = neighbours[j];

			// Some values only need to be calculated patch-wise.
			assemble_matrix_vector_block ( first_patch_ind,
			                               second_patch_ind,
			                               is,
			                               mvs, diag_quad_d );
			
			set_matrix_entries( first_patch_ind, second_patch_ind, mvs );
		}
	}
}

/**
	This function gives a highly efficient implementation of PUM.
*/
template <int dim>
void
pum_convergence<dim>::assemble_matrix_vector_block ( int                                patch_index_row,
                                                     int                                patch_index_column,
                                                     integration_scratch&               is,
                                                     matrix_vector_particulars&         mvs,
                                                     diagnostic_quad_details *          diag_quad_d )
{
	// Test functions correspond to the rows whereas
	// trial functions correspond to the columns.
	vector<global_basis_function<dim> > test;
	vector<global_basis_function<dim> > trial;

	global_approx_space_->get_global_basis_by_patch( patch_index_row,    test );
	global_approx_space_->get_global_basis_by_patch( patch_index_column, trial );
	
	integration_domain_decomposition<dim> tmp;

	int decomposition = pum_d_->integration_decomposition;
	
	while ( !global_basis_function<dim>::get_integration_domain_decomposition( test[0], trial[0], tmp, decomposition ) ) {
		--decomposition;
		assert( decomposition >= 0 );
	}
	
	tmp.intersection = cov_struct_->patches[patch_index_row];
	
	if ( patch_index_column != patch_index_row ) {
		tmp.intersection.clip_against( cov_struct_->patches[patch_index_column] );
	}
	
	bool boundary_check_needed = false;
	int    intersection_first_entry  = -1;
	int    intersection_num_inside   =  0;
	
	if ( cov_struct_->is_boundary_patch[patch_index_row] &&
	     cov_struct_->is_boundary_patch[patch_index_column] )
	{
		get_first_entry_num_inside ( pum_d_->domain,
		                             tmp.intersection,
		                             0, -1, intersection_first_entry,
		                             intersection_num_inside );
		
		if ( intersection_first_entry != -1 ) {
			boundary_check_needed = true;
		}
	}
	
#if 0
	clog << "MMM Debug integration points above and below domain for 1.3 1.3 level 3.\n";
	{
		const double * bext = static_cast<const box<dim>&>(tmp.intersection).get();
		
		if ( (bext[2]<1.0) && (1.0<bext[3]) ) {
			
			if ( cov_struct_->is_boundary_patch[patch_index_row] ) {
				clog << patch_index_row << " is a boundary patch\n";
			} else {
				clog << patch_index_row << " is not a boundary patch\n";
			}
			
			if ( cov_struct_->is_boundary_patch[patch_index_column] ) {
				clog << patch_index_column << " is a boundary patch\n";
			} else {
				clog << patch_index_column << " is not a boundary patch\n";
			}
			
			assert ( boundary_check_needed );
		}
	}
#endif

	//clog << "Debug force interior integrals.\n";
	//boundary_check_needed =false;
	
	// Check which of the tmp boxes require restriction of integration and
	// keep in a vector<bool>. This will be much simpler than anything else.
	// 2009-06-24 ML.
	
	int num_box = tmp.box_vector.size();
	
	int test_size  = test.size();
	int trial_size = trial.size();
	
	assert( test_size  > 0 );
	assert( trial_size > 0 );
	
	int block_size = test_size * trial_size;
	
	// Used to prevent too much copying when retrieving box_interior quadrature points.
	bool    last_was_box_interior = false;
	bool    this_is_box_interior  = true; // Will only be false when using more complicated cubature.
	
	mvs.pu_val_function_val_vec_vals          = 0.0;
	mvs.pu_val_function_val_squared_vec_vals  = 0.0;
	mvs.pu_val_function_grad_squared_vec_vals = 0.0;
	
	mvs.approx_val_rhs_fun_vec_vals.resize(test_size,0.0);
	
	mvs.approx_val_vec_vals.resize(test_size, 0.0);
	mvs.approx_val_function_val_vec_vals.resize(test_size, 0.0);
	mvs.approx_grad_function_grad_vec_vals.resize(test_size, 0.0);
	
	mvs.mass_matrix_vals.resize(block_size, 0.0);
	mvs.stiffness_matrix_vals.resize(block_size, 0.0);
	
	if ( bleach_matrix_created_ ) {
		mvs.bleach_matrix_vals.resize(block_size, 0.0);
	}
	
	if ( siggia_matrix_created_ ) {
		mvs.siggia_matrix_vals.resize( block_size, 0.0 );
	}
	
	//PetscPrintf ( PETSC_COMM_WORLD, "Num box is %i.\n", num_box );
	
	double    centre[2];
	
	for ( int b=0; b<num_box; ++b ) {		
		if ( !boundary_check_needed ) {
			
			// Even if a boundary check is not needed, the decomposed box
			// may still be entirely outside the domain. (e.g. two boundary boxes with
			// intersection outside the domain.) Unfortunately this increases the cost
			// in the interior. 2009-07-25 ML.
			tmp.box_vector[b].get_centre_point ( centre );
			// point_in_polygon is insensitive to boundary orientation. 2009-07-25 ML.
			if ( !point_in_polygon ( pum_d_->domain, centre[0], centre[1] ) ) {
					// Handle the case where box is entirely outside the domain.
				continue;
			}
		} else {
			// Boundary check is needed.

			int    tmpbox_first_entry = -1;
			int    tmpbox_num_inside  =  0;
			
			get_first_entry_num_inside ( pum_d_->domain,
			                             tmp.box_vector[b],
			                             intersection_first_entry,
			                             intersection_num_inside,
			                             tmpbox_first_entry,
			                             tmpbox_num_inside );
			
			if ( tmpbox_first_entry == -1 ) {
				tmp.box_vector[b].get_centre_point ( centre );
				
				/*
					Point in polygon is insensitive to boundary orientation.
					2009-07-25 ML.
				
				Delete - >
					in	anti	keep	continue
					Y	Y	Y	N
					Y	N	N	Y
					N	Y	N	Y
					N	N	Y	N
				<-delete
				*/
				// point_in_polygon is insensitive to boundary orientation. 2009-07-25 ML.	
				if ( !point_in_polygon ( pum_d_->domain, centre[0], centre[1] ) ) {
					// Handle the case where box is entirely outside the domain.
					continue;
				}
				
				this_is_box_interior = true;
			} else {
				this_is_box_interior = false;
				
				generate_gauss_legendre_rule_local ( pum_d_->boundary_base_rule,
				                                     tmp.box_vector[b],
				                                     pum_d_->domain,
				                                     tmpbox_first_entry,
				                                     pum_d_->anticlockwise_dom,
				                                     is.boundary_quad_rule );
				
				is.restricted_local_quad_pts.resize ( is.boundary_quad_rule.points.size() );
				is.restricted_local_quad_pts = is.boundary_quad_rule.points;
				
				is.restricted_local_quad_weights.resize ( is.boundary_quad_rule.weights.size() );
				is.restricted_local_quad_weights = is.boundary_quad_rule.weights;
				
				last_was_box_interior = false;
			}
		}
		
		if ( this_is_box_interior && !last_was_box_interior ) {
			// This is an interior box so we require the quadrature points.
			// The last box was not an interior box so we need to copy the point
			// in order to use them.
#if 1
			is.restricted_local_quad_pts.resize ( pum_d_->interior_box_cubature.points.size() );
			is.restricted_local_quad_pts = pum_d_->interior_box_cubature.points;
			
			is.restricted_local_quad_weights.resize ( pum_d_->interior_box_cubature.weights.size() );
			is.restricted_local_quad_weights = pum_d_->interior_box_cubature.weights;
#else
			is.restricted_local_quad_pts.resize ( pum_d_->box_cubature_points.size() );
			is.restricted_local_quad_pts = pum_d_->box_cubature_points;
			
			is.restricted_local_quad_weights.resize ( pum_d_->box_cubature_weights.size() );
			is.restricted_local_quad_weights = pum_d_->box_cubature_weights;
#endif
			last_was_box_interior = true;
		}
		
		const box<dim>& restricted_local_box = tmp.box_vector[b];

		
		// Random segfaulting will be observed if more than one thread
		// attempts to write into diag_quad_d memory.
		// 2009-06-24 ML.
		// Debugging using the integration points exported showed that
		// decomposed regions lying completely outside the domain were
		// not being detected!
		// 2009-07-06 ML.
#pragma omp critical
		if ( diag_quad_d ) {
			// Copy out location of points for checking.
			restricted_local_box.map_local_to_global( is.restricted_local_quad_pts,
			                                          diag_quad_d->interior_points );
#if 0
			clog << "MMM Debug: check integration points inside (0,1)x(0,1.0)\n";
			{
				int npts = diag_quad_d->interior_points.size()/2;

				for ( int i=0; i<npts; ++i ) {
					if ( diag_quad_d->interior_points[2*i]<0.0 ) {
						// Found some points below the domain.
						
						int    tmpbox_first_entry = -1;
						int    tmpbox_num_inside  =  0;
						
						get_first_entry_num_inside ( pum_d_->domain,
						                             tmp.box_vector[b],
						                             intersection_first_entry,
						                             intersection_num_inside,
						                             tmpbox_first_entry,
						                             tmpbox_num_inside );
						
						assert ( tmpbox_first_entry != -1 );
					}
					
					assert ( !(diag_quad_d->interior_points[2*i]<0.0) );
					assert ( !(diag_quad_d->interior_points[2*i]>1.0) );
					assert ( !(diag_quad_d->interior_points[2*i+1]<0.0) );
					if ( diag_quad_d->interior_points[2*i+1]>1.0 ) {
						clog << "bad number : "<< diag_quad_d->interior_points[2*i+1] << "\n";
					}
					assert ( !(diag_quad_d->interior_points[2*i+1]>1.0) );
				}

			}
#endif
				
			int vala_size = diag_quad_d->interior_points.size();
			
			copy ( &(diag_quad_d->interior_points[0]),
			       &(diag_quad_d->interior_points[0])+vala_size,
			       back_inserter ( diag_quad_d->all_interior_points ) );
		}
		
		if ( this_is_box_interior ) {
			is.mapped_quadrature_final_factor = restricted_local_box.measure()/(1<<dim);
			assert ( is.mapped_quadrature_final_factor > 0.0 );
		} else {
			is.mapped_quadrature_final_factor = 1.0;
		}
		
		calculate_entries ( restricted_local_box,
		                    test,
		                    trial,
		                    is,
		                    mvs,
		                    patch_index_row == patch_index_column/*,
		                    this_is_box_interior*/ );
	}
}

template <int dim>
void
pum_convergence<dim>::calculate_entries ( const box<dim>&                      restricted_local_box,
                                          vector<global_basis_function<dim> >& test,
                                          vector<global_basis_function<dim> >& trial,
                                          integration_scratch&                 is,
                                          matrix_vector_particulars&           mvs,
                                          bool                                 patch_indices_match/*,
                                          bool                                 this_is_box_interior*/ )
{
	
	int test_size  = test.size();
	int trial_size = trial.size();
	
	int num_coord = is.restricted_local_quad_pts.size();
	int num_pts   = num_coord/dim;
	
	assert ( pum_d_->fun_ptr );
	differentiable_function<dim> &   global_fun = *(pum_d_->fun_ptr);
	function<dim> &                  rhs_fun    = *(pum_d_->rhs_ptr);
	
	(restricted_local_box).map_local_to_global ( is.restricted_local_quad_pts,
	                                             is.global_pts_first );
	
	global_fun.global_evaluate ( restricted_local_box,
	                             is.global_pts_first,
	                             is.global_fun_val );
	
	global_fun.global_evaluate_grad ( restricted_local_box,
	                                  is.global_pts_first,
	                                  is.global_fun_grad );
	
	rhs_fun.global_evaluate ( restricted_local_box,
	                          is.global_pts_first,
	                          is.rhs_fun_val );
	
	const box<dim>& first_local_box      = test[0].access_box();
	const box<dim>& second_local_box     = trial[0].access_box();
	
	for ( int i=0; i<test_size; ++i ) {
		global_basis_function<dim> &    first = test[i];
		
		
		if ( 0==i ) {
			(first_local_box).map_restricted_local_to_local( restricted_local_box,
				is.restricted_local_quad_pts,
				is.local_pts_first );
			
			first.pu_local_evaluate_and_grad ( is.local_pts_first,
			                                   is.global_basis_fun_scratch,
			                                   is.pu_vals_first,
			                                   is.pu_grad_first );
			
			if ( bleach_matrix_created_ ) {
				pum_d_->bleach_indicator_fun->global_evaluate ( restricted_local_box, is.global_pts_first, is.indicator_vals );
				assert ( is.indicator_vals.min() > -1e-10 );
				assert ( is.indicator_vals.max() < 1 + 1e-10 );
			}
			
			if ( patch_indices_match ) {
					// Global function times pu to integrate the global
					// function.
				is.global_fun_val_pu_val.resize ( num_pts );
				is.global_fun_val_pu_val = is.pu_vals_first;
				assert ( is.global_fun_val.size()==num_pts );
				is.global_fun_val_pu_val *= is.global_fun_val;
				
					// Global function squared times pu to integrate
					// the square of the global function.
				is.global_fun_val_squared_pu_val.resize( num_pts );
				is.global_fun_val_squared_pu_val = is.global_fun_val_pu_val;
				is.global_fun_val_squared_pu_val *= is.global_fun_val;
				
					// Global function grad dot grad in preparation.
				is.global_fun_grad_global_fun_grad.resize( num_coord );
				assert ( is.global_fun_grad.size() == num_coord );
				is.global_fun_grad_global_fun_grad = is.global_fun_grad;
				is.global_fun_grad_global_fun_grad *= is.global_fun_grad;
				
					// Global function grad squared times pu to allow integration
					// of global function grad squared.
				is.global_fun_grad_squared_pu_val.resize ( num_pts, 0.0 );
				for ( int d=0; d<dim; ++d ) {
					is.global_fun_grad_squared_pu_val += is.global_fun_grad_global_fun_grad[ slice(d,num_pts,dim) ];
				}
				is.global_fun_grad_squared_pu_val *= is.pu_vals_first;
				
				if ( 1 /*this_is_box_interior*/ ) {
						// *= will invalidate the values so the values should not be reused.
					
					is.global_fun_val_pu_val *= is.restricted_local_quad_weights;
					mvs.pu_val_function_val_vec_vals
						+= is.global_fun_val_pu_val.sum() * is.mapped_quadrature_final_factor;
					
					is.global_fun_val_squared_pu_val *= is.restricted_local_quad_weights;
					mvs.pu_val_function_val_squared_vec_vals
						+= is.global_fun_val_squared_pu_val.sum() * is.mapped_quadrature_final_factor;

					is.global_fun_grad_squared_pu_val *= is.restricted_local_quad_weights;
					mvs.pu_val_function_grad_squared_vec_vals
						+= is.global_fun_grad_squared_pu_val.sum() * is.mapped_quadrature_final_factor; 
				}
			}
		}
		
		first.local_approx_local_evaluate_and_grad ( is.local_pts_first,
		                                             is.global_basis_fun_scratch,
		                                             is.local_approx_vals_first,
		                                             is.local_approx_grad_first );
		
		first.local_evaluate_and_grad( is.local_pts_first,
		                               is.pu_vals_first,
		                               is.pu_grad_first,
		                               is.local_approx_vals_first,
		                               is.local_approx_grad_first,
		                               is.global_basis_fun_scratch,
		                               is.vals_first,
		                               is.grad_first );
		
		if ( patch_indices_match ) {
				// For RHS vector and projection.
			is.approx_val_rhs_fun_val.resize(num_pts);
			is.approx_val_rhs_fun_val = is.vals_first;
			is.approx_val_rhs_fun_val *= is.rhs_fun_val;
			
				// For middle term of L2 norm evaluation.
			is.approx_val_global_fun_val.resize ( num_pts );
			is.approx_val_global_fun_val = is.vals_first;
			is.approx_val_global_fun_val *= is.global_fun_val;
			
				// For middle term of H1 norm evaluation.
			is.approx_grad_global_fun_grad.resize ( num_coord );
			is.approx_grad_global_fun_grad = is.grad_first;
			is.approx_grad_global_fun_grad *= is.global_fun_grad;
			
				// Seemed implausible that the W12 error was staying fairly
				// constant and it was due to *= here instead of +=.
				// 2009-02-25 ML.
			is.approx_grad_dot_global_fun_grad.resize(num_pts,0.0);
			for ( int d=0; d<dim; ++d ) {
				is.approx_grad_dot_global_fun_grad
					+= is.approx_grad_global_fun_grad[ slice(d,num_pts,dim) ];
			}
			
			if ( 1 /*this_is_box_interior*/ ) {
				is.approx_val_rhs_fun_val *= is.restricted_local_quad_weights;
				mvs.approx_val_rhs_fun_vec_vals[i]
					+= is.approx_val_rhs_fun_val.sum() * is.mapped_quadrature_final_factor;
				
					// is.vals_first is not to be over-written as is reused.
				is.vals_third.resize ( num_pts );
				is.vals_third = is.vals_first;
				is.vals_third *= is.restricted_local_quad_weights;
				mvs.approx_val_vec_vals[i]
					+= is.vals_third.sum() * is.mapped_quadrature_final_factor;
				
				is.approx_val_global_fun_val *= is.restricted_local_quad_weights;
				mvs.approx_val_function_val_vec_vals[i]
					+= is.approx_val_global_fun_val.sum() * is.mapped_quadrature_final_factor;
				
				is.approx_grad_dot_global_fun_grad *= is.restricted_local_quad_weights;
				mvs.approx_grad_function_grad_vec_vals[i]
					+= is.approx_grad_dot_global_fun_grad.sum() * is.mapped_quadrature_final_factor;
			}
		}

		if ( siggia_matrix_created_ )
		{
			pum_d_->equilibrium_fun->global_evaluate ( restricted_local_box,
			                                           is.global_pts_first,
			                                           is.equilibrium_vals );
			
			pum_d_->equilibrium_fun->global_evaluate_grad ( restricted_local_box,
			                                                is.global_pts_first,
			                                                is.equilibrium_grad );
		}
		
		for ( int j=0; j<trial_size; ++j ) {
			global_basis_function<dim> &    second = trial[j];
			
			if ( j==0 ) {
				(second_local_box).map_restricted_local_to_local( restricted_local_box,
					is.restricted_local_quad_pts,
					is.local_pts_second );
				second.pu_local_evaluate_and_grad ( is.local_pts_second,
				                                    is.global_basis_fun_scratch,
				                                    is.pu_vals_second,
				                                    is.pu_grad_second );
			}
			
			second.local_approx_local_evaluate_and_grad ( is.local_pts_second,
			                                              is.global_basis_fun_scratch,
			                                              is.local_approx_vals_second,
			                                              is.local_approx_grad_second );
			
			second.local_evaluate_and_grad( is.local_pts_second,
			                                is.pu_vals_second,
			                                is.pu_grad_second,
			                                is.local_approx_vals_second,
			                                is.local_approx_grad_second,
			                                is.global_basis_fun_scratch,
			                                is.vals_second,
			                                is.grad_second );
			
				// For the mass matrix.
			is.first_approx_val_second_approx_val.resize( num_pts );
			is.first_approx_val_second_approx_val = is.vals_first;
			is.first_approx_val_second_approx_val *= is.vals_second;
			
				// For the stiffness matrix.
			is.first_approx_grad_second_approx_grad.resize(num_coord);
			is.first_approx_grad_second_approx_grad = is.grad_first;
			is.first_approx_grad_second_approx_grad *= is.grad_second;
			is.first_approx_grad_dot_second_approx_grad.resize(num_pts, 0.0);
			for ( int d=0; d<dim; ++d ) {
				is.first_approx_grad_dot_second_approx_grad
					+= is.first_approx_grad_second_approx_grad[ slice(d,num_pts,dim) ];
			}
			
			if ( 1 /*this_is_box_interior*/ ) {
				
				is.first_approx_val_second_approx_val *= is.restricted_local_quad_weights;
				mvs.mass_matrix_vals[ i*trial_size + j ]
					+= is.first_approx_val_second_approx_val.sum() * is.mapped_quadrature_final_factor;
				
				if ( bleach_matrix_created_ ) {
					assert ( is.indicator_vals.size() == is.first_approx_val_second_approx_val.size() );
					is.first_approx_val_second_approx_val *= is.indicator_vals;
					mvs.bleach_matrix_vals[ i*trial_size + j ] += is.first_approx_val_second_approx_val.sum() * is.mapped_quadrature_final_factor;
				}
				
				is.first_approx_grad_dot_second_approx_grad *= is.restricted_local_quad_weights;
				mvs.stiffness_matrix_vals[ i*trial_size + j ]
					+= is.first_approx_grad_dot_second_approx_grad.sum() * is.mapped_quadrature_final_factor;
			}
			
			if ( siggia_matrix_created_ )
			{
				valarray<double> &    grad_product = is.grad_second;
				
				assert ( is.grad_first.size() == grad_product.size() );
				grad_product = is.grad_first;
				
				assert ( grad_product.size() == is.equilibrium_grad.size() );
				grad_product *= is.equilibrium_grad;
				
				assert ( is.grad_first.size() == is.equilibrium_grad.size() );
				
				is.grad_dot_grad.resize ( num_pts );
				is.grad_dot_grad = grad_product[ slice(0, num_pts, dim) ]; // Initialization.
				for ( int d=1; d<dim; ++d ) {
					is.grad_dot_grad += grad_product[ slice( d, num_pts, dim ) ];
				}
				
				valarray<double> &    equilibrium_integrand = is.grad_dot_grad;
				
				assert ( equilibrium_integrand.size() == is.vals_second.size() );
				equilibrium_integrand *= is.vals_second;
				
				assert ( equilibrium_integrand.size() == is.equilibrium_vals.size() );
				assert ( is.equilibrium_vals.min() > 0.0 );
				equilibrium_integrand /= is.equilibrium_vals;
				
				assert ( is.equilibrium_vals.min() > 0.0 );
				
				equilibrium_integrand *= is.restricted_local_quad_weights;
				mvs.siggia_matrix_vals[ i*trial_size + j ]
					+= equilibrium_integrand.sum() * is.mapped_quadrature_final_factor;
			}
		}
	}
}


template <int dim>
PetscErrorCode
pum_convergence<dim>::set_matrix_entries ( int                           first,
                                           int                           second,
                                           matrix_vector_particulars&    mvs )
{
	assert ( petsc_objects_created_ );
	
	int num_dof_first  = patch_to_num_dof_[first];
	int start_first  = patch_to_start_index_[first];

	valarray<int> idxm( num_dof_first );
	idxm[0] = start_first;
	for ( int i=1; i<num_dof_first; ++i ) {
		idxm[i] = start_first + i;
	}
	
// 	PetscPrintf ( PETSC_COMM_SELF, "Suspect last entries: handling %i pu_val_function_val_vec_ insert %f.\n",
// 	              first,
// 	              mvs.pu_val_function_val_vec_vals );
// 	
	PetscErrorCode ierr;

	/*
		Debugging helped find that we should be adding values to the patch-based
		vectors rather than inserting them as we did previously.
		
		When vectors were assembled separately from matrices, the vectors values
		were inserted patch-wise. As we now assemble matrices and vectors together,
		we necessarily assemble patch intersection-wise and so must accumulated
		the values.
		
		The incorrectness showed itself as monotone decreasing values with higher
		levels of refinement and a large spike in one corner in projecting constant
		functions. Worked on this until 6am last night with no success.
		2009-02-21 14:18 ML.
		
		Behaviour of pum_convergence-1 now matches behaviour of petsc_solver-3.
		Duplicate assembly was carried out where a patch is integrated when first==second
		but also for overlaps with various neighbour patches.
		
		Care has to be taken in the calculation to avoid unneeded work without missing
		out a required calculation. Care has to be taken in the setting of values to
		ensure correct addition without duplication in overlap regions.
		2009-02-21 18:14 ML.
	*/
	
#pragma omp critical
	if ( first==second ) {
		// Debugging helped find that start_first was being used in
		// place of `first' in these three statements so insertion was attempted outside
		// the correct range. ML.
		ierr = VecSetValue( pu_val_function_val_vec_,
				first,
				mvs.pu_val_function_val_vec_vals,
				INSERT_VALUES);
	
		ierr = VecSetValue( pu_val_function_val_squared_vec_,
				first,
				mvs.pu_val_function_val_squared_vec_vals,
				INSERT_VALUES);
	
		ierr = VecSetValue( pu_val_function_grad_squared_vec_,
				first,
				mvs.pu_val_function_grad_squared_vec_vals,
				INSERT_VALUES);

		// 	PetscPrintf ( PETSC_COMM_SELF, "Suspect last entries: handling %i ksp_rhs_vec_ insert %f.\n",
		// 	              first,
		// 	              mvs.approx_val_rhs_fun_vec_vals[0] );
		
		ierr = VecSetValues( ksp_rhs_vec_,
				num_dof_first,
				&idxm[0],
				&mvs.approx_val_rhs_fun_vec_vals[0],
				INSERT_VALUES ); //CHKERRQ(ierr);
	
		ierr = VecSetValues( approx_val_vec_,
				num_dof_first,
				&idxm[0],
				&mvs.approx_val_vec_vals[0],
				INSERT_VALUES ); //CHKERRQ(ierr);
	
		ierr = VecSetValues( approx_val_function_val_vec_,
				num_dof_first,
				&idxm[0],
				&mvs.approx_val_function_val_vec_vals[0],
				INSERT_VALUES ); //CHKERRQ(ierr);
		
		ierr = VecSetValues( approx_grad_function_grad_vec_,
				num_dof_first,
				&idxm[0],
				&mvs.approx_grad_function_grad_vec_vals[0],
				INSERT_VALUES ); //CHKERRQ(ierr);
	}	
	
	int num_dof_second = patch_to_num_dof_[second];
	int start_second = patch_to_start_index_[second];
	
	valarray<int> idxn( num_dof_second );
	idxn[0] = start_second;
	for ( int i=1; i<num_dof_second; ++i ) {
		idxn[i] = start_second + i;
	}

#pragma omp critical
	ierr = MatSetValues( mass_matrix_,
	                     num_dof_first,    &idxm[0],
	                     num_dof_second,   &idxn[0],
	                     &mvs.mass_matrix_vals[0],
	                     INSERT_VALUES ); CHKERRQ(ierr);

#pragma omp critical
	ierr = MatSetValues( stiffness_matrix_,
	                     num_dof_first,    &idxm[0],
	                     num_dof_second,   &idxn[0],
	                     &mvs.stiffness_matrix_vals[0],
	                     INSERT_VALUES ); CHKERRQ(ierr);
	
	if ( bleach_matrix_created_ ) {
		//		if ( (mvs.bleach_matrix_vals.max()>0.0) || (mvs.bleach_matrix_vals.min()<0.0) ) {
#pragma omp critical
			ierr = MatSetValues ( bleach_matrix_,
			                      num_dof_first,  &idxm[0],
			                      num_dof_second, &idxn[0],
			                      &mvs.bleach_matrix_vals[0],
			                      INSERT_VALUES ); CHKERRQ(ierr);
		//}
	}
	
	if ( siggia_matrix_created_ ) {
#pragma omp critical
		ierr = MatSetValues ( siggia_matrix_,
		                      num_dof_first,  &idxm[0],
		                      num_dof_second, &idxn[0],
		                      &mvs.siggia_matrix_vals[0],
		                      INSERT_VALUES ); CHKERRQ(ierr);
	}
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::petsc_objects_all_assemble ()
{	
	PetscErrorCode ierr;
	
	ierr = VecAssemblyBegin ( pu_val_function_val_vec_ ); CHKERRQ(ierr);
	ierr = VecAssemblyEnd (   pu_val_function_val_vec_ ); CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin ( pu_val_function_val_squared_vec_ ); CHKERRQ(ierr);
	ierr = VecAssemblyEnd (   pu_val_function_val_squared_vec_ ); CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin ( pu_val_function_grad_squared_vec_ ); CHKERRQ(ierr);
	ierr = VecAssemblyEnd (   pu_val_function_grad_squared_vec_ ); CHKERRQ(ierr);
	
	
	ierr = VecAssemblyBegin ( ksp_rhs_vec_ ); CHKERRQ(ierr);
	ierr = VecAssemblyEnd (   ksp_rhs_vec_ ); CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin ( approx_val_vec_ ); CHKERRQ(ierr);
	ierr = VecAssemblyEnd (   approx_val_vec_ ); CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin ( approx_val_function_val_vec_ ); CHKERRQ(ierr);
	ierr = VecAssemblyEnd (   approx_val_function_val_vec_ ); CHKERRQ(ierr);
	
	ierr = VecAssemblyBegin ( approx_grad_function_grad_vec_ ); CHKERRQ(ierr);
	ierr = VecAssemblyEnd (   approx_grad_function_grad_vec_ ); CHKERRQ(ierr);
	
	
	ierr = MatAssemblyBegin ( mass_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	ierr = MatAssemblyEnd (   mass_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin ( stiffness_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	ierr = MatAssemblyEnd (   stiffness_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	
	if ( bleach_matrix_created_ ) {
		ierr = MatAssemblyBegin ( bleach_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
		ierr = MatAssemblyEnd (   bleach_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	}
		
	if ( siggia_matrix_created_ ) {
		ierr = MatAssemblyBegin ( siggia_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
		ierr = MatAssemblyEnd (   siggia_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	}
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::load_mass_matrix ( double aa )
{
	assert ( petsc_objects_created_ );
	PetscErrorCode ierr;
	
// Keep for symmetry assert if wanted.
// 	if ( comm_size_==1 ) {
// 		PetscTruth test;
// 		MatIsSymmetric ( mass_matrix_, 1, &test );
// 		if ( PETSC_FALSE == test  ) {
// 			PetscPrintf ( PETSC_COMM_WORLD, "MASS MATRIX IS NOT SYMMETRIC!\n" );
// 		} else {
// 			PetscPrintf ( PETSC_COMM_WORLD, "Phew, mass matrix is at least symmetric.\n" );
// 		}
// 	}

#if 0 // Keep to check asymmetry.
	PetscPrintf ( PETSC_COMM_WORLD, "MASS MATRIX IS NOT SYMMETRIC?\n" );
	Mat    tmp_mat;
	MatTranspose(mass_matrix_,&tmp_mat);
	MatAXPY( tmp_mat, -1,mass_matrix_,SAME_NONZERO_PATTERN);
	MatView ( tmp_mat, PETSC_VIEWER_STDOUT_WORLD );
	MatDestroy ( tmp_mat );
#endif
	
	if ( ksp_matrix_populated_ ) {
		ierr = MatCopy ( mass_matrix_, ksp_matrix_, SAME_NONZERO_PATTERN ); CHKERRQ(ierr);
	} else {
		ierr = MatCopy ( mass_matrix_, ksp_matrix_, DIFFERENT_NONZERO_PATTERN ); CHKERRQ(ierr);
		ksp_matrix_populated_ = true;
	}
	
	if ( (aa>1.0 || aa<1.0) ) {
		ierr = MatScale ( ksp_matrix_, static_cast<PetscScalar>(aa) );
	}
	
	//PetscPrintf(PETSC_COMM_WORLD, "Hello ksp_matrix-.");
	//MatView ( ksp_matrix_, PETSC_VIEWER_STDOUT_WORLD );
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::load_stiffness_matrix ( double aa )
{
	assert ( petsc_objects_created_ );
	PetscErrorCode ierr;
	
	if ( ksp_matrix_populated_ ) {
		ierr = MatCopy ( stiffness_matrix_, ksp_matrix_, SAME_NONZERO_PATTERN ); CHKERRQ(ierr);
	} else {
		ierr = MatCopy ( stiffness_matrix_, ksp_matrix_, DIFFERENT_NONZERO_PATTERN ); CHKERRQ(ierr);
		ksp_matrix_populated_ = true;
	}
	
	if ( (aa>1.0 || aa<1.0) ) {
		ierr = MatScale ( ksp_matrix_, static_cast<PetscScalar>(aa) );
	}
	
	//PetscPrintf(PETSC_COMM_WORLD, "Hello ksp_matrix-.");
	//MatView ( ksp_matrix_, PETSC_VIEWER_STDOUT_WORLD );
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::add_scaled_stiffness_matrix ( double aa )
{
	assert ( petsc_objects_created_ );
	assert ( ksp_matrix_populated_  );
	PetscErrorCode ierr = MatAXPY ( ksp_matrix_,
	                                static_cast<PetscScalar>(aa),
	                                stiffness_matrix_,
	                                SAME_NONZERO_PATTERN );
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::add_scaled_bleach_matrix ( double aa )
{
	assert ( petsc_objects_created_ );
	assert ( bleach_matrix_created_ );
	PetscErrorCode ierr = MatAXPY ( ksp_matrix_,
	                                static_cast<PetscScalar>(aa),
	                                bleach_matrix_,
	                                DIFFERENT_NONZERO_PATTERN );
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::add_scaled_siggia_matrix ( double aa )
{
	assert ( siggia_matrix_created_ );
	PetscErrorCode ierr = MatAXPY ( ksp_matrix_,
	                                static_cast<PetscScalar>(aa),
	                                siggia_matrix_,
	                                SAME_NONZERO_PATTERN );
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::species_set_initial_proportions ( double* in )
{
	assert ( in );
	assert ( petsc_objects_created_ );
	assert ( num_species_ > 0 );
	
	PetscErrorCode ierr = 0;
	
	this->load_mass_matrix();
	this->solve(false);
	
	PetscReal a;
	
	for ( int i=0; i<num_species_; ++i ) {
		ierr = VecCopy ( ksp_solution_, species_previous_vec_[i] );CHKERRQ(ierr);
		a = in[i];
		ierr = VecScale ( species_previous_vec_[i], a );CHKERRQ(ierr);
	}
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::species_set_matrices ( species_system_data& species_sys_dat,
                                             bool                 spot_bleach_on,
                                             valarray<bool>*      species_mask )
{
	PetscErrorCode ierr;
	
	assert ( num_species_ > 0 );
	if ( num_species_ < 1 ) return 0;
	
	double time_step_dt = species_sys_dat.time_step_dt;
	
	double    leaving_species;
	for ( int i=0; i<num_species_; ++i ) {
		if ( species_mask && (false==(*species_mask)[i]) ) continue;
		
		leaving_species = 0.0;
		for ( int j=0; j<num_species_; ++j ) {
#if 0 // Removing volume adjust.
			if ( i>j ) {
				leaving_species += species_sys_dat.reaction_coefficients[i][j]*species_sys_dat.volume_adjust[i][j];
			} else {
				leaving_species += species_sys_dat.reaction_coefficients[i][j];
			}
#endif
			leaving_species += species_sys_dat.reaction_coefficients[i][j];
		}
		
		load_mass_matrix ( 1.0 +
		                   time_step_dt*(species_sys_dat.global_bleach_coefficient[i]
		                                 + leaving_species ) );
		
		ksp_operators_set_ = false;
		
		if ( species_sys_dat.diffusion_coefficient[i]>0.0 ) {
			add_scaled_stiffness_matrix ( time_step_dt*species_sys_dat.diffusion_coefficient[i] );
		}
		
		if ( spot_bleach_on && (species_sys_dat.spot_bleach_coefficient[i]>0.0) ) {
			add_scaled_bleach_matrix ( time_step_dt*species_sys_dat.spot_bleach_coefficient[i] );
		}
		
		assert ( species_matrices_vectors_created_ );
		if ( species_matrices_populated_[i] ) {
			ierr = MatCopy ( ksp_matrix_, species_matrices_[i], SAME_NONZERO_PATTERN ); CHKERRQ(ierr);
		} else {
			ierr = MatCopy ( ksp_matrix_, species_matrices_[i], DIFFERENT_NONZERO_PATTERN ); CHKERRQ(ierr);
			species_matrices_populated_[i] = true;
		}
	}
	
	return ierr;
}

template <int dim>
PetscErrorCode
pum_convergence<dim>::species_solve ( species_system_data& species_sys_dat )
{
	assert ( num_species_ > 0 );
	
	PetscErrorCode ierr = 0;
	
	PetscReal    a;
	
	for ( int i=0; i<num_species_; ++i ) {
#if 0 // Reduce to a single MatMult.
		ierr = MatMult ( mass_matrix_, species_previous_vec_[i], ksp_rhs_vec_ ); CHKERRQ(ierr);
		
		for ( int j=0; j<num_species_; ++j ) {
			if ( j==i ) continue;
			
			ierr = MatMult ( mass_matrix_, species_previous_vec_[j], calculation_vec_ );CHKERRQ(ierr);
			
			a = species_sys_dat.time_step_dt * species_sys_dat.reaction_coefficients[j][i];
			
			ierr = VecAXPY ( ksp_rhs_vec_, a, calculation_vec_ );CHKERRQ(ierr);
		}
#else
		ierr=VecCopy(species_previous_vec_[i],calculation_vec_);CHKERRQ(ierr);
		
		for ( int j=0; j<num_species_; ++j ) {
			if ( j==i ) continue;
			
			a = species_sys_dat.time_step_dt * species_sys_dat.reaction_coefficients[j][i];
			
			ierr=VecAXPY(calculation_vec_, a, species_previous_vec_[j]);
		}
		
		if ( 0==i ) {
			a = species_sys_dat.time_step_dt * species_sys_dat.coupling_param;
			ierr=VecAXPY(calculation_vec_, -a, species_previous_vec_[0]);CHKERRQ(ierr);
			ierr=VecAXPY(calculation_vec_,  a, species_previous_vec_[2]);CHKERRQ(ierr);
		} else if ( 2==i ) {
			a = species_sys_dat.time_step_dt * species_sys_dat.coupling_param;
			a *= species_sys_dat.focal_percentage/(100.0-species_sys_dat.focal_percentage);
			ierr=VecAXPY(calculation_vec_,  a, species_previous_vec_[0]);CHKERRQ(ierr);
			ierr=VecAXPY(calculation_vec_, -a, species_previous_vec_[2]);CHKERRQ(ierr);
		}
		
		ierr=MatMult(mass_matrix_,calculation_vec_,ksp_rhs_vec_);CHKERRQ(ierr);
#endif
		
		assert ( species_matrices_populated_[i] );
		
		ierr = MatCopy ( species_matrices_[i], ksp_matrix_, SAME_NONZERO_PATTERN ); CHKERRQ(ierr);
		
		if ( !ksp_operators_set_ ) {
			ierr = KSPSetOperators(ksp_, ksp_matrix_, ksp_matrix_, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
			ksp_operators_set_ = true;
		}
		
		ierr = KSPSolve ( ksp_, ksp_rhs_vec_, species_current_vec_[i] ); CHKERRQ(ierr);
	}
	
	for ( int i=0; i<num_species_; ++i ) {
		ierr = VecCopy ( species_current_vec_[i], species_previous_vec_[i] );CHKERRQ(ierr);
	}
	return ierr;
}

template <int dim>
void
pum_convergence<dim>::get_dof_structure ( dof_structure::ptr& ds ) const
{
	if ( !ds ) {
		ds.reset ( new dof_structure );
	}
	
	assert ( ds );
	
	ds->patch_to_num_dof.resize( patch_to_num_dof_.size() );
	ds->patch_to_num_dof = patch_to_num_dof_;
	
	ds->patch_to_start_index.resize( patch_to_start_index_.size() );
	ds->patch_to_start_index = patch_to_start_index_;
	
	ds->computenode_to_num_patches.resize( computenode_to_num_patches_.size() );
	ds->computenode_to_num_patches = computenode_to_num_patches_;
	
	ds->computenode_to_start_patch.resize( computenode_to_start_patch_.size() );
	ds->computenode_to_start_patch = computenode_to_start_patch_;
	
	ds->computenode_to_num_dof.resize( computenode_to_num_dof_.size() );
	ds->computenode_to_num_dof = computenode_to_num_dof_;
	
	ds->patch_indices_on_computenode.resize( patch_indices_on_computenode_.size() );
	ds->patch_indices_on_computenode = patch_indices_on_computenode_;
	
	ds->computenode_to_start_dof.resize( computenode_to_start_dof_.size() );
	ds->computenode_to_start_dof = computenode_to_start_dof_;
	
}


template <int dim>
PetscErrorCode
pum_convergence<dim>::solve ( solution<dim>& sol, bool output_iter )
{
	PetscErrorCode ierr;
	
	if ( !ksp_operators_set_ ) {
		ierr = KSPSetOperators(ksp_, ksp_matrix_, ksp_matrix_, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
		ksp_operators_set_ = true;
	}
	
	ierr = KSPSolve ( ksp_, ksp_rhs_vec_, ksp_solution_ ); CHKERRQ(ierr);

#if 0// Keep for debug.
	if ( computenode_to_num_dof_.sum() < 50){
		PetscPrintf( PETSC_COMM_WORLD, "Debug slight imperfections. Mass matrix:\n");
		MatView ( mass_matrix_, PETSC_VIEWER_STDOUT_WORLD );
		PetscPrintf( PETSC_COMM_WORLD, "Debug slight imperfections. RHS vector:\n");
		VecView ( ksp_rhs_vec_, PETSC_VIEWER_STDOUT_WORLD );
		PetscPrintf( PETSC_COMM_WORLD, "Debug slight imperfections. Solution:\n");
		VecView ( ksp_solution_, PETSC_VIEWER_STDOUT_WORLD );
		PetscPrintf( PETSC_COMM_WORLD, "Debug Helmholtz. Stiffness matrix:\n");
		MatView ( stiffness_matrix_, PETSC_VIEWER_STDOUT_WORLD );
	}
#endif

	if ( output_iter ) {
		int iter;
		
		ierr = KSPGetIterationNumber ( ksp_, &iter ); CHKERRQ(ierr);
		ierr = PetscPrintf ( PETSC_COMM_WORLD, "KSPSolve : %D iterations.\n", iter ); CHKERRQ(ierr);
	}
	
	int sendcount = computenode_to_num_dof_[comm_rank_];
	
#ifndef NDEBUG
	{
		int low, high;
		VecGetOwnershipRange ( ksp_solution_, &low, &high );
		assert ( high - low == sendcount );
	}
#endif
	
	// PetscScalar must be double.
	double* sendbuf = 0;
	
	ierr = VecGetArray ( ksp_solution_, &sendbuf ); CHKERRQ(ierr);
	
	if ( comm_size_ == 1 ) {
		assert ( comm_rank_ == 0 );
		assert ( sendbuf );
		assert ( sol.global_coefficients_.size() == computenode_to_num_dof_[0] );
		copy ( sendbuf, sendbuf + computenode_to_num_dof_[0], &sol.global_coefficients_[0] );
	} else {
		// void* sendbuf, int sendcount, MPI_Datatype sendtype,
		// void* recvbuf
		// int* recvcounts
		// int* displs

		MPI_Allgatherv ( sendbuf, sendcount, MPI_DOUBLE,
				&sol.global_coefficients_[0],
				&computenode_to_num_dof_[0],
				&computenode_to_start_dof_[0],
				MPI_DOUBLE,
				PETSC_COMM_WORLD );
	}

	ierr = VecRestoreArray ( ksp_solution_, &sendbuf ); CHKERRQ(ierr);
	
	{
		dof_structure::ptr tmp_dof_st;
		get_dof_structure( tmp_dof_st );
		assert ( tmp_dof_st );
		sol_->set_dof_structure( tmp_dof_st );
	}
	
	return ierr;
}


template <int dim>
void
pum_convergence<dim>::set_sum_all_pu_test ( int step )
{
	valarray<double>    tmp ( 0.0, computenode_to_num_dof_.sum() );
	int num_patches = computenode_to_num_patches_.sum();
	for ( int i=0; i<num_patches; ++i ) {
		tmp[step*i] = 1.0;
	}
	
	sol_->set_coefficients ( tmp );
}

template class pum_convergence<2>;

