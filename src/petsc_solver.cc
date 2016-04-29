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

#include "petsc_solver.hh"

#include <algorithm>
#include <iostream>
#include <cassert>

using std::clog;
using std::copy;
using std::cout;

template<int dim>
petsc_solver<dim>::petsc_solver()
:
petsc_objects_created_(false),
mass_matrix_created_(false),
bleach_matrix_created_(false),
stiffness_matrix_created_(false),
siggia_matrix_created_(false),
solution1_created_(false),
domain_integration_vec_created_(false),
bleach_integration_vec_created_(false),
ksp_operators_set_(false)
{
	MPI_Comm_rank( PETSC_COMM_WORLD, &comm_rank_ );
	MPI_Comm_size( PETSC_COMM_WORLD, &comm_size_ );
}

template<int dim>
petsc_solver<dim>::~petsc_solver()
{
	// Destructors must not throw or otherwise fail!
	(void) destroy_petsc_objects();
}

template<int dim>
PetscErrorCode
petsc_solver<dim>::create_petsc_objects( const valarray<int>& patch_to_num_dof,
                                         int                  row_nonzero )
{
	assert( !petsc_objects_created_ );

	if ( petsc_objects_created_ ) {
		destroy_petsc_objects();
	}
	
	/**
	 * Keep a copy of patch_to_num_dof. (1)
	 */
	patch_to_num_dof_.resize( patch_to_num_dof.size() );
	patch_to_num_dof_ = patch_to_num_dof;
	
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
		int global_size = computenode_to_num_dof_.sum();
		
		// cout << "Mike partition is " << computenode_to_num_dof_[comm_rank_] << "\n";
		
		PetscErrorCode ierr;
		
		ierr = VecCreate(PETSC_COMM_WORLD, &petsc_solution_); CHKERRQ(ierr);
		//ierr = VecSetSizes(petsc_solution_, computenode_to_num_dof_[comm_rank_], PETSC_DECIDEglobal_size); CHKERRQ(ierr);
		ierr = VecSetSizes(petsc_solution_, computenode_to_num_dof_[comm_rank_], global_size); CHKERRQ(ierr);
		ierr = VecSetFromOptions(petsc_solution_); CHKERRQ(ierr);
		
#if 0 // Profiling motivated. 2008-12-04 Michael LI.
		ierr = MatCreate(PETSC_COMM_WORLD, &petsc_matrix_); CHKERRQ(ierr);
		ierr = MatSetSizes( petsc_matrix_,
		                    computenode_to_num_dof_[comm_rank_],
		                    computenode_to_num_dof_[comm_rank_],
		                    global_size,
		                    global_size ); CHKERRQ(ierr);
#endif
		ierr = MatCreateMPIAIJ ( PETSC_COMM_WORLD,
		                         computenode_to_num_dof_[comm_rank_],
		                         computenode_to_num_dof_[comm_rank_],
		                         global_size,
		                         global_size,
		                         row_nonzero,
		                         PETSC_NULL,
		                         0,
		                         PETSC_NULL,
		                         &petsc_matrix_ );
		ierr = MatSetFromOptions(petsc_matrix_); CHKERRQ(ierr);

		// petsc-3.0.0 doesn't have this option any more.
		// 2009-02-18 ML.
		//		ierr = MatSetOption( petsc_matrix_, MAT_ROWS_SORTED ); CHKERRQ(ierr);

#if 0 // Profiling motivated. 2008-12-04 Michael LI.
		ierr = MatCreate(PETSC_COMM_WORLD, &mass_matrix_); CHKERRQ(ierr);
		ierr = MatSetSizes( mass_matrix_,
		                    computenode_to_num_dof_[comm_rank_],
		                    computenode_to_num_dof_[comm_rank_],
		                    global_size,
		                    global_size ); CHKERRQ(ierr);
#endif
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

		// petsc-3.0.0 doesn't have this option any more.
		// 2009-02-18 ML.
		//		ierr = MatSetOption( mass_matrix_, MAT_ROWS_SORTED ); CHKERRQ(ierr);
		mass_matrix_created_ = true;
	
#if 0 // Profiling motivated. 2008-12-04 Michael LI.
		ierr = MatCreate(PETSC_COMM_WORLD, &stiffness_matrix_); CHKERRQ(ierr);
		ierr = MatSetSizes( stiffness_matrix_,
		                    computenode_to_num_dof_[comm_rank_],
		                    computenode_to_num_dof_[comm_rank_],
		                    global_size,
		                    global_size ); CHKERRQ(ierr);
#endif
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

		// petsc-3.0.0 doesn't have this option any more.
		// 2009-02-18 ML.
		//		ierr = MatSetOption( stiffness_matrix_, MAT_ROWS_SORTED ); CHKERRQ(ierr);
		stiffness_matrix_created_ = true;
		
#if 0 // Profiling motivated. 2008-12-04 Michael LI.
		ierr = MatCreate(PETSC_COMM_WORLD, &bleach_matrix_); CHKERRQ(ierr);
		ierr = MatSetSizes( bleach_matrix_,
		                    computenode_to_num_dof_[comm_rank_],
		                    computenode_to_num_dof_[comm_rank_],
		                    global_size,
		                    global_size ); CHKERRQ(ierr);
#endif
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

		// petsc-3.0.0 doesn't have this option any more.
		// 2009-02-18 ML.
		//		ierr = MatSetOption( bleach_matrix_, MAT_ROWS_SORTED ); CHKERRQ(ierr);
		bleach_matrix_created_ = true;
		
#if 0 // Profiling motivated. 2008-12-04 Michael LI.
		ierr = MatCreate(PETSC_COMM_WORLD, &siggia_matrix_); CHKERRQ(ierr);
		ierr = MatSetSizes( siggia_matrix_,
		                    computenode_to_num_dof_[comm_rank_],
		                    computenode_to_num_dof_[comm_rank_],
		                    global_size,
		                    global_size ); CHKERRQ(ierr);
#endif
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

		// petsc-3.0.0 doesn't have this option any more.
		// 2009-02-18 ML.
		//		ierr = MatSetOption( siggia_matrix_, MAT_ROWS_SORTED ); CHKERRQ(ierr);
		siggia_matrix_created_ = true;

		
		ierr = VecDuplicate(petsc_solution_, &rhs_vec_); CHKERRQ(ierr);
		
		ierr = KSPCreate(PETSC_COMM_WORLD, &ksp_); CHKERRQ(ierr);
				
		// Attempt to reduce number of iterations increasing too much as system size increase.
		// e.g.
		// level 6: 12 iterations
		// level 7: 21 iterations
		// level 8: 33 iterations taking 24% run time.
		ierr = KSPSetInitialGuessNonzero ( ksp_, PETSC_TRUE ); CHKERRQ(ierr);
				
		// The user may use -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
		ierr = KSPSetFromOptions(ksp_); CHKERRQ(ierr);
		
		petsc_objects_created_ = true;

#if 1 // ndef NDEBUG
		{
			int first_row, last_row;
			MatGetOwnershipRange ( petsc_matrix_, &first_row, &last_row );
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
		
		return ierr;
	}

}

template<int dim>
PetscErrorCode
petsc_solver<dim>::destroy_petsc_objects()
{
    PetscErrorCode ierr;
	if ( petsc_objects_created_ ) {
		ierr = VecDestroy(petsc_solution_); CHKERRQ( ierr );
		ierr = MatDestroy(petsc_matrix_); CHKERRQ( ierr );
		ierr = VecDestroy(rhs_vec_); CHKERRQ( ierr );
		ierr = KSPDestroy(ksp_); CHKERRQ(ierr);
		ksp_operators_set_ = false;
		petsc_objects_created_ = false;
	}
	if ( mass_matrix_created_ ) {
		ierr = MatDestroy ( mass_matrix_ );  CHKERRQ(ierr);
		mass_matrix_created_ = false;
	}
	if ( bleach_matrix_created_ ) {
		ierr = MatDestroy ( bleach_matrix_ ); CHKERRQ(ierr);
		bleach_matrix_created_ = false;
	}
	if ( stiffness_matrix_created_ ) {
		ierr = MatDestroy ( stiffness_matrix_ ); CHKERRQ(ierr);
		stiffness_matrix_created_ = false;
	}
	if ( siggia_matrix_created_ ) {
		ierr = MatDestroy ( siggia_matrix_ ); CHKERRQ(ierr);
		siggia_matrix_created_ = false;
	}
	if ( solution1_created_ ) {
		ierr = VecDestroy ( solution1_ ); CHKERRQ(ierr);
		solution1_created_ = false;
	}
	if ( domain_integration_vec_created_ ) {
		ierr = VecDestroy ( domain_integration_vec_ ); CHKERRQ(ierr);
		domain_integration_vec_created_ = false;
	}
	if ( bleach_integration_vec_created_ ) {
		ierr = VecDestroy ( bleach_integration_vec_ ); CHKERRQ(ierr);
		bleach_integration_vec_created_ = false;
	}
	return ierr;
}

template <int dim>
void
petsc_solver<dim>::get_patch_indices_on_computenode ( valarray<int>& out ) const
{
	out.resize( patch_indices_on_computenode_.size() );
	out = patch_indices_on_computenode_;
}

// Const qualifier not present on vals due to use of pg compiler.
// If vals is const, it must be copied before accessing contents with pointer.
// 2008-08-20 Mike Li.
template <int dim>
PetscErrorCode
petsc_solver<dim>::set_matrix_entries ( int first, int second, valarray<double>& vals )
{
	
	int num_dof_first  = patch_to_num_dof_[first];
	int num_dof_second = patch_to_num_dof_[second];
	
	int start_first  = patch_to_start_index_[first];
	int start_second = patch_to_start_index_[second];

	valarray<int> idxm( num_dof_first );
	idxm[0] = start_first;
	for ( int i=1; i<num_dof_first; ++i ) {
		idxm[i] = start_first + i;
	}
	
	valarray<int> idxn( num_dof_second );
	idxn[0] = start_second;
	for ( int i=1; i<num_dof_second; ++i ) {
		idxn[i] = start_second + i;
	}
	//#ifndef NDEBUG
	//std::clog << "\t\t" << vals << "\n";
	//#endif

	assert ( vals.size() == num_dof_first*num_dof_second );

//	valarray<double> tmp ( vals.size() );
//	tmp = vals;

	PetscErrorCode ierr;
#pragma omp critical
	ierr = MatSetValues( petsc_matrix_,
	                                    num_dof_first,    &idxm[0],
	                                    num_dof_second,   &idxn[0],
	                                    &vals[0],
	                                    INSERT_VALUES );
	
	CHKERRQ(ierr);
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::set_matrix_entries ( int                   first,
                                        int                   second,
                                        valarray<double> &    mass_matrix_vals,
                                        valarray<double> &    bleach_matrix_vals,
                                        valarray<double> &    stiffness_matrix_vals,
                                        valarray<double> *    siggia_matrix_vals )
{
	assert ( mass_matrix_created_ );
	assert ( bleach_matrix_created_ );
	assert ( stiffness_matrix_created_ );
	assert ( siggia_matrix_created_ );
	
	int num_dof_first  = patch_to_num_dof_[first];
	int num_dof_second = patch_to_num_dof_[second];
	
	int start_first  = patch_to_start_index_[first];
	int start_second = patch_to_start_index_[second];

	valarray<int> idxm( num_dof_first );
	idxm[0] = start_first;
	for ( int i=1; i<num_dof_first; ++i ) {
		idxm[i] = start_first + i;
	}
	
	valarray<int> idxn( num_dof_second );
	idxn[0] = start_second;
	for ( int i=1; i<num_dof_second; ++i ) {
		idxn[i] = start_second + i;
	}

	assert ( mass_matrix_vals.size() == num_dof_first*num_dof_second );
	assert ( bleach_matrix_vals.size() == num_dof_first*num_dof_second );
	assert ( stiffness_matrix_vals.size() == num_dof_first*num_dof_second );
	
#ifndef NDEBUG
	if ( siggia_matrix_vals ) {
		assert ( siggia_matrix_vals->size() == num_dof_first*num_dof_second );
	}
#endif

	PetscErrorCode outer_ierr;


#pragma omp critical
	{
		PetscErrorCode ierr = MatSetValues( mass_matrix_,
	                                    num_dof_first,    &idxm[0],
	                                    num_dof_second,   &idxn[0],
	                                    &mass_matrix_vals[0],
	                                    INSERT_VALUES );
		outer_ierr = ierr;
	}
	
#pragma omp critical
	{
		PetscErrorCode ierr = MatSetValues( bleach_matrix_,
	                                    num_dof_first,    &idxm[0],
	                                    num_dof_second,   &idxn[0],
	                                    &bleach_matrix_vals[0],
	                                    INSERT_VALUES );
		outer_ierr = ierr;
	}
	
#pragma omp critical
	{
		PetscErrorCode ierr = MatSetValues( stiffness_matrix_,
	                                    num_dof_first,    &idxm[0],
	                                    num_dof_second,   &idxn[0],
	                                    &stiffness_matrix_vals[0],
	                                    INSERT_VALUES );
		outer_ierr = ierr;
	}
	
	if ( siggia_matrix_vals ) {
		// Valgrind helped to debug problem with blow up of the solution
		// using openmp but disappears without openmp.
		// 2008-12-05 Michael LI.
#pragma omp critical
		{
			PetscErrorCode ierr = MatSetValues( siggia_matrix_,
		                                    num_dof_first,    &idxm[0],
		                                    num_dof_second,   &idxn[0],
	        	                            &(siggia_matrix_vals->operator[](0)),
		                                    INSERT_VALUES );
			outer_ierr = ierr;
		}
	}
	
	return outer_ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::matrix_assembly_begin ()
{
	PetscErrorCode ierr = MatAssemblyBegin( petsc_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::matrix_assembly_end ()
{
	PetscErrorCode ierr = MatAssemblyEnd( petsc_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	//#ifndef NDEBUG
	if ( computenode_to_num_dof_.sum() < 50 ) {
		MatView ( petsc_matrix_, PETSC_VIEWER_STDOUT_WORLD );
	}
	//#endif
	return ierr;
}

#if 0 // Appears unused - will not update for siggia_matrix. 2008-12-04 Michael LI.
template <int dim>
PetscErrorCode
petsc_solver<dim>::quick_setup_matrices_after_first_assembly ()
{
	PetscErrorCode ierr;
	
	if ( !mass_matrix_created_ ) {
		ierr = MatDuplicate ( petsc_matrix_, MAT_DO_NOT_COPY_VALUES, &mass_matrix_ ); CHKERRQ(ierr);
		mass_matrix_created_ = true;
	}
	
	if ( !stiffness_matrix_created_ ) {
		ierr = MatDuplicate ( petsc_matrix_, MAT_DO_NOT_COPY_VALUES, &stiffness_matrix_ ); CHKERRQ(ierr);
		stiffness_matrix_created_ = true;
	}
	
	if ( !bleach_matrix_created_ ) {
		ierr = MatDuplicate ( petsc_matrix_, MAT_DO_NOT_COPY_VALUES, &bleach_matrix_ ); CHKERRQ(ierr);
		bleach_matrix_created_ = true;
	}
	
	return ierr;
}
#endif

template <int dim>
PetscErrorCode
petsc_solver<dim>::quick_assemble_matrices ()
{
	PetscErrorCode ierr;
	ierr = MatAssemblyBegin ( mass_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	ierr = MatAssemblyEnd (   mass_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin ( bleach_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	ierr = MatAssemblyEnd (   bleach_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin ( stiffness_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	ierr = MatAssemblyEnd (   stiffness_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	
	ierr = MatAssemblyBegin ( siggia_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	ierr = MatAssemblyEnd (   siggia_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::load_mass_matrix ( double aa )
{
	assert ( mass_matrix_created_ );
	PetscErrorCode ierr;
	static bool first_load = true;
	if ( first_load ) {
		ierr = MatCopy ( mass_matrix_, petsc_matrix_, DIFFERENT_NONZERO_PATTERN ); CHKERRQ(ierr);
		first_load = false;
	} else {
		ierr = MatCopy ( mass_matrix_, petsc_matrix_, SAME_NONZERO_PATTERN ); CHKERRQ(ierr);
	}
	
	if ( (aa>1.0 || aa<1.0) ) {
		ierr = MatScale ( petsc_matrix_, static_cast<PetscScalar>(aa) );
	}
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::add_scaled_stiffness_matrix ( double aa )
{
	assert ( stiffness_matrix_created_ );
	PetscErrorCode ierr = MatAXPY ( petsc_matrix_,
	                                static_cast<PetscScalar>(aa),
	                                stiffness_matrix_,
	                                SAME_NONZERO_PATTERN );
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::add_scaled_bleach_matrix ( double aa )
{
	assert ( bleach_matrix_created_ );
	PetscErrorCode ierr = MatAXPY ( petsc_matrix_,
	                                static_cast<PetscScalar>(aa),
	                                bleach_matrix_,
	                                SAME_NONZERO_PATTERN );
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::add_scaled_siggia_matrix ( double aa )
{
	assert ( siggia_matrix_created_ );
	PetscErrorCode ierr = MatAXPY ( petsc_matrix_,
	                                static_cast<PetscScalar>(aa),
	                                siggia_matrix_,
	                                SAME_NONZERO_PATTERN );
	
	ksp_operators_set_ = false;
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::set_rhs_vector_entries ( int patch_index, valarray<PetscScalar>& vals )
{
	int num_entries = patch_to_num_dof_[patch_index];
	valarray<int> indices( num_entries );
	int start_patch_index = patch_to_start_index_[patch_index];
	for ( int i=0; i<num_entries; ++i ) {
		indices[i] = start_patch_index + i;
	}
	
	PetscErrorCode ierr;
	
#pragma omp critical
	ierr = VecSetValues( rhs_vec_, num_entries, &indices[0], &vals[0], INSERT_VALUES );
	
	CHKERRQ(ierr);
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::rhs_vector_assembly_begin ()
{
	PetscErrorCode ierr = VecAssemblyBegin( rhs_vec_ ); CHKERRQ(ierr);
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::rhs_vector_assembly_end ()
{
	PetscErrorCode ierr = VecAssemblyEnd( rhs_vec_ ); CHKERRQ(ierr);
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::save_as_mass_matrix ()
{
	PetscErrorCode ierr;
	
	if ( !mass_matrix_created_ ) {
		ierr = MatDuplicate ( petsc_matrix_, MAT_COPY_VALUES, &mass_matrix_ ); CHKERRQ(ierr);
		mass_matrix_created_ = true;
	} else {
		ierr = MatCopy ( petsc_matrix_, mass_matrix_, DIFFERENT_NONZERO_PATTERN ); CHKERRQ(ierr);
		// ierr = MatAssemblyBegin ( mass_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
		// ierr = MatAssemblyEnd (   mass_matrix_, MAT_FINAL_ASSEMBLY ); CHKERRQ(ierr);
	}
	
	// PetscPrintf ( PETSC_COMM_WORLD, "save_as_mass_matrix MatView\n" );
	// MatView ( mass_matrix_, PETSC_VIEWER_STDOUT_WORLD );
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::save_as_stiffness_matrix ()
{
	PetscErrorCode ierr;
	
	if ( !stiffness_matrix_created_ ) {
		ierr = MatDuplicate ( petsc_matrix_, MAT_COPY_VALUES, &stiffness_matrix_ );
		CHKERRQ(ierr);
		stiffness_matrix_created_ = true;
		return ierr;
	} else {
		ierr = MatCopy ( petsc_matrix_, stiffness_matrix_, DIFFERENT_NONZERO_PATTERN );
		CHKERRQ(ierr);
	}
	
	// PetscPrintf ( PETSC_COMM_WORLD, "save_as_stiffness_matrix MatView.\n" );
	// MatView ( petsc_matrix_, PETSC_VIEWER_STDOUT_WORLD );
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::save_solution1 ()
{
	PetscErrorCode ierr;
	
	if ( !solution1_created_ ) {
		ierr = VecDuplicate ( petsc_solution_, &solution1_ ); CHKERRQ(ierr);
		solution1_created_ = true;
	}
	
	ierr = VecCopy ( petsc_solution_, solution1_ ); CHKERRQ(ierr);
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::save_domain_integration_vec ()
{
	PetscErrorCode ierr;
	
	if ( !domain_integration_vec_created_ ) {
		ierr = VecDuplicate ( rhs_vec_, &domain_integration_vec_ ); CHKERRQ(ierr);
		domain_integration_vec_created_ = true;
	}
	
	ierr = VecCopy ( rhs_vec_, domain_integration_vec_ ); CHKERRQ(ierr);
	
	return ierr;
}


template <int dim>
PetscErrorCode
petsc_solver<dim>::save_bleach_integration_vec ()
{
	PetscErrorCode ierr;
	
	if ( !bleach_integration_vec_created_ ) {
		ierr = VecDuplicate ( rhs_vec_, &bleach_integration_vec_ ); CHKERRQ(ierr);
		bleach_integration_vec_created_ = true;
	}
	
	ierr = VecCopy ( rhs_vec_, bleach_integration_vec_ ); CHKERRQ(ierr);
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::view_petsc_matrix ( const char* filename )
{
	PetscErrorCode ierr;
	if ( filename == 0 ) {
		ierr = MatView ( petsc_matrix_, PETSC_VIEWER_STDOUT_WORLD ); CHKERRQ(ierr);
	} else {
		PetscViewer ascii_file;
		PetscViewerASCIIOpen ( PETSC_COMM_WORLD, filename, &ascii_file );
		ierr = MatView ( petsc_matrix_, ascii_file ); CHKERRQ(ierr);
		PetscViewerDestroy ( ascii_file );
	}
	return ierr;
}


template <int dim>
PetscErrorCode
petsc_solver<dim>::view_mass_matrix ( const char* filename )
{
	PetscErrorCode ierr;
	if ( filename == 0 ) {
		ierr = MatView ( mass_matrix_, PETSC_VIEWER_STDOUT_WORLD ); CHKERRQ(ierr);
	} else {
		PetscViewer ascii_file;
		PetscViewerASCIIOpen ( PETSC_COMM_WORLD, filename, &ascii_file );
		ierr = MatView ( mass_matrix_, ascii_file ); CHKERRQ(ierr);
		PetscViewerDestroy ( ascii_file );
	}
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::view_stiffness_matrix ( const char* filename )
{
	PetscErrorCode ierr;
	if ( filename == 0 ) {
		ierr = MatView ( stiffness_matrix_, PETSC_VIEWER_STDOUT_WORLD ); CHKERRQ(ierr);
	} else {
		PetscViewer ascii_file;
		PetscViewerASCIIOpen ( PETSC_COMM_WORLD, filename, &ascii_file );
		ierr = MatView ( stiffness_matrix_, ascii_file ); CHKERRQ(ierr);
		PetscViewerDestroy ( ascii_file );
	}
	
	return ierr;
}


template <int dim>
PetscErrorCode
petsc_solver<dim>::view_bleach_matrix ( const char* filename )
{
	PetscErrorCode ierr;
	if ( filename == 0 ) {
		ierr = MatView ( bleach_matrix_, PETSC_VIEWER_STDOUT_WORLD ); CHKERRQ(ierr);
	} else {
		PetscViewer ascii_file;
		PetscViewerASCIIOpen ( PETSC_COMM_WORLD, filename, &ascii_file );
		ierr = MatView ( bleach_matrix_, ascii_file ); CHKERRQ(ierr);
		PetscViewerDestroy ( ascii_file );
	}
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::view_solution ( const char* filename )
{
	PetscErrorCode ierr;
	if ( filename == 0 ) {
		ierr = VecView ( petsc_solution_, PETSC_VIEWER_STDOUT_WORLD ); CHKERRQ(ierr);
	} else {
		PetscViewer ascii_file;
		PetscViewerASCIIOpen ( PETSC_COMM_WORLD, filename, &ascii_file );
		ierr = VecView ( petsc_solution_, ascii_file ); CHKERRQ(ierr);
		PetscViewerDestroy ( ascii_file );
	}
	return ierr;
}



template <int dim>
double
petsc_solver<dim>::integrate_solution_over_domain ()
{
	assert ( domain_integration_vec_created_ );
	double tmp;
	PetscErrorCode ierr = VecDot ( petsc_solution_, domain_integration_vec_, &tmp ); CHKERRQ(ierr);
	return tmp;
}

template <int dim>
double
petsc_solver<dim>::integrate_solution_over_bleach_region ()
{
	assert ( bleach_integration_vec_created_ );
	double tmp;
	PetscErrorCode ierr = VecDot ( petsc_solution_, bleach_integration_vec_, &tmp ); CHKERRQ(ierr);
	return tmp;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::scale_current_and_add_mass_matrix ( double aa )
{
	assert ( mass_matrix_created_ );
	// Y = a Y + X
	PetscPrintf ( PETSC_COMM_WORLD, "PetscScalar a is %f.\n", aa );
	PetscErrorCode ierr = MatAYPX ( petsc_matrix_,
	                                static_cast<PetscScalar>(aa),
	                                mass_matrix_,
	                                SAME_NONZERO_PATTERN );
	
	// PetscPrintf ( PETSC_COMM_WORLD, "scale_current_and_add_mass_matrix MatView\n" );
	// MatView ( petsc_matrix_, PETSC_VIEWER_STDOUT_WORLD );
	
	CHKERRQ(ierr);
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::zero_rhs_vec ()
{
	assert ( petsc_objects_created_ );
	
	PetscErrorCode ierr = VecZeroEntries ( rhs_vec_ ); CHKERRQ(ierr);
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::zero_petsc_solution ()
{
	assert ( petsc_objects_created_ );
	
	PetscErrorCode ierr = VecZeroEntries ( petsc_solution_ ); CHKERRQ(ierr);
	
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::rhs_vec_add_mass_matrix_solution1 ()
{
	assert ( petsc_objects_created_ );
	assert ( mass_matrix_created_ );
	assert ( solution1_created_ );
	
	PetscErrorCode ierr;
	ierr = MatMultAdd ( mass_matrix_, solution1_, rhs_vec_, rhs_vec_ ); CHKERRQ(ierr);
	return ierr;
}

template <int dim>
PetscErrorCode
petsc_solver<dim>::rhs_vec_mass_mult_previous_solution ()
{
	assert ( petsc_objects_created_ );
	assert ( mass_matrix_created_ );
	
	PetscErrorCode ierr;
	ierr = MatMult ( mass_matrix_, petsc_solution_, rhs_vec_ ); CHKERRQ(ierr);
	return ierr;
}

template <int dim>
void
petsc_solver<dim>::get_dof_structure ( dof_structure::ptr& ds ) const
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
petsc_solver<dim>::solve ( solution<dim>& sol, bool output_iter )
{
	assert ( sol.global_approx_space_ && "Have global approximation space." );
	assert ( sol.global_coefficients_.size() == computenode_to_num_dof_.sum() );
	
	PetscErrorCode ierr;
	
	// Moved from creation routine.
	if ( !ksp_operators_set_ ) {
		ierr = KSPSetOperators(ksp_, petsc_matrix_, petsc_matrix_, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
		// ierr = KSPSetOperators(ksp_, petsc_matrix_, petsc_matrix_, SAME_PRECONDITIONER); CHKERRQ(ierr);
		ksp_operators_set_ = true;
	}
	
	// Keep for future debugging distribution problem at level 8.
	// VecSet ( petsc_solution_, 1.0 );
	


#pragma omp critical
	ierr = KSPSolve ( ksp_, rhs_vec_, petsc_solution_ );
	
// Used for validating behaviour of pum_convergence-1. Keep.
// 	static bool once = false;
// 	if ( !once ) {
// 		PetscPrintf ( PETSC_COMM_WORLD, "Comparison with pum_convergence-1.\n" );
// 		MatView ( petsc_matrix_, PETSC_VIEWER_STDOUT_WORLD );
// 		VecView ( rhs_vec_, PETSC_VIEWER_STDOUT_WORLD );
// 		VecView ( petsc_solution_, PETSC_VIEWER_STDOUT_WORLD );
// 		once = true;
// 	}
	
	CHKERRQ(ierr);
	if ( output_iter ) {
		int iter;
		
		ierr = KSPGetIterationNumber ( ksp_, &iter ); CHKERRQ(ierr);
		ierr = PetscPrintf ( PETSC_COMM_WORLD, "KSPSolve : %D iterations.\n", iter ); CHKERRQ(ierr);
	}
	
	int sendcount = computenode_to_num_dof_[comm_rank_];
	
#ifndef NDEBUG
	{
		int low, high;
		VecGetOwnershipRange ( petsc_solution_, &low, &high );
		assert ( high - low == sendcount );
	}
#endif
	
	// PetscScalar must be double.
	double* sendbuf = 0;
	
	ierr = VecGetArray ( petsc_solution_, &sendbuf ); CHKERRQ(ierr);
	
	if ( comm_size_ == 1 ) {
		// cout << "Begin process of copying solution.\n";
		assert ( comm_rank_ == 0 );
		assert ( sendbuf );
		copy ( sendbuf, sendbuf + computenode_to_num_dof_[0], &sol.global_coefficients_[0] );
	} else {
		// cout << "Begin process of distributing solution.\n";
		
		// void* sendbuf, int sendcount, MPI_Datatype sendtype,
		// void* recvbuf
		// int* recvcounts
		// int* displs
#pragma omp critical
		MPI_Allgatherv ( sendbuf, sendcount, MPI_DOUBLE,
				&sol.global_coefficients_[0],
				&computenode_to_num_dof_[0],
				&computenode_to_start_dof_[0],
				MPI_DOUBLE,
				PETSC_COMM_WORLD );
	}

	ierr = VecRestoreArray ( petsc_solution_, &sendbuf ); CHKERRQ(ierr);
	
	{
		dof_structure::ptr tmp_dof_st;
		get_dof_structure( tmp_dof_st );
		assert ( tmp_dof_st );
		sol.set_dof_structure( tmp_dof_st );
	}
	
	return ierr;
}

//

template class petsc_solver<2>;

#if 0

template<int dim>
void
petsc_solver<dim>::set_pum_discretization( shared_ptr<pum_discretization<dim> >& d )
{
	assert( d );
	dis = d;
	assert( dis );
	have_discretization = true;
}

template<int dim>
void
petsc_solver<dim>::update()
{
	assert( have_discretization );
	if ( !petsc_objects_created_ ) {
		create_petsc_objects();
	}
	updated = true;
}
template<int dim>
void
petsc_solver<dim>::assemble_for_initial_projection()
{
	assert(dis);
	assert(dis->is_updated());
	assert(updated);

	dis->assemble_petsc_matrix_( petsc_matrix_ );
	dis->assemble_rhs_vec_tor(rhs_vec_);

	MPI_Barrier( PETSC_COMM_WORLD );
}

template<int dim>
void
petsc_solver<dim>::species_set_function( const string& species, const pfunction& species_func )
{
	assert( dis->is_species(species) );
	assert( species_func );
	assert( species_func->is_global_function() );

	species_functions[species] = species_func;
}

template<int dim>
void
petsc_solver<dim>::solve_initial_projection()
{
	assert( petsc_objects_created_ );

	KSPSolve( ksp_, rhs_vec_, petsc_solution_ );

// TODO: Consider removing.
#if 0
#ifndef NDEBUG
	KSPView(ksp_,PETSC_VIEWER_STDOUT_WORLD);
#endif
#endif

}

#endif
