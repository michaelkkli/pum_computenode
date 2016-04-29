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

#ifdef _OPENMP
#include <omp.h>
#endif

//#undef TRACE
//#define TRACE cout << __FILE__ << ":" << __LINE__

#include "box.hh"
#include "box_utils.hh"
#include "cover_structure.hh"
#include "basic_dbinary_tree.hh"
#include "basic_dbinary_tree_utils.hh"
#include "dbinary_tree.hh"
#include "dbinary_tree_utils.hh"
#include "extract_geometry_2d.hh"
#include "geometry.hh"
#include "geometry_utils.hh"
#include "global_approximation_space.hh"
#include "global_basis_function.hh"
#include "indicator_function.hh"
#include "integration_scheme.hh"
#include "integration_scratch.hh"
#include "line_segments.hh"
#include "petsc_solver.hh"
#include "polynomial.hh"
#include "quadrature_rule.hh"
#include "singleimage.hh"
#include "solution.hh"
#include "vtk_output.hh"

#include <petscts.h>

#include <boost/shared_ptr.hpp>

// For char* basename ( char* ).
#include <libgen.h>
#include <sys/resource.h>
#include <sys/stat.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <valarray>

using boost::shared_ptr;

using std::clog;
using std::copy;
using std::cout;
using std::endl;
using std::exp;
using std::max;
using std::ofstream;
using std::slice;
using std::istringstream;
using std::ostringstream;
using std::valarray;

#define TRACE std::clog << __FILE__ << ":" << __LINE__ << "\n"

class constant_fn : public function<2> {
public:
	constant_fn ( double val=1.0 ) : value_(val) {
		this->set_global_function();
	}
private:
	double evaluate ( const double* x ) const {
		return value_;
	}
private:
	double value_;
};

class my_func : public function<2> {
public:
	my_func( double A, double B ) : A_(A), B_(B) {
		this->set_global_function();
	}
private:
	double evaluate( const double* x ) const {
		assert( x );
		return A_ * exp ( -10.0 * ( x[0]*x[0] + x[1]*x[1] ) * B_ );
	}
private:
	double A_, B_;
};

class square_dis_func : public function<2> {
public:
	square_dis_func () {
		this->set_global_function();
	}
private:
	double evaluate( const double* x ) const {
		assert( x );

		if ( -1.0 < x[0] && x[0] < 1.0 &&
		     -1.0 < x[1] && x[1] < 1.0 ) {
			return 3.0;
		} else {
			return 1.0;
		}
	}
};

class circle_dis_func : public function<2> {
public:
	circle_dis_func () {
		this->set_global_function();
	}
private:
	double evaluate( const double* x ) const {
		assert( x );
		double x0_x0 = x[0]*x[0];
		double x1_x1 = x[1]*x[1];
		
		if ( sqrt ( x0_x0 + x1_x1 ) < 2.0 ) {
			return 0.0;
		} else {
			return 1.0;
		}
	}
};

class smooth_bump : public function<2> {
public:
	smooth_bump ( double inner_r = 2.9, double outer_r = 3,
	              double inner_val = 0.5, double outer_val = 5.0 )
		:
		inner_radius_ ( inner_r ), outer_radius_ ( outer_r ),
		inner_value_ ( inner_val ), outer_value_ ( outer_val )
	{
		this->set_global_function();
		
		// Found term outer_value_ incorrectly added here.
		// Became apparent when visualizing function to project.
		bump_scaling_ = (inner_value_-outer_value_) / exp(-1.0);
	}
private:
	double evaluate( const double* x ) const {
		assert( x );

		assert( x );
		double x0_x0 = x[0]*x[0];
		double x1_x1 = x[1]*x[1];
		double r     = sqrt ( x0_x0 + x1_x1 );
		
		assert ( r > 0.0 - 1e-10 );
		
		if (  r < inner_radius_ ) {
			return inner_value_;
		} else if ( r > outer_radius_ ) {
			return outer_value_;
		} else {
			double little_r = (r - inner_radius_)/(outer_radius_-inner_radius_);
			return outer_value_ + bump_scaling_ * exp ( -1.0/(1.0 - little_r*little_r) );
		}
	}
private:
	double inner_radius_;
	double outer_radius_;
	double inner_value_;
	double outer_value_;
	
	double bump_scaling_;
};

enum assembler_matrix_type { mass_matrix, stiffness_matrix, bleach_matrix };
enum assembler_vector_type { projection_vector, domain_integration_vector, bleach_integration_vector };

struct matrix_assembler {
	typedef shared_ptr<matrix_assembler> ptr;
	
	matrix_assembler () : integration_decomposition ( 0 ) { }
	
	void efficient_assemble_matrices ( petsc_solver<2> & solver, bool do_siggia_matrix=false ) {
		assert( cover_st );

		const vector<vector<int> >& index_to_neighbour_indices = cover_st->index_to_neighbour_indices;

		valarray<int> patch_indices_on_computenode;

		solver.get_patch_indices_on_computenode( patch_indices_on_computenode );

		assert( patch_indices_on_computenode.size() > 0 ); // Not strictly necessary.

		int num = patch_indices_on_computenode.size();

		valarray<double> mass_matrix_vals;
		valarray<double> bleach_matrix_vals;
		valarray<double> stiffness_matrix_vals;
		valarray<double> siggia_matrix_vals;
		
		valarray<double>* siggia_ptr=0;

		int first_patch_ind;
		int num_neigh;
		int second_patch_ind;
		
		integration_scratch    scratch;

		// Valgrind helped to debug a problem with blow up of the solution
		// when using the siggia matrix. Multiple threads were accessing
		// the same memory. Definitely some openmp problem. 2008-12-05 Michael LI.
#pragma omp parallel for private(mass_matrix_vals,stiffness_matrix_vals, bleach_matrix_vals,siggia_ptr, siggia_matrix_vals,first_patch_ind,num_neigh,second_patch_ind, scratch)
		for ( int i=0; i<num; ++i ) {
			if ( do_siggia_matrix ) {
				siggia_ptr = &siggia_matrix_vals;
			}
#ifdef _OPENMP
			if ( 0 == i ) {
				cout << "There are "
					<< omp_get_num_threads()
					<< " OpenMP threads and the maximum number has been set to "
					<< omp_get_max_threads()
					<< ".\n";
			}
#endif
			
			first_patch_ind = patch_indices_on_computenode[i];

			// Debugging found that `i' was used here incorrectly in place
			// of first_patch_ind which happened to be correct for serial execution
			// but not parallel execution. 2008-08-20 Mike Li.
			const vector<int>& neighbours = index_to_neighbour_indices[first_patch_ind];

			num_neigh = neighbours.size();

			assert( num_neigh > 0 ); // Not strictly necessary.
			
			for ( int j=0; j<num_neigh; ++j ) {

				second_patch_ind = neighbours[j];

#if 0 // ndef NDEBUG // Keep for future debug.
				std::clog << "Assembling matrix patch block " << first_patch_ind << ", " << second_patch_ind << ".\n";
#endif
				
			efficient_assemble_matrix_block ( first_patch_ind,
			                                  second_patch_ind,
			                                  mass_matrix_vals,
			                                  bleach_matrix_vals,
			                                  stiffness_matrix_vals,
			                                  scratch,
			                                  siggia_ptr );
				
				// PETSc is not thread safe.
				// A critical section is now marked at the point of
				// value insertion into the PETSc matrix.
				solver.set_matrix_entries( first_patch_ind,
				                           second_patch_ind,
				                           mass_matrix_vals,
				                           bleach_matrix_vals,
				                           stiffness_matrix_vals,
				                           siggia_ptr );
			}
		}
	}
	
	void efficient_assemble_matrix_block ( int                      patch_index_row,
	                             int                      patch_index_column,
	                             valarray<double> &       mass_matrix_vals,
	                             valarray<double> &       bleach_matrix_vals,
	                             valarray<double> &       stiffness_matrix_vals,
	                             integration_scratch &    scratch,
	                             valarray<double> *       siggia_matrix_vals ) {
#if 0 // def _OPENMP
		cout << "Patch row " << patch_index_row
			<< ", patch column " << patch_index_column
			<< " calculated by OpenMP thread "
			<< omp_get_thread_num() << "\n";
#endif
		vector<global_basis_function<2> > test;
		vector<global_basis_function<2> > trial;

		global_approx_sp->get_global_basis_by_patch( patch_index_row,    test );
		global_approx_sp->get_global_basis_by_patch( patch_index_column, trial );

		int test_size  = test.size();
		int trial_size = trial.size();

		assert( test_size  > 0 );
		assert( trial_size > 0 );

		int mat_size = test_size * trial_size;

		mass_matrix_vals.resize( mat_size );
		bleach_matrix_vals.resize( mat_size );
		stiffness_matrix_vals.resize( mat_size );
		if ( siggia_matrix_vals ) {
			siggia_matrix_vals->resize ( mat_size );
		}
		

		integrator->integrate_matrix_block ( test,
							trial,
							*bleach_indicator_fun,
							mass_matrix_vals,
							bleach_matrix_vals,
							stiffness_matrix_vals,
							scratch,
							integration_decomposition,
							(equilibrium_fun) ? &(*equilibrium_fun) : 0,
							siggia_matrix_vals );

	}

	void assemble_rhs_vector ( petsc_solver<2>& solver, assembler_vector_type vector_type ) {
		assert( cover_st );
		
		valarray<int> patch_indices_on_computenode;
		solver.get_patch_indices_on_computenode( patch_indices_on_computenode );
		
		assert( patch_indices_on_computenode.size() > 0 ); // Not strictly necessary.
		
		int num = patch_indices_on_computenode.size();
		
		valarray<double> tmp_vals;
		
		integration_scratch    scratch;
		
		
		// Debugging random behaviour when using a member of a member of scratch.
		// Caused by one thread overwriting memory that should be for the other thread.
		// 2008-09-15 Mike Li.
#pragma omp parallel for private(tmp_vals, scratch)
		for ( int i=0; i<num; ++i ) {
			assemble_rhs_vector_block( patch_indices_on_computenode[i], tmp_vals, vector_type, scratch );
			
			// PETSc is not thread safe.
			// A critical section is now marked at the point of
			// value insertion into the PETSc vector.
			solver.set_rhs_vector_entries ( patch_indices_on_computenode[i], tmp_vals );
		}
	}
	
	void assemble_rhs_vector_block ( int                      patch_index,
	                                 valarray<double> &       values,
	                                 assembler_vector_type    vector_type,
	                                 integration_scratch &    scratch ) {
		vector<global_basis_function<2> > test;
		
		global_approx_sp->get_global_basis_by_patch( patch_index, test );
		
		int size = test.size();
		
		values.resize(size);
		
		for ( int i=0; i<size; ++i ) {
			assemble_single_rhs_vector_entry( test[i], values[i], vector_type, scratch );
		}	
	}
	

	void assemble_single_rhs_vector_entry ( const global_basis_function<2>&    test,
	                                        double&                            rhs_vector_entry,
	                                        assembler_vector_type              vector_type,
	                                        integration_scratch &              scratch )
	{
		assert( global_function_to_project );
		if ( vector_type == projection_vector ) {
			rhs_vector_entry = integrator->integrate_product ( *global_function_to_project,
			                                                   test,
			                                                   scratch,
			                                                   integration_decomposition );
		} else if ( vector_type == domain_integration_vector ) {
			rhs_vector_entry = integrator->integrate_global_basis_function ( test,
			                                                                 scratch,
			                                                                 integration_decomposition );
		} else if ( vector_type == bleach_integration_vector ) {
			rhs_vector_entry = integrator->integrate_product ( *bleach_indicator_fun,
				test,
				scratch,
				integration_decomposition );
		} else {
			abort ();
		}
	}

	global_approximation_space<2>::ptr global_approx_sp;
	integration_scheme<2>::ptr         integrator;
	function<2>::ptr                   global_function_to_project;
	function<2>::ptr                   bleach_indicator_fun;
	differentiable_function<2>::ptr    equilibrium_fun;
	cover_structure<2>::ptr            cover_st;
	int                                integration_decomposition;
};

static char help[] = "petsc_solver-3";

int main (int argc, char* argv[]) {
	PetscInitialize( &argc, &argv, PETSC_NULL, help );
	
	// Allow calculation of execution run time.
	PetscLogDouble totalruntime_t0, totalruntime_t1;
	PetscGetTime ( &totalruntime_t0 );
	
#ifdef NDEBUG
	PetscPrintf ( PETSC_COMM_WORLD, "NDEBUG has been defined.\n" );
#else
	PetscPrintf ( PETSC_COMM_WORLD, "Preprocessor symbol NDEBUG not defined.\n" );
#endif
	
#ifdef _OPENMP
	PetscPrintf ( PETSC_COMM_WORLD, "_OPENMP has been defined.\n" );
#else
	PetscPrintf ( PETSC_COMM_WORLD, "Preprocessor symbol _OPENMP not defined.\n" );
#endif
	
	PetscMPIInt mpi_size, mpi_rank;
	MPI_Comm_size ( PETSC_COMM_WORLD, &mpi_size );
	MPI_Comm_rank ( PETSC_COMM_WORLD, &mpi_rank );
		
	PetscLogDouble log_t0, log_t1;
	
	// Used to divide work between MPI processes.
	int round_robin_counter = 0;

	PetscTruth usage_requested = PETSC_FALSE;
	PetscOptionsHasName ( PETSC_NULL, "-help", &usage_requested );
	if ( usage_requested ) {
		cout << "usage: " << basename(argv[0]) << " [options with default values follow]\n";
	}
	
	PetscTruth    min_iter_info = PETSC_TRUE;
	PetscOptionsGetTruth ( PETSC_NULL, "-min_iter_info", &min_iter_info, PETSC_NULL );

	PetscReal initial_value1 = 1.0;
	PetscOptionsGetReal ( PETSC_NULL, "-initial_value1", &initial_value1, PETSC_NULL );
	PetscPrintf ( PETSC_COMM_WORLD, "\t-initial_value1 is %f.\n", initial_value1 );

	PetscReal domain_size = 10.0;
	PetscOptionsGetReal ( PETSC_NULL, "-domain_size", &domain_size, PETSC_NULL );
	if ( domain_size < 0.0 ) {
		domain_size = 10.0;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-domain_size is %f.\n", domain_size );
	
	PetscReal global_bleach1 = 0.0;
	PetscOptionsGetReal ( PETSC_NULL, "-global_bleach1", &global_bleach1, PETSC_NULL );
	PetscPrintf ( PETSC_COMM_WORLD, "\t-global_bleach1 is %f.\n", global_bleach1 );
	
	PetscReal bleach_x_centre1 = 256.0;
	PetscOptionsGetReal ( PETSC_NULL, "-bleach_x_centre1", &bleach_x_centre1, PETSC_NULL );
	PetscPrintf ( PETSC_COMM_WORLD, "\t-bleach_x_centre1 is %f.\n", bleach_x_centre1 );
	
	PetscReal bleach_y_centre1 = 256.0;
	PetscOptionsGetReal ( PETSC_NULL, "-bleach_y_centre1", &bleach_y_centre1, PETSC_NULL );
	PetscPrintf ( PETSC_COMM_WORLD, "\t-bleach_y_centre1 is %f.\n", bleach_y_centre1 );
	
	PetscReal bleach_radius1 = 20.;
	PetscOptionsGetReal ( PETSC_NULL, "-bleach_radius1", &bleach_radius1, PETSC_NULL );
	if ( bleach_radius1 < 0.0 ) {
		bleach_radius1 = 20.;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-bleach_radius1 is %f.\n", bleach_radius1 );
	
	PetscReal bleach_param1 = 100.0;
	PetscOptionsGetReal ( PETSC_NULL, "-bleach_param1", &bleach_param1, PETSC_NULL );
	if ( bleach_param1 < 0.0 ) {
		bleach_param1 = 100.0;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-bleach_param1 is %f.\n", bleach_param1 );
	
	PetscInt begin_bleach1 = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-begin_bleach1", &begin_bleach1, PETSC_NULL );
	PetscPrintf ( PETSC_COMM_WORLD, "\t-begin_bleach1 is %i.\n", begin_bleach1 );
	
	PetscInt end_bleach1 = 51;
	PetscOptionsGetInt ( PETSC_NULL, "-end_bleach1", &end_bleach1, PETSC_NULL );
	PetscPrintf ( PETSC_COMM_WORLD, "\t-end_bleach1 is %i.\n", end_bleach1 );
	
	
	PetscInt vtk_out = 3;
	PetscOptionsGetInt ( PETSC_NULL, "-vtk_out", &vtk_out, PETSC_NULL );
	if ( vtk_out < 0 || vtk_out > 3 ) {
		PetscPrintf ( PETSC_COMM_WORLD, "Invalid value for vtk_out.\n" );
		vtk_out = 3;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-vtk_out is %i.\n", vtk_out );
	
	PetscInt iterations = 72;
	PetscOptionsGetInt(PETSC_NULL,"-iterations", &iterations, PETSC_NULL);
	if ( iterations<0 ) {
		iterations = 72;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-iterations is %i.\n", iterations );
	
	PetscReal diffusion1 = 33.0;
	PetscOptionsGetReal ( PETSC_NULL, "-diffusion1", &diffusion1, PETSC_NULL );
	if ( diffusion1 < 0.0 ) {
		diffusion1 = 33.0;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-diffusion1 is %f.\n", diffusion1 );
	
	PetscReal time_step_dt = 0.754;
	PetscOptionsGetReal ( PETSC_NULL, "-dt", &time_step_dt, PETSC_NULL );
	if ( !(time_step_dt > 0.0) ) {
		time_step_dt = 0.754;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-dt is %f.\n", time_step_dt );
	
	PetscReal parent_box_expansion_factor = 1.2;
	PetscOptionsGetReal(PETSC_NULL,"-parent_box_expansion_factor", &parent_box_expansion_factor, PETSC_NULL);
	if ( parent_box_expansion_factor<=0.0 ) {
		parent_box_expansion_factor = 1.2;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-parent_box_expansion_factor is %f.\n", parent_box_expansion_factor );
	
	PetscInt dtree_initial_refinement = 4;
	PetscOptionsGetInt(PETSC_NULL,"-dtree_initial_refinement", &dtree_initial_refinement, PETSC_NULL);
	if ( dtree_initial_refinement<0 ) {
		dtree_initial_refinement = 4;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-dtree_initial_refinement is %i.\n", dtree_initial_refinement );
	
	PetscInt viz_refinement = dtree_initial_refinement;
	PetscOptionsGetInt(PETSC_NULL,"-viz_refinement", &viz_refinement, PETSC_NULL);
	if ( viz_refinement < 0 ) {
		PetscPrintf ( PETSC_COMM_WORLD, "Invalid value for viz_refinement.\n" );
		viz_refinement = dtree_initial_refinement;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-viz_refinement is %i.\n", viz_refinement );
	
	PetscInt dtree_max_refinement = max ( dtree_initial_refinement, viz_refinement ) + 2;
	PetscOptionsGetInt(PETSC_NULL,"-dtree_max_refinement", &dtree_max_refinement, PETSC_NULL);
	if ( dtree_max_refinement < dtree_initial_refinement
	     || dtree_max_refinement < viz_refinement ) {
		     PetscPrintf ( PETSC_COMM_WORLD, "Invalid value for dtree_max_refinement.\n" );
		     dtree_max_refinement = max ( dtree_initial_refinement, viz_refinement ) + 2;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-dtree_max_refinement is %i.\n", dtree_max_refinement );
	
	PetscReal dtree_cover_factor = 1.3;
	PetscOptionsGetReal(PETSC_NULL,"-dtree_cover_factor", &dtree_cover_factor, PETSC_NULL);
	if ( dtree_cover_factor<=0.0 ) {
		dtree_cover_factor = 1.3;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-dtree_cover_factor is %f.\n", dtree_cover_factor );
	
	PetscInt initial_function_choice = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-initial_function_choice", &initial_function_choice, PETSC_NULL );
	if ( initial_function_choice < 1 || 5 < initial_function_choice ) {
		initial_function_choice = 1;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-initial_function_choice is %i.\n", initial_function_choice );
	
	PetscInt assembly_quadrature_rule = 4;
	PetscOptionsGetInt ( PETSC_NULL, "-assembly_quadrature_rule", &assembly_quadrature_rule, PETSC_NULL );
	if ( assembly_quadrature_rule < 1 || 5 < assembly_quadrature_rule ) {
		assembly_quadrature_rule = 4;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-assembly_quadrature_rule is %i.\n", assembly_quadrature_rule );
	
	PetscInt local_basis_order = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-local_basis_order", &local_basis_order, PETSC_NULL );
	if ( local_basis_order < 1 || 2 < local_basis_order ) {
		local_basis_order = 1;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-local_basis_order is %i.\n", local_basis_order );
	
	PetscInt integration_decomposition = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-integration_decomposition", &integration_decomposition, PETSC_NULL );
	if ( integration_decomposition < 0 || 1 < integration_decomposition ) {
		integration_decomposition = 1;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-integration_decomposition is %i.\n", integration_decomposition );
	
	PetscInt boundary_detail = 100;
	PetscOptionsGetInt ( PETSC_NULL, "-boundary_detail", &boundary_detail, PETSC_NULL );
	PetscPrintf ( PETSC_COMM_WORLD, "\t-boundary_detail is %i.\n", boundary_detail );
	
	string    dom_string;
	{
		char    tmp_dom[256];
		PetscTruth    given = PETSC_FALSE;
		PetscOptionsGetString ( PETSC_NULL, "-dom", tmp_dom, 255, &given );
		if ( given ) {
			dom_string.assign ( tmp_dom );
		}
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-dom is \"%s\".\n", dom_string.c_str() );
	
	string    image1_string;
	{
		char    tmp_dom[256];
		PetscTruth    given = PETSC_FALSE;
		PetscOptionsGetString ( PETSC_NULL, "-image1", tmp_dom, 255, &given );
		if ( given ) {
			image1_string.assign ( tmp_dom );
		}
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-image1 is \"%s\".\n", image1_string.c_str() );
	
	PetscReal rescale_image1 = 0.5;
	PetscOptionsGetReal(PETSC_NULL,"-rescale_image1", &rescale_image1, PETSC_NULL);
	PetscPrintf ( PETSC_COMM_WORLD, "\t-rescale_image1 is %f.\n", rescale_image1 );
	
	string    manifest1_filename;
	{
		char    tmp_arr[256];
		PetscTruth    given = PETSC_FALSE;
		PetscOptionsGetString ( PETSC_NULL, "-manifest1", tmp_arr, 255, &given );
		if ( given ) {
			manifest1_filename.assign ( tmp_arr );
		}
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-manifest1 is \"%s\".\n", manifest1_filename.c_str() );

	string    manifest1_dir;
	{
		char    tmp_arr[256];
		PetscTruth    given = PETSC_FALSE;
		PetscOptionsGetString ( PETSC_NULL, "-manifest1_dir", tmp_arr, 255, &given );
		if ( given ) {
			manifest1_dir.assign ( tmp_arr );
			if ( manifest1_dir[ manifest1_dir.size()-1 ] != '/' ) {
				manifest1_dir.append ( "/" );
			}
		}
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-manifest1_dir is \"%s\".\n", manifest1_dir.c_str() );

	enum simulation_mode { frap_flip, siggia_clflip, L2_projection };

	PetscInt    simulation_mode_num = 1;
	simulation_mode    sim_mode = frap_flip;
	
	PetscOptionsGetInt ( PETSC_NULL, "-simulation_mode", &simulation_mode_num, PETSC_NULL );
	if ( simulation_mode_num < 0 || simulation_mode_num > 2 ) {
		abort();
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-simulation_mode is %i.\n", simulation_mode_num );

	if ( simulation_mode_num == 0 ) {
		PetscPrintf ( PETSC_COMM_WORLD, "Simulation mode: time-stepping diffusion with optional bleaching.\n" );
		sim_mode = frap_flip;
	}
	if ( simulation_mode_num == 1 ) {
		PetscPrintf ( PETSC_COMM_WORLD, "Simulation mode: Siggia-modified time-stepping diffusion with optional bleaching.\n" );
		sim_mode = siggia_clflip;
	}
	if ( simulation_mode_num == 2 ) {
		PetscPrintf ( PETSC_COMM_WORLD, "Simulation mode: .\n" );
		sim_mode = L2_projection;
	}


	if ( usage_requested ) {
		return 0;
	}
	
	MPI_Barrier ( PETSC_COMM_WORLD );
	
	vector<string>    manifest1_content;
	// Load manifest data if the relevant simulation mode is requested.
	if ( sim_mode == L2_projection ) {
		if ( manifest1_filename == "" ) {
			cerr << "A text file listing tiff images must be provided (please give an argument to the option -manifest1).\n";
			abort();
		}
		ifstream    file ( manifest1_filename.c_str() );
		if ( !file ) {
			cerr << "There was a problem opening the manifest1 file \""
			     << manifest1_filename << "\".\n";
			abort();
		}
		string      ln;
		while ( getline ( file, ln ) ) {
			manifest1_content.push_back ( ln );
		}
		int size = manifest1_content.size();
		struct stat    ignore;
		ostringstream    ss;
		int ret;
		for ( int i=0; i<size; ++i ) {
			ss << manifest1_dir << manifest1_content[i];
			ret = stat ( (ss.str()).c_str(), &ignore );
			//cout << "(ss.str()).c_str() is " << (ss.str()).c_str() << "\n";
			if ( ret != 0 ) {
				cerr << "The file \"" << manifest1_content[i]
				     << "\" named in \"" << manifest1_filename
				     << "\" does not exist.\n";
				abort();
			}
			ss.str ( "" );
		}
		iterations = manifest1_content.size();
		cout << "Iterations has been set to " << iterations << ".\n";
	}
	
#if 1
	const int num_doubles = 12;
	double coords[] = { 0.0, -10.0, 1.2, -7.0, 5.4, 2.0, 3.0, 7.8, -3.0, 9.5, -0.7, 0.5 };
#else
	const int num_doubles = 8;
	double coords[num_doubles] = { 0.0, 0.0,
			0.0, 1.0,
			1.0, 1.0,
			1.0, 0.0 };
#endif

#if 0
	shared_ptr<geometry<2> > geom( new geometry<2> );
#endif

	valarray<double> vcoords( num_doubles );
	
	if ( dom_string == "" ) {
	
		for ( int i=0; i<num_doubles; ++i ) {
			vcoords[i] = coords[i];
		}
		
		generate_circular_boundary( boundary_detail, domain_size, 0, vcoords );
	} else {
		simple_off_input_geometry_2d ( dom_string, vcoords );
	}
	
	line_segments    domain_segs;
	{
		int tmp_size = vcoords.size();
		domain_segs.segments.resize ( tmp_size +2 );
		copy ( &vcoords[0], &(vcoords[0])+tmp_size, &(domain_segs.segments[0]) );
		
		// Debugging found the last point was being set incorrectly.
		// Using -2 and -1 to adjust indices. 2008-12-01 Michael LI.
		domain_segs.segments[ tmp_size ]   = vcoords[0];
		domain_segs.segments[ tmp_size+1 ] = vcoords[1];
		
		int num_segs = tmp_size/2;
		
		if ( round_robin_counter++ % mpi_size == mpi_rank )
		{
			PetscPrintf( PETSC_COMM_SELF, "Domain outline output.\n" );
			ofstream    file ( "petsc_solver-3_stepper_labelled_bdry.gnuplot" );
			gnuplot_output<2> ( domain_segs, file, (num_segs<1000) );
		}
	}
	
	box<2>    parent_box;
	get_bounding_box ( domain_segs, parent_box );
	parent_box.scale ( parent_box_expansion_factor );
	
	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		PetscPrintf( PETSC_COMM_SELF, "Parent box file output.\n" );
		ofstream file_15( "petsc_solver-3_parent_box" );
		gp_draw_single( parent_box, file_15 );
		file_15.close();
	}

	// Preparing to eliminate dbinary_tree from global_approx_space.
	// 2009-02-11
	basic_dbinary_tree<2>::ptr    sim_basic_dtree ( new basic_dbinary_tree<2> );
	sim_basic_dtree->initialize ( parent_box, dtree_max_refinement );
	
	refinement_structure<2>::ptr    ref_struct ( new refinement_structure<2> );
	ref_struct->dtree = sim_basic_dtree;
	ref_struct->levels["0"] = dtree_initial_refinement;
	ref_struct->cover_factor = dtree_cover_factor;

	
	// No longer needed here.
	// vtk_output    bdry_vtkout;
	
	cover_structure<2>::ptr    cov_struct ( new cover_structure<2> );
	cov_struct->ref_structure = ref_struct;
	
	{

		set<string>   bdry_keys;
		
		{
			vector<string>    vec_bdry_keys;
			vector<int>       entry_seg;
			
			get_intersect_keys_entry_indices ( *ref_struct, domain_segs, vec_bdry_keys, entry_seg );
			
			copy ( vec_bdry_keys.begin(), vec_bdry_keys.end(),
			inserter( bdry_keys, bdry_keys.begin() ) );
		}
#if 0 // No longer needed.	
		make_boundary_vtk_output ( domain_segs,
	                                   *ref_struct,
	                                   bdry_vtkout,
	                                   outside, &bdry_keys );
#endif
	
		map<string, set<string> > &    nbrs = cov_struct->patch_neighbour_keys;
	
		vector<string> tree_all_keys;
		
		PetscPrintf ( PETSC_COMM_WORLD, "Beginning generate_keys_and_neighbours.\n" ); PetscGetTime ( &log_t0 );
		
		//generate_keys_and_neighbours< 2 >( *ref_struct, tree_all_keys, nbrs );
		
		// Performance testing.
		{
			PetscPrintf( PETSC_COMM_WORLD, "Testing code is being executed.\n" );
			map<string,pair<int,int> >    unused;
			
			generate_keys_and_neighbours ( *ref_struct,
                                    tree_all_keys,
                                    nbrs,
                                    &domain_segs,
                                    &unused );

		}
		PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, "generate_keys_and_neighbours took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		set<string>    inside_poss_bdry_keys;
		get_inside_keys ( tree_all_keys, *sim_basic_dtree, domain_segs, inside_poss_bdry_keys );
		
		PetscPrintf ( PETSC_COMM_WORLD, "Beginning set_union.\n" ); PetscGetTime ( &log_t0 );
		set_union ( bdry_keys.begin(), bdry_keys.end(),
			inside_poss_bdry_keys.begin(), inside_poss_bdry_keys.end(),
			inserter( cov_struct->patch_keys, cov_struct->patch_keys.begin() )
		);
		PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, "set_union took %f minutes.\n", (log_t1-log_t0)/60.0 );
	}

	PetscPrintf ( PETSC_COMM_WORLD, "Beginning generate_index_information.\n" ); PetscGetTime ( &log_t0 );
	generate_index_information ( *cov_struct );
	PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, "generate_index_information took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscGetTime ( &log_t0 );
	generate_patch_information ( *cov_struct );
	PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, "generate_patch_information took %f minutes.\n", (log_t1-log_t0)/60.0 );

#if 0
	PetscPrintf ( PETSC_COMM_WORLD, "Deprecated geom use for boundary.\n" );
	// Set up the geometry.
	geom->create_boundary( "chloroplast envelope", vcoords );
	geom->create_region( "chloroplast", "+chloroplast envelope" );
	
	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		ofstream file_10( "petsc_solver-3_chloroplast" );
		geom->gp_draw( "chloroplast", file_10 );
		file_10.close();
	}

	geom->update();
#endif
	
	
#if 0	
	// Set up the d-binary tree used to form the
	// partition of unity method cover.
	shared_ptr<dbinary_tree<2> > dtree( new dbinary_tree<2> );
	dtree->set_max_level ( dtree_max_refinement );
	dtree->set_geometry( geom, parent_box_expansion_factor );
#endif
	

	
#if 0
	PetscGetTime ( &log_t0 );
	dtree->set_initial_refinement( dtree_initial_refinement ); // 1D 2^n, 2D 4^n, 3D 8^n so refinement can be expensive.
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "dbinary_tree initial refinement took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	dtree->set_cover_factor( dtree_cover_factor ); // Possibly use "length extension" to be descriptive.
#endif
	
	// Create covers for the regions and boundaries.
	// TODO: allow the option not to cover unwanted areas.
	
	// clog << "Begin update of dbinary_tree.\n";

#if 0
	PetscGetTime ( &log_t0 );
	dtree->update();
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "dbinary_tree update took %f minutes.\n", (log_t1-log_t0)/60.0 );
#endif

#if 0 // Keep. Better way of doing this now.
	clog << "Begin gnuplot file output.\n";

	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		ofstream file( "petsc_solver-3_boundary_cover" );
		dtree->gp_draw_boundary_cover( "chloroplast envelope", file, 0.3 );
		file.close();
	}

	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		ofstream file( "petsc_solver-3_region_boundary_cover" );
		dtree->gp_draw_region_boundary_cover( "chloroplast", file, 0.1 );
		file.close();
	}

	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		ofstream file( "petsc_solver-3_region_interior_cover" );
		dtree->gp_draw_region_interior_cover( "chloroplast", file, 0.2 );
		file.close();
	}
#endif
	
	// dtree->dump_cout();
	
	global_approximation_space<2>::ptr global_approx_space;

	// create_petsc_objects must be adjusted to allow for the correct number of expected
	// non-zero entries per row.
	switch ( local_basis_order ) {
	case 1:
		global_approx_space.reset ( new global_approximation_space<2> ( monomial1 ) );
		break;
	case 2:
		global_approx_space.reset ( new global_approximation_space<2> ( monomial2 ) );
		break;
	default:
		global_approx_space.reset ( new global_approximation_space<2> );
		break;
	}
	
	matrix_assembler::ptr              mat_assembler ( new matrix_assembler );
	
#if 0
	assert ( dtree );
	// Needed for point intersection search.
	global_approx_space->give_dbinary_tree ( dtree );
	global_approx_space->set_support ( "chloroplast" );
#endif
	
	{
		cover_structure<2>::ptr cs;
		
#if 0
		assert( dtree );
		PetscPrintf ( PETSC_COMM_WORLD, "Begin creation of cover_structure.\n" );
		PetscGetTime ( &log_t0 );
		dtree->get_cover_structure( "chloroplast", cs );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "get_cover_structure took %f minutes.\n", (log_t1-log_t0)/60.0 );
		cs->ref_structure = ref_struct;
#endif
		

	
		// Replacing dbinary_tree
		cs = cov_struct;
	
		mat_assembler->cover_st = cs;
		
		global_approx_space->give_cover_structure( cs );
	}
	
	mat_assembler->global_approx_sp = global_approx_space;
	
	switch ( assembly_quadrature_rule ) {
	case 1:
		(mat_assembler->integrator).reset( new integration_scheme<2>( gauss1 ) );
		break;
	case 2:
		(mat_assembler->integrator).reset( new integration_scheme<2>( gauss2 ) );
		break;
	case 3:
		(mat_assembler->integrator).reset( new integration_scheme<2>( gauss3 ) );
		break;
	case 4:
		(mat_assembler->integrator).reset( new integration_scheme<2>( gauss4 ) );
		break;
	case 5:
		(mat_assembler->integrator).reset( new integration_scheme<2>( gauss5 ) );
		break;
	default:
		(mat_assembler->integrator).reset( new integration_scheme<2>( gauss4 ) );
	}
	
	switch ( initial_function_choice ) {
	case 1:
		mat_assembler->global_function_to_project.reset( new constant_fn ( initial_value1 ) );
		break;
	case 2:
		mat_assembler->global_function_to_project.reset( new my_func ( domain_size, 1.0/(domain_size*domain_size) ) );
		break;
	case 3:
		mat_assembler->global_function_to_project.reset( new square_dis_func );
		break;
	case 4:
		mat_assembler->global_function_to_project.reset( new circle_dis_func );
		break;
	case 5:
		mat_assembler->global_function_to_project.reset( new smooth_bump );
		break;
	default:
		mat_assembler->global_function_to_project.reset( new circle_dis_func );
		break;
	}


	singleimage::ptr    image_initial;
	
	if ( image1_string != "" ) {
		image_initial.reset ( new singleimage );
		image_initial->load ( image1_string, 1 );
		
		image_initial->rescale_percentage ();
		
		mat_assembler->global_function_to_project = image_initial;
		mat_assembler->equilibrium_fun = image_initial;
	}
	
	
	// Prepare bleach indicator function.
	valarray<double>    bleach_region_coords;
	
	double bleach_centre[2];
	
	// (dtree->access_parent_box()).get_centre_point ( bleach_centre );
	bleach_centre[0] = bleach_x_centre1;
	bleach_centre[1] = bleach_y_centre1;
	
	generate_circular_boundary( 30, bleach_radius1, bleach_centre, bleach_region_coords );
	shared_ptr<geometry<2> >    bleach_geom ( new geometry<2> );
	bleach_geom->create_boundary ( "bleach boundary", bleach_region_coords );
	bleach_geom->create_region ( "bleach region", "+bleach boundary" );
	bleach_geom->update ();
	
#if 0
	{
		ofstream file ( "petsc_solver-3_geom_bleach_region" );
		bleach_geom->gp_draw( "bleach region", file );
	}
#endif
	
	indicator_function<2>::ptr    bleach_ind_fun ( new indicator_function<2> ( "bleach region", bleach_geom ) );
	
	// Set bleach indicator function.
	mat_assembler->bleach_indicator_fun = bleach_ind_fun;
	
	mat_assembler->integration_decomposition = integration_decomposition;
		
	valarray<int> patch_to_num_dof;
	global_approx_space->get_patch_to_num_dof( patch_to_num_dof );
	
	petsc_solver<2>::ptr pet_solver ( new petsc_solver<2> );

	// The number of non-zero entries per row of the global matrix must
	// be set correctly or else there will be a sudden jump in computation
	// time.
	
	// cout << "Begin creation of PETSc objects.\n";
	PetscGetTime ( &log_t0 );
	if ( local_basis_order == 1 ) {
		pet_solver->create_petsc_objects ( patch_to_num_dof, 27 );
	} else if ( local_basis_order == 2 ) {
		pet_solver->create_petsc_objects ( patch_to_num_dof, 54 );
	}
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "create_petsc_objects took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	solution<2>::ptr sol ( new solution<2> );
	
	PetscGetTime ( &log_t0 );
	sol->set_global_approximation_space ( global_approx_space, patch_to_num_dof.sum() );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "set_global_approximation_space took %f minutes.\n", (log_t1-log_t0)/60.0 );

#if 0
	{
		cout << "Mike early abort.\n";
		PetscErrorCode ierr = PetscFinalize();
		exit(0);
	}
#endif
	
	PetscPrintf ( PETSC_COMM_WORLD, "Beginning efficient_assemble_matrices. " );
	PetscGetTime ( &log_t0 );
	mat_assembler->efficient_assemble_matrices ( *pet_solver, sim_mode == siggia_clflip ); // Assemble Siggia matrix for simulation mode 1.
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "efficient_assemble_matrices took %f minutes.\n", (log_t1-log_t0)/60.0 );

	PetscGetTime ( &log_t0 );
	pet_solver->quick_assemble_matrices ();
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "quick_assemble_matrices took %f minutes.\n",  (log_t1-log_t0)/60.0 );

	PetscGetTime ( &log_t0 );
	pet_solver->load_mass_matrix ();
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "load_mass_matrix took %f minutes.\n",  (log_t1-log_t0)/60.0 );

	//	pet_solver->view_bleach_matrix ();
	
#if 0
	pet_solver->view_mass_matrix ();
	cout << "Hello";
	pet_solver->view_petsc_matrix ();
#endif
	
	// Vector to allow integration over the whole domain with a simple dot product.
	PetscGetTime ( &log_t0 );
	mat_assembler->assemble_rhs_vector ( *pet_solver, domain_integration_vector );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "assemble_rhs_vector took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	pet_solver->rhs_vector_assembly_begin ();
	pet_solver->rhs_vector_assembly_end ();
	pet_solver->save_domain_integration_vec ();
	
	// Vector to allow integration over the bleach region with a simple dot product.
	mat_assembler->assemble_rhs_vector ( *pet_solver, bleach_integration_vector );
	pet_solver->rhs_vector_assembly_begin ();
	pet_solver->rhs_vector_assembly_end ();
	pet_solver->save_bleach_integration_vec ();
	
	// cout << "Begin computenode local assembly of right-hand side vector.\n";
	PetscGetTime ( &log_t0 );
	mat_assembler->assemble_rhs_vector ( *pet_solver, projection_vector );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Local assembly of right-hand side vector took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	// cout << "Begin global assembly of right-hand side vector.\n";
	PetscGetTime ( &log_t0 );
	pet_solver->rhs_vector_assembly_begin ();
	pet_solver->rhs_vector_assembly_end ();
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Assembly of global right-hand side vector took %f minutes.\n", (log_t1-log_t0)/60.0 );


	
	// cout << "Begin global solve.\n";
	PetscGetTime ( &log_t0 );
	pet_solver->solve( *sol );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Solve and distribution of solution took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	if ( global_bleach1 > 0.0 ) {
		PetscGetTime ( &log_t0 );
		pet_solver->load_mass_matrix ( 1.0 + time_step_dt*global_bleach1 );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "load_mass_matrix (scaled) took %f minutes.\n",  (log_t1-log_t0)/60.0 );
	}
	
	// pet_solver->view_solution ( "petsc_solver-3_stepper_check_large_coeff" );
	

	if ( (sim_mode == frap_flip) || (sim_mode == siggia_clflip) ) {
		pet_solver->add_scaled_stiffness_matrix ( time_step_dt*diffusion1 );
	}
	
	if ( sim_mode == siggia_clflip ) {
		pet_solver->add_scaled_siggia_matrix ( -time_step_dt*diffusion1 );
	}

	
	vector<string> domain_keys;
	{
		vector<string> domain_all_keys;
	
		generate_keys< 2 >( viz_refinement, domain_all_keys );	
		
#if 0
		assert ( domain_all_keys.size() > 0 );
		PetscGetTime ( &log_t0 );
		dtree->get_intersect_region ( "chloroplast", domain_all_keys, domain_keys );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "get_intersect_region took %f minutes.\n", (log_t1-log_t0)/60.0 );
#endif
	}
	
	basic_dbinary_tree<2>::ptr basic_dtree ( new basic_dbinary_tree<2> );
	basic_dtree->initialize ( parent_box, viz_refinement );
	
	
	//vector<string>    viz_ref_keys;
	//generate_keys ( viz_ref_struct, viz_ref_keys );
	
	vtk_output::ptr vtk_st ( new vtk_output );
	
	basic_dtree->get_grid_points_3D ( vtk_st->points_3D );

#if 1 // Testing for new boundary handling.
	{
		refinement_structure<2>   viz_ref_struct;
		viz_ref_struct.dtree = basic_dtree;
		viz_ref_struct.levels["0"] = viz_refinement;
	
		vtk_output    bdry_vtkout;
		
		// This used to be set<string>. Hope it doesn't mess up the
		// union later. 2009-03-06 ML.
		vector<string>   viz_bdry_keys;
		
		// VTK seems to give a clockwise winding boundary.
		make_boundary_vtk_output ( domain_segs,
		                           viz_ref_struct,
		                           bdry_vtkout,
		                           outside, &viz_bdry_keys );

		vector<string>    tree_all_keys;
		generate_keys ( viz_ref_struct, tree_all_keys );

		set<string>    inside_poss_bdry_keys;
		
		get_inside_keys ( tree_all_keys, *basic_dtree, domain_segs, inside_poss_bdry_keys );
		
		// This replaces the dbinary_tree generation of keys.
		domain_keys.clear();
		
		PetscPrintf ( PETSC_COMM_WORLD, "This set union must work with unsorted vector.\n" );
		
		set<string>    sorted_viz_bdry_keys;
		
		copy ( viz_bdry_keys.begin(),
		       viz_bdry_keys.end(),
		       inserter ( sorted_viz_bdry_keys, sorted_viz_bdry_keys.end() )
		);
		
		set_union ( sorted_viz_bdry_keys.begin(), sorted_viz_bdry_keys.end(),
			inside_poss_bdry_keys.begin(), inside_poss_bdry_keys.end(),
			back_inserter( domain_keys )
		);

		bdry_vtkout.scalars.resize ( bdry_vtkout.points_3D.size() / 3, 0.0 );

		if ( vtk_out > 0 && (round_robin_counter++ % mpi_size == mpi_rank) ) {
			PetscPrintf ( PETSC_COMM_WORLD, "Boundary test output file.\n" );
			vtk_simple_output( bdry_vtkout,
		        	         "petsc_solver-3_bdry_test.vtp" );
//			vtk_append_raw ( bdry_vtkout,
//		        	         "petsc_solver-3_bdry_test.vtp" );
		}
	}
#endif
	
	int num_pts = (vtk_st->points_3D).size() / 3;
	
	valarray<double> points_two_D ( 2*num_pts );
	const valarray<double>& points_three_D = vtk_st->points_3D;
	
	points_two_D[ slice(0, num_pts, 2) ] = points_three_D[ slice(0, num_pts, 3) ];
	points_two_D[ slice(1, num_pts, 2) ] = points_three_D[ slice(1, num_pts, 3) ];
	
	basic_dtree->get_box_grid_connectivity_offsets ( domain_keys, vtk_st->connectivity, vtk_st->offsets );
#if 1
	if ( (domain_keys.size()<1000) && (vtk_out > 0) && (round_robin_counter++ % mpi_size == mpi_rank) ) {
		PetscPrintf ( PETSC_COMM_SELF, "Labelled cover output file.\n" );
		ofstream    file ( "petsc_solver-3_stepper_labelled_cover.gnuplot" );
		gnuplot_output<2> ( domain_keys, *basic_dtree, file, true );
	}
#endif

	if ( vtk_out > 0 && (round_robin_counter++ % mpi_size == mpi_rank) )
	{
		// PetscGetTime ( &log_t0 );
		mat_assembler->global_function_to_project->global_evaluate ( box<2>(), points_two_D, vtk_st->scalars );
		// PetscGetTime ( &log_t1 );
		// PetscPrintf ( PETSC_COMM_WORLD, "Evaluate of function to project took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		assert ( (vtk_st->scalars).size() == num_pts );
		
		
		vtk_st->points_3D[ slice(2, vtk_st->scalars.size(), 3) ] = vtk_st->scalars;
		
		if ( vtk_out == 3 ) {
			vtk_append_raw ( *vtk_st, "petsc_solver-3_function_to_project.vtp" );
		} else if ( vtk_out ==2 ) {
			vtk_simple_output( *vtk_st, "petsc_solver-3_function_to_project.vtp", true );
		} else {
			vtk_simple_output( *vtk_st, "petsc_solver-3_function_to_project.vtp" );
		}
	}
	
	solution_evaluation_structure    sol_eval_struct;
	
	PetscPrintf ( PETSC_COMM_WORLD, "Beginning generation of solution_evaluation_structure.\n" );
	PetscGetTime ( &log_t0 );
	global_approx_space->get_solution_evaluation_structure ( points_two_D, sol_eval_struct );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "get_solution_evaluation_structure took %f minutes.\n", (log_t1-log_t0)/60.0 );

	
	// Visualization of solution.
	if ( !(sim_mode==L2_projection) && (vtk_out>0) && (round_robin_counter++ % mpi_size == mpi_rank) )
	{
		// Solution evaluation structure was originally not being used.
		// This showed itself when during debugging a regression when
		// trying to eliminate usage of dbinary_tree. The first frame was sometimes
		// seen to be correct with all subsequent frames incorrect. It would
		// also happen that the first frame was incorrect with subsequent frames
		// being correct and this was due to use of sol_eval_struct elsewhere but not here.
		PetscGetTime ( &log_t0 );
		sol->global_evaluate ( points_two_D, sol_eval_struct, vtk_st->scalars );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "Solution evaluation took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		PetscGetTime ( &log_t0 );

		vtk_st->points_3D[ slice(2, vtk_st->scalars.size(), 3) ] = vtk_st->scalars;
		
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "Preparing vtk information took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		if ( vtk_out == 3 ) {
			vtk_append_raw ( *vtk_st, "petsc_solver-3_stepper_0.vtp" );
		} else if ( vtk_out ==2 ) {
			vtk_simple_output( *vtk_st, "petsc_solver-3_stepper_0.vtp", true );
		} else {
			vtk_simple_output( *vtk_st, "petsc_solver-3_stepper_0.vtp" );
		}
	}
		
	solution<2>::vec_ptr   solution_ring ( mpi_size );
	for ( int i=0; i<mpi_size; ++i ) {
		solution_ring[i].reset ( new solution<2> );
		solution_ring[i]->set_global_approximation_space ( global_approx_space, patch_to_num_dof.sum() );
	}
	
	solution<2>::vec_ptr   viz_solution_ring ( mpi_size );
	for ( int i=0; i<mpi_size; ++i ) {
		viz_solution_ring[i].reset ( new solution<2> );
		viz_solution_ring[i]->set_global_approximation_space ( global_approx_space, patch_to_num_dof.sum() );
	}
	
	{
		int num_to_do;
		if ( L2_projection ) {
			num_to_do = iterations;
		} else {
			num_to_do = iterations + 1;
		}
		pvd_output pvd;
		pvd.head = "petsc_solver-3_stepper_";
		pvd.middle.resize ( num_to_do );
		ostringstream oss;
		for ( int i=0; i<num_to_do; ++i ) {
			oss.str ( "" );
			oss << i;
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "petsc_solver-3_stepper.pvd" );
		pvd_simple_output ( pvd, file );
	}
	
	PetscLogDouble iterations_t0, iterations_t1;
	PetscGetTime ( &iterations_t0 );
	
#if 0 // Moving further up to use for first frame.
	solution_evaluation_structure    sol_eval_struct;
	
	PetscPrintf ( PETSC_COMM_WORLD, "Beginning generation of solution_evaluation_structure.\n" );
	PetscGetTime ( &log_t0 );
	global_approx_space->get_solution_evaluation_structure ( points_two_D, sol_eval_struct );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "get_solution_evaluation_structure took %f minutes.\n", (log_t1-log_t0)/60.0 );
#endif
	
	ofstream br_file;
	if ( mpi_rank == 0 ) {
		br_file.open ( "petsc_solver-3_stepper_bleach_region_total.gnuplot" );
		br_file.precision ( 16 );
		br_file << "plot '-' w l\n";
	}
	
	valarray<double>    br_total (iterations+1);
	PetscGetTime ( &log_t0 );
	br_total[0] = pet_solver->integrate_solution_over_bleach_region();
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "integrate_solution_over_bleach_region took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	if ( mpi_rank == 0 ) {
		br_file << "0.0 " << br_total[0] << "\n";
	}
	
	// wcs is whole chloroplast slice.
	ofstream wcs_file;
	if ( mpi_rank == 0 ) {
		wcs_file.open ( "petsc_solver-3_stepper_whole_chloroplast_slice_total.gnuplot" );
		wcs_file.precision ( 16 );
		wcs_file << "plot '-' w l\n";
	}
	
	valarray<double>    wcs_total (iterations+1);
	
	// integrate_solution_over_domain calls collective VecDot.
	PetscGetTime ( &log_t0 );
	wcs_total[0] = pet_solver->integrate_solution_over_domain();
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "integrate_solution_over_domain took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	if ( mpi_rank == 0 ) {
		wcs_file << "0.0 " << wcs_total[0] << "\n";
	}
	
	// Visualization output is generated using viz_solution_ring.
	// No output is possible for t == 0. The final solution is 
	// created at t == iterations - 1 and the final output is forced
	// when t == iterations without a solve. At all other times, output
	// occurs when t % mpi_size == 0 (except at t == 0 which is impossible).
	for ( int t=0; t<=iterations; ++t ) {		
		//#pragma omp parallel sections
		{
			//#pragma omp section
			if ( (vtk_out > 0) && (0!=t) && ( iterations == t || (t % mpi_size == 0) ) ) {
				// PetscPrintf ( PETSC_COMM_WORLD, "Entering visualization section on iteration %i.\n", t );
				PetscGetTime ( &log_t0 );
				
				// We have solved t iterations up to now and we must work out
				// when solution_ring and viz_solution_ring were swapped.
				int viz_start =  t % mpi_size ==0 ? t - mpi_size : t - ( t % mpi_size );
				
				// t is the number of time steps we have solved for and hence
				// either there are mpi_size outputs to be done, or there are
				// t % mpi_size remaining at the end of the for-loop.
	
				int remaining = t % mpi_size == 0 ? mpi_size : t % mpi_size;
					
				if ( (vtk_out > 0) && (mpi_rank < remaining) ) {
				
					ostringstream name;
					if ( sim_mode == L2_projection ) {
						name << "petsc_solver-3_stepper_" << viz_start + mpi_rank << ".vtp";
					} else {
						name << "petsc_solver-3_stepper_" << viz_start + mpi_rank + 1 << ".vtp";
					}
					
					// PetscPrintf ( PETSC_COMM_SELF, "Computenode %i is attempting to output viz %i.\n", mpi_rank, viz_start + mpi_rank + 1 );
					
					// viz_solution_ring[mpi_rank]->global_evaluate ( points_two_D, vtk_st->scalars );
					
					viz_solution_ring[mpi_rank]->global_evaluate ( points_two_D,                                                           sol_eval_struct,
					                                               vtk_st->scalars );
					
					vtk_st->points_3D[ slice(2, vtk_st->scalars.size(), 3) ] = vtk_st->scalars;
					
					
					if ( vtk_out == 3 ) {
						vtk_append_raw ( *vtk_st, (name.str()).c_str() );
					} else if ( vtk_out ==2 ) {
						vtk_simple_output( *vtk_st, (name.str()).c_str(), true );
					} else {
						vtk_simple_output( *vtk_st, (name.str()).c_str() );
					}
				}
				PetscGetTime ( &log_t1 );
				if ( !min_iter_info ) {
					PetscPrintf ( PETSC_COMM_WORLD,
					              "Visualization within iteration took %f minutes.\n",
					              (log_t1-log_t0)/60.0 );
				}
				// cout << mpi_rank << " is waiting in visualization.\n";
			}
			//#pragma omp section
			if ( iterations != t ) {
				if ( sim_mode == L2_projection ) {
					ostringstream    ss;
					ss << manifest1_dir << manifest1_content[t];
					
					image_initial.reset ( new singleimage );
					image_initial->load ( ss.str(), 1 );
					image_initial->rescale_percentage ();
					
					// Change to a more consistent way of rescaling that does not
					// vary between images. 2009-02-13 ML.
					// image_initial->rescale ( rescale_image1 );
					
					mat_assembler->global_function_to_project = image_initial;
					mat_assembler->assemble_rhs_vector ( *pet_solver, projection_vector );
					pet_solver->rhs_vector_assembly_begin ();
					pet_solver->rhs_vector_assembly_end ();
				} else {
					if ( begin_bleach1 == t ) {
						pet_solver->add_scaled_bleach_matrix ( time_step_dt*bleach_param1 );
					}
					
					if ( end_bleach1 == t ) {
						pet_solver->load_mass_matrix ( 1.0 + time_step_dt*global_bleach1 );
						
						pet_solver->add_scaled_stiffness_matrix ( time_step_dt*diffusion1 );
						
						if ( sim_mode == siggia_clflip ) {
							pet_solver->add_scaled_siggia_matrix ( -time_step_dt*diffusion1 );
						}
					}
					
					pet_solver->rhs_vec_mass_mult_previous_solution ();
				}
				pet_solver->solve ( *solution_ring[ t % mpi_size ], !min_iter_info );
				
				
				// integrate member functions call collective VecDot.
				br_total[t+1]  = pet_solver->integrate_solution_over_bleach_region();
				wcs_total[t+1] = pet_solver->integrate_solution_over_domain();
		
				if ( mpi_rank == 0 ) {
					br_file << (t+1)*time_step_dt << " "
						<< br_total[t+1] << "\n";
					wcs_file << (t+1)*time_step_dt << " "
						<< wcs_total[t+1] << "\n";
				}
				
				
				// PetscPrintf ( PETSC_COMM_SELF, "Computenode %i is attempting to solve into solution_ring %i.\n", mpi_rank, t%mpi_size );
				// cout << mpi_rank << " is waiting in solve.\n";
			}
		}

		//#pragma omp single
		if ( (t % mpi_size == mpi_size - 1) || (iterations - 1 == t) ) {
			// Make solutions available for visualization after every
			// mpi_size solutions are completed and after the very last solve.
			solution_ring.swap ( viz_solution_ring );
		}
		if ( min_iter_info ) {
			PetscPrintf ( PETSC_COMM_WORLD, "*" );
		} else {
			if ( iterations != t ) {
				PetscGetTime ( &iterations_t1 );
				PetscLogDouble elapsed = (iterations_t1-iterations_t0)/60.0;
				PetscPrintf ( PETSC_COMM_WORLD,
					"%i of %i iterations in %f minutes. Estimated %f minutes left (%i x %f seconds).\n",
					t+1,
					iterations,
					elapsed,
					(iterations-t-1)*elapsed/(t+1),
					(iterations-t-1),
					60.0*elapsed/(t+1)
					);
			}
		}
	}
	
	
	
	{
		// wcs is whole chloroplast slice.
		ofstream file;
		if ( mpi_rank == 0 ) {
			file.open ( "petsc_solver-3_stepper_whole_chloroplast_slice_normalized.gnuplot" );
			file.precision ( 16 );
			file << "plot [*:*] [0:1] '-' w l\n";
			
			wcs_total /= wcs_total.max();
			
			valarray<double> &    wcs_normalized = wcs_total;
			
			assert ( wcs_total.size() == iterations + 1 );
			
			for ( int t=0; t<iterations+1; ++t ) {
				file << t*time_step_dt << "\t" << wcs_normalized[t] << "\n";
			}
		}
	}
	
	{
		// wcs is whole chloroplast slice.
		ofstream file;
		if ( mpi_rank == 0 ) {
			file.open ( "petsc_solver-3_stepper_bleach_region_normalized.gnuplot" );
			file.precision ( 16 );
			file << "plot [*:*] [0:1] '-' w l\n";
			
			br_total /= br_total.max();
			
			valarray<double> &    br_normalized = br_total;
			
			assert ( br_total.size() == iterations + 1 );
			
			for ( int t=0; t<iterations+1; ++t ) {
				file << t*time_step_dt << "\t" << br_normalized[t] << "\n";
			}
		}
	}
	
	{
		rusage    r_st;
		int ret = getrusage ( RUSAGE_SELF, &r_st );
		if ( !ret ) {
			PetscPrintf ( PETSC_COMM_WORLD, "maximum resident set size was %i.\n", r_st.ru_maxrss );
		} else {
			PetscPrintf ( PETSC_COMM_WORLD, "There was a problem.\n" );
		}
	}
	
	// Destroy petsc objects before finalize!
	pet_solver.reset();
	
	// Hope to avoid incomplete output returned by PBS.
	MPI_Barrier ( PETSC_COMM_WORLD );
	
	PetscGetTime ( &totalruntime_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Total run time was %f minutes.\n", (totalruntime_t1-totalruntime_t0)/60.0 );
	PetscErrorCode ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}
