#ifdef _OPENMP
#include <omp.h>
#endif

#undef TRACE
#define TRACE cout << __FILE__ << ":" << __LINE__

#include "box.hh"
#include "cover_structure.hh"
#include "basic_dbinary_tree.hh"
#include "dbinary_tree.hh"
#include "dbinary_tree_utils.hh"
#include "geometry.hh"
#include "global_approximation_space.hh"
#include "global_basis_function.hh"
#include "integration_scheme.hh"
#include "petsc_solver.hh"
#include "polynomial.hh"
#include "quadrature_rule.hh"
#include "solution.hh"
#include "vtk_output.hh"

#include <petscts.h>

#include <boost/shared_ptr.hpp>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <valarray>

using boost::shared_ptr;
using std::clog;
using std::cout;
using std::endl;
using std::exp;
using std::ofstream;
using std::slice;
using std::valarray;

#define TRACE std::clog << __FILE__ << ":" << __LINE__ << "\n"

class unity : public function<2> {
public:
	unity () {
		this->set_global_function();
	}
private:
	double evaluate( const double* x ) const {
		return 1.0;
	}
};

class my_func : public function<2> {
public:
	my_func() {
		this->set_global_function();
	}
private:
	double evaluate( const double* x ) const {
		assert( x );
		return x[0]*x[1]/10.0;
	}
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
		
		if ( sqrt ( x0_x0 + x1_x1 ) < 1.0 ) {
			return 3.0;
		} else {
			return 1.0;
		}
	}
};

class smooth_bump : public function<2> {
public:
	smooth_bump ( double inner_r = 2.9, double outer_r = 3,
	              double inner_val = 2.0, double outer_val = 1.0 )
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

struct matrix_assembler {
	typedef shared_ptr<matrix_assembler> ptr;
	
	matrix_assembler () : integration_decomposition ( 0 ) { }
	
#if 0
	void simple_assemble_all_and_solve ( petsc_solver<2>& solver, solution<2>& sol ) {
		
		clog << "Begin computenode local assembly of matrix.\n";
		
		assemble_matrix ( solver );
		
		clog << "Begin global assembly of matrix.\n";
		
		solver.matrix_assembly_begin ();
		solver.matrix_assembly_end ();
		
		assemble_rhs_vector ( solver );
		
		solver.rhs_vector_assembly_begin ();
		solver.rhs_vector_assembly_end ();
			
		clog << "Begin global solve.\n";
		
		solver.solve( sol );
	}
#endif

	void assemble_matrix ( petsc_solver<2>& solver ) {
		assert( cover_st );

		const vector<vector<int> >& index_to_neighbour_indices = cover_st->index_to_neighbour_indices;
		
		valarray<int> patch_indices_on_computenode;
		
		solver.get_patch_indices_on_computenode( patch_indices_on_computenode );

		assert( patch_indices_on_computenode.size() > 0 ); // Not strictly necessary.

		int num = patch_indices_on_computenode.size();
		
#if 0 // ndef NDEBUG // Keep for future debug.
		std::clog << "Patch indices on computenode are : ";
		for ( int i=0; i<num; ++i ) {
			if ( i!=0 ) {
				std::clog << ", ";
			}
			std::clog << patch_indices_on_computenode[i];
		}
		std::clog << "\n";
#endif

		valarray<double> tmp_vals;
		
#ifdef _OPENMP
		if ( omp_get_max_threads() < omp_get_num_procs() ) {
			omp_set_num_threads( omp_get_num_procs() );
		}
#endif

		int first_patch_ind;
		int num_neigh;
		int second_patch_ind;
		
#pragma omp parallel for private(tmp_vals,first_patch_ind,num_neigh,second_patch_ind)
		for ( int i=0; i<num; ++i ) {
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

#if 0 //def _OPENMP // Keep for future debug.
			cout << omp_get_thread_num()
				<< " is doing patch "
				<< first_patch_ind
				<< " which has "
				<< num_neigh
				<< " neighbours.\n";
#endif
			
			assert( num_neigh > 0 ); // Not strictly necessary.
			
#if 0 // ndef NDEBUG // Keep for future debug.
			std::clog << "Starting matrix patch block row " << first_patch_ind << ".\n";
			std::clog << "Claimed neighbours are : ";
			for ( int j=0; j<num_neigh; ++j ) {
				if ( j!=0 ) {
					std::clog << ", ";
				}
				std::clog << neighbours[j];
			}
#endif

			for ( int j=0; j<num_neigh; ++j ) {

				second_patch_ind = neighbours[j];

#if 0 // ndef NDEBUG // Keep for future debug.
				std::clog << "Assembling matrix patch block " << first_patch_ind << ", " << second_patch_ind << ".\n";
#endif
				
				assemble_matrix_block( first_patch_ind,
				                       second_patch_ind,
				                       tmp_vals );
				
				// PETSc is not thread safe.
				// A critical section is now marked at the point of
				// value insertion into the PETSc matrix.
				solver.set_matrix_entries( first_patch_ind, second_patch_ind, tmp_vals );
			}
		}
	}

	void assemble_matrix_block ( int patch_index_row, int patch_index_column, valarray<double>& vals ) {
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

		vals.resize( test_size * trial_size );

		for ( int i=0; i<test_size; ++i ) {
			for ( int j=0; j<trial_size; ++j ) {
#if 0 //ndef NDEBUG // Keep for future debug.
				std::clog << "\tAssembling row " << i << ", column " << j << " of matrix patch block " << patch_index_row << ", " << patch_index_column << ".\n";
#endif
				assemble_single_matrix_entry( test[i], trial[j], vals[ i*trial_size + j ] );
			}
		}
#if 0 // ndef NDEBUG // Keep for future debug.
		std::clog << "\t\tRow finished for " << patch_index_row << ", " << patch_index_column << " of matrix patch blocks.\n";
		for ( int i=0; i<test_size; ++i ) {
			std::clog << "\t\t";
			for ( int j=0; j<trial_size; ++j ) {
				std::clog << vals[ i*trial_size + j ];
				if ( j != trial_size-1 ) {
					std::clog << "  ";
				}
			}
			std::clog << "\n";
		}
#endif
	}

	void assemble_single_matrix_entry ( const global_basis_function<2>& test,
	                                    const global_basis_function<2>& trial,
	                                    double&                         matrix_entry )
	{
		assert( integrator );
		
		matrix_entry     = integrator->integrate_product_global_basis_function ( test, trial, integration_decomposition );
	
#if 0 // def _OPENMP
		cout << "matrix_entry calculated by OpenMP thread "
			<< omp_get_thread_num() << " : " << matrix_entry << "\n";
#endif
		
#if 0 //ndef NDEBUG // Keep for future debug.
		std::clog << "\t\tMatrix entry is " << matrix_entry << "\n";
#endif
	}
	
	void assemble_rhs_vector ( petsc_solver<2>& solver ) {
		assert( cover_st );
		
		valarray<int> patch_indices_on_computenode;
		solver.get_patch_indices_on_computenode( patch_indices_on_computenode );
		
		assert( patch_indices_on_computenode.size() > 0 ); // Not strictly necessary.
		
		int num = patch_indices_on_computenode.size();
		
		valarray<double> tmp_vals;
		
#pragma omp parallel for private(tmp_vals)
		for ( int i=0; i<num; ++i ) {
			assemble_rhs_vector_block( patch_indices_on_computenode[i], tmp_vals );
			
			// PETSc is not thread safe.
			// A critical section is now marked at the point of
			// value insertion into the PETSc vector.
			solver.set_rhs_vector_entries ( patch_indices_on_computenode[i], tmp_vals );
		}
	}
	
	void assemble_rhs_vector_block ( int patch_index, valarray<double>& values ) {
		vector<global_basis_function<2> > test;
		
		global_approx_sp->get_global_basis_by_patch( patch_index, test );
		
		int size = test.size();
		
		values.resize(size);
		
		for ( int i=0; i<size; ++i ) {
			assemble_single_rhs_vector_entry( test[i], values[i] );
		}	
	}
	

	void assemble_single_rhs_vector_entry ( const global_basis_function<2>& test,
	                                        double&                         rhs_vector_entry )
	{
		assert( global_function_to_project );
		rhs_vector_entry = integrator->integrate_product ( *global_function_to_project, test, integration_decomposition );
	}

	global_approximation_space<2>::ptr global_approx_sp;
	integration_scheme<2>::ptr         integrator;
	function<2>::ptr                   global_function_to_project;
	cover_structure<2>::ptr            cover_st;
	int                                integration_decomposition;
};

static char help[] = "petsc_solver-4";

int main (int argc, char* argv[]) {
	PetscInitialize( &argc, &argv, static_cast<char*>(0), help );
	
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
	
	int petsc_log_stages[5];
	PetscLogStageRegister ( &petsc_log_stages[0], "Update dbinary_tree." );
	PetscLogStageRegister ( &petsc_log_stages[1], "Cover structure creation." );
	PetscLogStageRegister ( &petsc_log_stages[2], "Gnuplot file output." );
	PetscLogStageRegister ( &petsc_log_stages[3], "Assemble and solve." );
	PetscLogStageRegister ( &petsc_log_stages[4], "VTK output." );
	
	const int& log_stage_dbinary_tree       = petsc_log_stages[0];
	const int& log_cover_structure_creation = petsc_log_stages[1];
	const int& log_stage_gnuplot_output     = petsc_log_stages[2];
	const int& log_stage_assemble_solve     = petsc_log_stages[3];
	const int& log_stage_vtk_output         = petsc_log_stages[4];

	PetscLogDouble log_t0, log_t1;
	
	// Use to divide work between MPI processes.
	int round_robin_counter = 0;
	
	PetscReal parent_box_expansion_factor = 1.2;
	PetscOptionsGetReal(PETSC_NULL,"-parent_box_expansion_factor", &parent_box_expansion_factor, PETSC_NULL);
	if ( parent_box_expansion_factor<=0.0 ) {
		parent_box_expansion_factor = 1.2;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "parent_box_expansion_factor is %f.\n", parent_box_expansion_factor );
	
	PetscInt dtree_initial_refinement = 2;
	PetscOptionsGetInt(PETSC_NULL,"-dtree_initial_refinement", &dtree_initial_refinement, PETSC_NULL);
	if ( dtree_initial_refinement<0 ) {
		dtree_initial_refinement = 2;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "dtree_initial_refinement is %i.\n", dtree_initial_refinement );
	
	PetscInt dtree_max_refinement = dtree_initial_refinement + 2;
	PetscOptionsGetInt(PETSC_NULL,"-dtree_max_refinement", &dtree_max_refinement, PETSC_NULL);
	if ( dtree_max_refinement<0 ) {
		dtree_max_refinement = dtree_initial_refinement + 2;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "dtree_max_refinement is %i.\n", dtree_max_refinement );
	
	PetscReal dtree_cover_factor = 1.3;
	PetscOptionsGetReal(PETSC_NULL,"-dtree_cover_factor", &dtree_cover_factor, PETSC_NULL);
	if ( dtree_cover_factor<=0.0 ) {
		dtree_cover_factor = 1.3;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "dtree_cover_factor is %f.\n", dtree_cover_factor );
	
	PetscInt initial_function_choice = 5;
	PetscOptionsGetInt ( PETSC_NULL, "-initial_function_choice", &initial_function_choice, PETSC_NULL );
	if ( initial_function_choice < 1 || 5 < initial_function_choice ) {
		initial_function_choice = 5;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "initial_function_choice is %i.\n", initial_function_choice );
	
	PetscInt assembly_quadrature_rule = 5;
	PetscOptionsGetInt ( PETSC_NULL, "-assembly_quadrature_rule", &assembly_quadrature_rule, PETSC_NULL );
	if ( assembly_quadrature_rule < 1 || 5 < assembly_quadrature_rule ) {
		assembly_quadrature_rule = 5;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "assembly_quadrature_rule is %i.\n", assembly_quadrature_rule );
	
	PetscInt local_basis_order = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-local_basis_order", &local_basis_order, PETSC_NULL );
	if ( local_basis_order < 1 || 2 < local_basis_order ) {
		local_basis_order = 1;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "local_basis_order is %i.\n", local_basis_order );
	
	PetscInt integration_decomposition = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-integration_decomposition", &integration_decomposition, PETSC_NULL );
	if ( integration_decomposition < 0 || 1 < integration_decomposition ) {
		integration_decomposition = 1;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "integration_decomposition is %i.\n", integration_decomposition );
	
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
	
	valarray<double> vcoords( num_doubles );
	for ( int i=0; i<num_doubles; ++i ) {
		vcoords[i] = coords[i];
	}
	
	// Set up the geometry.
	shared_ptr<geometry<2> > geom( new geometry<2> );
	geom->create_boundary( "chloroplast envelope", vcoords );
	geom->create_region( "chloroplast", "+chloroplast envelope" );
	
	ofstream file_10( "petsc_solver-4_chloroplast" );
	
	geom->gp_draw( "chloroplast", file_10 );
	file_10.close();
	
	geom->update();
	
	
	// Set up the d-binary tree used to form the
	// partition of unity method cover.
	shared_ptr<dbinary_tree<2> > dtree( new dbinary_tree<2> );
	dtree->set_max_level ( dtree_max_refinement );
	dtree->set_geometry( geom, parent_box_expansion_factor );
	
	ofstream file_15( "petsc_solver-4_parent_box" );
	
	dtree->gp_draw_parent_box( file_15 );
	file_15.close();
	
	PetscGetTime ( &log_t0 );
	dtree->set_initial_refinement( dtree_initial_refinement ); // 1D 2^n, 2D 4^n, 3D 8^n so refinement can be expensive.
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "dbinary_tree initial refinement took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	dtree->set_cover_factor( dtree_cover_factor ); // Possibly use "length extension" to be descriptive.
	
	// Create covers for the regions and boundaries.
	// TODO: allow the option not to cover unwanted areas.
	
	clog << "Begin update of dbinary_tree.\n";


	PetscGetTime ( &log_t0 );
	
	PetscLogStagePush ( log_stage_dbinary_tree );
	dtree->update();
	PetscLogStagePop (  );
		
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "dbinary_tree update took %f minutes.\n", (log_t1-log_t0)/60.0 );

	
#if 0
	clog << "Begin gnuplot file output.\n";
	PetscLogStagePush ( log_stage_gnuplot_output );

	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		ofstream file( "petsc_solver-4_boundary_cover" );
		dtree->gp_draw_boundary_cover( "chloroplast envelope", file, 0.3 );
		file.close();
	}

	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		ofstream file( "petsc_solver-4_region_boundary_cover" );
		dtree->gp_draw_region_boundary_cover( "chloroplast", file, 0.1 );
		file.close();
	}

	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		ofstream file( "petsc_solver-4_region_interior_cover" );
		dtree->gp_draw_region_interior_cover( "chloroplast", file, 0.2 );
		file.close();
	}
	PetscLogStagePop (  );
#endif
	

	// dtree->dump_cout();
	
	global_approximation_space<2>::ptr global_approx_space;
	
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
	
	PetscLogStagePush ( log_cover_structure_creation );
	{
		cover_structure<2>::ptr cs;
		assert( dtree );
		cout << "Begin creation of cover_structure.\n";
		PetscGetTime ( &log_t0 );
		dtree->get_cover_structure( "chloroplast", cs );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "get_cover_structure took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
		mat_assembler->cover_st = cs;
		
		global_approx_space->give_cover_structure( cs );
		assert( !cs ); // Safety reset.
	}
	PetscLogStagePop ();
	
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
		mat_assembler->global_function_to_project.reset( new unity );
		break;
	case 2:
		mat_assembler->global_function_to_project.reset( new my_func );
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
	
	mat_assembler->integration_decomposition = integration_decomposition;
		
	valarray<int> patch_to_num_dof;
	
	global_approx_space->get_patch_to_num_dof( patch_to_num_dof );
	
	petsc_solver<2>::ptr pet_solver ( new petsc_solver<2> );
	
	cout << "Begin creation of PETSc objects.\n";
	PetscGetTime ( &log_t0 );
	pet_solver->create_petsc_objects ( patch_to_num_dof );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "create_petsc_objects took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	solution<2>::ptr sol ( new solution<2> );
	
	PetscGetTime ( &log_t0 );
	sol->set_global_approximation_space ( global_approx_space, patch_to_num_dof.sum() );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "set_global_approximation_space took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscLogStagePush ( log_stage_assemble_solve );

	PetscPrintf ( PETSC_COMM_WORLD, "Begin computenode local assembly of matrix.\n" );
	PetscGetTime ( &log_t0 );
	mat_assembler->assemble_matrix ( *pet_solver );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Local assembly of matrix took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	
	PetscPrintf ( PETSC_COMM_WORLD, "Begin global assembly of matrix.\n" );
	PetscGetTime ( &log_t0 );
	pet_solver->matrix_assembly_begin ();
	pet_solver->matrix_assembly_end ();
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Assembly of global matrix took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	cout << "Begin computenode local assembly of right-hand side vector.\n";
	PetscGetTime ( &log_t0 );
	mat_assembler->assemble_rhs_vector ( *pet_solver );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Local assembly of right-hand side vector took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	cout << "Begin global assembly of right-hand side vector.\n";
	PetscGetTime ( &log_t0 );
	pet_solver->rhs_vector_assembly_begin ();
	pet_solver->rhs_vector_assembly_end ();
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Assembly of global right-hand side vector took %f minutes.\n", (log_t1-log_t0)/60.0 );

	cout << "Begin global solve.\n";
	PetscGetTime ( &log_t0 );
	pet_solver->solve( *sol );
	PetscGetTime ( &log_t1 );
	PetscPrintf ( PETSC_COMM_WORLD, "Solve and distribution of solution took %f minutes.\n", (log_t1-log_t0)/60.0 );
	
	PetscLogStagePop ();

	PetscLogStagePush ( log_stage_vtk_output );

	basic_dbinary_tree<2>::ptr basic_dtree ( new basic_dbinary_tree<2> );
	basic_dtree->initialize ( dtree->access_parent_box(), dtree_max_refinement );
	
	// Visualization of solution.
	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		cout << "Entering visualization of solution part.\n";
		PetscGetTime ( &log_t0 );
		
		vtk_output::ptr vtk_st ( new vtk_output );
		
		basic_dtree->get_grid_points_3D ( vtk_st->points_3D );
		
		int num_pts = (vtk_st->points_3D).size() / 3;
		
		valarray<double> points_two_D ( 2*num_pts );
		const valarray<double>& points_three_D = vtk_st->points_3D;
		
		points_two_D[ slice(0, num_pts, 2) ] = points_three_D[ slice(0, num_pts, 3) ];
		points_two_D[ slice(1, num_pts, 2) ] = points_three_D[ slice(1, num_pts, 3) ];
		
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_SELF, "Generation of 3D points took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		PetscGetTime ( &log_t0 );
		sol->global_evaluate ( points_two_D, vtk_st->scalars );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_SELF, "Solution evaluation took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		PetscGetTime ( &log_t0 );
		
		vtk_st->points_3D[ slice(2, vtk_st->scalars.size(), 3) ] = vtk_st->scalars;
		vector<string> all_keys_to_max_level;
		
		generate_keys<2> ( dtree_max_refinement, all_keys_to_max_level );
		
		basic_dtree->get_box_grid_connectivity_offsets ( all_keys_to_max_level, vtk_st->connectivity, vtk_st->offsets );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_SELF, "Preparing vtk information took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		{
			ofstream file ( "petsc_solver-4_solution.vtp" );
			PetscGetTime ( &log_t0 );
			vtk_simple_output( *vtk_st, file );
			PetscGetTime ( &log_t1 );
			PetscPrintf ( PETSC_COMM_SELF, "vtk_simple_output of solution took %f minutes.\n", (log_t1-log_t0)/60.0 );
		}
	}
	PetscLogStagePop (  );
	
	if ( round_robin_counter++ % mpi_size == mpi_rank )
	{
		cout << "Visualization of function to project.\n";
		PetscGetTime ( &log_t0 );
		
		vtk_output::ptr vtk_st ( new vtk_output );
		
		basic_dtree->get_grid_points_3D ( vtk_st->points_3D );
		
		int num_pts = (vtk_st->points_3D).size() / 3;
		
		valarray<double> points_two_D ( 2*num_pts );
		const valarray<double>& points_three_D = vtk_st->points_3D;
		
		points_two_D[ slice(0, num_pts, 2) ] = points_three_D[ slice(0, num_pts, 3) ];
		points_two_D[ slice(1, num_pts, 2) ] = points_three_D[ slice(1, num_pts, 3) ];
		
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_SELF, "Generation of 3D points took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		PetscGetTime ( &log_t0 );
		mat_assembler->global_function_to_project->global_evaluate ( box<2>(), points_two_D, vtk_st->scalars );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_SELF, "Evaluate of function to project took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		assert ( (vtk_st->scalars).size() == num_pts );
		
		PetscGetTime ( &log_t0 );
		
		vtk_st->points_3D[ slice(2, vtk_st->scalars.size(), 3) ] = vtk_st->scalars;
		vector<string> all_keys_to_max_level;
		
		generate_keys<2> ( dtree_max_refinement, all_keys_to_max_level );
		
		basic_dtree->get_box_grid_connectivity_offsets ( all_keys_to_max_level, vtk_st->connectivity, vtk_st->offsets );
		PetscGetTime ( &log_t1 );
		PetscPrintf ( PETSC_COMM_SELF, "Preparing vtk information took %f minutes.\n", (log_t1-log_t0)/60.0 );
		
		{
			ofstream file ( "petsc_solver-4_function_to_project.vtp" );
			PetscGetTime ( &log_t0 );
			vtk_simple_output( *vtk_st, file );
			PetscGetTime ( &log_t1 );
			PetscPrintf ( PETSC_COMM_SELF, "vtk_simple_output of function to project took %f minutes.\n", (log_t1-log_t0)/60.0 );
		}
	}
	
	// Destroy petsc objects before finalize!
	pet_solver.reset();
	
	PetscErrorCode ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}
