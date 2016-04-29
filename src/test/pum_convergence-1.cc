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

#include "extract_geometry_2d.hh"
#include "global_approximation_space.hh"
#include "pum_convergence.hh"
#include "refinement_structure.hh"
#include "vtk_output.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

using std::copy;
using std::cout;
using std::fabs;
using std::inserter;
using std::log10;
using std::ostringstream;
using std::set_union;
using std::string;
using std::vector;

/*
	Rewrite 'a' as (a)(10^b).
	Used to tabulate convergence data. The interface is based on the
	base-2 frexp function.
*/
double frexp_dec ( double a, int* b ) {	
	assert ( b );
	*b = 0;
	
	if ( 0.0 == fabs(a) ) {
		return a;
	}
	
	while ( fabs(a) < 1 ) {
		a *= 10.0;
		--(*b);
	}
	
	while ( !(fabs(a)<10.0) ) {
		a /= 10.0;
		++(*b);
	}
	
	return a;
}

class simple_fun : public differentiable_function<2> {
public:
	simple_fun ( ) { this->set_global_function(); }
private:
	double evaluate    ( const double * co ) const {
		// L2 projection
		//return 1.0;
		
		// Unit box Helmholtz. (1)
		return (2.0*co[0]-3.0)*co[0]*co[0];
	}
	void evaluate_grad ( const double * co, double *   grad ) const {
		// L2 projection
		//grad[0]=0.0; grad[1]=0.0;
	
		// Unit box Helmholtz. (1)
		grad[0] = 6.0*(co[0]-1.0)*co[0]; grad[1] = 0.0;
	}
};

class rhs_fun : public function<2> {
public:
	rhs_fun ( ) { this->set_global_function(); }
private:
	double evaluate    ( const double * co ) const {
		// L2 projection
		//return 1.0;

		// Unit box Helmholtz. (1)
		return co[0]*(co[0]*(2.0*co[0]-3.0) - 12.0)+6.0;
	
		// Poisson
		//return -6.0*co[0];
	}
};

static const char* help = "pum_convergence-1";

int main ( int argc, char* argv[] ) {
	const int num_decimal_digits = std::numeric_limits<double>::digits10;

	PetscInitialize( &argc, &argv, PETSC_NULL, help );

	PetscMPIInt                        comm_rank;
	PetscMPIInt                        comm_size;

	MPI_Comm_rank( PETSC_COMM_WORLD, &comm_rank );
	MPI_Comm_size( PETSC_COMM_WORLD, &comm_size );

	int round_robin_counter = 0;

	PetscLogDouble log_t0, log_t1;
	
	PetscPrintf ( PETSC_COMM_WORLD,
	              "====== Beginning execution of pum_convergence-1 ======\n" );
	
#ifdef NDEBUG
	PetscPrintf ( PETSC_COMM_WORLD, "NDEBUG has been defined.\n" );
#else
	PetscPrintf ( PETSC_COMM_WORLD, "Preprocessor symbol NDEBUG not defined.\n" );
#endif
	
#ifdef _OPENMP
	PetscPrintf ( PETSC_COMM_WORLD, "_OPENMP has been defined and omp_get_max_threads returned %i.\n", omp_get_max_threads() );
#else
	PetscPrintf ( PETSC_COMM_WORLD, "Preprocessor symbol _OPENMP not defined.\n" );
#endif

	PetscPrintf ( PETSC_COMM_WORLD, "std::numerical_limits<double> is %i.\n", num_decimal_digits );

	if ( comm_size > 1 ) {
		PetscPrintf ( PETSC_COMM_WORLD, "There are %i MPI processes.\n", comm_size );
	}

	string    dom_string;
	{
		char    tmp_dom[256];
		PetscTruth    given = PETSC_FALSE;
		PetscOptionsGetString ( PETSC_NULL, "-dom", tmp_dom, 255, &given );
		if ( given ) {
			dom_string.assign ( tmp_dom );
		} else {
			dom_string = "unit_box.off";
		}
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-dom is \"%s\".\n", dom_string.c_str() );
	
	PetscTruth    anticlockwise_dom = PETSC_TRUE;
	PetscOptionsGetTruth ( PETSC_NULL, "-anticlockwise_dom", &anticlockwise_dom, PETSC_NULL );
	
	PetscInt initial_refinement = 4;
	PetscOptionsGetInt(PETSC_NULL,"-initial_refinement", &initial_refinement, PETSC_NULL);
	if ( initial_refinement<0 ) {
		initial_refinement = 4;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-initial_refinement is %i.\n", initial_refinement );
	
	PetscInt max_refinement = initial_refinement;
	PetscOptionsGetInt(PETSC_NULL,"-max_refinement", &max_refinement, PETSC_NULL);
	PetscPrintf ( PETSC_COMM_WORLD, "\t-max_refinement is %i.\n", max_refinement );

	PetscReal parent_box_expansion_factor = 1.01;
	PetscOptionsGetReal(PETSC_NULL,"-parent_box_expansion_factor", &parent_box_expansion_factor, PETSC_NULL);
	if ( !(parent_box_expansion_factor>1.0) ) {
		PetscPrintf( PETSC_COMM_WORLD, "\tInvalid -parent_box_expansion_factor so using a default value.\n" );
		parent_box_expansion_factor = 1.01;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-parent_box_expansion_factor is %f.\n", parent_box_expansion_factor );
	
	PetscReal cover_factor = 1.3;
	PetscOptionsGetReal(PETSC_NULL,"-cover_factor", &cover_factor, PETSC_NULL);
	if ( cover_factor<=0.0 ) {
		cover_factor = 1.3;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-cover_factor is %f.\n", cover_factor );

	PetscInt local_basis_order = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-local_basis_order", &local_basis_order, PETSC_NULL );
	if ( local_basis_order < 1 || 2 < local_basis_order ) {
		local_basis_order = 1;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-local_basis_order is %i.\n", local_basis_order );
	
	PetscInt box_interior_quad_rule = 4;
	PetscOptionsGetInt ( PETSC_NULL, "-box_interior_quad_rule", &box_interior_quad_rule, PETSC_NULL );
	if ( box_interior_quad_rule < 1 ) {
		box_interior_quad_rule = 4;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-box_interior_quad_rule is %i.\n", box_interior_quad_rule );

	PetscInt integration_decomposition = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-integration_decomposition", &integration_decomposition, PETSC_NULL );
	if ( integration_decomposition < 0 || 1 < integration_decomposition ) {
		integration_decomposition = 1;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-integration_decomposition is %i.\n", integration_decomposition );

	valarray<double> vcoords;
	simple_off_input_geometry_2d ( dom_string, vcoords );
	assert ( vcoords.size()<1e7 );
	
	pum_details<2>::ptr    details ( new pum_details<2> );
	
	line_segments &    domain_segs = details->domain;
	{
		int tmp_size = vcoords.size();
		domain_segs.segments.resize ( tmp_size +2 );
		copy ( &vcoords[0], &(vcoords[0])+tmp_size, &(domain_segs.segments[0]) );
		
		// Debugging found the last point was being set incorrectly.
		// Using -2 and -1 to adjust indices. 2008-12-01 Michael LI.
		domain_segs.segments[ tmp_size ]   = vcoords[0];
		domain_segs.segments[ tmp_size+1 ] = vcoords[1];
		
		// Keep in case of need to inspect
		 int num_segs = tmp_size/2;
		if ( ++round_robin_counter % comm_size == comm_rank ) {
			ofstream    file ( "pum_convergence-1_labelled_bdry.gnuplot" );
			gnuplot_output<2> ( domain_segs, file, (num_segs<1000) );
		}
	}
	
	box<2>    parent_box;
	get_bounding_box ( domain_segs, parent_box );
	
#if 0
	cout << "Debug: set unit box parent_box.\n";
	double parboxmod[4] = { 0.0, 1.0, 0.0, 1.0 };
	parent_box.set( parboxmod );
#endif
	
	parent_box.scale ( parent_box_expansion_factor );
	
	basic_dbinary_tree<2>::ptr &    basic_dtree = details->dtree;
	basic_dtree.reset ( new basic_dbinary_tree<2> );
	basic_dtree->initialize ( parent_box, max_refinement );
	
	refinement_structure<2>::ptr &    ref_struct = details->ref_struct;
	ref_struct.reset ( new refinement_structure<2> );
	ref_struct->dtree = basic_dtree;
	ref_struct->levels["0"] = initial_refinement;
	ref_struct->cover_factor = cover_factor;

	details->viz_ref_struct.reset ( new refinement_structure<2> );
	details->viz_ref_struct->dtree = basic_dtree;
	
	// If any larger than max_refinement, basic_dtree must be considered
	// or there will be an infinite loop during the processing. 2009-06-27 ML.
	details->viz_ref_struct->levels["0"] = max_refinement;
	
	details->viz_ref_struct->cover_factor = cover_factor;

	valarray<double>                               box_cubature_points;
	valarray<double>                               box_cubature_weights;

	// create_petsc_objects must be adjusted to allow for the correct number of expected
	// non-zero entries per row.
	switch ( local_basis_order ) {
	case 1:
		details->local_approx_space = monomial1;
		break;
	case 2:
		details->local_approx_space = monomial2;
		break;
	default:
		details->local_approx_space = monomial1;
		break;
	}
	
	generate_gauss_legendre_rule( box_interior_quad_rule, details->interior_box_cubature );
	generate_gauss_legendre_rule( box_interior_quad_rule, details->boundary_base_rule );
	
	details->fun_ptr.reset ( new simple_fun );
	details->rhs_ptr.reset ( new rhs_fun );
	
	details->integration_decomposition = integration_decomposition;
	details->anticlockwise_dom = anticlockwise_dom;
	
	pum_convergence<2>::ptr    pum_simulation ( new pum_convergence<2> );
	convergence_results    results;

	const int max_iterations = max_refinement - initial_refinement + 1;

	// Modelled after struct convergence_results.
	// Using vectors to allow arbitrary iterations in
	// the future.
	vector<double>    cover_factor_vec;
	vector<int>       quadrature_order_vec;
	vector<int>       refinement_level;
	vector<int>       degrees_of_freedom;

	vector<double>    fun_val_integral;
	vector<double>    approx_val_integral;
	
	vector<double>    approx_val_squared;
	vector<double>    approx_val_fun_val;
	vector<double>    fun_val_squared;
	
	vector<double>    approx_grad_squared;
	vector<double>    approx_grad_fun_grad;
	vector<double>    fun_grad_squared;
	
	vector<double>    L2_error;
	vector<double>    W12_seminorm;
	vector<double>    W12_error;
	
	vector<double>    rel_L2_error;
	vector<double>    rel_W12_seminorm;
	vector<double>    rel_W12_error;
	
	if ( ++round_robin_counter % comm_size == comm_rank ) {
		int num_to_do = max_iterations;
		pvd_output pvd;
		pvd.head = "pum_convergence-1_solution_";
		pvd.middle.resize ( num_to_do );
		ostringstream oss;
		for ( int i=0; i<num_to_do; ++i ) {
			oss.str ( "" );
			oss << i;
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "pum_convergence-1_solution.pvd" );
		pvd_simple_output ( pvd, file );
	}
	
	if ( ++round_robin_counter % comm_size == comm_rank ) {
		int num_to_do = max_iterations;
		pvd_output pvd;
		pvd.head = "pum_convergence-1_L2_error_integrand_";
		pvd.middle.resize ( num_to_do );
		ostringstream oss;
		for ( int i=0; i<num_to_do; ++i ) {
			oss.str ( "" );
			oss << i;
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "pum_convergence-1_L2_error_integrand.pvd" );
		pvd_simple_output ( pvd, file );
	}
	
	PetscLogDouble single_iteration_t0, single_iteration_t1;
	
	// Errors are important enough that we allow drawing of the graph before completion of simulations.
	ofstream    errors_file;
	if ( 0 == comm_rank ) {
		errors_file.open ( "pum_convergence-1_errors.gnuplot" );
		if ( !errors_file ) {
			abort();
		}
		errors_file.precision(num_decimal_digits);
		errors_file << "set logscale xy\n";
		errors_file << "plot 'pum_convergence-1_errors.gnuplot' using 2:3 w linesp t 'L2', '' using 2:4 w linesp t 'W^{1,2} semi-norm', '' using 2:5 w linesp t 'W^{1,2}'\n";
		
		// As we use the same file to hold our data we don't want our data being
		// interpreted as commands.
		errors_file << "exit\n";
	}
	
	shared_ptr<diagnostic_quad_details> diag_quad_d ( new diagnostic_quad_details );
	
	pum_return_information    pum_return_info;
	
	for ( int iter=0; iter<max_iterations; ++iter ) {
		
		cover_factor_vec.push_back ( cover_factor );
		quadrature_order_vec.push_back ( box_interior_quad_rule );
		
		PetscGetTime ( &single_iteration_t0 );
		
		PetscPrintf ( PETSC_COMM_WORLD, "====== ITERATION %i with REFINEMENT LEVEL %i. =====\n", iter, initial_refinement + iter );
		ref_struct->levels["0"] = initial_refinement + iter;
		
		// This allows for arbirary ordering of refinement levels.
		refinement_level.push_back ( initial_refinement + iter );
		
		bool diag_output = false;
		
		PetscPrintf ( PETSC_COMM_WORLD, "Beginning initialize.\n" ); PetscGetTime ( &log_t0 );
		if ( ref_struct->levels["0"] < 7 ) {
			PetscPrintf ( PETSC_COMM_WORLD, "Integration points will be collected and the cover will be drawn for this iteration.\n" );
			pum_simulation->initialize ( details, &pum_return_info, &(*diag_quad_d) );
			
			ostringstream    name;
			name << "pum_convergence-1_cover_" << iter << ".gnuplot";
			
			pum_simulation->draw_cover ( (name.str()).c_str(), true );
			diag_output = true;
		} else {
			pum_simulation->initialize ( details, &pum_return_info );
			diag_output = false;
		}
		PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, "initialize took %.2f minutes.\n", (log_t1-log_t0)/60.0 );

		PetscPrintf ( PETSC_COMM_WORLD,
		              "\tThere are %i degrees of freedom across %i MPI processes.\n",
		              pum_return_info.num_dof,
		              comm_size );
		
		degrees_of_freedom.push_back ( pum_return_info.num_dof );
		
		if ( diag_quad_d && diag_output ) {
			ostringstream    name;
			name << "pum_convergence-1_diag_quad_d_" << iter << ".gnuplot";

			int num_d_pts = diag_quad_d->all_interior_points.size()/2;
			
			cout << "Num pts is " << num_d_pts << "\n";
			
			ofstream    file ( (name.str()).c_str() );
			
			file << "set size square\nplot '-' w p pt 0\n";
			for ( int i=0; i<num_d_pts; ++i ) {
				file << diag_quad_d->all_interior_points[2*i] << " "
					<< diag_quad_d->all_interior_points[2*i+1] << "\n";
			}
			diag_quad_d->all_interior_points.clear();
		}

		//pum_simulation->set_problem_L2_projection();
		pum_simulation->set_problem_Helmholtz ( -1.0, 1.0 );
		//pum_simulation->set_problem_Poisson ( -1.0 );
		


		PetscPrintf ( PETSC_COMM_WORLD, "Beginning solve_and_get_results.\n" ); PetscGetTime ( &log_t0 );
		pum_simulation->solve_and_get_results ( results );
		PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, "solve_and_get_results took %.2f minutes.\n", (log_t1-log_t0)/60.0 );
		
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Function                     : %.15f\n", results.function_val_integral );
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Approx                       : %.15f\n\n", results.approx_val_integral );

		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Function squared             : %.15f\n", results.function_val_squared_integral );
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Function approx product      : %.15f\n", results.approx_val_function_val_integral );
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Approx squared               : %.15f\n\n", results.approx_val_squared_integral );

		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Function grad squared        : %.15f\n", results.function_grad_squared_integral );
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Function grad dot approx grad: %.15f\n", results.approx_grad_function_grad_integral );
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Approx grad squared          : %.15f\n\n", results.approx_grad_squared_integral );

		fun_val_integral.push_back    ( results.function_val_integral );
		approx_val_integral.push_back ( results.approx_val_integral );
		
		approx_val_squared.push_back ( results.approx_val_squared_integral );
		approx_val_fun_val.push_back ( results.approx_val_function_val_integral );
		fun_val_squared.push_back    ( results.function_val_squared_integral );
		
		approx_grad_squared.push_back  ( results.approx_grad_squared_integral );
		approx_grad_fun_grad.push_back ( results.approx_grad_function_grad_integral );
		fun_grad_squared.push_back     ( results.function_grad_squared_integral );
		
		L2_error.push_back (
		        results.approx_val_squared_integral
		    - 2*results.approx_val_function_val_integral
		    +   results.function_val_squared_integral );
		if ( !(L2_error.back() > 0) ) {
			PetscPrintf ( PETSC_COMM_WORLD, "L2 norm squared should not be negative: %f should be zero or minus zero.\n",
			              L2_error.back() );
		}
		L2_error.back() = L2_error.back()<0 ? -L2_error.back() : L2_error.back();
		
		W12_seminorm.push_back (
		        results.approx_grad_squared_integral
		    - 2*results.approx_grad_function_grad_integral
		    + results.function_grad_squared_integral );
		
		if ( !(W12_seminorm.back() > 0) ) {
			PetscPrintf ( PETSC_COMM_WORLD, "W12 semi-norm squared should not be negative: %f should be zero or minus zero.\n",
			              W12_seminorm.back() );
		}
		W12_seminorm.back() = W12_seminorm.back()<0 ? -W12_seminorm.back() : W12_seminorm.back();
		
		W12_error.push_back( L2_error.back() + W12_seminorm.back() );
#if 0	
		// Redundant for checking gnuplot output. 2009-06-24 ML.
		W12_error.back() = results.approx_val_squared_integral
		                      - 2*results.approx_val_function_val_integral
		                      +   results.function_val_squared_integral
		                      + results.approx_grad_squared_integral
		                      - 2*results.approx_grad_function_grad_integral
		                      + results.function_grad_squared_integral;
#endif
		

		
		// Got some nan - from -0.0 I believe.
		
		// Relative error norms.
		{
			double tmp_nonzero;
			
			tmp_nonzero = sqrt( results.function_val_squared_integral );
			L2_error.back()     = sqrt( L2_error.back() );
			rel_L2_error.push_back ( L2_error.back() );
			if ( tmp_nonzero>0 ) {
				rel_L2_error.back() /= tmp_nonzero;
			}
			
			tmp_nonzero = sqrt( results.function_grad_squared_integral );
			W12_seminorm.back() = sqrt( W12_seminorm.back() );
			rel_W12_seminorm.push_back ( W12_seminorm.back() );
			if ( tmp_nonzero>0 ) {
				rel_W12_seminorm.back() /= tmp_nonzero;
			}
		
			tmp_nonzero = sqrt( results.function_val_squared_integral + results.function_grad_squared_integral );
			W12_error.back() = sqrt( W12_error.back() );
			rel_W12_error.push_back ( W12_error.back() );
			if ( tmp_nonzero>0 ) {
				rel_W12_error.back() /= tmp_nonzero;
			}
		}
		

		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Rel. L2 error norm           : %.15f\n", rel_L2_error.back() );
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Rel W^{1,2} seminorm         : %.15f\n", rel_W12_seminorm.back() );
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\t * Rel. W^{1,2} error norm      : %.15f\n", rel_W12_error.back() );
		
		if ( iter%comm_size == comm_rank ) {
			//pum_simulation->set_sum_all_pu_test ( 3 );
			{
				ostringstream    name;
				name << "pum_convergence-1_solution_" << iter << ".vtp";
				ostringstream    error_name;
				error_name << "pum_convergence-1_L2_error_integrand_" << iter << ".vtp";
				
				pum_simulation->output_vtk( name.str(), error_name.str() );
			}
		}
		
		
		if ( 0 == comm_rank ) {
			errors_file << refinement_level.back()   << " "
				<< degrees_of_freedom.back() << " "
				<< rel_L2_error.back()           << " "
				<< rel_W12_seminorm.back()       << " "
				<< rel_W12_error.back()          << "\n";
			errors_file.flush();
		}
		
		PetscGetTime ( &single_iteration_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "\tIteration took %.2f minutes.\n", (single_iteration_t1-single_iteration_t0)/60.0 );

	}
	
	if ( ++round_robin_counter % comm_size == comm_rank ) {
		ofstream    file ( "pum_convergence-1_fun_val_integral.gnuplot" );
		file.precision(num_decimal_digits);
		
		file << "plot '-' using 1:3 w l\n";
		
		for ( int i=0; i<fun_val_integral.size(); ++i ) {
			file << refinement_level[i]   << " "
			     << degrees_of_freedom[i] << " "
			     << fun_val_integral[i]   << "\n";
		}
	}
	
	if ( ++round_robin_counter % comm_size == comm_rank ) {
		ofstream    file ( "pum_convergence-1_approx_val_integral.gnuplot" );
		file.precision(num_decimal_digits);
		
		file << "plot '-' using 1:3 w l\n";
		
		for ( int i=0; i<approx_val_integral.size(); ++i ) {
			file << refinement_level[i]   << " "
			     << degrees_of_freedom[i] << " "
			     << approx_val_integral[i]   << "\n";
		}
	}
	
	{
		assert( max_iterations == cover_factor_vec.size() );
		assert( max_iterations == quadrature_order_vec.size() );
		assert( max_iterations == refinement_level.size() );
		assert( max_iterations == degrees_of_freedom.size() );
		
		assert( max_iterations == fun_val_integral.size() );
		assert( max_iterations == approx_val_integral.size() );
		
		assert( max_iterations == approx_val_squared.size() );
		assert( max_iterations == approx_val_fun_val.size() );
		assert( max_iterations == fun_val_squared.size() );
		
		assert( max_iterations == approx_grad_squared.size() );
		assert( max_iterations == approx_grad_fun_grad.size() );
		assert( max_iterations == fun_grad_squared.size() );
		
		assert( max_iterations == L2_error.size() );
		assert( max_iterations == W12_seminorm.size() );
		assert( max_iterations == W12_error.size() );
		
		assert( max_iterations == rel_L2_error.size() );
		assert( max_iterations == rel_W12_seminorm.size() );
		assert( max_iterations == rel_W12_error.size() );

		pair<double,int>    rel_L2_error_p;
		pair<double,int>    rel_W12_seminorm_p;
		
		for ( int iter=0; iter<max_iterations; ++iter ) {
			rel_L2_error_p.first = frexp_dec( rel_L2_error[iter], &rel_L2_error_p.second );
			rel_W12_seminorm_p.first = frexp_dec ( rel_W12_seminorm[iter], &rel_W12_seminorm_p.second );
			
			if ( 0==iter ) {
				PetscPrintf ( PETSC_COMM_WORLD, "$%2.1f$ & $%i$ & $%i$ & $%i$ & $%i$ & $%4.3f_{%+i}$ & $%4.3f_{%+i}$ & & \\\\\n",
				              cover_factor_vec[iter],
				              quadrature_order_vec[iter],
				              local_basis_order,
				              refinement_level[iter],
				              degrees_of_freedom[iter],
				              rel_L2_error_p.first,
				              rel_L2_error_p.second,
				              rel_W12_seminorm_p.first,
				              rel_W12_seminorm_p.second );
				
			} else {
				//PetscPrintf ( PETSC_COMM_WORLD, "$%2.1f$ & $%i$ & $%i$ & $%i$ & $%i$ & $%4.3f_{%+i}$ & $%4.3f_{%+i}$ & $%+4.3f$ & $%+4.3f$ \\\\\n",
				//	cover_factor_vec[iter],
				//	quadrature_order_vec[iter],
				//	local_basis_order,
				PetscPrintf ( PETSC_COMM_WORLD, "      &     &     & $%i$ & $%i$ & $%4.3f_{%+i}$ & $%4.3f_{%+i}$ & $%+4.3f$ & $%+4.3f$ \\\\\n",
				              refinement_level[iter],
				              degrees_of_freedom[iter],
				              rel_L2_error_p.first,
				              rel_L2_error_p.second,
				              rel_W12_seminorm_p.first,
				              rel_W12_seminorm_p.second,
				              log10(L2_error[iter]/L2_error[iter-1])/log10(static_cast<double>(degrees_of_freedom[iter])/degrees_of_freedom[iter-1]),
				              log10(W12_seminorm[iter]/W12_seminorm[iter-1])/log10(static_cast<double>(degrees_of_freedom[iter])/degrees_of_freedom[iter-1]));
			}
		}
	}
	
	pum_simulation.reset(); // Destroys petsc objects before finalize.
	
	PetscFinalize();
}
