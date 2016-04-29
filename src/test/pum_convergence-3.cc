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

// Four species test.

#include "extract_geometry_2d.hh"
#include "global_approximation_space.hh"
#include "pum_convergence.hh"
#include "refinement_structure.hh"
#include "singleimage.hh"
#include "vtk_output.hh"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
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
		return 100.0;
		
		// Unit box Helmholtz. (1)
		//return (2.0*co[0]-3.0)*co[0]*co[0];
	}
	void evaluate_grad ( const double * co, double *   grad ) const {
		// L2 projection
		grad[0]=0.0; grad[1]=0.0;
	
		// Unit box Helmholtz. (1)
		// grad[0] = 6.0*(co[0]-1.0)*co[0]; grad[1] = 0.0;
	}
};

class rhs_fun : public function<2> {
public:
	rhs_fun ( ) { this->set_global_function(); }
private:
	double evaluate    ( const double * co ) const {
		// L2 projection
		return 100.0;

		// Unit box Helmholtz. (1)
		//return co[0]*(co[0]*(2.0*co[0]-3.0) - 12.0)+6.0;
	
		// Poisson
		//return -6.0*co[0];
	}
};

#ifndef itoa

using std::stringstream;

const char* itoa ( int in ) {
	ostringstream  ss ("");
	ss << static_cast<int>(in);
	return (ss.str()).c_str();
}
#endif

static const char* help = "pum_convergence-3";

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
	              "====== Beginning execution of pum_convergence-3 ======\n" );
	
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
	
	PetscTruth    image1_siggia = PETSC_FALSE;
	if ( "" == image1_string ) {
		image1_siggia = PETSC_FALSE;
	}
	PetscOptionsGetTruth ( PETSC_NULL, "-image1_siggia", &image1_siggia, PETSC_NULL );
	
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
	
	PetscInt time_step_subdivision = 1;
	PetscOptionsGetInt ( PETSC_NULL, "-time_step_subdivision", &time_step_subdivision, PETSC_NULL );
	if ( time_step_subdivision < 1 ) {
		time_step_subdivision = 1;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-time_step_subdivision is %i.\n", time_step_subdivision );


	PetscReal global_bleach1 = 0.005;
	
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
	
	
	PetscInt num_species = 4;
	PetscOptionsGetInt ( PETSC_NULL, "-num_species", &num_species, PETSC_NULL );
	if ( num_species < 0 ) {
		time_step_subdivision = 0;
	}
	PetscPrintf ( PETSC_COMM_WORLD, "\t-num_species is %i.\n", num_species );
	
	species_system_data    species_sys_dat;
	
	vector<PetscReal>    initial_proportions ( num_species, 1.0/num_species );
	if ( num_species>=4 ) {
		initial_proportions[0] = 0.005;
		initial_proportions[1] = 0.005;
		initial_proportions[2] = 0.495;
		initial_proportions[3] = 0.495;
	}
	
	// Legacy support for scalar equation. 2009-11-08 ML.
	PetscReal diffusion1 = 33.0;
	PetscReal bleach_param1 = 100.0;
	
	// PetscReal is assumed to be double for the addresses/pointers
	// below to be correct. 2009-11-08 ML.
	
	if ( num_species > 0 ) {
		valarray<PetscReal>&    diffusion = species_sys_dat.diffusion_coefficient;
		valarray<PetscReal>&    bleach_param = species_sys_dat.spot_bleach_coefficient;
		valarray<PetscReal>&    global_bleach = species_sys_dat.global_bleach_coefficient;
		diffusion.resize    ( num_species, 33.0 );
		bleach_param.resize ( num_species, 0.0 );
		global_bleach.resize ( num_species, 0.0 );
		
		species_sys_dat.spot_bleach_coefficient.resize       ( num_species, 0.0 );
		species_sys_dat.reaction_coefficients.resize         ( num_species      );
		for ( int i=0; i<num_species; ++i ) {
			species_sys_dat.reaction_coefficients[i].resize ( num_species, 0.0 );
		}
		
		string    diffusion_name;
		for ( int i=0; i<num_species; ++i ) {
			diffusion_name = "-diffusion";
			diffusion_name.append ( itoa(i) );
			PetscOptionsGetReal ( PETSC_NULL, diffusion_name.c_str(), &diffusion[i], PETSC_NULL );
			PetscPrintf ( PETSC_COMM_WORLD, "\t%s is %f.\n", diffusion_name.c_str(), diffusion[i] );
		}
		
		string    bleach_param_name;
		for ( int i=0; i<num_species; ++i ) {
			bleach_param_name = "-bleach_param";
			bleach_param_name.append ( i ? itoa(i) : "0" );
			PetscOptionsGetReal ( PETSC_NULL, bleach_param_name.c_str(), &bleach_param[i], PETSC_NULL );
			PetscPrintf ( PETSC_COMM_WORLD, "\t%s is %f.\n", bleach_param_name.c_str(), bleach_param[i] );
		}
		
		string    global_bleach_name;
		for ( int i=0; i<num_species; ++i ) {
			global_bleach_name = "-global_bleach";
			global_bleach_name.append ( i ? itoa(i) : "0" );
			PetscOptionsGetReal ( PETSC_NULL, global_bleach_name.c_str(), &global_bleach[i], PETSC_NULL );
			PetscPrintf ( PETSC_COMM_WORLD, "\t%s is %f.\n", global_bleach_name.c_str(), global_bleach[i] );
		}
		
		diffusion1     = diffusion[1];
		bleach_param1  = bleach_param[1];
		global_bleach1 = global_bleach[1];
		
		vector<valarray<double> >& reaction_coefficients = species_sys_dat.reaction_coefficients;

		if ( num_species>=4 ) {
			reaction_coefficients[0][1] = 0.5;
			reaction_coefficients[1][0] = 0.5;
			reaction_coefficients[2][3] = 0.5;
			reaction_coefficients[3][2] = 0.5;
			reaction_coefficients[0][2] = 0.0;      // Mobile from focal plane out.
			reaction_coefficients[2][0] = 0.0; // Into focal plane.
		}
		
		string    initial_proportion_name;
		for ( int i=0; i<num_species; ++i ) {
			initial_proportion_name = "-initial_proportion";
			initial_proportion_name.append ( itoa(i) );
			PetscOptionsGetReal ( PETSC_NULL, initial_proportion_name.c_str(), &initial_proportions[i], PETSC_NULL );
			PetscPrintf ( PETSC_COMM_WORLD, "\t%s is %f.\n", initial_proportion_name.c_str(), initial_proportions[i] );
		}
		
		string    reaction_name;
		for ( int i=0; i<num_species; ++i ) {
			for ( int j=0; j<num_species; ++j ) {
				if ( i==j ) continue;
				
				reaction_name = "-a";
				if ( 0==i ) {
					reaction_name.append ( "0" );
				} else {
					reaction_name.append ( itoa(i) );
				}
				reaction_name.append ( "t" );
				if ( 1==j ) {
					reaction_name.append ( "1" );
				} else {
					reaction_name.append ( itoa(j) );
				}
				
				string tmp_str = reaction_name;
				
				PetscOptionsGetReal ( PETSC_NULL, tmp_str.c_str(), &(reaction_coefficients[i][j]), PETSC_NULL );
				PetscPrintf ( PETSC_COMM_WORLD, "\t%s is %f.\n", tmp_str.c_str(), reaction_coefficients[i][j] );
			}
		}
		
		species_sys_dat.coupling_param = 0.5;
		PetscOptionsGetReal ( PETSC_NULL, "-ac", &(species_sys_dat.coupling_param), PETSC_NULL );
		
		species_sys_dat.focal_percentage = 1.0;
		PetscOptionsGetReal ( PETSC_NULL, "-focal_percentage", &(species_sys_dat.focal_percentage), PETSC_NULL );
#if 0 // Removing volume adjust.
		vector<valarray<double> >& volume_adjust = species_sys_dat.volume_adjust;
		volume_adjust.resize(num_species);
		// volume_adjust[i][j] must have i>j.
		for ( int i=1; i<num_species; ++i ) {
			volume_adjust[i].resize ( i, 1.0 );
		}
		volume_adjust[2][0] = 0.01/0.99;
	
		string    volume_adjust_name;
		for ( int i=1; i<num_species; ++i ) {
			for ( int j=0; j<i; ++j ) {
				volume_adjust_name = "-c";
				volume_adjust_name.append ( itoa(i) );
				volume_adjust_name.append ( "f" ); // "from"
				volume_adjust_name.append ( itoa(j) );
				
				string tmp_str = volume_adjust_name;
				
				PetscOptionsGetReal ( PETSC_NULL, tmp_str.c_str(), &(volume_adjust[i][j]), PETSC_NULL );
				PetscPrintf ( PETSC_COMM_WORLD, "\t%s is %f.\n", tmp_str.c_str(), volume_adjust[i][j] );
			}
		}
#endif
		
	} else {
		PetscOptionsGetReal ( PETSC_NULL, "-diffusion1", &diffusion1, PETSC_NULL );
		PetscPrintf ( PETSC_COMM_WORLD, "\t-diffusion1 is %f.\n", diffusion1 );
		PetscOptionsGetReal ( PETSC_NULL, "-bleach_param1", &bleach_param1, PETSC_NULL );
		PetscPrintf ( PETSC_COMM_WORLD, "\t-bleach_param1 is %f.\n", bleach_param1 );
		PetscOptionsGetReal ( PETSC_NULL, "-global_bleach1", &global_bleach1, PETSC_NULL );
		PetscPrintf ( PETSC_COMM_WORLD, "\t-global_bleach1 is %f.\n", global_bleach1 );
	}
		
	
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
			ofstream    file ( "pum_convergence-3_labelled_bdry.gnuplot" );
			gnuplot_output<2> ( domain_segs, file, (num_segs<1000) );
		}
	}
	
	box<2>    parent_box;
	get_bounding_box ( domain_segs, parent_box );
	
	// Set the aspect ratio to be 1.0.
	// 2009-08-10 ML.
	regularize_bounding_box( parent_box );
	
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
	
	
	
	// Use image1 for projection, if given, otherwise
	// us the constant function. 2009-07-26 ML.
	if ( image1_string != "" ) {
		singleimage::ptr    image_initial;
		image_initial.reset ( new singleimage );
		image_initial->load ( image1_string, 1 );
		
		image_initial->rescale_percentage ();
		
		details->rhs_ptr = image_initial;
		details->fun_ptr = image_initial;
		
		if ( image1_siggia ) {
			details->equilibrium_fun = image_initial;
		}
	} else {
		details->fun_ptr.reset ( new simple_fun );
		details->rhs_ptr.reset ( new rhs_fun );
	}
	
	details->integration_decomposition = integration_decomposition;
	details->anticlockwise_dom = anticlockwise_dom;
	

	details->bleach_indicator_fun.reset ( new circular_indicator_function ( bleach_x_centre1, bleach_y_centre1, bleach_radius1 ) );
	assert ( details->bleach_indicator_fun );
	
	pum_convergence<2>::ptr    pum_simulation ( new pum_convergence<2> );
	convergence_results    results;

	const int max_iterations = 72;

	vector<double>    fun_val_integral;
	vector<double>    approx_val_integral;
	

	
	if ( ++round_robin_counter % comm_size == comm_rank ) {
		int num_to_do = max_iterations;
		pvd_output pvd;
		pvd.head = "pum_convergence-3_solution_";
		pvd.middle.resize ( num_to_do );
		ostringstream oss;
		for ( int i=0; i<num_to_do; ++i ) {
			oss.str ( "" );
			oss << i;
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "pum_convergence-3_solution.pvd" );
		pvd_simple_output ( pvd, file );
	}
	
	if ( ++round_robin_counter % comm_size == comm_rank ) {
		int num_to_do = max_iterations;
		pvd_output pvd;
		pvd.head = "pum_convergence-3_solution_";
		ostringstream oss;
		
		pvd.middle.resize ( 5 );
		int numbers[] = { 0, 2, 51, 52, 71 };
		
		for ( int i=0; i<5; ++i ) {
			oss.str ( "" );
			oss << numbers[i];
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "pum_convergence-3_selection.pvd" );
		pvd_simple_output ( pvd, file );
	}
	
	PetscLogDouble single_iteration_t0, single_iteration_t1;
	
	
	shared_ptr<diagnostic_quad_details> diag_quad_d ( new diagnostic_quad_details );
	
	pum_return_information    pum_return_info;
	
	{
		ref_struct->levels["0"] = initial_refinement;
		
		bool diag_output = false;
		
		PetscPrintf ( PETSC_COMM_WORLD, "Beginning initialize.\n" ); PetscGetTime ( &log_t0 );
		if ( ref_struct->levels["0"] < 7 ) {
			PetscPrintf ( PETSC_COMM_WORLD, "Integration points will be collected and the cover will be drawn for this iteration.\n" );
			pum_simulation->initialize ( details, &pum_return_info, &(*diag_quad_d), true, image1_siggia, num_species );
			
			ostringstream    name;
			name << "pum_convergence-3_cover.gnuplot";
			
			pum_simulation->draw_cover ( (name.str()).c_str(), true );
			diag_output = true;
		} else {
			pum_simulation->initialize ( details, &pum_return_info, 0, true, image1_siggia, num_species );
			diag_output = false;
		}
		PetscGetTime ( &log_t1 ); PetscPrintf ( PETSC_COMM_WORLD, "initialize took %.2f minutes.\n", (log_t1-log_t0)/60.0 );
		
		PetscPrintf ( PETSC_COMM_WORLD,
		              "\tThere are %i degrees of freedom across %i MPI processes.\n",
		              pum_return_info.num_dof,
		              comm_size );
		
		if ( diag_quad_d && diag_output ) {
			ostringstream    name;
			name << "pum_convergence-3_diag_quad_d.gnuplot";
			
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
		
		
	}
	
	double    full_step    = 0.764;
	double    time_step_dt = full_step / time_step_subdivision;
	
	species_sys_dat.time_step_dt = time_step_dt;
	
	ofstream    species_data_file;
	vector<double>    times;
	vector<double>    clflip_integral;
	
	if ( 0==comm_rank ) {
		species_data_file.open( "pum_convergence-3_species_data.gnuplot" );
		species_data_file.precision(16);
		
		{
			ofstream   file ( "pum_convergence-3_species.gnuplot" );
			file << "plot [*:*] 'pum_convergence-3_species_data.gnuplot' u 1:2 w l, '' u 1:3 w l, '' u 1:4 w l, '' u 1:5 w l\n";
		}
		
		{
			ofstream   file ( "pum_convergence-3_species01.gnuplot" );
			file << "plot [*:*] 'pum_convergence-3_species_data.gnuplot' u 1:2 w l, '' u 1:3 w l\n";
		}
		{
			ofstream   file ( "pum_convergence-3_species23.gnuplot" );
			file << "plot [*:*] 'pum_convergence-3_species_data.gnuplot' u 1:4 w l, '' u 1:5 w l\n";
		}
		
		{
			ofstream   file ( "pum_convergence-3_species0p1.gnuplot" );
			file << "plot [*:*] 'pum_convergence-3_species_data.gnuplot' u 1:($2+$3) w l\n";
		}
	}
	
	
	double    initial_integral;
	
	valarray<bool>    species_mask (4);
	species_mask[0] = true;
	species_mask[1] = true;
	species_mask[2] = false;
	species_mask[3] = false;
	
	vector<double>    species_integrals ( max_iterations*num_species, 0.0 );
	valarray<double>  tmp_species_integrals(num_species);
	
	// 2 prebleach, 50 bleach, 20 post bleach
	for ( int iter=0; iter<max_iterations; ++iter ) {
		PetscGetTime ( &single_iteration_t0 );

		if ( 0==iter ) {
			pum_simulation->species_set_initial_proportions ( &initial_proportions[0] );
		} else if ( 1==iter ) {
			// Bleach off.
			pum_simulation->species_set_matrices ( species_sys_dat, false );
		} else if ( 2==iter ) {
			// Bleach on.
			pum_simulation->species_set_matrices ( species_sys_dat, true, &species_mask );
		} else if ( 52==iter ) {
			// Bleach off.
			pum_simulation->species_set_matrices ( species_sys_dat, false, &species_mask );
		}
		
		if ( 0==iter ) {
			// Solved carried out in initial proportion setting.
		} else {
			for ( int j=0; j<time_step_subdivision; ++j ) {
				pum_simulation->species_solve ( species_sys_dat );
			}
		}
		
		pum_simulation->species_integrate_over_domain ( tmp_species_integrals );
		
		for ( int i=0; i<num_species; ++i ) {
			species_integrals[iter*num_species + i] = tmp_species_integrals[i];
		}
		
		times.push_back ( iter*full_step );
		
		if ( 0==comm_rank ) {
			species_data_file << times.back() << "\t";
			for ( int i=0; i<num_species; ++i ) {
				if ( 0!=i ) species_data_file << "\t";
				
				species_data_file << species_integrals[iter*num_species + i];
			}
			species_data_file << "\n";
		}
		
#if 0 // Nothing for species yet.
		if ( 0==iter ) {
			initial_integral = pum_simulation->integrate_solution_over_domain();
			clflip_integral.push_back( 100.0 );
		} else {
			clflip_integral.push_back( 100.0*pum_simulation->integrate_solution_over_domain()/initial_integral );
		}
		
		if ( 0==comm_rank ) {
			species_data_file << times.back() << "\t" << clflip_integral.back() << "\n";
		}
#endif
		
		if ( iter%comm_size == comm_rank ) {
			
			if ( true ) {
				ostringstream    name;
				name << "pum_convergence-3_solution_" << iter << ".vtp";
				if ( ref_struct->levels["0"] < 7 ) {
					// Second argument is to deactivate error output.
					// Third is to give cell data.
					pum_simulation->species_sum_output_vtk( name.str(), "", 2 );
				} else {
					pum_simulation->species_sum_output_vtk( name.str(), "", 2  );
				}
			}
		}
		

		PetscGetTime ( &single_iteration_t1 );
		PetscPrintf ( PETSC_COMM_WORLD, "\tIteration took %.2f s.\n", (single_iteration_t1-single_iteration_t0) );

	}

	if ( 0==comm_rank ) {
		ofstream   file ( "pum_convergence-3_species0p1.gnuplot" );
		file << "plot [*:*] [50:100] 'pum_convergence-3_species_data.gnuplot' u 1:(100.0*($2+$3)/"
			<< species_integrals[0]+species_integrals[1] << ") w l\n";
	}
#if 0 // Keep for integral info.
	if ( ++round_robin_counter % comm_size == comm_rank ) {
		ofstream    file ( "pum_convergence-3_fun_val_integral.gnuplot" );
		file.precision(num_decimal_digits);
		
		file << "plot '-' using 1:3 w l\n";
		
		for ( int i=0; i<fun_val_integral.size(); ++i ) {
			file << refinement_level[i]   << " "
			     << degrees_of_freedom[i] << " "
			     << fun_val_integral[i]   << "\n";
		}
	}
	
	if ( ++round_robin_counter % comm_size == comm_rank ) {
		ofstream    file ( "pum_convergence-3_approx_val_integral.gnuplot" );
		file.precision(num_decimal_digits);
		
		file << "plot '-' using 1:3 w l\n";
		
		for ( int i=0; i<approx_val_integral.size(); ++i ) {
			file << refinement_level[i]   << " "
			     << degrees_of_freedom[i] << " "
			     << approx_val_integral[i]   << "\n";
		}
	}
#endif
	
	pum_simulation.reset(); // Destroys petsc objects before finalize.
	
	PetscFinalize();
}
