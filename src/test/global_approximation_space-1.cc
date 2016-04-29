#include "global_approximation_space.hh"
#include "dbinary_tree.hh"
#include "geometry.hh"

#include "polynomial.hh"

#include <petscts.h>

#include <fstream>
#include <iostream>
#include <valarray>

using std::cout;
using std::endl;
using std::ofstream;
using std::valarray;

#define dim 2

static char help[] = "global_approximation_space-1";

int main (int argc, char* argv[]) {
	
	PetscInitialize( &argc, &argv, static_cast<char*>(0), help );
	
	PetscReal parent_box_expansion_factor = 1.0;
	PetscOptionsGetReal(PETSC_NULL,"-parent_box_expansion_factor", &parent_box_expansion_factor, PETSC_NULL);
	if ( parent_box_expansion_factor<=0.0 ) {
		parent_box_expansion_factor = 1.0;
	}
	
	PetscInt dtree_initial_refinement = 4;
	PetscOptionsGetInt(PETSC_NULL,"-dtree_initial_refinement", &dtree_initial_refinement, PETSC_NULL);
	if ( dtree_initial_refinement<0 ) {
		dtree_initial_refinement = 4;
	}
	
	PetscReal dtree_cover_factor = 1.3;
	PetscOptionsGetReal(PETSC_NULL,"-dtree_cover_factor", &dtree_cover_factor, PETSC_NULL);
	if ( dtree_cover_factor<=0.0 ) {
		dtree_cover_factor = 1.3;
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

	valarray<double> vcoords( num_doubles );
	for ( int i=0; i<num_doubles; ++i ) {
		vcoords[i] = coords[i];
	}
	
	// Set up the geometry.
	geometry<dim>::ptr geom( new geometry<dim> );
	geom->create_boundary( "chloroplast envelope", vcoords );
	geom->create_region( "chloroplast", "+chloroplast envelope" );
	
	ofstream file_10( "global_approximation_space-1_chloroplast" );
	
	geom->gp_draw( "chloroplast", file_10 );
	file_10.close();
	
	geom->update();

	
	// Set up the d-binary tree used to form the
	// partition of unity method cover.
	dbinary_tree<dim>::ptr dtree( new dbinary_tree<dim> );
	dtree->set_geometry( geom, parent_box_expansion_factor );
	
	ofstream file_15( "global_approximation_space-1_parent_box" );
	
	dtree->gp_draw_parent_box( file_15 );
	file_15.close();
	
	dtree->set_initial_refinement( dtree_initial_refinement ); // 1D 2^n, 2D 4^n, 3D 8^n so refinement can be expensive.
	
	dtree->set_cover_factor( dtree_cover_factor ); // Possibly use "length extension" to be descriptive.
	
	// Create covers for the regions and boundaries.
	// TODO: allow the option not to cover unwanted areas.
	dtree->update();
	
	{
		ofstream file( "global_approximation_space-1_boundary_cover" );
		dtree->gp_draw_boundary_cover( "chloroplast envelope", file, 0.3 );
	}
	{
		ofstream file( "global_approximation_space-1_region_boundary_cover" );
		dtree->gp_draw_region_boundary_cover( "chloroplast", file, 0.1 );
	}
	{
		ofstream file( "global_approximation_space-1_region_interior_cover" );
		dtree->gp_draw_region_interior_cover( "chloroplast", file, 0.2 );
	}
	
	// dtree->dump_cout();
	
	// Set up the partition of unity approximation.
	// The cover should not need to be adjusted
	// after the above finalization step.
	
	global_approximation_space<dim>::ptr global_approx_space( new global_approximation_space<dim> );

	
	PetscErrorCode ierr = PetscFinalize(); CHKERRQ(ierr);
}
