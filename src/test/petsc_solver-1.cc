#include "polynomial.hh"

#include "petsc_solver.hh"

#include <petscts.h>
#include <boost/shared_ptr.hpp>
#include <map>
using boost::shared_ptr;
using std::map;

static char help[] = "petsc_solver-1";

#define dim 2

int main ( int argc, char* argv[] ) {
	
	PetscInitialize( &argc, &argv, static_cast<char*>(0), help );
	
	int stages[3] = {0,1,2};
	
	PetscLogStageRegister(&stages[0],"dbinary_tree update");
	PetscLogStageRegister(&stages[1],"pum_discretization update");
	
	PetscReal parent_box_expansion_factor = 1.0;
	PetscOptionsGetReal(PETSC_NULL,"-parent_box_expansion_factor", &parent_box_expansion_factor, PETSC_NULL);
	if ( parent_box_expansion_factor<=0.0 ) {
		parent_box_expansion_factor = 1.0;
	}
	
	PetscInt dtree_initial_refinement = 4;
	PetscOptionsGetInt(PETSC_NULL,"-dtree_initial_refinement", &dtree_initial_refinement, PETSC_NULL);
	if ( dtree_initial_refinement<=0 ) {
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
	shared_ptr<geometry<dim> > geom( new geometry<dim> );
	geom->create_boundary( "chloroplast envelope", vcoords );
	geom->create_region( "chloroplast", "+chloroplast envelope" );
	geom->update();
	
	PetscPrintf(PETSC_COMM_WORLD, "parent_box_expansion_factor is %f\n", parent_box_expansion_factor);
	PetscPrintf(PETSC_COMM_WORLD, "dtree_initial_refinement is %i\n", dtree_initial_refinement);
	PetscPrintf(PETSC_COMM_WORLD, "dtree_cover_factor is %f\n", dtree_cover_factor);
	// Set up the d-binary tree used to form the
	// partition of unity method cover.
	shared_ptr<dbinary_tree<dim> > dtree( new dbinary_tree<dim> );
	dtree->set_geometry( geom, parent_box_expansion_factor );
	dtree->set_initial_refinement( dtree_initial_refinement ); // 1D 2^n, 2D 4^n, 3D 8^n so refinement can be expensive.
	
	dtree->set_cover_factor( dtree_cover_factor ); // Possibly use "length extension" to be descriptive.
	
	PetscLogStagePush(stages[0]);
	
	
	// Create covers for the regions and boundaries.
	// TODO: allow the option not to cover unwanted areas.
	dtree->update();
	
	PetscLogStagePop();
	
	shared_ptr<integration_scheme<dim> > integrator( new integration_scheme<dim> );
	
	shared_ptr<pum_discretization<dim> > sim( new pum_discretization<dim> );
	sim->set_dbinary_tree( dtree );
	sim->create_species( "gfp", "chloroplast" );
	
	// (xy)^2 in 2D; (xyz)^2 in 3D
	shared_ptr<function<dim> >
			initial_condition ( new polynomial<dim>(1.0, 2, 2 ) );
	initial_condition->set_global_function();
	assert( initial_condition->is_global_function() );
	
	#if 0 // TODO: remove
	sim->species_set_function( "gfp",
				   initial_condition
				 );
				 #endif
				 
	sim->set_integration_scheme( integrator );
	
	PetscLogStagePush(stages[1]);
	
	sim->update();
	
	PetscLogStagePop();
#if 0
	shared_ptr<petsc_solver<dim> > solver( new petsc_solver<dim> );
	solver->set_pum_discretization( sim );
	solver->update();
	
	MPI_Barrier( PETSC_COMM_WORLD );
	
	solver->assemble_for_initial_projection();

	
	MPI_Barrier( PETSC_COMM_WORLD );

	solver->solve_initial_projection();
	
	
	solver.reset(); // Destroy PETSc objects.
#endif
	
	PetscErrorCode ierr = PetscFinalize(); CHKERRQ(ierr);
}
