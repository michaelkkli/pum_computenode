#include "vtk_output.hh"
#include "basic_dbinary_tree.hh"
#include "box.hh"
#include "cover_structure.hh"
#include "dbinary_tree.hh"
#include "dbinary_tree_utils.hh"
#include "geometry.hh"
//#include "global_approximation_space.hh"
//#include "global_basis_function.hh"
#include "integration_scheme.hh"
#include "petsc_solver.hh"
#include "polynomial.hh"
//#include "solution.hh"

#include <petscts.h>

#include <boost/shared_ptr.hpp>

#include <cassert>
#include <fstream>
#include <iostream>
#include <valarray>

using boost::shared_ptr;
using std::cout;
using std::endl;
using std::ofstream;
using std::slice;
using std::valarray;

class my_func : public function<3> {
public:
	my_func() {
		this->set_global_function();
	}
private:
	double evaluate( const double* x ) const {
		assert( x );
		return -x[0]*x[1];
	}
};

static char help[] = "vtk_output-2";

int main (int argc, char* argv[]) {
	
	PetscInitialize( &argc, &argv, static_cast<char*>(0), help );
	
	PetscReal parent_box_expansion_factor = 1.0;
	PetscOptionsGetReal(PETSC_NULL,"-parent_box_expansion_factor", &parent_box_expansion_factor, PETSC_NULL);
	if ( parent_box_expansion_factor<=0.0 ) {
		parent_box_expansion_factor = 1.0;
	}
	
	PetscInt dtree_initial_refinement = 2;
	PetscOptionsGetInt(PETSC_NULL,"-dtree_initial_refinement", &dtree_initial_refinement, PETSC_NULL);
	if ( dtree_initial_refinement<0 ) {
		dtree_initial_refinement = 2;
	}

	PetscInt dtree_max_refinement = dtree_initial_refinement;
	PetscOptionsGetInt(PETSC_NULL,"-dtree_max_refinement", &dtree_max_refinement, PETSC_NULL);
	if ( dtree_max_refinement<0 ) {
		dtree_max_refinement = dtree_initial_refinement;
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
	shared_ptr<geometry<2> > geom( new geometry<2> );
	geom->create_boundary( "chloroplast envelope", vcoords );
	geom->create_region( "chloroplast", "+chloroplast envelope" );
	
	ofstream file_10( "vtk_output-2_chloroplast" );
	
	geom->gp_draw( "chloroplast", file_10 );
	file_10.close();
	
	geom->update();
	
	
	// Set up the d-binary tree used to form the
	// partition of unity method cover.
	shared_ptr<dbinary_tree<2> > dtree( new dbinary_tree<2> );
	dtree->set_geometry( geom, parent_box_expansion_factor );
	
	ofstream file_15( "vtk_output-2_parent_box" );
	
	dtree->gp_draw_parent_box( file_15 );
	file_15.close();
	
	dtree->set_initial_refinement( dtree_initial_refinement ); // 1D 2^n, 2D 4^n, 3D 8^n so refinement can be expensive.
	
	dtree->set_cover_factor( dtree_cover_factor ); // Possibly use "length extension" to be descriptive.
	
	// Create covers for the regions and boundaries.
	// TODO: allow the option not to cover unwanted areas.
	dtree->update();
	
	{
		ofstream file( "vtk_output-2_boundary_cover" );
		dtree->gp_draw_boundary_cover( "chloroplast envelope", file, 0.3 );
		file.close();
	}
	{
		ofstream file( "vtk_output-2_region_boundary_cover" );
		dtree->gp_draw_region_boundary_cover( "chloroplast", file, 0.1 );
		file.close();
	}
	{
		ofstream file( "vtk_output-2_region_interior_cover" );
		dtree->gp_draw_region_interior_cover( "chloroplast", file, 0.2 );
		file.close();
	}
	
	//dtree->dump_cout();
	
	cover_structure<2>::ptr cs;
	assert( dtree );
	dtree->get_cover_structure( "chloroplast", cs );
	
	basic_dbinary_tree<2>::ptr basic_dtree ( new basic_dbinary_tree<2> );
	basic_dtree->initialize ( dtree->access_parent_box(), dtree_max_refinement );
	
	vtk_output::ptr vtk_st ( new vtk_output );
	
	basic_dtree->get_grid_points_3D ( vtk_st->points_3D );

#if 0
	std::clog << "Check points\n";
	for ( int i=0; i<vtk_st->points_3D.size(); ++i ) {
		std::clog << vtk_st->points_3D[i] << " ";
	}
#endif
	
	int num_scalars = vtk_st->points_3D.size()/3;

#if 0
	vtk_st->scalars.resize( num_scalars );
	
	for ( int i=0; i<num_scalars; ++i ) {
		vtk_st->scalars[i] = static_cast<double>(i);//num_scalars;
	}
#endif

	function<3>::ptr gen_scalars( new my_func );
	
	gen_scalars->global_evaluate ( box<3>(), vtk_st->points_3D, vtk_st->scalars );
	
	vtk_st->scalars /= 10;
	
	vtk_st->points_3D[ slice(2, vtk_st->scalars.size(), 3) ] = vtk_st->scalars;
	
	vector<string> all_keys_to_max_level;
	
	generate_keys<2> ( dtree_max_refinement, all_keys_to_max_level );
	
	basic_dtree->get_box_grid_connectivity_offsets ( all_keys_to_max_level, vtk_st->connectivity, vtk_st->offsets );
	
	{
		ofstream file ( "vtk_output-2_chloroplast_cover.vtp" );
		vtk_simple_output( *vtk_st, file );
	}
	
	PetscErrorCode ierr = PetscFinalize(); CHKERRQ(ierr);
	return 0;
}
