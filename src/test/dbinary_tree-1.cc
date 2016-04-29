#include "dbinary_tree.hh"

#include <fstream>
using std::ofstream;

#define DIM 2

int main() {

	// This test will fail as no bounding box is set.
	assert( false );

	dbinary_tree<DIM> first_tree;
	
	double ext[6] = { 10.0, 100.0, 5.0, 30.0, 0.0, 1.0 };
	box<DIM> parent;
	parent.set( &ext[0] );
	
	first_tree.set_parent_box( parent );
	
	first_tree.update();
	
	first_tree.set_initial_refinement( 3 );
	
	first_tree.dump_cout();
	
	ofstream fvec_box2( "dbinary_tree-1_dump_gp.dat" );
	
	first_tree.dump_gp( fvec_box2 );
	
	ofstream template_boxes_f( "dbinary_tree-1_dump_gp_template_boxes.dat" );
	
	first_tree.dump_gp_template_boxes( template_boxes_f );
}
