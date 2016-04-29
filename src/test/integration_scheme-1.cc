#include "box.hh"
#include "global_basis_function.hh"
#include "integration_scheme.hh"
#include "cubic_spline_weight.hh"
#include "partition_of_unity_function.hh"
#include "polynomial.hh"


#include <algorithm>
#include <iostream>
using std::copy;
using std::cout;
using std::endl;



int main() {

	/*
	   A B
	   C D
	*/

	box<2/*dim*/>::vec_ptr boxes(4);

	for ( int i=0; i<4; ++i ) {
		boxes[i].reset( new box<2/*dim*/> );
	}

	double Apts[] = { -1.0,  0.0,  0.0, 1.0 };
	double Bpts[] = {  0.0,  1.0,  0.0, 1.0 };
	double Cpts[] = { -1.0,  0.0, -1.0, 0.0 };
	double Dpts[] = {  0.0,  1.0, -1.0, 0.0 };

	boxes[0]->set( Apts );
	boxes[1]->set( Bpts );
	boxes[2]->set( Cpts );
	boxes[3]->set( Dpts );

	for( int i=0; i<4; ++i ) {
		assert( !boxes[i]->empty() );
		boxes[i]->scale(1.2);
		assert( !boxes[i]->empty() );
	}

	differentiable_function<2/*dim*/>::ptr patch_weight( new cubic_spline_weight<2> );
	differentiable_function<2/*dim*/>::ptr local_approx( new polynomial<2>( 1.0, 3, 3 ) );

	global_basis_function<2/*dim*/> global_basis_fn;
	global_basis_fn.set( boxes, 0, patch_weight, local_approx );

	integration_scheme<2> int_sch;
	cout.precision(16);
	cout << "Default decomposition level of integration : "
	     << int_sch.integrate_global_basis_function( global_basis_fn ) << "\n"
	     << "Level one                                  : "
	     << int_sch.integrate_global_basis_function( global_basis_fn, 1 ) << "\n";

	cout << endl;
}
