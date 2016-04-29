#include "cubic_spline_weight.hh"
#include "quadrature_rule.hh"

#include <fstream>
#include <iostream>
#include <valarray>
#include <vector>
using std::cout;
using std::ofstream;
using std::valarray;
using std::vector;

class constant : public function<1> {
private:
	double evaluate ( const double* x ) const
	{
		return x[0]*x[0]*x[0]*x[0];
	}
};


/**
	Convergence test for Gauss-Legendre quadrature applied to
	the cubic spline weight.
*/
int main ( int argc, char* argv[] ) {
	cout << "Begin\n";
	
	box<1>    bx;
	{
		double vals[] = {-1.0, 1.0};
		bx.set(vals);
	}

	cubic_spline_weight<1>    weight_fn;
	//constant    weight_fn;

	// Allow evaluation on subset of its support.
	weight_fn.set_global_function();

	{
		ofstream    file ( "quadrature_rule-2_quick_test_weight" );
		const int num = 100;
		valarray<double> local_x_vals(num);

		double    delta = 2.0/(num-1);
		
		for ( int i=0; i<num; ++i ) {
			local_x_vals[i] = -1.0 + i*delta;
		}

		valarray<double>    global_x_vals;
		bx.map_local_to_global ( local_x_vals, global_x_vals );
		
		valarray<double>    results;
		weight_fn.local_evaluate( bx, local_x_vals, results );


		
		file << "plot '-' w l\n";
		for ( int i=0; i<num; ++i ) {
			file << global_x_vals[i] << "\t" << results[i] << "\n";
		}
		file << "\n";
	}

	valarray<double>    tmp;

	quadrature_rule_type    qrt[5];
	qrt[0] = gauss1;    qrt[1] = gauss2;    qrt[2] = gauss3;
	qrt[3] = gauss4;    qrt[4] = gauss5;

	vector<double>    integrals;

	for ( int i=0; i<5; ++i ) {
		quadrature_rule<1>    quad_rule ( qrt[i] );

		weight_fn.local_evaluate( bx, quad_rule.access_quadrature_points(), tmp );

		cout << "Before mult weights. Sum is " << tmp.sum() << "\n";
		
		tmp *= quad_rule.access_quadrature_weights();

		cout << "After mult weights. Sum is " << tmp.sum() << "\n";

		// Multiple by bx.measure()*0.5 for adjustment factor.
		integrals.push_back( tmp.sum()*bx.measure()*0.5 );
	}

	double    correct_val = 0.5;
	{
		ofstream    file ( "quadrature_rule-2_error_graph.gnuplot" );

		file << "plot '-' w l";
		file << ", " << correct_val << " t \"Correct value.\"";
		file << "\n";
		for ( int i=0; i<5; ++i ) {
			file << i << "\t" << integrals[i] << "\n";
		}
		file << "\n";

		cout << "Hello\n";
	}
}
