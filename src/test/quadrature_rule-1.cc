#include "box.hh"
#include "box_utils.hh"
#include "function.hh"
#include "line_segments.hh"
#include "point_utils.hh"
#include "polynomial.hh"
#include "quadrature_rule.hh"

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <valarray>

using std::cout;
using std::fstream;
using std::memcpy;
using std::ofstream;
using std::slice;
using std::stringstream;

class integrand : public function<2> {
public:
	integrand () {
		this->set_global_function ();
	}
private:
	double evaluate ( const double* co ) const {
		// return 1.0;
		// return co[0];
		return co[0]*co[1];
	}
};

int main ( int argc, char* argv[] ) {
	box<2>    box1;
	double box_ext[] = { 0.0, 1.0, 0.0, 1.0 };
	box1.set ( box_ext );


	line_segments    segs;
	segs.segments.resize ( 4 );

	int num_segs = (segs.segments.size()/2) - 1;


	double seg_endpts[] = { 0.0, -1.0, 1.0, 2.0 };
	memcpy ( &segs.segments[0], seg_endpts, 4*sizeof(double) );


	int num_quadpts = 4;

	if ( argc == 2 ) {
		stringstream ss;
		ss << argv[1];
		ss >> num_quadpts;
		if ( !ss ) {
			cout << "Illegal argument given!\n";
			abort();
		}
	}

	simple_quadrature_rule<1>    base_rule;
	generate_gauss_legendre_rule ( num_quadpts, base_rule );

	simple_quadrature_rule<2>    out_rule;
	out_rule.points.resize  ( 2*num_quadpts*num_quadpts, 3.45 );
	out_rule.weights.resize ( num_quadpts*num_quadpts, 4.56   );

	int num_sections = num_quadpts;

	// double final_quad_factor = box1.measure()/4;
	// final_quad_factor *= (box_ext[3]-box_ext[2])/2;
	double final_quad_factor = (box_ext[3]-box_ext[2])/2;
	
	for ( int i=0; i<num_sections; ++i ) {
		
		out_rule.points[ slice(2*num_quadpts*i + 1, num_quadpts, 2 ) ] = base_rule.points[i];

		valarray<double>    local_horiz_line (4);
		local_horiz_line[0] = -1.0;
		local_horiz_line[1] = base_rule.points[i];
		local_horiz_line[2] = 1.0;
		local_horiz_line[3] = base_rule.points[i];

		valarray<double>    global_horiz_line (4);

		box1.map_local_to_global ( local_horiz_line, global_horiz_line );

		double param1 = -3.0, param2 = -3.0;
		double pt[2];

		bool cross = intersect_lines_2d ( &global_horiz_line[0], &segs.segments[0], param1, param2, pt );

		// Notice no indexing on second weights only.
		out_rule.weights[ slice(num_quadpts*i, num_quadpts, 1) ] = base_rule.weights*base_rule.weights[i]*(pt[0] - box_ext[0])/2;

		if ( out_rule.weights.max() < 0 ) {
			cout << "Weights not strictly positive.\n";
			abort();
		}
		
		if ( cross ) {
			//cout << "They cross.\n";
		} else {
			cout << "Don't cross.\n";
			abort();
		}
		
		//cout << "param1 is " << param1 << "\n";
		
		box<1>    local_compress_box;
		{
			double compress[2] = { -1.0, param1 };
			local_compress_box.set ( compress );
		}

		valarray<double>    compressed_section;
		local_compress_box.map_local_to_global ( base_rule.points, compressed_section );

		out_rule.points[ slice(2*num_quadpts*i, num_quadpts, 2) ] = compressed_section;
	}



	valarray<double>    global_outpoints;
	box1.map_local_to_global ( out_rule.points, global_outpoints );

	{
		// Integrate on part of standard box.
		valarray<double>    temp_val ( out_rule.weights.size() );

		integrand    my_fun;
		my_fun.global_evaluate ( box1, global_outpoints, temp_val );

		assert ( temp_val.size () == out_rule.weights.size() );
		temp_val *= out_rule.weights;

		cout << "Max weights is " << out_rule.weights.max() << "\n";

		cout.precision ( 15 );
		cout << "Result of integral is " << temp_val.sum()*final_quad_factor << "\n";
	}
	
	{
		ofstream    file ( "quadrature_rule-1_draw.gnuplot" );
		file << "#!/usr/bin/gnuplot\n";
		file << "set size square\n";
		file << "plot [-1:2] [-1:2] '-' w l, '-' w l, '-' w dots\n";
		gp_draw_single ( box1, file );
		file << "e\n";
		gp_draw<2> ( segs, file );
		file << "e\n";
		gp_draw_points<2> ( global_outpoints, file );
		file << "e\n";
	}



#if 0 // Changing box constructor.
  cout << "Creating 1D quadrature rule of type gauss3\n";
  quadrature_rule<1> quadrule1( gauss3 );

  cout << "Creating zero polynomial\n";
  polynomial<1> poly1( 0.0 );

  cout << "Integral on (-1,1) is "
       << quadrule1.integrate( box<1>(-1,1), poly1 )
       << "\n";

  cout << "Integral on (-100,100) is "
       << quadrule1.integrate( box<1>(-100,100), poly1 )
       << "\n";

  cout << "Creating polynomial constant 1.0\n";
  polynomial<1> poly2( 1.0 );

  cout << "Integral on (-1,1) is "
       << quadrule1.integrate( box<1>(-1,1), poly2 )
       << "\n";

  cout << "Integral on (-100,100) is "
       << quadrule1.integrate( box<1>(-100,100), poly2 )
       << "\n";

  cout << "Creating polynomial x\n";
  polynomial<1> poly3( 1.0, 1 );

  cout << "Integral on (0,1) is "
       << quadrule1.integrate( box<1>(0,1), poly3 )
       << "\n";

  cout << "Integral on (-1,1) is "
       << quadrule1.integrate( box<1>(-1,1), poly3 )
       << "\n";
#endif
}

