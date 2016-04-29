

#include "function_utils.hh"

#include "point_utils.hh"

#include "box.hh"

#include <fstream>

#include <valarray>
using std::valarray;

#include <cassert>

class user_function : public function<1> {
public:
  user_function ()
  {
	  this->set_global_function();
  }
  ~user_function () { }
private:
  double evaluate ( const double* co ) const {
    return ( 1 - co[0] ) * ( 1 + co[0] );
  }
};

int main () {

  const double A = -2.0, B = 3.0;
  double ext[2] = { A, B };

  const int N = 10; // Number of points in subdivision.

  box<1> domain;
  domain.set( &ext[0] );

  user_function global_function;

  vector<valarray<double> > points_values(2);
  // 0 is domain_points
  // is userfunction values

  valarray<double>& domain_points = points_values[0];

  make_subdivisions ( A, B, N, domain_points );

  assert ( domain_points.size() == N );

  global_function.global_evaluate( domain, domain_points, points_values[1] );

  std::ofstream file1 ( "single_patch_approx-1_out.txt" );

  output_values ( points_values, file1 );

}
