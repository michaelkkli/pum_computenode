//#include <petsc.h>

#include <valarray>
using std::valarray;
#include <cmath>
using std::rand;
#include <fstream>
using std::ofstream;
#include <iostream>
using std::cout;
using std::endl;

int main( int argc, char* argv[] ) {
  const int dim = 2;
  const int num = 400;
  valarray<double> pts( dim*num );
  for ( int i=0; i<dim*num; ++i ) {
    pts[i] = static_cast<double>( rand() )/RAND_MAX;
  }
  ofstream gnuplot0 ( "points_data" );
  for ( int p=0; p<num; ++p ) {
    for ( int d=0; d<dim; ++d ) {
      if ( d == 0 ) {
	gnuplot0 << pts[p*dim+d];
      } else {
	gnuplot0 << " " << pts[p*dim+d];
      }
    }
    gnuplot0 << "\n";
  }
  gnuplot0.close();
}
