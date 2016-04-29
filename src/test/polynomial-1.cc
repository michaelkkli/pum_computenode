#include "polynomial.hh"

#include <iostream>
using std::cout;

int main(){
  polynomial<1> p1( 3.0, 2 );
  polynomial<2> p2( 3.0, 1, 1 );
  polynomial<3> p3( 3.0, 1, 1, 1 );
}
