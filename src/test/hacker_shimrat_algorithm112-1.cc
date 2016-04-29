#include "hacker_shimrat_algorithm112.hh"

#include <cassert>

int main ()
{
	{
		// Triangle.
		double n = 3;
		double triangle[] = { 0.0, 0.0,
				0.0, 1.0,
				1.0, 0.0
		};
		double x[] = { 0.0, 0.0, 1.0 };
		double y[] = { 0.0, 1.0, 0.0 };
		assert( !point_in_polygon( n, x, y, 1.0, 1.0 ) );
		assert( point_in_polygon( n, x, y, 0.001, 0.001 ) );
	}
	{
		double coords[] = { 0.0, -10.0,
		1.2, -7.0,
		5.4, 2.0,
		3.0, 7.8,
		-3.0, 9.5,
		-0.7, 0.5 };
		
		double x[] = {   0.0,  1.2, 5.4, 3.0, -3.0, -0.7 };
		double y[] = { -10.0, -7.0, 2.0, 7.8,  9.5,  0.5 };
		
		assert( point_in_polygon( 6, x, y, 0.0, 0.0 ) );
		assert( !point_in_polygon( 6, x, y, -3.0, -10.0 ) );
		assert( !point_in_polygon( 6, x, y, 3.0, -5.0 ) );
	}
}
