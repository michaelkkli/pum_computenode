/*
	Copyright (C) 2009 Michael Li
	This file is part of the Computenode Library.

	The Computenode Library is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "boundary_element_utils.hh"

#include <algorithm>
#include <valarray>
#include <cassert>
#include <cstdlib>

using std::abs;
using std::fabs;
using std::copy;
using std::valarray;

template<int dim>
void line_segment_get_point ( const double* line, double t, double* co )
{
	assert( line && co );
	
	for ( int d=0; d<dim; ++d ) {
		co[d] = 0.5*( (line[d] + line[dim + d]) + t*(line[dim+d] - line[d]) );
	}
}

/**
	Line segment is parametrized by parameter `t'.
		r = 0.5*(a+b) + 0.5*t*(b-a)
	The plane is given by
		r.n = p.n
	for `p' a point on the plane and `n' a normal to the plane.
	Substituting the first equation into the second and rearranging.
		t = (2p.n - (a+b).n)/(b-a).n = (2p -a -b).n/(b-a).n
*/
template<int dim>
bool intersect_line_with_plane ( const double* line,
				 const double* outward_normal,
				 const double* point_on_plane,
			         double* intersection_parameter )
{
	assert( line && outward_normal && point_on_plane );
	valarray<double> a( dim );
	valarray<double> b( dim );
	valarray<double> n( dim );
	valarray<double> p( dim );
	
	copy( line,       line + dim,   &a[0] );
	copy( line + dim, line + 2*dim, &b[0] );
	copy( outward_normal, outward_normal+dim, &n[0] );
	copy( point_on_plane, point_on_plane+dim, &p[0] );
	
	// (2p-a-b)
	valarray<double> tmp( p	);
	tmp *= 2.0;
	tmp -= a;
	tmp -= b;
	
	// (2p-a-b).n
	tmp *= n;
	double value = tmp.sum();
	
	// (b-a)
	tmp =  b;
	tmp -= a;
	
	// (2p-a-b).n/(b-a).n
	tmp *= n;
	double denom = tmp.sum();
	
	// Regard the line as parallel to the plane if
	// the angle between the line segment and the norm
	// is close to pi/2.
	// We take 2*pi/10000 to be a very small angle and if
	// the angle between (b-a) and `n' differs by less than
	// this small angle from pi/2, we regard the angle as pi/2.
	// cos( pi/2 + 2*pi/10000 ) = -0.0006283184893761124389
	if ( fabs( static_cast<long double>(denom) ) < 0.0001 ) {
		valarray<double> ang( tmp );
		ang *= tmp;
		double len = sqrt( ang.sum() );
		double normalized = denom / len;
		if ( fabs( static_cast<long double>(normalized) ) < 0.0006283184893761124389 ) {
			return false;
		}
	}
	
	value /= denom;
	
	if ( intersection_parameter ) {
		*intersection_parameter = value;
	}
	
	if ( fabs( static_cast<long double>(value) ) <= 1.0 ) {
		return true;
	} else {
		return false;
	} 
}

#if 0
void gp_draw_single( const boundary_element<1>& be, ostream& str, double level )
{
	const double* pt = be.get_point( 0 );
	str << pt[0] << " " << level << "\n";
}
#endif

void gp_draw_single( const boundary_element<2>& be, ostream& str, double level )
{
	const double* pt = be.get_point( 0 );
	str << pt[0] << " " << pt[1] << " " << level << "\n";
	pt = be.get_point( 1 );
	str << pt[0] << " " << pt[1] << " " << level << "\n";
	str << "\n"; // Leave line to end boundary element.
}

#if 0
void gp_draw_single( const boundary_element<3>& be, ostream& str, double ignored )
{
	const double* pt = be.get_point( 0 );
	str << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
	pt = be.get_point( 1 );
	str << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
	pt = be.get_point( 2 );
	str << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
	pt = be.get_point( 0 );
	str << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
	str << "\n"; // Leave line to end boundary element.
}
#endif

template<int dim>
void
gp_draw( const vector<boundary_element<dim> >& vec, ostream& str )
{
  for ( size_t i=0; i<vec.size(); ++i ) {
    gp_draw_single( vec[i], str );
    str << "\n";
  }
}

//

template bool
intersect_line_with_plane<1> ( const double* line, const double*, const double*, double* );
template bool
intersect_line_with_plane<2> ( const double* line, const double*, const double*, double* );
template bool
intersect_line_with_plane<3> ( const double* line, const double*, const double*, double* );

// template void gp_draw( const vector<boundary_element<1> >&, ostream& );
template void gp_draw( const vector<boundary_element<2> >&, ostream& );
// template void gp_draw( const vector<boundary_element<3> >&, ostream& );
