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

#include "point_utils.hh"

void make_points_1d ( double a,
		      double b,
		      int n,
		      valarray<double>& out,
		      bool on_boundary=true )
{
  assert( a < b );
  assert( n > 0 );
  out.resize( n );
  if ( on_boundary ) {
    double d = (b-a)/(n-1);
    for ( int i=0; i<n; ++i ) {
      out[i] = a + d*i;
    }
  } else {
    double d = (b-a)/(n+1);
    for ( int i=0; i<n; ++i ) {
      out[i] = a + d*(i+1);
    }
  }
}

void make_subdivisions ( double a, double b, int n, valarray<double>& out )
{
  make_points_1d ( a, b, n, out, true );
}

void make_point_vec ( const valarray<double>& x_vec,
		      const valarray<double>& y_vec,
		      valarray<double>& points2D )
{
  int xsize = x_vec.size();
  int ysize = y_vec.size();
  if ( xsize == 0 || ysize == 0 || xsize != ysize ) {
    return;
  }
  points2D.resize( 2*xsize*ysize );

  for ( int y=0; y<ysize; ++y ){
    for ( int x=0; x<xsize; ++x ){

      // dim == 2
      int pos = 2*xsize*y + 2*x;

      points2D[ pos    ] = x_vec[x];
      points2D[ pos +1 ] = y_vec[y];

      //     std::cout << "(x,y) "<<x <<", " << y << "pos "<< pos << "\n";

    }
  }
}

void make_point_vec ( const valarray<double>& x_vec,
		      const valarray<double>& y_vec,
		      const valarray<double>& z_vec,
		      valarray<double>& points3D )
{
  int xsize = x_vec.size();
  int ysize = y_vec.size();
  int zsize = z_vec.size();

  if ( xsize == 0 || ysize == 0 || zsize == 0 ) {
    return;
  }

  points3D.resize( 3*xsize*ysize*zsize );

  for ( int z=0; z<zsize; ++z ) {
    for ( int y=0; y<ysize; ++y ){
      for ( int x=0; x<xsize; ++x ){

	// dim == 3
	int pos = 3*xsize*ysize*z + 3*xsize*y + 3*x;

	points3D[ pos    ] = x_vec[x];
	points3D[ pos +1 ] = y_vec[y];
	points3D[ pos +2 ] = z_vec[z];
      }
    }
  }
}

template<int dim>
void make_point_vec_on_box ( const box<dim>& bx,
			     valarray<double>& out,
			     int num_x,
			     int num_y,
			     int num_z )
{
  if ( bx.empty() ) {
    return;
  }
  const double* vals = bx.get();
  if ( dim == 1 ) {
    make_subdivisions( vals[0], vals[1], num_x, out );
    return;
  } else if ( dim == 2 ) {

    valarray<double> vecx;
    valarray<double> vecy;

    make_subdivisions( vals[0], vals[1], num_x, vecx );
    make_subdivisions( vals[2], vals[3], num_y, vecy );

    make_point_vec( vecx, vecy, out );

    return;

  } else if ( dim == 3 ) {

    valarray<double> vecx;
    valarray<double> vecy;
    valarray<double> vecz;

    make_subdivisions( vals[0], vals[1], num_x, vecx );
    make_subdivisions( vals[2], vals[3], num_y, vecy );
    make_subdivisions( vals[4], vals[5], num_z, vecz );

    make_point_vec( vecx, vecy, vecz, out );

    return;

  } else {
    out.resize(0);
    return;
  }
}

template<int dim>
void make_point_vec_interior_box ( const box<dim>& bx,
				   valarray<double>& out,
				   int num_x,
				   int num_y,
				   int num_z )
{
  if ( bx.empty() || num_x<=0 || num_y<=0 || num_z<=0 ) {
    out.resize(0);
    return;
  }
  const double* vals = bx.get();
  if ( dim == 1 ) {

    double sep = (vals[1]-vals[0])/(num_x+1);
    
    double ext[2];
    ext[0] = vals[0]+sep;
    ext[1] = vals[1]-sep;

    box<1> tbox;
    tbox.set( &ext[0] );

    make_point_vec_on_box ( tbox, out, num_x );

    return;

  } else if ( dim == 2 ) {

    std::cout << "In\n";

    double sepx = (vals[1]-vals[0])/(num_x+1);
    double sepy = (vals[3]-vals[2])/(num_y+1);
    
    double ext[4] = { vals[0]+sepx, vals[1]-sepx, vals[2]+sepy, vals[3]-sepy  };

    box<2> tbox;
    tbox.set( &ext[0] );

    std::cout << "Before\n";

    make_point_vec_on_box ( tbox, out, num_x, num_y );

    std::cout << "After\n";

    return;

  } else if ( dim == 3 ) {

    double sepx = (vals[1]-vals[0])/(num_x+1);
    double sepy = (vals[3]-vals[2])/(num_y+1);
    double sepz = (vals[5]-vals[4])/(num_z+1);
    
    double ext[6] = { vals[0]+sepx, vals[1]-sepx,
		 vals[2]+sepy, vals[3]-sepy,
		 vals[4]+sepz, vals[5]-sepz };

    box<3> tbox;
    tbox.set( &ext[0] );

    make_point_vec_on_box ( tbox, out, num_x, num_y, num_z );

    return;

  }
}

void output_values ( const vector<valarray<double> >& v_vals,
		     std::ostream& out )
{
  const int vec_size = v_vals.size();
  if ( vec_size == 0 ) {
    return;
  }
  const int val_size = v_vals[0].size();
  if ( val_size == 0 ) {
    return;
  }

#ifndef NDEBUG
  for ( int i=0; i<vec_size; ++i ) {
    assert ( v_vals[i].size() == val_size );
  }
#endif // NDEBUG

  for ( int val=0; val<val_size; ++val ) {
    for ( int vec=0; vec<vec_size; ++vec ) {
      if ( vec == vec_size - 1 ) {
	out << v_vals[vec][val];
      } else {
	out << v_vals[vec][val] << " ";
      }
    }
    out << "\n";
  }
}

template<int dim>
void
gp_draw_points ( const valarray<double>& pts, std::ostream& out, double z_offset )
{
  assert( pts.size() % dim == 0 );
  if ( pts.size() % dim != 0 ) {
    return;
  }
  int num_points = pts.size() / dim;
  if ( 1 == dim ){
    for ( int i=0; i<num_points; ++i ) {
      std::cout << "Check " << i << "\n";
      out << pts[i] << " " << 0.0 << "\n";
    }
    return;
  } else if ( 2 == dim ) {
    for ( int i=0; i<num_points; ++i ) {
      int pos = 2*i;
      out << pts[ pos ] << " " << pts[ pos +1 ] << " " << 0.0 << "\n"; 
    }
    return;
  } else if ( 3 == dim ) {
    for ( int i=0; i<num_points; ++i ) {
      int pos = 3*i;
      out << pts[ pos ] << " "
	  << pts[ pos +1 ] << " "
	  << pts[ pos +2 ]+z_offset << "\n"; 
    }
    return;
  }
}

template<int dim>
void gp_draw_points ( const valarray<double>& pts, const valarray<double>& vals, ostream& out,  int block, double z_offset )
{
	assert( pts.size() % dim == 0 );
	int num_pts = pts.size()/dim;
	assert( vals.size() == num_pts );
	int count = 0;
	for ( int i=0; i<num_pts; ++i ) {
		if ( block != -1 ) {
			// We make blocks to allow gnuplot to draw surfaces "with lines".
			if ( count == block ) {
				out << "\n";
				count = 0;
			}
		}
		for ( int d=0; d<dim; ++d ) {
			out << pts[ dim*i +d ] << " ";
		}
		out << vals[i] + z_offset << "\n";
		++count;
	}
	out << "\n";
}

template<int dim>
void
gp_draw_labels ( string* prefix,
		 const valarray<double>& pts,
		 const vector<string>* labels,
		 std::ostream& out,
		 double offx,
		 double offy,
		 double offz )
{
  assert( pts.size() % dim == 0 );
  if ( pts.size() % dim != 0 ) {
    return;
  }
  int num_points = pts.size() / dim;

  string pre( "" );
  if ( prefix ) {
    pre = *prefix;
  }

  const vector<string>* plabels;
  vector<string> temp_labels;

  if ( labels ) {
    plabels = labels;
  } else {
    plabels = &temp_labels;
    temp_labels.resize( num_points );
    for ( int i=0; i<num_points; ++i ) {
      std::string str;
      std::stringstream ss;
      ss << i;
      ss >> str;
      temp_labels[i] = str;
    }
  }

  const vector<string>& lab = *plabels;

  if ( 1 == dim ){
    for ( int i=0; i<num_points; ++i ) {
      out << "set label \""
	  << lab[i] << "\" at "
	  << offx+pts[i] << ","
	  << offy <<","
	  << offz <<"\n";
    }
  } else if ( 2 == dim ) {
    for ( int i=0; i<num_points; ++i ) {
      int pos = 2*i;
      out << "set label \""
	  << lab[i]
	  << "\" at "
	  << offx+pts[ pos ] << ","
	  << offy+pts[ pos +1 ] << ","
	  << offz <<"\n"; 
    }
  } else if ( 3 == dim ) {
    for ( int i=0; i<num_points; ++i ) {
      int pos = 3*i;
      out <<"set label \""
	  << lab[i]
	  << "\" at "
	  << offx+pts[ pos ] << ","
	  << offy+pts[ pos +1 ] << ","
	  << offz+pts[ pos +2 ] << "\n"; 
    }
  }
  
}

#if 0
template<int dim>
void not_predicate_point_zero_values ( const interior_predicate<dim>& pred, const valarray<double>& points, valarray<double>& values )
{
	assert( points.size()%dim == 0 );
	
	int num_pts = points.size()/dim;
	assert( num_pts == values.size() );
	
	for ( int i=0; i<num_pts; ++i ) {
		if ( !pred(&points[dim*i]) ) {
			values[i] = 0.0;
		}
	}
}
#endif

//

template void make_point_vec_on_box<1> ( const box<1>&,
					 valarray<double>&,
					 int,
					 int, 
					 int );

template void make_point_vec_on_box<2> ( const box<2>&,
					 valarray<double>&,
					 int,
					 int,
					 int );

template void make_point_vec_on_box<3> ( const box<3>&,
					 valarray<double>&,
					 int,
					 int,
					 int );

template void make_point_vec_interior_box<1> ( const box<1>& bx,
					       valarray<double>& out,
					       int num_x,
					       int num_y,
					       int num_z );

template void make_point_vec_interior_box<2> ( const box<2>& bx,
					       valarray<double>& out,
					       int num_x,
					       int num_y,
					       int num_z );

template void make_point_vec_interior_box<3> ( const box<3>& bx,
					       valarray<double>& out,
					       int num_x,
					       int num_y,
					       int num_z );


template void gp_draw_points<1> ( const valarray<double>&, std::ostream&, double );
template void gp_draw_points<2> ( const valarray<double>&, std::ostream&, double );
template void gp_draw_points<3> ( const valarray<double>&, std::ostream&, double );

template void gp_draw_points<1> ( const valarray<double>&, const valarray<double>&, ostream&, int, double );
template void gp_draw_points<2> ( const valarray<double>&, const valarray<double>&, ostream&, int, double );
template void gp_draw_points<3> ( const valarray<double>&, const valarray<double>&, ostream&, int, double );

template void gp_draw_labels<1> ( string* prefix,
				  const valarray<double>& pts,
				  const vector<string>* labels,
				  std::ostream& out,
				  double offset_x,
				  double offset_y,
				  double offset_z );

template void gp_draw_labels<2> ( string* prefix,
				  const valarray<double>& pts,
				  const vector<string>* labels,
				  std::ostream& out,
				  double offset_x,
				  double offset_y,
				  double offset_z );

template void gp_draw_labels<3> ( string* prefix,
				  const valarray<double>& pts,
				  const vector<string>* labels,
				  std::ostream& out,
				  double offset_x,
				  double offset_y,
				  double offset_z );
#if 0
template void not_predicate_point_zero_values<2> ( const interior_predicate<2>&, const valarray<double>&, valarray<double>& );
template void not_predicate_point_zero_values<3> ( const interior_predicate<3>&, const valarray<double>&, valarray<double>& );
#endif
