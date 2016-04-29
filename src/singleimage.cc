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

#include "singleimage.hh"

#include <tiffio.h>

#include <cmath>
#include <iostream>
using std::ceil;
using std::cout;
using std::endl;
using std::floor;

singleimage::singleimage() :
//	genericfunction( 2, 1 ),
	have_data( false ),
	roi_set( false )
{
	this->set_global_function();
}

singleimage::~singleimage()
{
}

singleimage::singleimage( const singleimage& other )
{
	deep_copy(other);
}

singleimage& singleimage::operator= ( const singleimage& other )
{
	deep_copy(other);
	return *this;
}

void singleimage::deep_copy( const singleimage& other )
{
	have_data = other.have_data;
	if ( other.have_data ) {
		x_size = other.x_size;
		y_size = other.y_size;
		data.resize ( x_size*y_size );
		copy( &other.data[0], &other.data[0]+x_size*y_size, &data[0] );
	}
	roi_set = other.roi_set;
	if ( other.roi_set ) {
		copy( other.roi_box, other.roi_box+4, roi_box );
		copy( other.centre,  other.centre+2,  centre );
		copy( other.half_extent, other.half_extent+2, half_extent );
	}
}

#if 0
// channel -1 for white
void
singleimage::draw_rectangle_on_image( int channel,
                                         const string& filein,
	                                 int x0,
	                                 int x1,
	                                 int y0,
	                                 int y1,
					 const string& fileout )
{
	CImg<unsigned char> base;
	base.load( filein.c_str() );
	
	if ( !base ) {
		cerr << "Could not open file.\n";
		abort();
	}
	
	if ( x0 > x1 || y0 > y1 || x0 < 0 || y0 < 0 || base.dimx() < x1 || base.dimy() < y1 ) {
		return;
	}
	
	unsigned char color[ 4 ] = { 0, 0, 0, 255 };
	assert( -1<= channel && channel <= base.dimv() );
	if ( channel < -1 || base.dimv() < channel ) {
		return;
	}
	if ( channel == -1 ) {
		color[0] = color[1] = color[2] = color[3] = static_cast<unsigned char>(255);
	} else {
		color[ channel ] = static_cast<unsigned char>(255);
	}
	base.draw_rectangle( x0, base.dimy()-y1, x1, base.dimy()-y0, &color[0], 1.0f, ~0U ); // Notice different order of the coordinates.
	base.save( fileout.c_str() );
}
#endif

#if 0
void
singleimage::draw_polygon_on_image( int channel,
	                               const string& filein,
	                               const valarray<double>& points_2d,
	                               const string& fileout )
{
	CImg<unsigned char> base;
	base.load( filein.c_str() );
	
	int base_y = base.dimy();
	
	
	unsigned char color[ 4 ] = { 0, 0, 0, 255 };
	assert( -1<= channel && channel <= base.dimv() );
	if ( channel < -1 || base.dimv() < channel ) {
		return;
	}
	if ( channel == -1 ) {
		color[0] = color[1] = color[2] = static_cast<unsigned char>(255);
	} else {
		color[ channel ] = static_cast<unsigned char>(255);
	}
	assert( points_2d.size() % 2 == 0 );
	int num_pts = points_2d.size() / 2;
	// Start at second point: point 1.
	for ( int i=1; i<num_pts; ++i ) {
		base.draw_line( points_2d[ 2*(i-1) ], base_y - points_2d[ 2*i + 1 ],
		                points_2d[ 2*i     ], base_y - points_2d[ 2*(i-1)     + 1 ],
		                &color[0], 1.0f, ~0U );
	}
	base.draw_line( points_2d[ 2*(num_pts-1) ], base_y - points_2d[ 1 ],
	                points_2d[ 0 ],             base_y - points_2d[ 2*(num_pts-1) + 1 ],
	                &color[0], 1.0f, ~0U );
	base.save( fileout.c_str() );
}
#endif

void
singleimage::load ( const string& filename, int channel )
{
	TIFF* tif_file = TIFFOpen ( filename.c_str(), "r" );
	if ( !tif_file ) {
		abort ();
	}
	
	uint32    width, height;
	
	TIFFGetField ( tif_file, TIFFTAG_IMAGEWIDTH, &width   );
	TIFFGetField ( tif_file, TIFFTAG_IMAGELENGTH, &height );
	
#if 0 // Marking for delete. 2008-12-05 Michael LI.
	cout << "Width of image is " << width
	     << ". Height of image is " << height << ".\n";
#endif
	
	size_t    npixels = width * height;
	
	uint32*    raster = (uint32*)_TIFFmalloc ( npixels * sizeof(uint32) );
	
	if ( raster == NULL ) {
		cerr << "Unable to allocate memory." << endl;
		abort ();
	}
	
	if ( !TIFFReadRGBAImage ( tif_file, width, height, raster, 0 ) ) {
		cerr << "Unable to read from tif file." << endl;
		abort ();
	}

	// Attempting -1 here to disallow attempts to evaluate off the ends.
	// The sizes are used to determine the parent box to use.
	// It may be that the subtraction of one should be done outside.
	// 2009-01-20 ML.
	// CANNOT subtract one from here as results in skewed picture.
	// Needs another think. 2009-01-20 ML.
	x_size = width;
	y_size = height;
	set_roi( 0, width, 0, height );
	data.resize ( npixels );
	
	int idx;
	
	switch ( channel ) {
	case 0:
	{
		// Debugging apparent inability to get red channel. It was due to
		// two uses of load. The second was being left unmodified to use the
		// red channel. 2009-01-26 ML.
		// cout << "Filename is " << filename <<". Getting red channel.\n";
		for ( int y=0; y<height; ++y ) {
			for ( int x=0; x<width; ++x ) {
				idx = y*width + x;
				data[idx] = TIFFGetR(raster[idx]);
			}
		}
		break;
	}
	case 1:
	{
		for ( int y=0; y<height; ++y ) {
			for ( int x=0; x<width; ++x ) {
				idx = y*width + x;
				data[idx] = TIFFGetG(raster[idx]);
			}
		}
		break;
	}
	case 2:
	{
		for ( int y=0; y<height; ++y ) {
			for ( int x=0; x<width; ++x ) {
				idx = y*width + x;
				data[idx] = static_cast<short>( TIFFGetB(raster[idx]) );
			}
		}
		break;
	}
	case 3:
	{
		for ( int y=0; y<height; ++y ) {
			for ( int x=0; x<width; ++x ) {
				idx = y*width + x;
				data[idx] = TIFFGetA(raster[idx]);
			}
		}
		break;
	}
	default:
		abort();
	}
	
	have_data = true;
	
	_TIFFfree ( raster );

	TIFFClose ( tif_file );

#if 0	
	{
		string tmp_name = "green_channel";
		tmp_name += out_id;
		tmp_name += ".gpdat";
		ofstream file ( tmp_name.c_str() );
		file << "splot '-' matrix every 10:10 w l\n";
		for ( int i=0; i<height; ++i ) {
			for ( int j=0; j<width; ++j ) {
				file << TIFFGetG(raster[ i*width + j ]) << " ";
			}
			file << "\n";
		}
	}
#endif
}

void
singleimage::rescale ( double scaling )
{
	assert ( have_data );
	data *= scaling;
}

double
singleimage::rescale_percentage ( double select_max )
{
	assert ( have_data );
	double max;
	if ( select_max > 0.0 ) {
		max = select_max;
	} else {
		max = data.max();
	}

	double divisor = 0.01*max;
//	cout << "Max is " << divisor << "\n";
	if ( divisor > 0.0 ) {
		data /= divisor;
	}
	return max;
}

#if 0
void
singleimage::load( const string& filename, int channel )
{
	assert( 0<= channel && channel <= 4 );
	if ( channel < 0 || 4 < channel ) {
		return;
	}
	CImg<unsigned char> tmp_cimg;
	tmp_cimg.load( filename.c_str() );
	if ( tmp_cimg ) {
		x_size = tmp_cimg.dimx();
		assert( x_size > 0 );
		y_size = tmp_cimg.dimy();
		assert( y_size > 0 );
	}
	if ( data ) {
		assert( have_data );
		delete[] data;
		data = 0;
	}
	try {
		data = new int[ y_size*x_size ];
	} catch ( ... ) {
		cerr << "Unable to allocate memory.\n";
		abort();
	}
	int v = tmp_cimg.dimv();
	assert( v==3 || v==4 );
	//unsigned char* puc = tmp_cimg.ptr();
	for ( int y=0; y<y_size; ++y ) {
		for ( int x=0; x<x_size; ++x ) {
			// Must reverse CImg's y direction.
			// The ``How pixel data are stored with CImg'' section of the CImg
			// documentation is empty. CImg seems to store the RGB channels partitioned
			// within its data region rather than with the colours interlaced as expected.
			// data[y*x_size + x]= static_cast<int>( puc[ (y_size-y-1)*x_size*v + x + channel ] );
			data[y*x_size + x] = static_cast<int>( tmp_cimg( x, y_size-y-1, 0, channel ) );
		}
	}
#ifndef NDEBUG // Dump image as it seems to be incorrect.
	ofstream img_dump( "dump_image.gnuplot" );
	assert( img_dump );
	for ( int y=0; y<y_size; ++y ) {
		for ( int x=0; x<x_size; ++x ) {
			img_dump << x << " " << y << " " << data[y*x_size+x] << "\n";
		}
		img_dump << "\n";
	}
#endif
	
	have_data = true;
}
#endif

int singleimage::get_x_size() const
{
	return x_size;
}
int singleimage::get_y_size() const
{
	return y_size;
}

void
singleimage::set_roi( double x0, double x1, double y0, double y1 )
{
	assert( x0 < x1 && y0 < y1 );
	
	roi_box[0]=x0;
	roi_box[1]=x1;
	roi_box[2]=y0;
	roi_box[3]=y1;
	
	centre[0] = 0.5*(x0+x1);
	centre[1] = 0.5*(y0+y1);
	
	half_extent[0] = 0.5*(x1-x0);
	half_extent[1] = 0.5*(y1-y0);
	
	roi_set = true;
}

#if 0
// Bilinear interpolation.
double
singleimage::evaluate_global( double x, double y ) const
{
	assert( have_data );
	assert( 0 <= x && x <= x_size && 0 <= y && y <= y_size );
	if ( x < 0 || x_size < x || y < 0 || y_size < y ) {
		return 100.0;
	}
	int across_x = 0;
	int down_y   = 0;
	
	while ( static_cast<double>(++across_x) <= x );
	while ( static_cast<double>(++down_y  ) <= y );
	
	assert( across_x <= x_size );
	assert( down_y   <= y_size );
	
	int across_x_minus = across_x - 1;
	int down_y_minus   = down_y   - 1;
	
	int NW = data[ down_y_minus*x_size + across_x_minus ], NE = data[ down_y_minus*x_size + across_x ];
	int SW = data[ down_y*x_size       + across_x_minus ], SE = data[ down_y*x_size       + across_x ];
	
	double along_top    = NW*(across_x - x) + NE*(x - across_x_minus);
	double along_bottom = SW*(across_x - x) + SE*(x - across_x_minus);
	
	return along_top*( down_y - y ) + along_bottom*( y - down_y_minus );
}
#endif


double
singleimage::evaluate_global( double x, double y ) const
{
	assert( have_data );
	assert( 0 <= x && x <= x_size && 0 <= y && y <= y_size );
	if ( x < -1e-10 || x_size+1e-10 < x || y < -1e-10 || y_size+1e-10 < y ) {
		return 777.0; // TODO: return fool-proof 0.0 when testing is over.
	}
	
	// Find top left corner of relevant box.
	int across_x = static_cast<int>( ceil(x) )-1;
	int down_y   = static_cast<int>( ceil(y) )-1;

	// Don't drop off the lower end.
	if ( across_x < 0 ) {
		across_x = 0;
	}
	if ( down_y < 0 ) {
		down_y = 0;
	}
	
	
#if 0
	int across_x = 0;
	int down_y   = 0;
	
	while ( ++across_x <= x ) { };
	while ( ++down_y <= y ) { };
	
	--across_x;
	--down_y;

	if ( (across_x != new_across_x)||(down_y!=new_down_y) ) {
		cout << "WRONG! - x " << x << ", y " << y << "\n";
		cout << "across_x " << across_x << ", new_across_x " << new_across_x << "\n";
		cout << "down_y " << down_y << ", new_down_y " << new_down_y << "\n";
	}
#endif
	
	assert( across_x < x_size );
	assert( down_y   < y_size );
	
	int across_x_plus = across_x + 1;
	int down_y_plus   = down_y   + 1;
	
	int NW = data[ down_y*x_size      + across_x ], NE = data[ down_y*x_size      + across_x_plus ];
	int SW = data[ down_y_plus*x_size + across_x ], SE = data[ down_y_plus*x_size + across_x_plus ];
	
	double along_top    = NW*(across_x_plus - x) + NE*(x - across_x);
	double along_bottom = SW*(across_x_plus - x) + SE*(x - across_x);
	
	return along_top*( down_y_plus - y ) + along_bottom*( y - down_y );
}

double
singleimage::evaluate ( const double* x ) const
{
	return evaluate_global ( x[0], x[1] );
}

void
singleimage::evaluate_grad ( const double* co, double* grad ) const
{
	const double &    x = co[0];
	const double &    y = co[1];
	
	assert( have_data );
	assert( 0 <= x && x <= x_size && 0 <= y && y <= y_size );

	// Find top left corner of relevant box.
	int across_x = static_cast<int>( ceil(x) )-1;
	int down_y   = static_cast<int>( ceil(y) )-1;

	// Don't drop off the lower end.
	if ( across_x < 0 ) {
		across_x = 0;
	}
	if ( down_y < 0 ) {
		down_y = 0;
	}
#if 0
	int across_x = 0;
	int down_y   = 0;
	
	// Find top left corner of relevant box.
	while ( ++across_x <= x ) { };
	while ( ++down_y <= y ) { };
	
	--across_x;
	--down_y;
#endif
	assert( across_x < x_size );
	assert( down_y   < y_size );
	
	int across_x_plus = across_x + 1;
	int down_y_plus   = down_y   + 1;
	
	int NW = data[ down_y*x_size      + across_x ], NE = data[ down_y*x_size      + across_x_plus ];
	int SW = data[ down_y_plus*x_size + across_x ], SE = data[ down_y_plus*x_size + across_x_plus ];
	
	double along_top    = NW*(across_x_plus - x) + NE*(x - across_x);
	double dx_along_top = -NW + NE;
	
	double along_bottom    = SW*(across_x_plus - x) + SE*(x - across_x);
	double dx_along_bottom = -SW + SE;
	
	// double value = along_top*( down_y_plus - y ) + along_bottom*( y - down_y );
	// dx
	grad[0] = dx_along_top*( down_y_plus - y ) + dx_along_bottom*( y - down_y );
	// dy
	grad[1] = -along_top + along_bottom;
}

double
singleimage::evaluate_local( double e0, double e1 ) const
{
	assert( have_data && roi_set );
	assert( -1. <= e0 && e0 <= 1. && -1. <= e1 && e1 <= 1. );
	if ( e0 < -1. || 1. < e0 || e1 < -1. || 1. < e1 ) {
		return 100.0; // TODO: return fool-proof 0.0 when testing is over.
	}
	return evaluate_global( e0*half_extent[0] + centre[0], e1*half_extent[1] + centre[1] );
}

void
singleimage::evaluate_global_3D_pts ( const double* in, int num_pts, double* out ) const
{
	assert ( in && out );
	
	size_t    idx = 0;
	
	for ( int i=0; i<num_pts; ++i ) {
		idx = 3*i;
		out[i] = evaluate_global ( in[idx], in[idx+1] );
		//cout << "x " << in[idx] << ", y " << in[idx+1] << ", out " << out[i] << ", evaluate_global "
		//     << evaluate_global ( in[idx], in[idx+1] ) << "\n";
	}
}

void
singleimage::save_local_gnuplot( const string& fileout,
	                            int points ) const
{
	assert( have_data && roi_set );
	ofstream out( fileout.c_str() );
	if ( !out ) {
		cerr << "Problem opening file.\n";
		abort();
	}
	double div[points];
	double sep = 2.0/(points-1);
	for ( int i=0; i<points; ++i ) {
		div[i] = -1.0 + i*sep;
	}
	for ( int y=0; y<points; ++y ) {
		for ( int x=0; x<points; ++x ) {
			out << div[x] << " " << div[y] << " " << evaluate_local( div[x], div[y] ) << "\n";
		}
		out << "\n";
	}
}
