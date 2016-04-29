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

#include "extract_geometry_2d.hh"

//#include "gf_singleimage.hh"

#include <vtkTIFFReader.h>
#include <vtkTIFFWriter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkImageThreshold.h>
#include <vtkImageOpenClose3D.h>
#include <vtkImageIslandRemoval2D.h>
#include <vtkMarchingSquares.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageMapper.h>
#include <vtkContourFilter.h>
#include <vtkGaussianSplatter.h>
#include <vtkContourGrid.h>
#include <vtkImageToStructuredPoints.h>
#include <vtkImageToPolyDataFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkStructuredPointsWriter.h>
#include <vtkImageMagnitude.h>
#include <vtkStripper.h>
#include <vtkImageMedian3D.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkHull.h>
#include <vtkImageQuantizeRGBToIndex.h>
#include <vtkLookupTable.h>
#include <vtkScalarsToColors.h>
#include <vtkStructuredPointsToPolyDataFilter.h>
#include <vtkDataSetToPolyDataFilter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkCellArray.h>
#include <vtkCell.h>

//#include <CImg.h>
//using namespace cimg_library;
//using cimg_library::CImg;

#include <petsc.h>

//#include <gdf/imageseries.hh>
//using namespace gdf;

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

using std::map;
using std::cout;
using std::cerr;
using std::getenv;
using std::min;
using std::max;
using std::floor;
using std::fstream;
using std::string;
using std::stringstream;
using std::getline;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::valarray;

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

template <typename Container>
void extract_geometry_2d ( const string& filename,
                           Container& out,
#if 0
                           int gaussian_kernel_size,
#endif
                           double binarize_threshold,

                           int close_kernel_size,
                           int island_keep_size )
{
	assert( binarize_threshold > 0. );
	assert( close_kernel_size > 0 );
	assert( island_keep_size > 0. );
	
	string prefix;
	{
		string::size_type dot_position = filename.find_last_of( "." );
		if ( dot_position == string::npos ) {
			prefix = filename;
		} else {
			prefix.assign( filename, 0, dot_position );
		}
	}
	
	vtkTIFFWriter* tiff_writer = vtkTIFFWriter::New();
	
	vtkTIFFReader* tiff_reader = vtkTIFFReader::New();
	tiff_reader->SetFileName( filename.c_str() );
	
	vtkImageGaussianSmooth* gaussian = vtkImageGaussianSmooth::New();
	
#if 0
	gaussian->SetRadiusFactor( gaussian_kernel_size );
	gaussian->SetInputConnection( tiff_reader->GetOutputPort() );
	
	{
		string tmp_name = prefix;
		tmp_name += " 05post_gaussian.tiff";
		tiff_writer->SetFileName( tmp_name.c_str() );
		tiff_writer->SetInputConnection( gaussian->GetOutputPort() );
		tiff_writer->Write();
	}
#endif
	
	vtkXMLImageDataWriter* xmlpd_writer = vtkXMLImageDataWriter::New();
	
#if 0
	xmlpd_writer->SetInputConnection( gaussian->GetOutputPort() );
	
	xmlpd_writer->SetFileName( "07post_gaussian.vti" );
	xmlpd_writer->Write();
#endif

	
	vtkImageThreshold* threshold = vtkImageThreshold::New();
	threshold->SetInValue( 255 );
	threshold->SetOutValue( 0 );
	
	cout << "Binarize_threshold is " << binarize_threshold << "\n";
	
	threshold->ThresholdByUpper( floor(0.5+binarize_threshold) );
	
	// threshold->SetInputConnection( gaussian->GetOutputPort() );
	threshold->SetInputConnection( tiff_reader->GetOutputPort() );
	
	tiff_writer->SetFileName( "10post_threshold.tiff" );
	tiff_writer->SetInputConnection( threshold->GetOutputPort() );
	tiff_writer->Write();

	vtkImageOpenClose3D* close_filt = vtkImageOpenClose3D::New();
	close_filt->SetOpenValue( 0 );
	close_filt->SetCloseValue( 255 );
	close_filt->SetKernelSize( close_kernel_size, close_kernel_size, 1 );
	close_filt->SetInputConnection( threshold->GetOutputPort() );

	tiff_writer->SetFileName( "20post_close.tiff" );
	tiff_writer->SetInputConnection( close_filt->GetOutputPort() );
	tiff_writer->Write();

	vtkImageIslandRemoval2D* island = vtkImageIslandRemoval2D::New();
	island->SetIslandValue( 255 );
	island->SetAreaThreshold( island_keep_size );
	island->SetInputConnection( close_filt->GetOutputPort() );
	island->SetReplaceValue( 0.0 );

	tiff_writer->SetFileName( "30post_island.tiff" );
	tiff_writer->SetInputConnection( island->GetOutputPort() );
	tiff_writer->Write();

	vtkImageMagnitude* magnitude = vtkImageMagnitude::New();
	magnitude->SetInputConnection( island->GetOutputPort() );

	vtkMarchingSquares* msquares = vtkMarchingSquares::New();
	msquares->SetInputConnection( magnitude->GetOutputPort() );
	msquares->SetNumberOfContours(1);
	msquares->SetValue( 0, 255 );
	
	vtkPolyDataWriter* polydata_writer = vtkPolyDataWriter::New();
	polydata_writer->SetFileName( "40post_marching.vtp" );
	polydata_writer->SetInputConnection( msquares->GetOutputPort() );
	polydata_writer->Write();
	polydata_writer->Delete();
	polydata_writer = vtkPolyDataWriter::New();
	
	vtkPolyData* final_polydata = vtkPolyData::New();
	
	vtkStripper* stripper = vtkStripper::New();
	
	
	stripper->SetInputConnection( msquares->GetOutputPort() );
	
	//Choose one or other.
	//stripper->SetOutput( final_polydata );
	msquares->SetOutput( final_polydata );
	
	//final_polydata->SetInputConnection( stripper->GetOutputPort() );
	
	polydata_writer->SetFileName( "50post_stripping.vtp" );
	polydata_writer->SetInputConnection( stripper->GetOutputPort() );
	polydata_writer->Write();
	
	vtkPoints* points = final_polydata->GetPoints();

	vtkCellArray* boundary_lines = final_polydata->GetLines();
	
	{
		int num_lines = boundary_lines->GetNumberOfCells();
		cout << "There are " << num_lines << " lines in the boundary.\n";
		
		// GetNextCell returns by reference the number of points in the cell and
		// a pointer to the points.
		vtkIdType npts;
		vtkIdType* pts;
		
		vtkIdType start_traversal;
		map<vtkIdType,vtkIdType> traversal;
		
		int count;
		boundary_lines->InitTraversal();
		for ( count=0; count < num_lines; ++count )
		{
			boundary_lines->GetNextCell(npts,pts);
			
			assert( pts );
			
			//cout << "npts was " << npts << "\n";
			assert( npts == 2 );
			// cout << "Traversal "<< pts[0] << " and " << pts[1] << "\n";
			
			if( count == 0 ) {
				start_traversal = pts[0];
			}
			
			traversal[pts[0]] = pts[1];
		}


		// Handle domains with holes.
		Container    prepare_out;
		prepare_out.resize ( 2*num_lines );

		{
			ofstream geometry_gnuplot( "geometry.gnuplot" );
			
			if ( !geometry_gnuplot ) {
				cout << "Unable to open geometry file.\n";
				abort();
			}
			
			double* tmp_d;
			
			boundary_lines->InitTraversal();
			
			double keep_first_point[3];
			double* tmp = 0;
			
			vtkIdType travel = start_traversal;

			// Delete. Handle domains with holes using prepare_out.
			// out.resize( 2*num_lines );
			
			int pt = 0;
			do {
				tmp = points->GetPoint( travel );
				assert( tmp!=0 );
				
				prepare_out[ 2*pt   ] = tmp[0];
				prepare_out[ 2*pt+1 ] = tmp[1];
				
				geometry_gnuplot << tmp[0] << " " << tmp[1] << " " << tmp[2] << "\n";
				assert( geometry_gnuplot );
				
				travel = traversal[travel];
				// cout << "travel is " << travel << ", traversal[travel] is " << traversal[travel] << "\n";
				
				++pt;
				assert( pt < 2000000 );
			} while ( travel!= start_traversal );

			out.resize ( 2*pt );

			for ( int i=0; i<2*pt; ++i ) {
				out[i] = prepare_out[i];
			}
			
			geometry_gnuplot << "\n";
			geometry_gnuplot.close();
		}
	}

#if 0	
	double bounds[6]; // x0, x1, y0, y1, z0, z1
	final_polydata->GetBounds( bounds );
	
	GF_SingleImage single;
	single.load( filename, 1 );
	single.set_roi( bounds[0], bounds[1], bounds[2], bounds[3] );
	single.save_local_gnuplot( "values.gnuplot", 100 );
	
	GF_SingleImage::draw_rectangle_on_image( -1,
	                                         filename,
	                                         bounds[0], bounds[1],
	                                         bounds[2], bounds[3],
	                                         "rectangle_on_image.tiff" );
#endif

	polydata_writer->Delete();
	final_polydata->Delete();
	stripper->Delete();
	msquares->Delete();
	magnitude->Delete();
	island->Delete();
	close_filt->Delete();
	threshold->Delete();
	xmlpd_writer->Delete();
	gaussian->Delete();
	tiff_reader->Delete();
	tiff_writer->Delete();
}

void apply_gaussian ( const string& filename, int gaussian_kernel_size )
{
	string prefix;
	{
		string::size_type dot_position = filename.find_last_of( "." );
		if ( dot_position == string::npos ) {
			prefix = filename;
		} else {
			prefix.assign( filename, 0, dot_position );
		}
	}
	
	vtkTIFFWriter* tiff_writer = vtkTIFFWriter::New();
	
	vtkTIFFReader* tiff_reader = vtkTIFFReader::New();
	tiff_reader->SetFileName( filename.c_str() );
	
	vtkImageGaussianSmooth* gaussian = vtkImageGaussianSmooth::New();
	gaussian->SetRadiusFactor( gaussian_kernel_size );
	gaussian->SetInputConnection( tiff_reader->GetOutputPort() );
	
	{
		string tmp_name = prefix;
		tmp_name += " post_gaussian.tiff";
		tiff_writer->SetFileName( tmp_name.c_str() );
		tiff_writer->SetInputConnection( gaussian->GetOutputPort() );
		tiff_writer->Write();
	}
	
	gaussian->Delete();
	tiff_reader->Delete();
	tiff_writer->Delete();
}

template <typename value_t>
void simple_off_input_geometry_2d ( const string &    filename, valarray<value_t> &    geom )
{
	ifstream    file ( filename.c_str() );
	if ( ! file ) {
		cerr << "Failed to open file.\n";
		return;
	}

	string   ln;

	getline ( file, ln ); // This should be "OFF" but we ignore first line regardless.

	getline ( file, ln );
	
	stringstream    ss(ln);

	int    num_pts;

	ss >> num_pts;
	if ( !ss.good() ) {
		cerr << "Unable to get the number of points.\n";
		return;
	}

	cout << "Number of points is " << num_pts << "\n";
	
	// Debugging found missing factor of two here. 2008-12-03 Michael LI.
	geom.resize(2*num_pts);
	value_t    tmp;
	for ( int i=0; i<num_pts; ++i ) {
		getline ( file, ln );
		if ( !file.good() ) {
			cerr << "Problem reading from file. It may be that you are not giving 3D points.\n";
			abort();
			return;
		}

		ss.str(ln);
		ss >> tmp;
		geom[2*i] = tmp;
		ss >> tmp;
		geom[2*i+1] = tmp;
	}
}

//

template void extract_geometry_2d<valarray<int> > ( const string&, valarray<int>&, double, int, int );
template void extract_geometry_2d<valarray<double> > ( const string&, valarray<double>&, double, int, int );

template void simple_off_input_geometry_2d<int> ( const string &, valarray<int> &    geom );
template void simple_off_input_geometry_2d<double> ( const string &, valarray<double> &    geom );
