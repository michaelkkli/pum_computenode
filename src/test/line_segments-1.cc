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

#include "basic_dbinary_tree.hh"
#include "extract_geometry_2d.hh"
#include "line_segments.hh"
#include "refinement_structure.hh"
#include "singleimage.hh"
#include "vtk_output.hh"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <valarray>

#include <sys/stat.h>

using std::cerr;
using std::copy;
using std::cout;
using std::inserter;
using std::ostringstream;
using std::set_difference;
using std::slice;
using std::string;
using std::system; // Run command in a shell using system ( "" )
using std::valarray;

struct averaging_values {
	vector<string>               labels;
	vector<valarray<double> >    values;
};

int main ( int argc, char* argv[] ) {

	string    manifest1_dir;
	
	if ( argc == 1 ) {
		manifest1_dir = "input/F2 Rep 1 blur/";
	} else {
		manifest1_dir.assign ( argv[1] );
		if ( manifest1_dir.at( manifest1_dir.size()-1 ) != '/' ) {
			manifest1_dir.append ( "/" );
		}
	}

	string    manifest1_filename = manifest1_dir + "manifest";
	
	string    dom_string = "extract_geometry_2d-1_contour.off";
	
	valarray<double> vcoords;
	simple_off_input_geometry_2d ( dom_string, vcoords );
	
// Force whole domain. This gives an anticlockwise-winding boundary
// whereas VTK give a clockwise-winding boundary.
#if 0
	vcoords.resize(8);
	vcoords[0] = 1.0;    vcoords[1] = 1.0;
	vcoords[2] = 510.0;    vcoords[3] = 1.0;
	vcoords[4] = 510.0;    vcoords[5] = 510.0;
	vcoords[6] = 1.0;    vcoords[7] = 510.0;
#endif

// Force whole domain. Clockwise-winding boundary.
#if 0
	vcoords.resize(8);
	vcoords[6] = 1.0;    vcoords[7] = 1.0;
	vcoords[4] = 510.0;    vcoords[5] = 1.0;
	vcoords[2] = 510.0;    vcoords[3] = 510.0;
	vcoords[0] = 1.0;    vcoords[1] = 510.0;
#endif
	
	line_segments    domain_segs;
	{
		int tmp_size = vcoords.size();
		domain_segs.segments.resize ( tmp_size +2 );
		copy ( &vcoords[0], &(vcoords[0])+tmp_size, &(domain_segs.segments[0]) );
		
		// Debugging found the last point was being set incorrectly.
		// Using -2 and -1 to adjust indices. 2008-12-01 Michael LI.
		domain_segs.segments[ tmp_size ]   = vcoords[0];
		domain_segs.segments[ tmp_size+1 ] = vcoords[1];
		
		// Keep in case of need to inspect
		// int num_segs = tmp_size/2;
		// ofstream    file ( "petsc_solver-3_stepper_labelled_bdry.gnuplot" );
		// gnuplot_output<2> ( domain_segs, file, (num_segs<1000) );
	}
	
	singleimage::ptr    image_initial;
	
	string    image1_string  = "initial_condition.tif";
	double    rescale_image1 = 0.3;
	
	// Used further below too.
	int col_channel = 1;
	
	image_initial.reset ( new singleimage );
	image_initial->load ( image1_string, col_channel );
#if 1
	image_initial->rescale_percentage ();
#else
	image_initial->rescale ( rescale_image1 );
#endif
	
	box<2>    parent_box;
	{
		double    ext[4];
		ext[0] = 0.0;
		ext[1] = image_initial->get_x_size()-1;
		ext[2] = 0.0;
		ext[3] = image_initial->get_y_size()-1;
		parent_box.set ( ext );
	}
	
	int    viz_refinement = 9;
	
	basic_dbinary_tree<2>::ptr basic_dtree ( new basic_dbinary_tree<2> );
	basic_dtree->initialize ( parent_box, viz_refinement );
	
	refinement_structure<2>   viz_ref_struct;
	viz_ref_struct.dtree = basic_dtree;
	viz_ref_struct.levels["0"] = viz_refinement;

	vtk_output    bdry_vtkout;
	
	vector<string>    vec_bdry_keys;
	set<string>   bdry_keys;
	
	// Set to outside for CLFLIP data. 2009-08-16 ML.
	make_boundary_vtk_output ( domain_segs,
	                           viz_ref_struct,
	                           bdry_vtkout,
	                           outside, &vec_bdry_keys );

	copy ( vec_bdry_keys.begin(),
	       vec_bdry_keys.end(),
	       inserter ( bdry_keys, bdry_keys.end() )
	);
	
	bdry_vtkout.scalars.resize ( bdry_vtkout.points_3D.size() / 3, 0.0 );

	image_initial->evaluate_global_3D_pts ( &(bdry_vtkout.points_3D[0]),
	                                        bdry_vtkout.points_3D.size()/3,
	                                        &(bdry_vtkout.scalars[0]) );
	
	bdry_vtkout.points_3D[ slice(2, bdry_vtkout.scalars.size(), 3) ] = bdry_vtkout.scalars;
	
#if 1
	vtk_append_raw ( bdry_vtkout,
	                 "line_segments-1_bdry_test.vtp" );
#else
	// Aid to debug.
	vtk_simple_output( bdry_vtkout,
	                   "line_segments-1_bdry_test.vtp" );
#endif

#if 0
	// Keep a copy.
	valarray<int>    bdry_connectivity ( bdry_vtkout.connectivity.size() );
	bdry_connectivity = bdry_vtkout.connectivity;
	
	valarray<int>    bdry_offsets ( bdry_vtkout.offsets.size() );
	bdry_offsets = bdry_vtkout.offsets;
#endif

	vector<string> tree_all_keys;
	generate_keys< 2 >( viz_ref_struct, tree_all_keys );
	
	set<string>    inside_poss_bdry_keys;
	get_inside_keys ( tree_all_keys, *basic_dtree, domain_segs, inside_poss_bdry_keys );
	
	// Let's get the inside-only keys by finding all inside keys and subtracting
	// boundary keys.
	
	set<string>    strict_inside_keys;
	set_difference ( inside_poss_bdry_keys.begin(),
	                 inside_poss_bdry_keys.end(),
	                 bdry_keys.begin(),
	                 bdry_keys.end(),
	                 inserter( strict_inside_keys,strict_inside_keys.begin() ) );
	
	cout << "There are " << bdry_keys.size()
	     << " boundary patches, " << inside_poss_bdry_keys.size()
	     << " patches with their centre inside the domain and "
	     << strict_inside_keys.size() << " patches lying strictly inside the domain.\n";
	
	vtk_output    strict_inside_vtkout;

	vector<string>    strict_inside_keys_vec ( strict_inside_keys.size() );
	copy ( strict_inside_keys.begin(),
	       strict_inside_keys.end(),
	       strict_inside_keys_vec.begin() );

	basic_dtree->get_grid_points_3D ( strict_inside_vtkout.points_3D );
	basic_dtree->get_box_grid_connectivity_offsets ( strict_inside_keys_vec,
	                                                 strict_inside_vtkout.connectivity,
	                                                 strict_inside_vtkout.offsets );

	strict_inside_vtkout.scalars.resize ( strict_inside_vtkout.points_3D.size()/3, 0.0 );

	image_initial->evaluate_global_3D_pts ( &(strict_inside_vtkout.points_3D[0]),
	                                        strict_inside_vtkout.points_3D.size()/3,
	                                        &(strict_inside_vtkout.scalars[0]) );

	strict_inside_vtkout.points_3D[ slice(2, strict_inside_vtkout.scalars.size(), 3) ] = strict_inside_vtkout.scalars;
	// For transfering height info.
	// vtk_st->points_3D[ slice(2, vtk_st->scalars.size(), 3) ] = vtk_st->scalars;
	
#if 1
	vtk_append_raw ( strict_inside_vtkout,
	                 "line_segments-1_strict_inside_test.vtp" );
#else
	// Present to aid debugging when vtk_append_raw
	// causes paraview segfault.
	vtk_simple_output ( strict_inside_vtkout,
	                 "line_segments-1_strict_inside_test.vtp" );
#endif
	
	// Let's make a unified file. The interior connectivity will come before the boundary
	// connectivity to reduce the shifting needed with the boundary offsets.
	
	vtk_output    unified_vtkout;
	unified_vtkout.points_3D.resize ( bdry_vtkout.points_3D.size() );
	unified_vtkout.points_3D = bdry_vtkout.points_3D;
	
	// No modifications need to be made to connectivity as numbering is consistent.
	unified_vtkout.connectivity.resize (
	    strict_inside_vtkout.connectivity.size() + bdry_vtkout.connectivity.size() );
	unified_vtkout.connectivity[ slice(0,strict_inside_vtkout.connectivity.size(),1) ] =
	    strict_inside_vtkout.connectivity;
	unified_vtkout.connectivity[ slice(strict_inside_vtkout.connectivity.size(), bdry_vtkout.connectivity.size(),1) ] =
	    bdry_vtkout.connectivity;

	// An adjustment will need to be made to offsets.
	size_t    strict_inside_connect_new_offset = strict_inside_vtkout.offsets[strict_inside_vtkout.offsets.size()-1];
	
	unified_vtkout.offsets.resize (
	    strict_inside_vtkout.offsets.size() + bdry_vtkout.offsets.size() );
	unified_vtkout.offsets[ slice(0,strict_inside_vtkout.offsets.size(),1) ] =
	    strict_inside_vtkout.offsets;
	unified_vtkout.offsets[ slice(strict_inside_vtkout.offsets.size(),bdry_vtkout.offsets.size(),1) ] =
	    bdry_vtkout.offsets + valarray<int> ( strict_inside_connect_new_offset, bdry_vtkout.offsets.size() );
	    
	unified_vtkout.scalars.resize ( bdry_vtkout.scalars.size() );
	unified_vtkout.scalars = bdry_vtkout.scalars;
	
#if 1
	vtk_append_raw ( unified_vtkout,
	                 "line_segments-1_unifed_test.vtp" );
#else
	// Present to aid debugging when vtk_append_raw
	// causes paraview segfault.
	vtk_simple_output ( unified_vtkout,
	                 "line_segments-1_unified_test.vtp" );
#endif



	vector<string>    manifest1_content;
	// Load manifest data if the relevant simulation mode is requested.
	{
		if ( manifest1_filename == "" ) {
			cerr << "A text file listing tiff images must be provided (please give an argument to the option -manifest1).\n";
			abort();
		}
		ifstream    file ( manifest1_filename.c_str() );
		if ( !file ) {
			cerr << "There was a problem opening the manifest1 file \""
				<< manifest1_filename << "\".\n";
			abort();
		}
		string      ln;
		while ( getline ( file, ln ) ) {
			manifest1_content.push_back ( ln );
		}
		int size = manifest1_content.size();
		struct stat    ignore;
		ostringstream    ss;
		int ret;
		for ( int i=0; i<size; ++i ) {
			ss << manifest1_dir << manifest1_content[i];
			ret = stat ( (ss.str()).c_str(), &ignore );
			//cout << "(ss.str()).c_str() is " << (ss.str()).c_str() << "\n";
			if ( ret != 0 ) {
				cerr << "The file \"" << manifest1_content[i]
					<< "\" named in \"" << manifest1_filename
					<< "\" does not exist.\n";
				abort();
			}
			ss.str ( "" );
		}
	}
	
	int   iterations = manifest1_content.size();
	
	cout << "Iterations has been set to " << iterations << ".\n";

	int   num_grid_pts = basic_dtree->get_num_grid_points();
	
	valarray<bool>    inside_mask ( false, unified_vtkout.points_3D.size()/3 );
	
	get_inside_points_3D ( &(unified_vtkout.points_3D[0]), num_grid_pts, domain_segs, &(inside_mask[0]) );
	
	int num_inside = 0;
	for ( int i=0; i<num_grid_pts; ++i ) {
		if ( inside_mask[i] ) {
			++num_inside;
		}
	}
	cout << "There are " << num_inside
	     << " grid points inside the domain.\n";

	averaging_values    avg_vals;
	avg_vals.labels.resize(3);
	avg_vals.labels[0] = "Total";
	avg_vals.labels[1] = "Above mean";
	avg_vals.labels[2] = "Below mean";
	
	avg_vals.values.resize(3);
	avg_vals.values[0].resize(iterations);
	avg_vals.values[1].resize(iterations);
	avg_vals.values[2].resize(iterations);

	double max_for_series;

	for ( int t=0; t<iterations; ++t ) {
		ostringstream    ss;
		ss << manifest1_dir << manifest1_content[t];
		
		// image_initial.reset ( new singleimage );
		image_initial->load ( ss.str(), col_channel );
		
		if ( 0 == t ) {
			max_for_series = image_initial->rescale_percentage ();
		} else {
			image_initial->rescale_percentage ( max_for_series );
		}

		
		image_initial->evaluate_global_3D_pts ( &(unified_vtkout.points_3D[0]),
	                                        unified_vtkout.points_3D.size()/3,
	                                        &(unified_vtkout.scalars[0]) );

		unified_vtkout.points_3D[ slice(2, unified_vtkout.scalars.size(), 3) ] = unified_vtkout.scalars;
		
		{
			ostringstream name;
			name << "line_segments-1_manifest_" << t << ".vtp";
			
			vtk_append_raw ( unified_vtkout, (name.str()).c_str() );
		}
				
		{
			const valarray<double> &    inside_values = (unified_vtkout.scalars)[ inside_mask ];
			(avg_vals.values[0])[t] = inside_values.sum();
			
			double    mean = inside_values.sum()/num_inside;
			{
				const valarray<bool> &    above_mean_mask = ( inside_values > mean );
				const valarray<double> &    above_mean_vals = inside_values[ above_mean_mask ];
				(avg_vals.values[1])[t] = above_mean_vals.sum();
			}
			{
				const valarray<bool> &    below_mean_mask = ( inside_values < mean );
				const valarray<double> &    below_mean_vals = inside_values[ below_mean_mask ];
				(avg_vals.values[2])[t] = below_mean_vals.sum();
			}
			
			valarray<double>    original_scalars ( unified_vtkout.scalars );
			
			// Set to zero if below the mean.
			for ( int i=0; i<unified_vtkout.scalars.size(); ++i ) {
				if ( unified_vtkout.scalars[i] < mean ) {
					unified_vtkout.scalars[i] = 0.0;
				}
			}
			unified_vtkout.points_3D[ slice(2, unified_vtkout.scalars.size(), 3) ] = unified_vtkout.scalars;
			
			{
				ostringstream name;
				name << "line_segments-1_over_mean_" << t << ".vtp";
				
				vtk_append_raw ( unified_vtkout, (name.str()).c_str() );
			}
			
			unified_vtkout.scalars = original_scalars;
			
			// Set to zero if over the mean.
			for ( int i=0; i<unified_vtkout.scalars.size(); ++i ) {
				if ( unified_vtkout.scalars[i] > mean ) {
					unified_vtkout.scalars[i] = 0.0;
				}
			}
			unified_vtkout.points_3D[ slice(2, unified_vtkout.scalars.size(), 3) ] = unified_vtkout.scalars;
			
			{
				ostringstream name;
				name << "line_segments-1_under_mean_" << t << ".vtp";
				
				vtk_append_raw ( unified_vtkout, (name.str()).c_str() );
			}
		}

	}
	
#if 0
	{
		pvd_output pvd;
		pvd.head = "line_segments-1_manifest_";
		ostringstream oss;
		
		pvd.middle.resize ( 5 );
		int numbers[] = { 0, 2, 51, 52, 71 };
		
		for ( int i=0; i<5; ++i ) {
			oss.str ( "" );
			oss << numbers[i];
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "line_segments-1_selection.pvd" );
		pvd_simple_output ( pvd, file );
	}
#endif
	
	{
		pvd_output pvd;
		pvd.head = "line_segments-1_manifest_";
		pvd.middle.resize ( iterations );
		ostringstream oss;
		for ( int i=0; i<iterations; ++i ) {
			oss.str ( "" );
			oss << i;
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "line_segments-1_manifest.pvd" );
		pvd_simple_output ( pvd, file );
	}
	
	{
		pvd_output pvd;
		pvd.head = "line_segments-1_over_mean_";
		pvd.middle.resize ( iterations );
		ostringstream oss;
		for ( int i=0; i<iterations; ++i ) {
			oss.str ( "" );
			oss << i;
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "line_segments-1_over_mean.pvd" );
		pvd_simple_output ( pvd, file );
	}
	
	{
		pvd_output pvd;
		pvd.head = "line_segments-1_under_mean_";
		pvd.middle.resize ( iterations );
		ostringstream oss;
		for ( int i=0; i<iterations; ++i ) {
			oss.str ( "" );
			oss << i;
			pvd.middle[i] = oss.str ();
		}
		pvd.tail = ".vtp";
		
		ofstream file ( "line_segments-1_under_mean.pvd" );
		pvd_simple_output ( pvd, file );
	}

	
	double    time_step_dt = 0.754;
	
	{
		for ( int i=0; i<avg_vals.values.size(); ++i ) {
			// The average would be found by dividing num_inside
			// but the normalized average would cancel out this factor.
			avg_vals.values[i] /= avg_vals.values[i].max();
		}
		ofstream file ( "line_segments-1_normalized_average_values.gnuplot" );
		
		file.precision(16);
		
		file << "plot [*:*] [0:1] 'line_segments-1_normalized_average_values.gnuplot' using 1:2 w l"
		     << " t \"" << avg_vals.labels[0] <<"\"";
		
		for ( int i=1; i<avg_vals.labels.size(); ++i ) {
			file << ", '' using 1:" << i+2 << " w l"
			     << " t \"" << avg_vals.labels[i] << "\"";
		}
		
		file << "\n";
		for ( int i=0; i<iterations; ++i ) {
			file << time_step_dt*i;
			for ( int j=0; j<avg_vals.values.size(); ++j ) {
				file << "\t" << (avg_vals.values[j])[i];
			}
			file << "\n";
		}
		file << "e\n";
	}
}
