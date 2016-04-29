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

#include "refinement_structure.hh"
#include "vtk_output.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <valarray>

// C99 int32_t, int63_t, uint32_t
// cstdint is experimental.
#include <stdint.h>

using std::copy;
using std::cout;
using std::ios;
using std::inserter;
using std::ofstream;
using std::ostringstream;
using std::slice;

vtk_output& vtk_output::operator= ( const vtk_output& other )
{
	this->points_3D.resize ( other.points_3D.size() );
	this->scalars.resize ( other.scalars.size() );
	this->connectivity.resize ( other.connectivity.size() );
	this->offsets.resize ( other.offsets.size() );
	
	this->points_3D    = other.points_3D;
	this->scalars      = other.scalars;
	this->connectivity = other.connectivity;
	this->offsets      = other.offsets;
}

bool little_endian ()
{
	int one = 1;
	char* p = reinterpret_cast<char*> ( &one );
	if ( p[0] ) {
		// Least significant bit in lowest address.
		return true;
	} else {
		return false;
	}
}

void base64 ( unsigned char * in, ostream & out, int len )
{
	if ( !little_endian() ) {
		abort ();
	}
#if 0
	static int count = 0;
	
	if ( count <1 ) {
		++count;
		char me[] = "Balloon." ;
		cout << "\"Self testing.\" is ";
		base64 ( reinterpret_cast<unsigned char*>(&me[0]), out, 8 );
		
		
	}
#endif
	
	char enc64[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	
	int len_incomplete_block = len % 3;
	int num_full_blocks      = (len - len_incomplete_block)/3;
	
	char w[4]; // Write block.
	
	unsigned char * s; // Start of block.
	for ( int i=0; i<num_full_blocks; ++i ) {
		if ( 0!=i && i%18==0 ) {
			out.put ( '\n' );
		}
		s = &in[3*i];
		
		// Right-shift 2 to extract first six
		// bits from 8-bit char.
		w[0] = enc64[ s[0] >> 2 ];
		
		//	cout << "First index is " << int(s[0] >> 2) << "\n";
		
		// Take highest two bits from first char
		// and four bits from second char.
		w[1] = enc64[ ( (0x03 & s[0]) << 4) | ((0xf0 & s[1]) >> 4 ) ];
		
		// Take four bits from second char and two bits
		// from third char.
		w[2] = enc64[ ((0x0f & s[1]) << 2) | ((0xc0 & s[2]) >> 6) ];

		w[3] = enc64[ (0x3f & s[2]) ];
		
		out.write ( w, 4 );
	}
	
	if ( len_incomplete_block > 0 ) {
		s += 3;
		if ( num_full_blocks % 18 == 0 ) {
			out.put ('\n');
		}
		w[0] = enc64[ s[0] >> 2 ];
		w[1] = enc64[ ( (0x03 & s[0]) << 4) | ((0xf0 & s[1]) >> 4 ) ];
		if ( 1 == len_incomplete_block ) {
			w[2] = '=';
			w[3] = '=';
		} else if ( 2 == len_incomplete_block ) {
			w[2] = enc64[ ((0x0f & s[1]) << 2) | ((0xc0 & s[2]) >> 6) ];
			w[3] = '=';
		}
		out.write ( w, 4 );
	}
}

void vtk_simple_output ( vtk_output& in, const string & filename, bool binary )
{
	ofstream out;
	if ( binary ) {
		out.open ( filename.c_str(), ios::binary | ios::out | ios::trunc );
	} else {
		out.open ( filename.c_str() );
	}
	
	int size_pts = in.points_3D.size();
	assert ( size_pts % 3 == 0 );
	
	int num_pts  = size_pts/3;
	assert ( num_pts == in.scalars.size() );
	assert ( 0<=in.connectivity.min() );
	assert ( in.connectivity.max() < num_pts );
	assert ( in.offsets.max() <= in.connectivity.size()+1 );
	
	int num_polys = in.offsets.size();

	int num_verts = 0;
#if 0
	if ( output_verts ) {
		num_verts = num_pts;
	}
#endif
	
	// HEAD
	
	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "\t<PolyData>\n"
		<< "\t\t<Piece NumberOfPoints=\"" << num_pts << "\" NumberOfVerts=\"" << num_verts << "\" NumberOfLines=\"0\"\n"
		<< "\t\t       NumberOfStrips=\"0\" NumberOfPolys=\"" << num_polys << "\">\n";
	
	// POINTS
	
	out << "\t\t\t<Points>\n";

	if ( binary ) {
		out << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">\n";
		uint32_t length = size_pts*sizeof( double );
		
		base64 ( reinterpret_cast<unsigned char*>( &length ), out, sizeof ( uint32_t ) );
		base64 ( reinterpret_cast<unsigned char*>(&(in.points_3D[0])), out, size_pts * sizeof( double ) );
		
	} else {
		out << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		
		// Using a ostringstream is not faster.
		for ( int i=0; i<size_pts; ++i ) {
			out << in.points_3D[i] << " ";
		}
	}


	
	out << "\n";
	
	out << "\t\t\t\t</DataArray>\n";
	out <<"\t\t\t</Points>\n";
	
#if 0
	if ( output_verts ) {
		out << "\t\t\t<Verts>\n";
		out << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
		for ( int i=0; i<num_pts; ++i ) {
			out << i << " ";
		}
		out << "\n";
		out << "\t\t\t\t</DataArray>\n";
		
		out << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
		for ( int i=0; i<num_pts; ++i ) {
			out << i+1 << " ";
		}
		out << "\n";
		out << "\t\t\t\t</DataArray>\n";
		out << "\t\t\t</Verts>\n";
	}
#endif
	
	
	// POINTDATA
	
	out << "\t\t\t<PointData Scalars=\"scalars\">\n";
	if (  binary ) {
		out << "\t\t\t\t<DataArray type=\"Float64\" Name=\"scalars\" format=\"binary\">\n";
		uint32_t length = in.scalars.size()* sizeof ( double );
		
		base64 ( reinterpret_cast<unsigned char*>( &length ), out, sizeof ( uint32_t ) );
		base64 ( reinterpret_cast<unsigned char*>(&(in.scalars[0])), out, in.scalars.size() * sizeof ( double ) );
		
#if 0
		out.write ( reinterpret_cast<char*>(&length), sizeof(uint32_t) );
		out.write ( reinterpret_cast<char*>( &(in.scalars[0]) ),
		                                       num_pts * sizeof ( double ) );
#endif
	} else {
		out << "\t\t\t\t<DataArray type=\"Float32\" Name=\"scalars\" format=\"ascii\">\n";
		
		// Using a ostringstream is not faster.
		for ( int i=0; i<num_pts; ++i ) {
			out << in.scalars[i] << " ";
		}
	}
	
	out << "\n";
	out << "\t\t\t\t</DataArray>\n";
	
	out << "\t\t\t</PointData>\n";
	
	// POLYS
	
	out << "\t\t\t<Polys>\n";
	
	
	int num_connectivity = in.connectivity.size();
	if ( binary ) {
		out << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n";
		uint32_t length = num_connectivity* sizeof(int);
		base64 ( reinterpret_cast<unsigned char*>(&length), out, sizeof(uint32_t) );
		base64 ( reinterpret_cast<unsigned char*>(&(in.connectivity[0])),
		         out, num_connectivity * sizeof(int) );
	} else {
		out << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
		for ( int i=0; i<num_connectivity; ++i ) {
			out << in.connectivity[i] << " ";
		}
	}
	out << "\n";
	out << "\t\t\t\t</DataArray>\n";

	int num_offsets = in.offsets.size();
	if ( binary ) {
		out << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n";
		uint32_t length = num_offsets*sizeof(int);

		base64 ( reinterpret_cast<unsigned char*>(&length), out, sizeof ( uint32_t ) );
		base64 ( reinterpret_cast<unsigned char*>(&(in.offsets[0]) ),
		         out,
		         num_offsets*sizeof(int) );
	} else {
		out << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
		for ( int i=0; i<num_offsets; ++i ) {
			out << in.offsets[i] << " ";
		}
	}
	out << "\n";
	out << "\t\t\t\t</DataArray>\n";
	
	out << "\t\t\t</Polys>\n";
	
	out << "\t\t</Piece>\n"
		<< "\t</PolyData>\n"
		<< "</VTKFile>";
}

// double is size 8 hence 64 bit.
// int is size 4 hence 32 bit.
// uint32_t number of occurrences of datatype
void vtk_append_raw ( vtk_output& in, const string & filename, bool do_cell_scalars )
{
	ofstream out ( filename.c_str(), ios::binary | ios::out | ios::trunc );

	
	int size_pts = in.points_3D.size();
	assert ( size_pts % 3 == 0 );
	
	int num_pts  = size_pts/3;
	assert ( num_pts == in.scalars.size() );
	assert ( 0<=in.connectivity.min() );
	assert ( in.connectivity.max() < num_pts );
	assert ( in.offsets.max() <= in.connectivity.size()+1 );
	
	int num_polys = in.offsets.size();
	
	vector<int> append_offsets;
	
	if ( do_cell_scalars ) {
		append_offsets.resize (5);
	} else {
		// Append scalars, points, connectivity, offsets.
		append_offsets.resize (4);
	}
	
#if 0 // Assuming number of occurrence of datatype.
	append_offsets[0] = 0;
	append_offsets[1] = num_pts + sizeof ( uint32_t );
	append_offsets[2] = append_offsets[1] + size_pts + sizeof ( uint32_t );
	append_offsets[3] = append_offsets[2] + in.connectivity.size() + sizeof ( uint32_t );
#else // Assuming number of chars
	append_offsets[0] = 0;
	append_offsets[1] = num_pts*sizeof(double) + sizeof ( uint32_t );
	append_offsets[2] = append_offsets[1] + size_pts*sizeof(double) + sizeof ( uint32_t );
	append_offsets[3] = append_offsets[2] + in.connectivity.size()*sizeof(int) + sizeof ( uint32_t );
	
	if ( do_cell_scalars ) {
		// Remember number of chars is of the prevous data piece.
		append_offsets[4] = append_offsets[3] +in.offsets.size()*sizeof(int) + sizeof(uint32_t);
	}
#endif
	
	// HEAD
	
	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	out << "\t<PolyData>\n"
		<< "\t\t<Piece NumberOfPoints=\"" << num_pts << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\"\n"
		<< "\t\t       NumberOfStrips=\"0\" NumberOfPolys=\"" << num_polys << "\">\n";
	
	out << "\t\t\t<PointData Scalars=\"scalars\">\n";
	out << "\t\t\t\t<DataArray type=\"Float64\" Name=\"scalars\" format=\"appended\" offset=\""
		<< append_offsets[0] << "\"/>\n";
	out << "\t\t\t</PointData>\n";
	
	out << "\t\t\t<Points>\n";
	out << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\""
		<< append_offsets[1] << "\"/>\n";
	out << "\t\t\t</Points>\n";
	
	out << "\t\t\t<Polys>\n";
	out << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""
		<< append_offsets[2] << "\"/>\n";
	out << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""
		<< append_offsets[3] << "\"/>\n";
	out << "\t\t\t</Polys>\n";
	
	if ( do_cell_scalars ) {
		assert ( in.cell_scalars.size() == num_polys );
		
		out << "\t\t\t<CellData Scalars=\"cell_scalars\">\n";
		out << "\t\t\t\t<DataArray type=\"Float64\" Name=\"cell_scalars\" format=\"appended\" offset=\""
			<< append_offsets[4] << "\"/>\n";
		out << "\t\t\t</CellData>\n";
	}
	
	out << "\t\t</Piece>\n";
	out << "\t</PolyData>\n";
	
	out << "\t<AppendedData encoding=\"raw\">\n";
	
	out << "_";
	
	uint32_t length;
	
	// Write out scalars.
	length  = num_pts * sizeof ( double );
	out.write ( reinterpret_cast<char*>( &length ), sizeof( uint32_t ) );
	out.write ( reinterpret_cast<char*>( &(in.scalars[0]) ),
	            num_pts * sizeof ( double ) );
	
	// Write out points.
	length = size_pts * sizeof( double );
	out.write ( reinterpret_cast<char*>( &length ), sizeof( uint32_t ) );
	out.write ( reinterpret_cast<char*>(&(in.points_3D[0])),
	            size_pts * sizeof( double ) );
	
	length = in.connectivity.size() * sizeof ( int );
	out.write ( reinterpret_cast<char*>( &length ), sizeof( uint32_t ) );
	out.write ( reinterpret_cast<char*>( &(in.connectivity[0]) ),
	            in.connectivity.size() * sizeof ( int ) );
	
	length = in.offsets.size() * sizeof ( int );
	out.write ( reinterpret_cast<char*>( &length ), sizeof( uint32_t ) );
	out.write ( reinterpret_cast<char*>( &(in.offsets[0]) ),
	            in.offsets.size() * sizeof ( int ) );
	
	if ( do_cell_scalars ) {
		length = in.cell_scalars.size() * sizeof ( double );
		out.write ( reinterpret_cast<char*>( &length ), sizeof ( uint32_t ) );
		out.write ( reinterpret_cast<char*>( &(in.cell_scalars[0]) ),
		            in.cell_scalars.size() * sizeof ( double ) );
	}
	
	
	out << "\t</AppendedData>\n";
	out << "</VTKFile>";
}

void pvd_simple_output ( const pvd_output & in, ofstream & out )
{
	out << "<?xml version=\"1.0\"?>\n";
	out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
	out << "\t<Collection>\n";
	
	int num = in.middle.size();
	
	for ( int i=0; i<num; ++i ) {
		out << "\t\t<DataSet timestep=\""
			<< i << "\" group=\"\" part=\"0\" file=\""
			<< in.head << in.middle[i] << in.tail
			<< "\"/>\n";
	}
	out << "\t</Collection>\n";
	out << "</VTKFile>";
}

void get_entry_connectivity_offsets ( vector<string> &         bdry_keys,
                                            vector<box<2> > &        boxes,
                                            vector<int> &            entry_seg,
                                            basic_dbinary_tree<2> &  dtree,
                                            line_segments &          ls,
                                            vector<double> &         first_entry,
                                            vector<int> &            connectivity,
                                            vector<int> &            offsets,
                                            int                      rebase_bdry,
                                            int                      rebase_ent,
                                            visualize_choice         viz_choice )
{
	int num = bdry_keys.size();
	assert ( boxes.size() == num );
	assert ( entry_seg.size() == num );
	
	int num_bdry_pts_inside;
	double entry_pt[2];
	double exit_pt[2];
	vector<int>    in_c;  // Inside corners.
	vector<int>    out_c; // Outside corners.
	
	// Grid indices.
	int    grid_ind[4];
	
	assert ( ls.segments.size() % 2 == 0 );
	int num_segs = (ls.segments.size()/2)-1;
	int tmp_next_seg; // To allow wrap around of segment indices.
	
	typedef vector<int>::iterator    vi_t;
	vi_t                             vi_it, vi_end;
	
	first_entry.clear();
	offsets.clear();
	
	// Must erase this at the end.
	// Used to allow taking a running sum.
	offsets.push_back(0);
	
	// Debugging found offsets[i] was incorrectly assumed
	// to be the back element but multiple push_backs are
	// used at each iteration of the for-loop so back()
	// is required.
	// 2008-12-07 Michael LI.
	
	int in_c_size, out_c_size;
	for ( int i=0; i<num; ++i ) {
	
		split_by_line_segments( boxes[i],
		                        ls,
		                        entry_seg[i],
		                        num_bdry_pts_inside,
		                        entry_pt,
		                        exit_pt,
		                        in_c,
		                        out_c );

		in_c_size  = in_c.size();
		out_c_size = out_c.size();

		first_entry.push_back (entry_pt[0]);
		first_entry.push_back (entry_pt[1]);
		
		dtree.get_box_grid_indices_winding( bdry_keys[i], grid_ind );
		
		tmp_next_seg = entry_seg[i] + 1;
		if ( tmp_next_seg >= num_segs ) {
			tmp_next_seg -= num_segs;
		}
		
#if 0 // Mark for delete. 2008-12-08 ML.
		cout << "num_bdry_pts_inside is " << num_bdry_pts_inside << "\n";
		cout << "in_c.size () is " << in_c.size()
		     << ", and out_c.size() is " << out_c.size() << "\n";

#endif
	//	cout << "Isolation\n";
	//	if ( i!=0 ) continue;

		
		
		
		if ( num_bdry_pts_inside == 0 ) {
			
			if ( (viz_choice==inside_outside)||(viz_choice==inside) ) {
				// Do inside.
				connectivity.push_back( i + rebase_ent );
				
				// Index wrap around.
				connectivity.push_back( ((i+1 < num)? i+1 : 0) + rebase_ent );
				
				vi_it  = in_c.begin();
				vi_end = in_c.end();
				for ( ; vi_it!=vi_end; ++vi_it ) {
					// No rebasing necessary for grid points.
					// (Until we have multiple sets of grid pts.)
					connectivity.push_back( grid_ind[ *vi_it ] );
				}
				
				assert ( in_c_size >= 0 );
				offsets.push_back ( offsets.back() + 2 + in_c_size );
			}
			
			if ( (viz_choice==inside_outside)||(viz_choice==outside) ) {
				// Do outside.
				connectivity.push_back( ((i+1 < num)? i+1 : 0) + rebase_ent );
				connectivity.push_back( i + rebase_ent );
				vi_it  = out_c.begin();
				vi_end = out_c.end();
				for ( ; vi_it!=vi_end; ++vi_it ) {
					// No rebasing necessary for grid points.
					// (Until we have multiple sets of grid pts.)
					connectivity.push_back( grid_ind[ *vi_it ] );
				}
				
				assert ( out_c_size >= 0 );
				offsets.push_back ( offsets.back() + 2 + out_c_size );
			}
			// Unnecessary continue statement found here.
		} else if ( (in_c_size==0) || (out_c_size==0) ) {
			// This may only happen when entry and exit are on the same edge.
			
			vector<int>    inside_bdry_ind;
			int bdry_ind = entry_seg[i];
			for ( int in_b=0; in_b<num_bdry_pts_inside; ++in_b ) {
				++bdry_ind;
				if ( bdry_ind > num_segs-1 ) {
					bdry_ind -= num_segs;
				}
				inside_bdry_ind.push_back( bdry_ind + rebase_bdry );
			}
			
			if ( (viz_choice==inside_outside)||(viz_choice==inside) ) {
				// Do inside.
				connectivity.push_back( i + rebase_ent );
				copy ( inside_bdry_ind.begin(),
				inside_bdry_ind.end(),
				back_inserter( connectivity ) );
				connectivity.push_back( ((i+1 < num)? i+1 : 0) + rebase_ent );
				vi_it  = in_c.begin();
				vi_end = in_c.end();
				for ( ; vi_it!=vi_end; ++vi_it ) {
					// No rebasing necessary for grid points.
					// (Until we have multiple sets of grid pts.)
					connectivity.push_back( grid_ind[ *vi_it ] );
				}
				
				assert ( in_c_size >= 0 ); assert ( num_bdry_pts_inside >= 0 );
				offsets.push_back ( offsets.back() + 2 + in_c_size + num_bdry_pts_inside );
			}
			
			if ( (viz_choice==inside_outside)||(viz_choice==outside) ) {
				// Do outside.
				connectivity.push_back( ((i+1 < num)? i+1 : 0) + rebase_ent );
				// Reverse inside boundary points.
				copy ( inside_bdry_ind.rbegin(),
				inside_bdry_ind.rend(),
				back_inserter( connectivity ) );
				connectivity.push_back( i + rebase_ent );
				vi_it  = out_c.begin();
				vi_end = out_c.end();
				for ( ; vi_it!=vi_end; ++vi_it ) {
					// No rebasing necessary for grid points.
					// (Until we have multiple sets of grid pts.)
					connectivity.push_back( grid_ind[ *vi_it ] );
				}
				
				assert ( out_c_size >= 0 ); assert ( num_bdry_pts_inside >= 0 );
				offsets.push_back ( offsets.back() + 2 + out_c_size + num_bdry_pts_inside );
			}
		} else {
			// This is the general case of entry and exit on distinct edges
			vector<int>    inside_bdry_ind;
			int bdry_ind = entry_seg[i];
			for ( int in_b=0; in_b<num_bdry_pts_inside; ++in_b ) {
				++bdry_ind;
				// Debugging found missing -1 here.
				// 2008-12-08 Michael LI.
				if ( bdry_ind > num_segs-1 ) {
					bdry_ind -= num_segs;
				}
				inside_bdry_ind.push_back( bdry_ind + rebase_bdry );
			}
			if ( (viz_choice==inside_outside)||(viz_choice==inside) ) {
				// Inside.
				connectivity.push_back( i + rebase_ent );
				copy ( inside_bdry_ind.begin(),
				inside_bdry_ind.end(),
				back_inserter( connectivity ) );
				connectivity.push_back( ((i+1 < num)? i+1 : 0) + rebase_ent );
				vi_it  = in_c.begin();
				vi_end = in_c.end();
				for ( ; vi_it!=vi_end; ++vi_it ) {
					// No rebasing necessary for grid points.
					// (Until we have multiple sets of grid pts.)
					connectivity.push_back( grid_ind[ *vi_it ] );
				}
				
				assert ( in_c_size >= 0 ); assert ( num_bdry_pts_inside >= 0 );
				offsets.push_back ( offsets.back() + 2 + in_c_size + num_bdry_pts_inside );
			}
			if ( (viz_choice==inside_outside)||(viz_choice==outside) ) {
				// Outside.
				connectivity.push_back( ((i+1 < num)? i+1 : 0) + rebase_ent );
				copy ( inside_bdry_ind.rbegin(),
				inside_bdry_ind.rend(),
				back_inserter( connectivity ) );
				connectivity.push_back( i + rebase_ent );
				vi_it  = out_c.begin();
				vi_end = out_c.end();
				for ( ; vi_it!=vi_end; ++vi_it ) {
					// No rebasing necessary for grid points.
					// (Until we have multiple sets of grid pts.)
					connectivity.push_back( grid_ind[ *vi_it ] );
				}
				
				assert ( out_c_size >= 0 ); assert ( num_bdry_pts_inside >= 0 );
				offsets.push_back ( offsets.back() + 2 + out_c_size + num_bdry_pts_inside );
			}
		}
	}

	offsets.erase ( offsets.begin() );

#ifndef NDEBUG
//	cout << "offsets is ";
	for ( int i=0; i<offsets.size(); ++i ) {
		if ( i!=0 ) {
			assert ( offsets[i-1] < offsets[i] );
//			cout << ", ";
		}
//		cout << offsets[i];
	}
#endif
}

void make_boundary_vtk_output ( line_segments &               ls,
                                refinement_structure<2> &     rs,
                                vtk_output &                  vtkout,
                                visualize_choice              viz_choice,
                                vector<string> *                 opt_bdry_keys_out,
                                vector<int> *                 opt_entry_set_out )
{

	
	vector<string>    bdry_keys;
	vector<int>       entry_seg;
	
	get_intersect_keys_entry_indices ( rs, ls, bdry_keys, entry_seg );
	
	if ( opt_bdry_keys_out ) {
		opt_bdry_keys_out->clear();
		opt_bdry_keys_out->reserve ( bdry_keys.size() );
		*opt_bdry_keys_out = bdry_keys;
/*		copy ( bdry_keys.begin(), bdry_keys.end(),
		       inserter( *opt_bdry_keys_out, opt_bdry_keys_out->begin() ) );*/
	}
	
	if ( opt_entry_set_out ) {
		opt_entry_set_out->clear();
		opt_entry_set_out->reserve ( entry_seg.size() );
		*opt_entry_set_out = entry_seg;
	}
	
	vector<box<2> > boxes;
	(rs.dtree)->get_box ( bdry_keys, boxes );
	
	valarray<double>  grid_pts;
	vector<double>    entry_pts;
	vector<int>       bdry_connectivity;
	vector<int>       bdry_offsets;
	
	assert ( rs.dtree );
	(rs.dtree)->get_grid_points_3D ( grid_pts );
	
	int num_grid_pts   = grid_pts.size()/3;
	int num_bdry_pts   = (ls.segments.size()/2)-1;
	
	get_entry_connectivity_offsets ( bdry_keys,
	                                 boxes,
	                                 entry_seg,
	                                 *(rs.dtree),
	                                 ls,
	                                 entry_pts,
	                                 bdry_connectivity,
	                                 bdry_offsets,
	                                 num_grid_pts,
	                                 num_grid_pts+num_bdry_pts,
	                                 viz_choice );

	int num_entry_pts  = entry_pts.size()/2;
	
	int tnum_pts = num_grid_pts + num_bdry_pts + num_entry_pts;
	
	vtkout.points_3D.resize( 3*tnum_pts, 0.0 );
	
	copy ( &grid_pts[0], &grid_pts[0] + grid_pts.size(), &vtkout.points_3D[0] );
	
	// x coord of boundary points.
	vtkout.points_3D[ slice( 3*num_grid_pts    , num_bdry_pts, 3 ) ]
		= ls.segments[ slice ( 0, num_bdry_pts, 2 ) ];
	// y coord of boundary points.
	vtkout.points_3D[ slice( 3*num_grid_pts + 1, num_bdry_pts, 3 ) ]
		= ls.segments[ slice ( 1, num_bdry_pts, 2 ) ];
	
	int    grid_bdry_offset = 3*(num_grid_pts+num_bdry_pts);
	for ( int i=0; i<num_entry_pts; ++i ) {
		vtkout.points_3D[ grid_bdry_offset+3*i  ] = entry_pts[ 2*i   ];
		vtkout.points_3D[ grid_bdry_offset+3*i+1] = entry_pts[ 2*i+1 ];
	}
	
	vtkout.scalars.resize ( tnum_pts, -16.0 ); // Test value.

	vtkout.connectivity.resize( bdry_connectivity.size() );
	copy ( bdry_connectivity.begin(),
	       bdry_connectivity.end(), &(vtkout.connectivity[0]) );

	vtkout.offsets.resize( bdry_offsets.size() );
	copy ( bdry_offsets.begin(),
	       bdry_offsets.end(), &(vtkout.offsets[0]) );
}
