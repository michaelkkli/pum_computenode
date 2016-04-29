#include "vtk_output.hh"

#include <fstream>

using std::ofstream;

int main () {
	
	vtk_output m;
	m.points_3D.resize( 6 * 3 );
	
	// x
	// y
	// z
	
	//	0,1	1,1	2,1
	//	0,0	1,0	2,0
	
	m.points_3D[0] = 0;     m.points_3D[3] = 1;     m.points_3D[6] = 2;
	m.points_3D[1] = 1;     m.points_3D[4] = 1;     m.points_3D[7] = 1;
	m.points_3D[2] = 0;     m.points_3D[5] = 0;     m.points_3D[8] = 0;
	
	m.points_3D[ 9] = 0;     m.points_3D[12] = 1;     m.points_3D[15] = 2;
	m.points_3D[10] = 0;     m.points_3D[13] = 0;     m.points_3D[16] = 0;
	m.points_3D[11] = 0;     m.points_3D[14] = 0;     m.points_3D[17] = 0;
	
	m.scalars.resize( 6 );
	m.scalars[0] = 0;
	m.scalars[1] = 1;
	m.scalars[2] = 2;
	m.scalars[3] = 3;
	m.scalars[4] = 4;
	m.scalars[5] = 5;
	
	// 0 1 2
	// 3 4 5
	
	m.connectivity.resize( 8 );
	m.connectivity[0] = 3;
	m.connectivity[1] = 4;
	m.connectivity[2] = 1;
	m.connectivity[3] = 0;
	
	m.connectivity[4] = 4;
	m.connectivity[5] = 5;
	m.connectivity[6] = 2;
	m.connectivity[7] = 1;
	
	m.offsets.resize(2);
	m.offsets[0] = 4;
	m.offsets[1] = 8;
	
	{
		ofstream file ( "first_vtk.vtp" );
		vtk_simple_output( m, file );
	}
}
