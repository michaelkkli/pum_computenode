#include <mpi.h>

#if 0
#include <boost/version.hpp>

#if (BOOST_VERSION % 100 >= 1) && (BOOST_VERSION / 100 % 1000 >= 34)
#include <boost/filesystem.hpp>
#else
#include <boost/filesystem/operations.hpp>   // includes path.hpp
#include <boost/filesystem/convenience.hpp>
#endif
#endif

#include <sstream>
#include <iostream>
#include <cstdlib>

#if 0
using namespace boost::filesystem;
#endif

int comm_size, comm_rank;

int main( int argc, char* argv[] ) {

  MPI_Init( &argc, &argv );

  MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
  MPI_Comm_rank( MPI_COMM_WORLD, &comm_rank );

  char* tmpdir = std::getenv( "TMPDIR" );

  if ( !tmpdir ) {
    std::cout << "Rank " << comm_rank << ": no environment variable TMPDIR\n";
  } else {
    std::cout << "Rank " << comm_rank << ": TMPDIR is " << tmpdir << "\n";
  }

#if 0
  path tmpdir_path;
  if( exists( tmpdir ) && is_directory( tmpdir ) ) {
    tmpdir_path = tmpdir;
    tmpdir_path /= "computenode";
  } else {
    tmpdir_path = "/tmp/computenode";
  }

  if( !exists( tmpdir_path ) ) {
    create_directory( tmpdir_path );
  }
#endif

  std::stringstream ss;
  ss << comm_rank;

#if 0
  create_directory( tmpdir_path / ss.str() );
#endif

  MPI_Finalize();

}
