#ifndef _COMPUTENODE_HH_
#define _COMPUTENODE_HH_

#include <mpi.h>

#include "box.hh"

#include <mpi.h>

#include <vector>
using std::vector;

template<int dim>
class computenode {
public:
  computenode( MPI_Comm = MPI_COMM_WORLD );
  ~computenode();
public:

private:
  MPI_Comm comm;
  int rank;
};

#endif // _COMPUTENODE_HH_
