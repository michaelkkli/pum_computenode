#include <petscksp.h>

#ifndef NDEBUG
#define petsc_assert(x) if(!(x)) SETERRQ(1,1,"petsc assert failed\n";
#else
#define petsc_assert(x)
#endif

static char help[] = "petsc-2";

#undef __FUNCT__
#define __FUNCT__ "main"
int main (int argc, char* argv[]) {
  PetscMPIInt comm_size, comm_rank;
  PetscErrorCode ierr;

  PetscInitialize( &argc, &argv, static_cast<char*>(0), help );

  MPI_Comm_rank( PETSC_COMM_WORLD, &comm_rank );
  MPI_Comm_size( PETSC_COMM_WORLD, &comm_size );

  PetscPrintf( PETSC_COMM_WORLD,
	       "There are %i compute nodes.\n", comm_size );

  PetscPrintf( PETSC_COMM_SELF,
	       "There are %i compute nodes and this is rank %i.\n",
	       comm_size, comm_rank ); // Repeat of size for paranoid check.

  PetscInt m1 = comm_rank;

  ierr = PetscOptionsGetInt(PETSC_NULL,"-m", &m1, PETSC_NULL); CHKERRQ(ierr);

  // Until we figure out how to abort a petsc program or `assert' properly.
  if ( m1 < 0 ) {
    m1 = comm_rank;
  }

  PetscPrintf(PETSC_COMM_SELF,
	      "Rank %i: my m1 is %i.\n",comm_rank, m1 );

  Vec vec1;

  ierr = VecCreate(PETSC_COMM_WORLD, &vec1); CHKERRQ(ierr);
  ierr = VecSetSizes(vec1,m1, PETSC_DECIDE); CHKERRQ(ierr);
  ierr = VecSetFromOptions(vec1); CHKERRQ(ierr);

  PetscInt global_size1;
  VecGetSize( vec1, &global_size1 );
  PetscPrintf(PETSC_COMM_WORLD,
	      "Global size of vec1 is %i.\n",
	      global_size1 );


  PetscInt size1;
  VecGetLocalSize( vec1, &size1 );
  PetscPrintf(PETSC_COMM_SELF,
	      "Rank %i: local size of vec1 is %i.\n",
	      comm_rank,
	      size1 );

  Mat mat1;
  MatCreate(PETSC_COMM_WORLD, &mat1);
  MatSetSizes(mat1, size1, global_size1,
	      PETSC_DECIDE, PETSC_DECIDE);
  ierr = MatSetFromOptions(mat1); CHKERRQ(ierr);

  PetscInt low1, high1;

  VecGetOwnershipRange( vec1, &low1, &high1 );
  PetscPrintf(PETSC_COMM_SELF,
	      "Rank %i: ownership range [%i,%i).\n",
	      comm_rank,
	      low1,
	      high1 );

  PetscScalar* array1;

  VecGetArray( vec1, &array1 );

  for (int i=0; i<size1; ++i ) {
    array1[i] = comm_rank;
  }

  VecRestoreArray( vec1, &array1 );

  VecAssemblyBegin(vec1);
  VecAssemblyEnd(vec1);

  MatAssemblyBegin(mat1,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat1,MAT_FINAL_ASSEMBLY);

  ierr = MatDestroy(mat1); CHKERRQ(ierr);
  ierr = VecDestroy(vec1); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);

}
