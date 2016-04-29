#include <petscksp.h>

#include <iostream>
using std::cout;

static char help[] = "solution-1";

int main( int argc, char* argv[] ) {
	PetscInitialize( &argc, &argv, static_cast<char*>(0), help );
	
	PetscErrorCode ierr;
	PetscTruth flg;
	
	PetscMPIInt comm_size, comm_rank;

	ierr = MPI_Comm_size( PETSC_COMM_WORLD, &comm_size ); CHKERRQ(ierr);
	ierr = MPI_Comm_rank( PETSC_COMM_WORLD, &comm_size ); CHKERRQ(ierr);
	
	int matrix_size = 5;
	PetscOptionsGetInt( PETSC_NULL, "-matrix_size", &matrix_size, &flg );

	Vec global_vec;

	ierr = VecCreate(PETSC_COMM_WORLD,&global_vec); CHKERRQ(ierr);
	ierr = VecSetSizes(global_vec,PETSC_DECIDE, matrix_size);CHKERRQ(ierr);
	ierr = VecSetFromOptions(global_vec);CHKERRQ(ierr);
	
	Vec local_vec;
	VecScatter ctx;

	// Creates the target vector.
	VecScatterCreateToAll( global_vec, &ctx, &local_vec );

	VecScatterBegin( ctx, global_vec, local_vec, INSERT_VALUES, SCATTER_FORWARD );
	VecScatterEnd( ctx, global_vec, local_vec, INSERT_VALUES, SCATTER_FORWARD );

	VecScatterDestroy(ctx);
	VecDestroy(local_vec);
	VecDestroy(global_vec);
	
	ierr = PetscFinalize(); CHKERRQ(ierr);
}
