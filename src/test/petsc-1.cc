#include <petscksp.h>

static char help[] = "petsc-1";

#undef __FUNCT__
#define __FUNCT__ "main"
int main ( int argc, char* argv[] ) {

	double abserr = 1e-4;
	int num_levels = 10;

	PetscErrorCode ierr;
	PetscTruth flg;

	PetscInitialize( &argc, &argv, static_cast<char*>(0), help );
	
	PetscOptionsGetReal( PETSC_NULL, "-abserr", &abserr, &flg );
	PetscOptionsGetInt( PETSC_NULL, "-num_levels", &num_levels, &flg );
	
	PetscPrintf( PETSC_COMM_WORLD, "abserr is %e\n", abserr );
	PetscPrintf( PETSC_COMM_WORLD, "num_levels is %i\n", num_levels );
	
	double tmp = abserr;
	for ( int i=0; i<num_levels; ++i ) {
		tmp = 0.25*tmp;
		PetscPrintf( PETSC_COMM_WORLD, "Level %i: %e\n", i, tmp );
	}
	
	ierr = PetscFinalize(); CHKERRQ(ierr);
}
