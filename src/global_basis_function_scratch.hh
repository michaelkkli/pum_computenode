#ifndef _GLOBAL_BASIS_FUNCTION_SCRATCH_HH_
#define _GLOBAL_BASIS_FUNCTION_SCRATCH_HH_

#include "pu_function_scratch.hh"

#include <valarray>
using std::valarray;

struct global_basis_function_scratch {
	valarray<double>    pu_grad;
	valarray<double>    local_approx_grad;
	
	valarray<double>    pu_vals;
	valarray<double>    local_approx_vals;
	
	pu_function_scratch    pu_fun_scratch;
};

#endif // _GLOBAL_BASIS_FUNCTION_SCRATCH_HH_
