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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <valarray>

#ifndef eF77
#define eF77(s) s##_
#endif

using std::cout;
using std::exp;
using std::for_each;
using std::stringstream;
using std::valarray;

extern "C" void eF77(dstevd) ( char* jobz,
                    int*  n,
                    double* d,
                    double* e,
                    double* z,
                    int*    ldz,
                    double* work,
                    int*    lwork,
                    int*    iwork,
                    int*    liwork,
                    int*    info );

void square ( double& x ) {
	x *= x;
}

void natural_exp ( double& x ) {
	x = exp(x);
}

int main ( int argc, char* argv[] ) {
	char jobz = 'V';

	int N = 10;

	if ( argc == 2 ) {
		stringstream ss;
		ss << argv[1];
		ss >> N;
		if ( !ss ) {
			cout << "Illegal argument given!\n";
			abort();
		}
	}

	int n = N; // Allow taking of non-const pointer.

	valarray<double>    d(0.0,N);
	valarray<double>    e(0.5,N-1);
	// Notice the 0.5 initialization so only division remains.
	for ( int i=0; i<N-1; ++i ) {
		e[i] /= sqrt( 1. - 1./(4*(i+1)*(i+1)) ); // Notice the +1 to account for indices.
	}

	int ldz = N;
	valarray<double>    z( ldz*N );

	int lwork = 1 + 4*N + N*N;
	valarray<double>    work( lwork );

	int liwork = 3 + 5*N;
	valarray<int>    iwork( liwork );

	int info;

	eF77(dstevd) ( &jobz, &n, &d[0], &e[0], &z[0], &ldz, &work[0], &lwork, &iwork[0], &liwork, &info );

	if ( info == 0 ) {
		cout << "Successful solve.\n";
	} else if ( info < 0 ) {
		cout << "Argument " << info << " had an illegal value.\n";
	} else {
		cout << info << " off-diagonal elements of E failed to converge to zero.\n";
	}
	cout << "Eigenvalues/points are:\n";
	for ( int i=0; i<ldz; ++i ) {
		cout << d[i] << "\n";
	}

	valarray<double>    quad_weights(N);
	for ( int i=0; i<N; ++i ) {
		quad_weights[i] = z[i*N];
	}
	quad_weights *= quad_weights;
	quad_weights *= 2.0;

	cout << "Weights are :\n";
        for ( int i=0; i<N; ++i ) {
                cout << quad_weights[i] << "\n";
        }

	{
		cout.precision ( 15 );

		cout << "Test integration.\n";
		valarray<double>    eval ( d );
		for_each ( &eval[0], &eval[0]+eval.size(), &natural_exp );
		eval *= quad_weights;
		cout << "Integral is " << eval.sum() << "\n";
	}
	
}

