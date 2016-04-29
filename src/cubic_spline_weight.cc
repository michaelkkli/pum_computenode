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

#include "cubic_spline_weight.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>

#undef TRACE
#define TRACE cout << __FILE__ << ":" << __LINE__ << "\n";

using std::abs;
using std::cout;
using std::fabs;
using std::slice;

template<int dim>
cubic_spline_weight<dim>::cubic_spline_weight()
{
	this->set_local_function();
}

template<int dim>
cubic_spline_weight<dim>::~cubic_spline_weight()
{
}

template<int dim>
double
cubic_spline_weight<dim>::evaluate ( const double* co ) const
{
	assert( co );
	
	double factor;
	double res = 1.0;
	
	for ( int d=0; d<dim; ++d ) {
		
		double v = co[d]<0.0 ? -co[d] : co[d];
		
		assert( v < 1.0 + 1e-10 );
		
		if ( v < 0.5 + 1e-10 ) {
			// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
			factor = 2.0/3.0 + 4.0*(v-1.0)*v*v;
			assert( factor > -1e-10 );
		} else if ( v < 1.0 + 1e-10 ) {
			// 4/3-4v+4v^2-(4/3)v^3
			//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
			//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
			//  = (4/3)( 1. -3v( 1. -v( 1. - v/3.) ) )
			//  = (4/3)( 1. -v( 3. -v( 3. - v) ) )
			factor = (4.0/3.0)*(1. - v*(3. - v*(3.-v)));
			assert( factor > -1e-10 );
		} else if ( v > 1.0 + 1e-10 ){
			factor = 0.0;
			assert( factor > -1e-10 );
		}
		res *= factor;
	}
	return res;
}

template<int dim>
void
cubic_spline_weight<dim>::evaluate ( valarray<double>& in,
                                     valarray<double>& out ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	out.resize( num_pts );
	
	for ( int i=0; i<num_pts; ++i ) {
		double* co = &in[dim*i];

		out[i] = 1.0;
		for ( int d=0; d<dim; ++d ) {
			
			// GDB found negative weights being returned and
			// discovered absolute values were needed so put here
			// for efficiency.
			double v = co[d]<0.0 ? -co[d] : co[d];
			
			assert( v <= 1.0 );
			
			if ( v < 0.5 + 1e-10 ) {
				// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
				out[i] *= 2.0/3.0 + 4.0*(v-1.0)*v*v;
			} else if ( v < 1.0 + 1e-10 ) {
				// 4/3-4v+4v^2-(4/3)v^3
				//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
				//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
				//  = (4/3)( 1. -3v( 1. -v( 1. - v/3.) ) )
				//  = (4/3)( 1. -v( 3. -v( 3. - v) ) )
				out[i] *= (4.0/3.0)*(1. - v*(3. - v*(3.-v)));
			} else if ( v > 1.0 ) {
				out[i] = 0.0;
				break;
			}
		}
		assert ( out[i] < pow(2./3., dim) + 1e-10 );
	}
	
	// Internal consistency testing.
#if 0//ndef NDEBUG
	for ( int i=0; i<num_pts; ++i ) {
		assert ( fabs ( out[i] - this->evaluate ( &in[dim*i] ) ) < 1e-10 );
	}
#endif
}


template<int dim>
void
cubic_spline_weight<dim>::evaluate ( valarray<double>& in, valarray<bool>& pred, valarray<double>& out ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	out.resize( num_pts );
	
	assert ( pred.size() == num_pts );
	
	for ( int i=0; i<num_pts; ++i ) {
		double* co = &in[dim*i];
		
		out[i] = 1.0;
		if ( pred[i] ) {
			for ( int d=0; d<dim; ++d ) {
				
				// GDB found negative weights being returned and
				// discovered absolute values were needed so put here
				// for efficiency.
				double v = co[d]<0.0 ? -co[d] : co[d];
				
				assert( v <= 1.0 );
				
				if ( v < 0.5 + 1e-10 ) {
					// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
					out[i] *= 2.0/3.0 + 4.0*(v-1.0)*v*v;
				} else if ( v < 1.0 + 1e-10 ) {
					// 4/3-4v+4v^2-(4/3)v^3
					//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
					//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
					//  = (4/3)( 1. -3v( 1. -v( 1. - v/3.) ) )
					//  = (4/3)( 1. -v( 3. -v( 3. - v) ) )
					out[i] *= (4.0/3.0)*(1. - v*(3. - v*(3.-v)));
				} else if ( v > 1.0 ) {
					out[i] = 0.0;
					break;
				}
			}
			assert ( out[i] < pow(2./3., dim) + 1e-10 );
		} else {
			out[i] = 0.0;
		}
	}
	
	// Internal consistency testing.
#if 0//ndef NDEBUG
	for ( int i=0; i<num_pts; ++i ) {
		if ( pred[i] ) {
			assert ( fabs ( out[i] - this->evaluate ( &in[dim*i] ) ) < 1e-10 );
		} else {
			assert ( fabs ( out[i] ) < 1e-10 );
		}
	}
#endif
}


template<int dim>
void
cubic_spline_weight<dim>::evaluate_grad ( const double* co, double* grad ) const
{
	assert( co && grad );
	double res = 0.0;
	for ( int d=0; d<dim; ++d ) {
		res = 1.0;
		for ( int i=0; i<dim; ++ i ) {
			double v = co[i]<0.0 ? -co[i] : co[i];
			if ( i == d ) {
				assert( v <= 1.0 );
				
				// Cubic spline derivative.
				if ( v < 0.5 + 1e-10 ) {
				//	// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
				
					// Derivative -8v + 12v^2 = 4(3v-2)v
					res *= 4.0*(3.0*v-2.0)*v;
					if ( co[i] < 0.0 ) {
						res *= -1.0;
					}
#if 0 // Would be bad if both branches are executed. 2008-12-04 Michael LI.
					if ( co[i] < 0.0 ) {
						res *= -4.0*(3.0*v-2.0)*v;
					} else {
						res *=  4.0*(3.0*v-2.0)*v;
					}
#endif
				} else if ( v < 1.0 +1e-10 ) {
					// 4/3-4v+4v^2-(4/3)v^3
					//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
					//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
					
					// Derivative -4 + 8v -4v^2 = -4(1-2v+v^2)
					//                          = -4(1+(-2+v)v)
					res *= -4.0*(1.0+(-2.0+v)*v);
					if ( co[i] < 0.0 ) {
						res *= -1.0;
					}
#if 0  // Would be bad if both branches are executed. 2008-12-04 Michael LI.
					if ( co[i] < 0.0 ) {
						res *=  4.0*(1.0+(-2.0+v)*v);
					} else {
						res *= -4.0*(1.0+(-2.0+v)*v);
					}
#endif
				}
#if 0//ndef NDEBUG
				else {
					assert ( false );
					res = 0.0;
					break;
				}
#endif
			} else {
				// Cubic spline function.
				if ( v < 0.5 + 1e-10 ) {
					// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
					res *= 2.0/3.0 + 4.0*(v-1.0)*v*v;
				} else if ( v < 1.0 + 1e-10 ) {
					// 4/3-4v+4v^2-(4/3)v^3
					//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
					//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
					//  = (4/3)( 1. -3v( 1. -v( 1. - v/3.) ) )
					//  = (4/3)( 1. -v( 3. -v( 3. - v) ) )
					res *= (4.0/3.0)*(1. - v*(3. - v*(3.-v)));
				}
#if 0//ndef NDEBUG
				else {
					assert ( false );
					res = 0.0;
					break;
				}
#endif
			}
		}
		grad[d] = res;
		assert ( fabs(res) > -1e-10 );
	}
}

template<int dim>
void
cubic_spline_weight<dim>::evaluate_and_grad ( const double *    co,
                                              double *          value,
                                              double *          grad ) const
{
	// Cubic spline evaluated in each dimension.
	double    co_vals[dim];
	
	// Partial derivative in each dimension.
	double    co_partial[dim];
	
	for ( int d=0; d<dim; ++d ) {
		double v = fabs( co[d] );
		if ( v < 0.5 + 1e-10 ) {
			co_vals[d] = 2.0/3.0 + 4.0*(v-1.0)*v*v;
		} else if ( v < 1.0 + 1e-10 ) {
			co_vals[d] = (4.0/3.0)*(1. - v*(3. - v*(3.-v)));
		} else {
			co_vals[d] = 0.0;
		}
		
		// Work out partial derivatives.
		if ( v < 0.5 + 1e-10 ) {
			co_partial[d] = 4.0*(3.0*v-2.0)*v;
		} else if ( v < 1.0 +1e-10 ) {
			co_partial[d] = -4.0*(1.0+(-2.0+v)*v);
		} else {
			co_partial[d] = 0.0;
		}
		
		// Adjust for sign of coord.
		if ( co[d] < 0.0 ) {
			co_partial[d] *= -1.0;
		}
	}
	
	*value = 1.0;
	for ( int d=0; d<dim; ++d ) {
		*value *= co_vals[d];
	}
	
	for ( int d=0; d<dim; ++d ) {
		grad[d] = co_partial[d];
		for ( int i=0; i<dim; ++i ) {
			if ( i != d ) {
				grad[d] *= co_vals[i];
			}
		}
	}
}

template<int dim>
void
cubic_spline_weight<dim>::evaluate_grad ( valarray<double>& in, valarray<double>& out ) const
{
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	double res = 0.0;
	
	for ( int pt=0; pt<num_pts; ++pt ) {
		double* co   = &in[dim*pt];
		double* grad = &out[dim*pt];
		
		for ( int d=0; d<dim; ++d ) {
			res = 1.0;
			for ( int i=0; i<dim; ++ i ) {
				double v = co[i]<0.0 ? -co[i] : co[i];
				if ( i == d ) {
					assert( v <= 1.0 );
					
				// Cubic spline derivative.
					if ( v < 0.5 + 1e-10 ) {
				//	// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
						
					// Derivative -8v + 12v^2 = 4(3v-2)v
						res *= 4.0*(3.0*v-2.0)*v;
						if ( co[i] < 0.0 ) {
							res *= -1.0;
						}
#if 0 // 2008-12-04 Michael LI.
						if ( co[i] < 0.0 ) {
							res *= -4.0*(3.0*v-2.0)*v;
						} else {
							res *=  4.0*(3.0*v-2.0)*v;
						}
#endif
					} else if ( v < 1.0 +1e-10 ) {
					// 4/3-4v+4v^2-(4/3)v^3
					//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
					//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
						
					// Derivative -4 + 8v -4v^2 = -4(1-2v+v^2)
					//                          = -4(1+(-2+v)v)
						res *= -4.0*(1.0+(-2.0+v)*v);
						if ( co[i] < 0.0 ) {
							res *= -1.0;
						}
#if 0 // 2008-12-04 Michael LI.
						if ( co[i] < 0.0 ) {
							res *=  4.0*(1.0+(-2.0+v)*v);
						} else {
							res *= -4.0*(1.0+(-2.0+v)*v);
						}
#endif
					}
#if 0//ndef NDEBUG
					else if ( v > 1.0 ){
						assert ( false );
						res = 0.0;
						break;
					}
#endif
				} else {
				// Cubic spline function.
					if ( v < 0.5 + 1e-10 ) {
					// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
						res *= 2.0/3.0 + 4.0*(v-1.0)*v*v;
					} else if ( v < 1.0 + 1e-10 ) {
					// 4/3-4v+4v^2-(4/3)v^3
					//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
					//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
					//  = (4/3)( 1. -3v( 1. -v( 1. - v/3.) ) )
					//  = (4/3)( 1. -v( 3. -v( 3. - v) ) )
						res *= (4.0/3.0)*(1. - v*(3. - v*(3.-v)));
					}
#if 0//ndef NDEBUG
					else if ( v > 1.0 ) {
						assert ( false );
						res = 0.0;
						break;
					}
#endif
				}
			}
			grad[d] = res;
			assert ( fabs(res) > 0.0 );
		}
	}
	// Internal consistency testing.
#if 0//ndef NDEBUG
	{
		double grad[dim];
		for ( int i=0; i<num_pts; ++i ) {
			evaluate_grad ( &in[dim*i], grad );
			for ( int d=0; d<dim; ++d ) {
				assert ( fabs ( out[dim*i+d] - grad[d] ) < 1e-10 );
			}
		}
	}
#endif
}

template<int dim>
void
cubic_spline_weight<dim>::evaluate_grad ( valarray<double>& in, valarray<bool>& pred, valarray<double>& out ) const
{
	/*
		This is the most expensive function in the assembly so we'll optimize aggressively.
		2008-12-04 Michael LI.
	*/
	int in_size = in.size();
	assert ( in_size % dim == 0 );
	int num_pts = in_size/dim;
	
	double res = 0.0;
	
	for ( int pt=0; pt<num_pts; ++pt ) {
		if ( !pred[pt] ) {
			for ( int i=0; i<dim; ++i ) {
				out[ dim*pt + i ] = 0.0;
			}
			continue;
		}
		double* co   = &in[dim*pt];
		double* grad = &out[dim*pt];
		
		for ( int d=0; d<dim; ++d ) {
			res = 1.0;
			for ( int i=0; i<dim; ++ i ) {
				double v = co[i]<0.0 ? -co[i] : co[i];
				if ( i == d ) {
					assert( v <= 1.0 );
					
				// Cubic spline derivative.
					if ( v < 0.5 + 1e-10 ) {
				//	// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
						
					// Derivative -8v + 12v^2 = 4(3v-2)v
					res *=  4.0*(3.0*v-2.0)*v;
					if ( co[i] < 0.0 ) {
						res *= -1.0;
					}
#if 0 // 2008-12-04 Michael LI.
						if ( co[i] < 0.0 ) {
							res *= -4.0*(3.0*v-2.0)*v;
						} else {
							res *=  4.0*(3.0*v-2.0)*v;
						}
#endif
#if 0//ndef NDEBUG
					} else if ( v < 1.0 +1e-10 ) {
#else
					} else {
#endif
					// 4/3-4v+4v^2-(4/3)v^3
					//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
					//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
						
					// Derivative -4 + 8v -4v^2 = -4(1-2v+v^2)
					//                          = -4(1+(-2+v)v)
						res *= -4.0*(1.0+(-2.0+v)*v);
						if ( co[i] < 0.0 ) {
							res *= -1.0;
						}
#if 0 // 2008-12-04 Michael LI.
						if ( co[i] < 0.0 ) {
							res *=  4.0*(1.0+(-2.0+v)*v);
						} else {
							res *= -4.0*(1.0+(-2.0+v)*v);
						}
#endif
					}
#if 0//ndef NDEBUG
					else if ( v > 1.0 ) {
						assert ( false );
						abort();
						res = 0.0;
						break;
					}
#endif
				} else {
				// Cubic spline function.
					if ( v < 0.5 + 1e-10 ) {
					// 2/3 - 4v^2 + 4v^3 = 2/3 +4(v-1)v^2
						res *= 2.0/3.0 + 4.0*(v-1.0)*v*v;
#ifndef NDEBUG
					} else if ( v < 1.0 + 1e-10 ) {
#else
					} else {
#endif
					// 4/3-4v+4v^2-(4/3)v^3
					//  = 4( 1./3. -v( 1. -v+ (v^2)/3.) )
					//  = 4( 1./3. -v( 1. -v( 1. - v/3.) ) )
					//  = (4/3)( 1. -3v( 1. -v( 1. - v/3.) ) )
					//  = (4/3)( 1. -v( 3. -v( 3. - v) ) )
						res *= (4.0/3.0)*(1. - v*(3. - v*(3.-v)));
					}
#if 0//ndef NDEBUG
					else if ( v > 1.0 ) {
						assert ( false );
						abort();
						res = 0.0;
						break;
					}
#endif
				}
			}
			grad[d] = res;
			assert ( fabs(res) > 0.0 );
		}
	}
	// Internal consistency testing.
#if 0//ndef NDEBUG
	{
		double grad[dim];
		for ( int i=0; i<num_pts; ++i ) {
			if ( pred[i] ) {
				evaluate_grad ( &in[dim*i], grad );
				for ( int d=0; d<dim; ++d ) {
					assert ( fabs ( out[dim*i+d] - grad[d] ) < 1e-10 );
				}
			} else {
				for ( int d=0; d<dim; ++d ) {
					assert ( fabs ( out[dim*i+d] ) < 1e-10 );
				}
			}
		}
	}
#endif
}


//

template class cubic_spline_weight<1>;
template class cubic_spline_weight<2>;
template class cubic_spline_weight<3>;
