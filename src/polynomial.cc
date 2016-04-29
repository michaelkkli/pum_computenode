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

#include "polynomial.hh"

#include <algorithm>
#include <iostream>
#include <cassert>
#include <cmath>

using std::cout;

template<int dim>
polynomial<dim>::polynomial() :
  number_of_terms(0), polynomial_order(0), common_multiplier( 1.0 )
{
	this->set_global_function();
}

template<int dim>
polynomial<dim>::polynomial ( double factor, int xpow, int ypow, int zpow ) :
  common_multiplier( 1.0 )
{
	// Only allowing low powers to help debugging.
	assert( xpow < 100 && ypow < 100 && zpow < 100 );
	
	this->set_global_function();
  if( xpow>=0 && ypow>=0 && zpow>=0 ) {
	  
	// GDB helped find that the line below is incorrect
	// as add_term also bumps the number of terms.
	// number_of_terms = 1;
	
	number_of_terms = 0; // This will be incremented by add_term but must be initialized.
	polynomial_order = xpow + ypow + zpow;
	this->add_term( factor, xpow, ypow, zpow );
  } else {
	  number_of_terms  = 0;
	  polynomial_order = 0;
  }
}

template<int dim>
polynomial<dim>::polynomial ( const polynomial& other )
{
  this->deep_copy( other );
}

template<int dim>
polynomial<dim>&
polynomial<dim>::operator= ( const polynomial& other )
{
  this->deep_copy( other );
  return *this;
}

template<int dim>
polynomial<dim>&
polynomial<dim>::operator+= ( const polynomial& other )
{
  if ( other.number_of_terms == 0 ) {
    return *this;
  } else {
    if ( this->common_multiplier == other.common_multiplier ) {
      for ( int i=0; i<other.number_of_terms; ++i ) {
	multipliers.push_back( other.multipliers[i] );
      }
    } else {
      for ( int i=0; i<number_of_terms; ++i ) {
	multipliers[i] *= common_multiplier;
      }
      common_multiplier = 1.0;
      for ( int i=0; i<other.number_of_terms; ++i ) {
	multipliers.push_back( other.common_multiplier * other.multipliers[i] );
      }
    }
    int psize = other.powers.size();
    assert( psize % dim == 0 );
    for ( int i=0; i<psize; ++i ) {
      powers.push_back( other.powers[i] );
    }
    this->number_of_terms += other.number_of_terms;
    this->polynomial_order = std::max( this->polynomial_order, other.polynomial_order );
  }
	return *this;
}

template<int dim>
polynomial<dim>::~polynomial()
{
}

template<int dim>
int polynomial<dim>::get_num_terms() const
{
  return number_of_terms;
}

template<int dim>
int polynomial<dim>::get_order() const
{
  return polynomial_order;
}

template<int dim>
void polynomial<dim>::clear()
{
  multipliers.clear();
  powers.clear();

  number_of_terms   = 0;
  polynomial_order  = 0;
  common_multiplier = 1.0;
}

template<int dim>
void
polynomial<dim>::add_term ( double factor,
			    int xpow,
			    int ypow,
			    int zpow )
{
    multipliers.push_back ( factor );
    int temp = xpow;
    powers.push_back ( xpow );
    if ( dim > 1 ) {
      temp += ypow;
      powers.push_back ( ypow );
    }
    if ( dim > 2 ) {
      temp += zpow;
      powers.push_back ( zpow );
    }
    ++number_of_terms;
    polynomial_order = std::max ( polynomial_order, temp );
}

template<int dim>
double polynomial<dim>::evaluate ( const double* co ) const
{
  if ( number_of_terms == 0 ) {
    return 0.0;
  }
  double poly_val = 0.0;
  double term_val;
  
  assert( number_of_terms == powers.size()/dim );
  
  for ( int i=0; i<number_of_terms; ++i ) {
	term_val = multipliers[i];
	for ( int d=0; d<dim; ++d ) {
		assert( dim*i+d < powers.size() );
		assert( powers[ dim*i + d ] >= 0 );
		term_val *= std::pow( co[d], powers[ dim*i + d ] );
	}
	poly_val += term_val;
  }
  
  // TODO: remove. Debug tiny values seen. Point sent in seen
  // to be tiny - of order 1e-310
  // assert( fabs( common_multiplier * poly_val ) > 1e-20 );
	
#if 0 //ndef NDEBUG // Keep for future debug.
	// Unusually large values are coming out of monomials. Debug.
	if ( number_of_terms == 1 && this->is_local_function() ) {
		std::clog << "Polynomial assumed monomial common mult " << common_multiplier << ".\n";
		for ( int i=0; i<number_of_terms; ++i ) {
			if ( i!=0 ) {
				std::clog << "+ ";
			}
			std::clog << multipliers[i] << " * ";
			for ( int d=0; d<dim; ++d ) {
				if ( d!=0 ) {
					std::clog << " ";
				}
				assert( dim*i+d < powers.size() );
				assert( powers[ dim*i + d ] >= 0 );
				std::clog << "(" << co[d] << ")^" << powers[ dim*i + d ];
			}
		}
		std::clog << " = " << common_multiplier * poly_val <<  "\n";
		assert ( fabs( static_cast<long double>(poly_val)<=2.0 ) );
	}
#endif
  
  return common_multiplier * poly_val;
}

template<int dim>
void
polynomial<dim>::evaluate_grad ( const double* co, double* out ) const
{
	assert( co && out );
  for ( int partial=0; partial<dim; ++partial ) {
	if ( number_of_terms == 0 ) {
		out[partial] = 0.0;
		continue;
	}
  
	double poly_val = 0.0;
	double term_val;
	
	for ( int i=0; i<number_of_terms; ++i ) {
		term_val = multipliers[i];
		if ( powers[ dim*i + partial ] == 0 ) {
			// This term is zero.
			out[partial]=0.0;
			continue;
		}
		for ( int d=0; d<dim; ++d ) {
			assert( dim*i+d < powers.size() );
			assert( powers[ dim*i + d ] >= 0 );
			if ( partial != d ) {
				// Ordinary factor.
				term_val *= std::pow( co[d], powers[ dim*i + d ] );
			} else {
				// Factor requires differentiation.
				if ( powers[dim*i+d]!=1 ) {
					// Factor makes a contribution to the term.
					// If the power is equal to one, we would needlessly multiply by 1.0 twice.
					term_val *= powers[dim*i + d]*std::pow( co[d], powers[dim*i + d]-1 );
				}
			}
		}
		poly_val += term_val;
	}
	out[partial] = common_multiplier * poly_val;
  }
  
  // TODO: remove. Debug tiny values seen. Point sent in seen
  // to be tiny - of order 1e-310
  // assert( std::abs( common_multiplier * poly_val ) > 1e-20 );
}

template<int dim>
void polynomial<dim>::deep_copy ( const polynomial& other )
{
	// TODO: remove after getting disappearing bug.
	// cout << "Copy constructor unexpectedly called.\n";
	
	
  this->number_of_terms   = other.number_of_terms;
  this->polynomial_order  = other.polynomial_order;
  this->common_multiplier = other.common_multiplier;
  this->multipliers       = other.multipliers;
  this->powers            = other.powers;
}

//

template class polynomial<1>;
template class polynomial<2>;
template class polynomial<3>;
