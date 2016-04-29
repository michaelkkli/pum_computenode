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

#ifndef _FUNCTION_HH_
#define _FUNCTION_HH_

#include "box.hh"
#include <boost/shared_ptr.hpp>
#include <valarray>
using boost::shared_ptr;
using std::valarray;

template<int dim> class differentiable_function;

/**
	The user needs to define a single member function to be able
	to make use of this abstract class. The user should also
	call one of set_global_function and set_local_function.
*/
template<int dim>
class function {
public:
  friend class differentiable_function<dim>;
  typedef shared_ptr<const function<dim> > const_ptr;
  typedef shared_ptr<function<dim> >       ptr;

  typedef shared_ptr<const function<dim> > const_pfunction;
  typedef shared_ptr<function<dim> >       pfunction;
private:
  enum function_type { local_function, global_function };
public:
  function();
  virtual ~function();
private:
  function( const function<dim>& );             // Not implemented.
  function<dim>& operator= ( const function<dim>& ); // Not implemented.
public:
  void set_global_function();
  void set_local_function();

  bool is_global_function() const;
  bool is_local_function() const;

private:
  /**
     The single member function that must be provided in a subclass
     inheriting from the abstract class function.
     N.B. It is private and const.
  */
  // Evaluation at single points.
  virtual
  double evaluate ( const double* ) const=0;

  // The user may optionally override this one.
  /**
     Evaluate a number of different points with a single member function
     call.
     @param valarray<double>& of size npts.
     @param valarray<double>& of size npts/dim.
  */
  virtual
  void evaluate ( valarray<double>& points,
                  valarray<double>& result ) const;
	
	virtual
		void evaluate ( valarray<double>& in, valarray<bool>& pred, valarray<double>& out ) const;
public:
  double global_evaluate ( const box<dim>&,
			   const double*  ) const;
  double local_evaluate ( const box<dim>&,
			  const double* ) const;
  double restricted_local_evaluate ( const box<dim>& restricted_local,
				     const box<dim>& local,
				     const double* ) const;

  // Evaluation at a number of point with a single call.
public:
  void global_evaluate ( const box<dim>&,
			 valarray<double>& points,
			 valarray<double>&       results ) const;
	
  void local_evaluate ( const box<dim>&,
			valarray<double>& points,
			valarray<double>&       results ) const;
	
	void local_evaluate ( const box<dim>&,
	                      valarray<double>& in,
	                      valarray<bool>&   pred,
	                      valarray<double>& out ) const;
	
  void restricted_local_evaluate ( const box<dim>&         restricted_local,
				   const box<dim>&         local,
				   valarray<double>& points,
				   valarray<double>&       results ) const;

private:
  function_type ftype;
};

// Abstract class so need the definitions.
// Use #include "function.cc" when inheriting from function.
// Reconsidering above comment as polynomial seems to do ok with just the header.

#endif // _FUNCTION_HH_
