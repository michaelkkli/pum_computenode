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

#ifndef _PUM_DISCRETIZATION_HH_
#define _PUM_DISCRETIZATION_HH_

#include "integration_scheme.hh"
#include "dbinary_tree.hh"
#include "box.hh"
#include "function.hh"
#include "geometry.hh"

#include <boost/shared_ptr.hpp>

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <valarray>

#include <petscts.h>

// Definition suggested by MPICH2 FAQ to get around the fact
// "both stdio.h and the MPI C++ interface use SEEK_SET, SEEK_CUR, and SEEK_END"
#define MPICH_IGNORE_CXX_SEEK
#include <mpi.h>

using boost::shared_ptr;

using std::map;
using std::ofstream;
using std::set;
using std::string;
using std::valarray;

template<int dim>
class pum_discretization {
public:
	typedef typename function<dim>::const_pfunction const_pfunction;
	typedef typename function<dim>::pfunction       pfunction;
	typedef typename integration_scheme<dim>::pintegration_scheme pintegration_scheme;
public:
	pum_discretization();
	~pum_discretization();
public:
	void set_dbinary_tree( shared_ptr<dbinary_tree<dim> >& );
	// Must have geometry before creating a cover.
	void create_species ( const string& species_name, const string& region_or_boundary );
	
	// TODO: remove.
	// void species_set_function( const string& species, const pfunction& ic );
	
	const const_pfunction& get_species_function( const string& species ) const;
	void set_integration_scheme( const pintegration_scheme int_sch );
	// int get_num_patches( const string& species ) const;
	bool is_species ( const string& species_name ) const;
	bool is_boundary_species( const string& species ) const;
	bool is_interior_species( const string& species ) const;
	int get_num_dof ( const string& species ) const;
	int get_num_patches ( const string& species ) const;
	
	// TODO: remove.
	// void set_coefficients( const string& species, const valarray<double>& );
	// const valarray<double>& access_species_coefficients ( const string& species ) const;
	
	void set_default_weight_and_local_basis(); // Keep for reset purposes
	void update_global_index_start ();
	void update_parallel_partitioning();
	void update_local_basis_sizes ();
	int get_local_basis_size( const string& species, const string& key ) const;
	void update();
	bool is_updated() const;
	int get_global_index_start( const string& species, const string& key ) const;
public:
	const const_pfunction& access_patch_weight( const string& species ) const;
	const map<string, vector<const_pfunction> >&
		access_local_basis( const string& species ) const;
public: // For the benefit of integration_scheme.
	void get_species_patch_neighbours( const string& species, const string& key, set<string>& out ) const;
	void get_patch( const string& species, const string& key, box<dim>& ) const;
	void get_patches_and_weights ( const string& species,
				       const vector<string>&,
				       vector<box<dim> >&,
				       vector<const_pfunction>& ) const;
	void get_weight( const string& species, const_pfunction& weight ) const;
	void get_local_basis( const string& species,
			      const string& key,
			      vector<const_pfunction>& res ) const;
	const vector<const_pfunction>& access_local_basis( const string& species,
								const string& key ) const;
	bool species_open_intersect_point( const string& species, const double* co ) const;
public: // For the benefit of petsc_solver
	int get_num_dof() const;
	const valarray<int>& access_dof_allocation( int comm_size ) const;
	void assemble_mass_matrix( Mat mass_matrix ) const;
	void assemble_rhs_vector( Vec rhs_vec) const;
private:
	string find_patch_containing_point( const string& species, const double* co ) const;
	string get_species_support( const string& species ) const;
	void get_patch_keys( const string& species, set<string>& ) const;
	const set<string>& access_patch_keys( const string& species ) const;
	bool is_species_patch_key( const string& species, const string& key ) const;
	void get_patches( const string& species, vector<box<dim> >& cover ) const;
private: // For evaluation.
	void get_local_coefficients( const string& species,
					const valarray<double>& coeff,
				     const string& key,
				     valarray<double>& out ) const;
	void evaluate_partition_of_unity( const string& species,
					  const double* co,
       					  valarray<double>& ) const;
	double evaluate_local_approximation( const string& species,
					const valarray<double>& coeff,
					     const string& key,
					     const double* co ) const;
	void evaluate_local_approximation( const string& species,
					const valarray<double>& coeff,
					   const double* co,
					   valarray<double>& ) const;
	double evaluate_pum_discretization_approximation( const string& species,
					const valarray<double>& coeff,
					 const double* co ) const;
	void evaluate_pum_discretization_approximation( const string& species,
					const valarray<double>& coeff,
					 const valarray<double>& in,
					 valarray<double>& res ) const;
public:
	void gp_draw_approximation( const string&           species,
	                            const valarray<double>& coeff,
	                            ofstream&,
	                            int                     points_per_axis=6/dim,
	                            double                  val_offset=0.0 );
private:
	bool have_coefficients;
	bool have_weights;
	bool have_local_basis;
	bool have_integration_scheme;
	shared_ptr<dbinary_tree<dim> >       dtree;
	shared_ptr<integration_scheme<dim> > integrator;
	// Map species names to region names and boundary names.
	map<string, string>            species_support;
private:
	// TODO: REMOVE.
	//map<string, valarray<double> >     present_coefficients;
	
	
	const_pfunction                    default_weight;
	vector<const_pfunction>            default_local_basis;
	map<string, const_pfunction>       present_patch_weight; // Species name as key.
	map<string, map<string, vector<const_pfunction> > > present_local_basis; // Species then patch key.
private:
	int                                  comm_size, comm_rank;
	bool                                 updated;
	map<string, map<string, int> >       present_global_index_start; // Species, patch key.
	map<string, map<string, int> >       present_local_basis_size;
	valarray<int>                        present_patch_allocation;
	map<int, map<string, set<string> > > present_parallel_partition;
	valarray<int>                        present_allocation_of_dof;
	#if 0 // TODO: remove.
private:
	map<string, const_pfunction>               species_initial_conditions;
	#endif
private:
	pum_discretization( const pum_discretization<dim>& );                  // Not implemented.
	pum_discretization<dim>& operator= ( const pum_discretization<dim>& ); // Not implemented.
};

#endif // _PUM_DISCRETIZATION_HH_
