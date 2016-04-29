#! /usr/bin/python
##
##	Copyright (C) 2009 Michael Li
##	This file is part of the Computenode Library.
##
##	The Computenode Library is free software: you can redistribute it and/or modify
##	it under the terms of the GNU General Public License as published by
##	the Free Software Foundation, either version 3 of the License, or
##	(at your option) any later version.
##
##	This program is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
##	GNU General Public License for more details.
##
##	You should have received a copy of the GNU General Public License
##	along with this program. If not, see <http://www.gnu.org/licenses/>.
##

import os

Help("""
\nComputenode\n
""")

# The directories should not be PathOptions as this would disallow use
# of symbolic links.
opts = Options()
opts.AddOptions(
        ('config',        'Enable configuration step',         True),
	('local_lib',     'Directory for local libraries',     '/usr/local/lib'),
	('fortran_lib',   'Directory for fortran libraries',   '' ),
	('boost_include', 'Directory for boost header files',  ''),
	('boost_lib',     'Directory for boost library files', '/usr/lib' ),
	('X11_lib',       'Directory for X11 library files',   '/usr/X11R6/lib64' ),
	('lapack_lib',    'Directory for LAPACK library',      '/usr/lib/atlas/sse2'),
        ('mpi_dir',       'Directory for MPI implementation',  '/usr/lib/mpich'),
        ('mpi_include',   'Directory for MPI header files',    '${mpi_dir}/include'),
        ('mpi_lib',       'Directory for MPI libraries',       '${mpi_dir}/lib'),
        ('petsc_dir',     'Directory for PETSc',               os.environ['PETSC_DIR']),
        ('petsc_arch',    'PETSc architecture',                os.environ['PETSC_ARCH']),
        ('petsc_include', 'Directory for PETSc header files', '${petsc_dir}/include'),
        ('petsc_lib',     'Directory for PETSc libraries',    '${petsc_dir}/lib/${petsc_arch}'),
	('vtk_include',   'Directory for VTK header files',   '/usr/include/vtk-5.0'),
	('vtk_lib',       'Directory for VTK libraries',      '/usr/lib'),
	('tiff_include',  'Directory for libtiff headers',    '/usr/include' ),
	('tiff_lib',      'Directory for libtiff libraries',  '/usr/lib'),
	BoolOption( 'openmp',        'Compile with OpenMP support'       , False ),
	BoolOption( 'fatalerrors',   'First error terminates compilation', True ),
        BoolOption( 'debug',         'Compile with debugging',             False ),
	BoolOption( 'define_ndebug', 'Define NDEBUG',                      True ),
        BoolOption( 'warn',          'Compile with warnings',              False ),
        BoolOption( 'profile',       'Compile with profiling',             False ),
	BoolOption( 'fastmath',      'Compile with fast-math',             False ),
	BoolOption( 'icc_fast',      'Intel Compiler flag fast',           False ),
	BoolOption( 'woodcrest',      'Intel Compiler flag fast',           False ),
	BoolOption( 'pgo_gen',       '(GCC) generate info for PGO',        False ),
	BoolOption( 'pgo_use',       '(GCC) profiling guided optimization', False ),
        EnumOption( 'optimize',      'Compile with optimization', '2', allowed_values=('s','0','1','2\
','3') ),
	EnumOption( 'mtune',       'Compile for architecture', 'generic', allowed_values=('generic', 'native', 'k8') ),
	EnumOption( 'compiler',    'Select compiler collection'         , 'GCC', allowed_values=('GCC', 'intel', 'pg') )
)

env = Environment( options = opts, CC='mpicc', CXX='mpicxx')
			  
# PATH, TERM, HOME needed for colorgcc.
# PGROUPD_LICENSE_FILE needed for pgi compilers on bluebird.
# LD_LIBRARY_PATH needed for mpicxx on francesca
environment_variables_to_append = [ 'PATH', 'TERM', 'HOME', 'LD_LIBRARY_PATH', 'PGROUPD_LICENSE_FILE' ]
	
for i in environment_variables_to_append:
	if i in os.environ:
		env.Append( ENV = { i : os.environ[ i ] } )

	  
platform      = ARGUMENTS.get( 'OS', Platform() )

Help( opts.GenerateHelpText( env ) )


	 
if env['debug']==1    :
        env.Append( CCFLAGS='-g', CXXFLAGS='-g', LINKFLAGS='-g' )

if env['define_ndebug']==1:
        env.Append( CPPDEFINES=['NDEBUG'] )
	if env['compiler']=='pg' :
		env.Append( CCFLAGS='-DNDEBUG', CXXFLAGS='-DNDEBUG' )

if env['warn']==1 :
	if env['compiler']=='GCC' :
        	env.Append( CCFLAGS='-Wfloat-equal -Wshadow -Wno-inline -Wall',
                            CXXFLAGS='-Woverloaded-virtual -Wfloat-equal -Wshadow -Wno-inline -Wall' )
# -Werror to force stop of compilation on error!! (don't use -Weffc++ with -Werror).

if env['profile']==1  :
        env.Append( CCFLAGS='-pg', CXXFLAGS='-pg', LINKFLAGS='-pg' )
	

if env['compiler']=='GCC' :
	env.Append( CCFLAGS='-O${optimize} -mtune=${mtune}', CXXFLAGS='-O${optimize} -mtune=${mtune}', LINKFLAGS='-O${optimize}' )
	if env['pgo_gen'] :
		env.Append( CCFLAGS='-fprofile-generate', CXXFLAGS='-fprofile-generate', LINKFLAGS='-fprofile-generate' )
	if env['pgo_use'] :
		env.Append( CCFLAGS='-fprofile-use', CXXFLAGS='-fprofile-use', LINKFLAGS='-fprofile-use' )
	if env['openmp']==1 :
		env.Append( CCFLAGS='-fopenmp', CXXFLAGS='-fopenmp', LINKFLAGS='-fopenmp' )
	if env['fastmath']==1 :
		env.Append( CCFLAGS='-ffast-math', CXXFLAGS='-ffast-math' )
	if env['fatalerrors']==1 :
		env.Append( CCFLAGS='-Wfatal-errors', CXXFLAGS='-Wfatal-errors' )

if env['compiler']=='pg' :
	if env['openmp']==1 :
		env.Append( CCFLAGS='-DBOOST_NO_EXCEPTIONS' )
		env.Append( CPPDEFINES=['BOOST_NO_EXCEPTIONS'], CCFLAGS='-mp', CXXFLAGS='-mp', LINKFLAGS='-mp' )

if env['compiler']=='intel' :
	env.Append( CCFLAGS='-O${optimize} ', CXXFLAGS='-O${optimize}', LINKFLAGS='-O${optimize}' )
	if env['icc_fast'] :
		env.Append( CCFLAGS='-fast', CXXFLAGS='-fast', LINKFLAGS='-fast' )
	if env['openmp']==1 :
		env.Append( CCFLAGS='-openmp', CXXFLAGS='-openmp', LINKFLAGS='-openmp' )
	if env['woodcrest']==1 :
		env.Append( CCFLAGS='-xT', CXXFLAGS='-xT', LINKFLAGS='-xT' )


# -ftrapping-math -ftrapv
# env.Append(CXXFLAGS   = ' -pg -g -O3 -pipe -ftrapping-math -ftrapv')

env.Append( CPPPATH=['#src',
                     '${boost_include}',
                     '${mpi_include}',
                     '${petsc_dir}/bmake/${petsc_arch}',
                     '${petsc_include}',
		     '${vtk_include}',
		     '${tiff_include}'],
            CXXFLAGS=[''],
            LIBPATH=['#src',
                     '${local_lib}',
                     '${fortran_lib}',
                     '${boost_lib}',
                     '${lapack_lib}',
                     '${mpi_lib}',
                     '${petsc_lib}',
                     '${X11_lib}',
		     '${vtk_lib}',
		     '${tiff_lib}'],
            LINKFLAGS = [''] )
#            LINKFLAGS = ['-pipe', '--as-needed'] )
	 
env.Append( LIBS=['petscts',
                  'petscsnes',
                  'petscksp',
                  'petscdm',
                  'petscmat',
                  'petscvec',
                  'petsccontrib',
                  'petsc',
		  'parmetis'] )

env.Append( LIBS = ['vtkImaging',
	'vtkGraphics',
	'vtkHybrid',
	'vtkParallel',
	'vtkRendering',
	'vtkIO',
	'vtkFiltering',
	'vtkGenericFiltering',
	'vtkCommon',
	'vtksys',
	'vtkWidgets'] )

env.Append( LIBS='tiff' )

if env['compiler']=='GCC' :
	env.Append( LIBS=['lapack','atlas','gfortran','X11'] )

if env['compiler']=='intel' :
	# Test
	#env.Append( LIBS=['flapack','fblas','ifcore'] )
	#env.Append( CCFLAGS='-xP', CXXFLAGS='-xP' )
	env.Append( LIBPATH='/opt/intel/mkl/10.0.011/lib/em64t' )
	env.Append( LIBS=['mkl_lapack','mkl','guide','pthread'] )

if env['compiler']=='pg' :
	env.Append( CXXFLAGS=['-fast'] )
	env.Append( LIBS=['lapack','blas','pgftnrtl'] )

# X11 lib not needed on francesca
#env.Append( LIBS=['computenode','boost_filesystem','mpich'])

# if env.GetOption('config')==1 and not env.GetOption('clean')==1 and not env.GetOption('help')==1 :
# 	conf = Configure( env )
# 	if not conf.CheckCHeader( 'mpi.h' ):
#     		print '\tmpi.h must be present.'
#     		Exit(1)
# 	if not conf.CheckCHeader( 'petsc.h' ):
#     		print '\tpetsc.h must be present.'
# 		Exit(1)
# 	env = conf.Finish()

# BuildDir specifies the target build directory and the location of the source code.
# The expression duplicate=0 prevents copying of the source code into the target directory
# before compilation.  The compilation work is delegated to the SConscript files located
# under the src tree.  Note the target build directory is specified before SConscript
# eg Major/SConscript

MainAlias='major'
BuildEnv=env.Copy()
#BuildPath='build/%s'%(MainAlias)
BuildPath='.'
#BuildDir( BuildPath, '.', duplicate=0)
SConscript( '%s/SConscript'%BuildPath, exports=['BuildEnv', 'MainAlias'])
#Alias('all',MainAlias)
#Default(MainAlias)

