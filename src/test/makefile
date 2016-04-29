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

# imagemagick convert -gaussian-blur radiusxsigma
# Setting radius zero give automatic setting of a suitable radius.
# for i in *.tif ; do convert "$i" -gaussian-blur 3x3 "$i" ; done


# ffmpeg -f image2 -i filename%04d.jpeg video.mpeg
# Adjust codec:
# ffmpeg -f image2 -i lev7_F2Rep1_.%04d.jpg  -sameq -vcodec mpeg2video -r 20 lev7_F2Rep1.mpeg2.mpeg 
# Limit file size to 6 MB :
# ffmpeg -f image2 -i lev7_F2Rep1_.%04d.jpg  -sameq -fs 6000000 lev7_F2Rep1.mov
# Slow down video:
# ffmpeg -i *.mov -f yuv4mpegpipe - | yuvfps -s 5:1 -r 25:1 | ffmpeg -f yuv4mpegpipe -i - -b 28800k -y output.avi
# Slow video:
# ffmpeg -y -f image2 -r 5 -i lev7_F2Rep1_.%04d.jpg  -sameq -r 25 test.mov

#in_folder=input/F2 Rep 1 blur
#domain=F2Rep 1 FRAP Series11_t00_ch01.tif.off

#in_folder=input/F2 Null 1
#domain=F2Null 1 FRAP Series31_t00_ch01.tif.off

in_folder=input/H Null 2
domain=HNull 2 FRAP Series25_t00_ch01.tif.off

all:

bdry:
	./petsc_solver-3  -dtree_initial_refinement 1 -boundary_detail 3

oldproj:
	mpirun -np 2 ./petsc_solver-3 -dom "F2Rep 1 FRAP Series11_t00_ch01.tif.off" \
                                      -dtree_initial_refinement 7 \
                                      -simulation_mode 2 \
                                      -manifest1_dir "input/F2 Rep 1 blur" \
                                      -manifest1 "input/F2 Rep 1 blur/manifest" \
                                      -min_iter_info

proj:
	mpirun -np 2 ./petsc_solver-3 -dom "${domain}" \
                                      -dtree_initial_refinement 7 \
                                      -simulation_mode 2 \
                                      -manifest1_dir "${in_folder}" \
                                      -manifest1 "${in_folder}/manifest" \
                                      -min_iter_info

# local_basis_order 1, assembly_quadrature_rule 2, lev 6, viz_refinement 8  shows stepped edges in bleach region
# local_basis_order 1, assembly_quadrature_rule 3, lev 6, viz_refinement 8  loses some detail on running
# local_basis_order 1, assembly_quadrature_rule 5, lev 6, viz_refinement 8  loses some detail on running
helsinkiclflip:
# global_bleach1 0.005 ok
# y is 511-
# diff 300 good for mobile
# 518 is 0.3 um/s
# 173 is 0.1
# 17.3 is 0.01
# diff param 3 good
	/usr/bin/time -v mpirun -np 2  \
	./petsc_solver-3 -dom extract_geometry_2d-1_contour.off \
	                              -dtree_initial_refinement 9 \
	                              -image1 initial_condition.tif \
	                              -parent_box_expansion_factor 1.000001 \
	                              -diffusion1 173 \
	                              -global_bleach1 0.0 \
	                              -bleach_x_centre1 246.0 \
	                              -bleach_y_centre1 211.0 \
	                              -bleach_param1 9.0 \
	                              -bleach_radius1 8.31 \
                                      -begin_bleach1 1 \
                                      -end_bleach1 51 \
                                      -local_basis_order 1 \
                                      -assembly_quadrature_rule 4 -vtk_out 3 -simulation_mode 1 \
	                              -min_iter_info false  \
	                              > petsc_solver-3_log.txt

nightclflip:
	/usr/bin/time -v mpiexec -np 1 ./petsc_solver-3 -dom extract_geometry_2d-1_contour_full.off \
                                      -dt 0.754      \
                                      -image1 initial_condition.tif \
                                      -dtree_initial_refinement 6 \
                                      -iterations 71 \
                                      -simulation_mode 1 \
                                      -diffusion1 120.0 \
	                              -global_bleach1 0.003 \
	                              -bleach_x_centre1 246.0 \
	                              -bleach_y_centre1 211.0 \
                                      -bleach_radius1 8.31 \
	                              -bleach_param1 1.0 \
                                      -begin_bleach1 1 \
                                      -end_bleach1 51 \
                                      -local_basis_order 1 \
				      -assembly_quadrature_rule 6 \
				      -min_iter_info false > petsc_solver-3_log.txt 2>&1

helsinkiconvergencebasic:
	/usr/bin/time -v mpirun -np 1 ./pum_convergence-1 -dom unit_box.off -parent_box_expansion_factor 1.00001 -initial_refinement 8 -max_refinement 8 \
		-ksp_view  2>&1 > pum_convergence-1_log.txt &

## Some flags that are probably only useful for very large problem sizes :
## -openmp_spawn_size 2 -pc_type openmp -openmp_pc_type asm -ksp_type preonly -ksp_gmres_restart 200 -ksp_view
## Aim for greater than 10 000 or 20 000 unknowns per process to justify communication cost.

helsinkiconvergenceone:
	/usr/bin/time -v mpiexec -np 1 ./pum_convergence-1 -dom unit_box.off -parent_box_expansion_factor 1.00000001 -initial_refinement 2 -max_refinement 9 \
		-pc_type asm -ksp_gmres_restart 200 -ksp_view -log_summary 2>&1 > pum_convergence-1_log.txt &

helsinkiconvergencetwo:
	/usr/bin/time -v mpiexec -np 2 ./pum_convergence-1 -dom unit_box.off -parent_box_expansion_factor 1.00001 -initial_refinement 2 -max_refinement 9 \
		-openmp_merge_size 2 -pc_type asm -ksp_gmres_restart 200 -ksp_view -log_summary 2>&1 > pum_convergence-1_log.txt &

vilniusconvergence:
	/usr/bin/time -v mpiexec -np 2 ./pum_convergence-1 \
		-initial_refinement 2 -max_refinement 9 \
		-dom unit_box.off -parent_box_expansion_factor 1.00001 \
		-local_basis_order 1 \
		-openmp_merge_size 2 -pc_type asm -ksp_gmres_restart 200 \
		-ksp_view -log_summary 2>&1 > pum_convergence-1_log.txt &

helsinkiconvergencesymm:
	/usr/bin/time -v mpiexec -np 1 ./pum_convergence-1 -dom unit_box.off -parent_box_expansion_factor 1.00001 -info -initial_refinement 9 -max_refinement 9 \
		-mat_type seqsbaij -mat_ignore_lower_triangular -pc_type icc

## Level 9 is first time that matrix memory usage is a problem.
nightconvergence:
	/usr/bin/time -v mpiexec -np 1 ./pum_convergence-1 \
		-dom unit_box.off \
		-parent_box_expansion_factor 1.01 \
		-cover_factor 1.3 \
		-initial_refinement 7 \
		-max_refinement 8 \
		-local_basis_order 1 \
		-box_interior_quad_rule 9 \
		-ksp_view \
		-ksp_rtol 1e-6 \
		-ksp_gmres_restart 500  > pum_convergence-1_log.txt 2>&1 &


		#-image1 initial_condition.tif             \
		#
		# diffusion1 173
newclflip:
	/usr/bin/time -v mpiexec -np 1 ./pum_convergence-2 \
		-dom extract_geometry_2d-1_contour.off     \
		-image1 ""           \
		-image1_siggia 0 \
		-anticlockwise_dom 0                       \
		-parent_box_expansion_factor 1.01          \
		-cover_factor 1.1                          \
		-box_interior_quad_rule 9                  \
		-initial_refinement 4                      \
		-diffusion1 17.3273                          \
		-time_step_subdivision 1                   \
	        -global_bleach1 0.006662587                      \
	        -bleach_x_centre1 246.0                    \
	        -bleach_y_centre1 418.0                    \
                -bleach_radius1 8.176                       \
	        -bleach_param1 0.6662587

species:
	/usr/bin/time -v mpiexec -np 1 ./pum_convergence-3 \
		-dom extract_geometry_2d-1_contour.off     \
		-image1 ""           \
		-image1_siggia 0 \
		-anticlockwise_dom 0                       \
		-parent_box_expansion_factor 1.01          \
		-cover_factor 1.1                          \
		-box_interior_quad_rule 9                  \
		-initial_refinement 5                      \
		-initial_proportion0 0.0                   \
		-initial_proportion1 1.0                   \
		-initial_proportion2 0.0                   \
		-initial_proportion3 0.0                   \
		-diffusion0 1732.73                        \
		-diffusion1 1.73273                        \
		-diffusion2 1732.73                        \
		-diffusion3 1.73273                        \
		-a0t1 0.0                                  \
		-a0t2 0.0 \
		-a0t3 0.0 \
		-a1t0 0.0 \
		-a1t2 0.0 \
		-a1t3 0.0 \
		-a2t0 0.0 \
		-a2t1 0.0 \
		-a2t3 0.0 \
		-a3t0 0.0 \
		-a3t1 0.0 \
		-a3t2 0.0 \
		-time_step_subdivision 1                   \
	        -global_bleach1 0.006662587                      \
	        -bleach_x_centre1 246.0                    \
	        -bleach_y_centre1 418.0                    \
                -bleach_radius1 8.176                       \
	        -bleach_param1 0.6662587                   \
		-num_species 4

petsc_solver-3_tarbz2:
	tar -jvcf out.tar.bz2 petsc_solver-3_*

line_segments-1_tarbz2:
	/usr/bin/time -v tar -jcf out.tar.bz2 line_segments-1_*


ccg:
	rm -f callgrind.out.*

cvtp:
	rm -f *.vtp

cmpi:
	rm -f PI[0-9]*

clean: ccg cvtp cmpi
