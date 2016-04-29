#!/bin/bash

## sp is species
## cl is class (e.g. FLIP or obs)

## Less than 0.01
coupling_param=0.01

focal_percentage=1.0

for sp in tp hc f1 f2; do
    for cl in white black; do
        for rep in `seq 0 4`; do

	case ${sp} in
	tp)
	bleach_x=254.0
	bleach_y=298.0
	conc_mobile=0.5
	conc_immobile=0.0
	multistep=1
	bind_unbind_param=0.0
	;;
	hc)
	bleach_x=246.0
	bleach_y=418.0
	conc_mobile=0.0
	conc_immobile=0.5
	multistep=1
	bind_unbind_param=0.0
	;;
	f1)
	bleach_x=280.0
	bleach_y=388.0
	conc_mobile=0.2
	conc_immobile=0.3
	multistep=1
	bind_unbind_param=0.05
	;;
	f2)
	bleach_x=246.0
	bleach_y=211.0
	conc_mobile=0.3
	conc_immobile=0.2
	multistep=1
	bind_unbind_param=0.05
	;;
	esac

	case ${cl} in
	'white')
	spot_bleach_param=6.662587
	;;
	'black')
	spot_bleach_param=0.0
	;;
	esac

        fill_other_reps="true"

	case ${rep} in
	0)
	;;
	1)
	continue
	;;
	2)
	continue
	;;
	3)
	continue
	;;
	4)
	continue
	;;
	esac

	/usr/bin/time -v mpiexec -np 1 ./pum_convergence-3 \
	        -time_step_subdivision ${multistep} \
		-dom ${sp}_sim_dom.off     \
		-image1 ""           \
		-image1_siggia 0 \
		-anticlockwise_dom 0                       \
		-parent_box_expansion_factor 1.01          \
		-cover_factor 1.1                          \
		-box_interior_quad_rule 9                  \
		-initial_refinement 4                      \
		-initial_proportion0 ${conc_mobile}                   \
		-initial_proportion1 ${conc_immobile}         \
		-initial_proportion2 ${conc_mobile}                   \
		-initial_proportion3 ${conc_immobile}                   \
		-diffusion0 1732.73                        \
		-diffusion1 1.73273                        \
		-diffusion2 1732.73                        \
		-diffusion3 1.73273                        \
		-a0t1 ${bind_unbind_param}                                  \
		-a0t2 0.0 \
		-a0t3 0.0 \
		-a1t0 ${bind_unbind_param} \
		-a1t2 0.0 \
		-a1t3 0.0 \
		-a2t0 0.0 \
		-a2t1 0.0 \
		-a2t3 ${bind_unbind_param} \
		-a3t0 0.0 \
		-a3t1 0.0 \
		-a3t2 ${bind_unbind_param} \
		-ac   ${coupling_param} \
		-focal_percentage ${focal_percentage}      \
	        -global_bleach0 0.006662587                      \
	        -global_bleach1 0.006662587                      \
	        -global_bleach2 0.0                      \
	        -global_bleach3 0.0                      \
	        -bleach_x_centre1 ${bleach_x}                    \
	        -bleach_y_centre1 ${bleach_y}                    \
                -bleach_radius1 8.176                       \
		-bleach_param0 ${spot_bleach_param} \
	        -bleach_param1 ${spot_bleach_param}                   \
		-bleach_param2 0.0 \
		-bleach_param3 0.0 \
		-num_species 4 && \
		mv pum_convergence-3_species_data.gnuplot pum_convergence-3_species_compare_${sp}_${cl}_${rep}.gnuplot

            if [ ${fill_other_reps} = "true" ]; then
		for i in `seq 1 4`; do
		    cp pum_convergence-3_species_compare_${sp}_${cl}_0.gnuplot \
		       pum_convergence-3_species_compare_${sp}_${cl}_${i}.gnuplot
	        done
	    fi
	done
    done
done

./species_compare.sh

R --no-save < clflip3_species_compare.R

