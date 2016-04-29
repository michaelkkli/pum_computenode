#!/bin/bash

if [ -f pum_convergence-3_species_compare_all.txt ]; then
    mv pum_convergence-3_species_compare_all.txt pum_convergence-3_species_compare_all.txt.bak
fi

awk '{print $1}' pum_convergence-3_species_compare_hc_0.gnuplot > pum_convergence-3_species_compare_all.txt

for i in tp hc f1 f2 ; do
    for class in white black; do
        for rep in `seq 0 4`; do
            awk '{print $2+$3}' pum_convergence-3_species_compare_${i}_${class}_${rep}.gnuplot > tmp_col.txt
	    mv pum_convergence-3_species_compare_all.txt tmp_old_all.txt
            paste tmp_old_all.txt tmp_col.txt > pum_convergence-3_species_compare_all.txt
	done
    done
done
