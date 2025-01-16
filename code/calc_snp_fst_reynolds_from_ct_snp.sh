#! /bin/bash
# Estimate Fst and diversity

# fixed parameters
pop1="ZI"
pop2="FR"

# demographic parameters that vary across chromosomal arms or combinations of populations
demo_models=("chrX_6_pop" "chr2R_5_pop" "chr3L_6_pop")


summary="../data/estimate_fst_snp.report"
for model in "${demo_models[@]}"; do
    in_ct_pop1="../data/ms_simulation_${model}_${pop1}.ct"
    in_ct_pop2="../data/ms_simulation_${model}_${pop2}.ct"
    out_fst="../data/fst_${model}_snp.txt"
    # estimate Fst
    python calc_snp_fst_reynolds_from_ct_snp.py --in_ct_pop1 $in_ct_pop1 --in_ct_pop2 $in_ct_pop2 --out_fst $out_fst

    # calculate the mean Fst
    mean_fst=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $out_fst)
    # calculate the median Fst (used '-g' to sort in numerical order while correctly handling scientific notation)
    median_fst=$(awk '{print $1}' $out_fst | sort -g | awk ' { a[i++]=$1; } END { if (i % 2 == 0) print (a[i/2-1] + a[i/2]) / 2; else print a[int(i/2)]; }')
    # print a header to the summary file
    if [ ! -f $summary ]; then
        echo -e "model\tmean_fst\tmedian_fst" > $summary
    fi
    echo -e "$model\t$mean_fst\t$median_fst" >> $summary
done