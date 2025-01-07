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
    # estimate Fst and diversity
    python calc_snp_fst_reynolds_from_ct_snp.py --in_ct_pop1 $in_ct_pop1 --in_ct_pop2 $in_ct_pop2 --out_fst $out_fst

    # calculate the mean Fst and diversity across windows
    mean_fst=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $out_fst)
    # calculate the median Fst and diversity across windows
    median_fst=$(awk '{print $1}' $out_fst | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
    # print a header to the summary file
    # print a header to the summary file
    if [ ! -f $summary ]; then
        echo -e "model\tmean_fst\tmedian_fst" > $summary
    fi
    echo -e "$model\t$mean_fst\t$median_fst" >> $summary
done