#! /bin/bash
# Estimate Fst and diversity

# fixed parameters
pop1="ZI"
pop2="FR"
win_size=5000

# demographic parameters that vary across chromosomal arms or combinations of populations
demo_models=("chrX_6_pop" "chr2R_5_pop" "chr3L_6_pop")


summary="../data/estimate_fst_diversity.report"
for model in "${demo_models[@]}"; do
    in_ct_pop1="../data/ms_simulation_${model}_${pop1}.ct"
    in_ct_pop2="../data/ms_simulation_${model}_${pop2}.ct"
    in_win="../data/ms_simulation_${model}.win"
    out_fst="../data/fst_${model}.txt"
    out_diversity_pop1="../data/diversity_${model}_${pop1}.txt"
    out_diversity_pop2="../data/diversity_${model}_${pop2}.txt"
    # estimate Fst and diversity
    python calc_snp_fst_reynolds_from_ct_window.py --win_size $win_size --in_win $in_win --in_ct_pop1 $in_ct_pop1 --in_ct_pop2 $in_ct_pop2 --out_fst $out_fst --out_diversity_pop1 $out_diversity_pop1 --out_diversity_pop2 $out_diversity_pop2

    # calculate the mean Fst and diversity across windows
    mean_fst=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $out_fst)
    mean_diversity_pop1=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $out_diversity_pop1)
    mean_diversity_pop2=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' $out_diversity_pop2)

    # calculate the median Fst and diversity across windows
    median_fst=$(awk '{print $1}' $out_fst | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }')
    median_diversity_pop1=$(awk '{print $1}' $out_diversity_pop1 | sort -g | awk ' { a[i++]=$1; } END { if (i % 2 == 0) print (a[i/2-1] + a[i/2]) / 2; else print a[int(i/2)]; }')
    median_diversity_pop2=$(awk '{print $1}' $out_diversity_pop2 | sort -g | awk ' { a[i++]=$1; } END { if (i % 2 == 0) print (a[i/2-1] + a[i/2]) / 2; else print a[int(i/2)]; }')

    # print a header to the summary file
    if [ ! -f $summary ]; then
        echo -e "model\tmean_fst\tmedian_fst\tmean_diversity_${pop1}\tmedian_diversity_${pop1}\tmean_diversity_${pop2}\tmedian_diversity_${pop2}" > $summary
    fi
    echo -e "$model\t$mean_fst\t$median_fst\t$mean_diversity_pop1\t$median_diversity_pop1\t$mean_diversity_pop2\t$median_diversity_pop2" >> $summary
done