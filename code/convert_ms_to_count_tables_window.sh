#! /bin/bash
# Convert the ms output to a genotype count table

# demographic parameters that vary across chromosomal arms or combinations of populations
demo_models=("chrX_6_pop" "chr2R_5_pop" "chr3L_6_pop")
pop1="ZI"
pop2="FR"
sample_size=18

for model in "${demo_models[@]}"; do
    in_sim="../data/ms_simulation_${model}.txt"
    out_ct_pop1="../data/ms_simulation_${model}_${pop1}.ct"
    out_ct_pop2="../data/ms_simulation_${model}_${pop2}.ct"
    out_win="../data/ms_simulation_${model}.win"
    python convert_ms_to_count_tables_window.py --in_ms_result $in_sim --out_ct_pop1 $out_ct_pop1 --out_ct_pop2 $out_ct_pop2 --out_win $out_win --sample_size $sample_size
done