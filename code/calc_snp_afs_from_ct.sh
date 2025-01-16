#! /bin/bash
# Estimate minor allele frequency spectrum (AFS) for each population

# fixed parameters
pop1="ZI"
pop2="FR"

# demographic parameters that vary across chromosomal arms or combinations of populations
demo_models=("chrX_6_pop" "chr2R_5_pop" "chr3L_6_pop")

for model in "${demo_models[@]}"; do
    in_ct_pop1="../data/ms_simulation_${model}_${pop1}.ct"
    in_ct_pop2="../data/ms_simulation_${model}_${pop2}.ct"
    out_afs_pop1="../data/afs_${pop1}_${model}.txt"
    out_afs_pop2="../data/afs_${pop2}_${model}.txt"
    
    # estimate minor allele frequency spectrum (AFS)
    python calc_snp_afs_from_ct.py --in_ct_pop1 $in_ct_pop1 --in_ct_pop2 $in_ct_pop2 --out_afs_pop1 $out_afs_pop1 --out_afs_pop2 $out_afs_pop2
done