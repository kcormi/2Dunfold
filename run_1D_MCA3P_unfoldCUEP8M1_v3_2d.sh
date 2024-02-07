#!/bin/bash
#SBATCH -p standard
#SBATCH --account=t3
#SBATCH --job-name=job-name
#SBATCH --mem=12G                     # memory 3GB (per job)
#SBATCH --time 10:00:00
#SBATCH --array=0-6
#SBATCH -n 10
obs=("spherocity" "thrust" "transverse_spherocity" "transverse_thrust" "broadening" "isotropy" "mass")

    python evaluate_iterations.py --config config/dataset_v7_MCA3P_unfoldCUEP8M1_v3_2d.json   --method multifold --migiter 0 1 --obs nparticle:nparticle_dim1_2D,${obs[${SLURM_ARRAY_TASK_ID}]} --eff-acc
#python merge_bs.py --input config/merge_MCA3P_unfoldCUEP8M1_v3_2d.txt --output results_finebin_v7_MCA3P_unfoldCUETP8M1_v3_2d_1/merge.csv
#python plot_iter.py --obs nparticle:nparticle_dim1_2D,spherocity nparticle:nparticle_dim1_2D,thrust nparticle:nparticle_dim1_2D,transverse_spherocity nparticle:nparticle_dim1_2D,transverse_thrust nparticle:nparticle_dim1_2D,broadening nparticle:nparticle_dim1_2D,isotropy nparticle:nparticle_dim1_2D,mass --config config/plot_1d_v7_MCA3P_unfoldCUEP8M1_v3_2d.json
