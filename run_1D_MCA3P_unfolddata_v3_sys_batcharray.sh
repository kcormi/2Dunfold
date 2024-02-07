#!/bin/bash
#SBATCH -p long
#SBATCH --account=t3
#SBATCH --job-name=job-name
#SBATCH --mem=60G                     # memory 3GB (per job)
#SBATCH --time 2-12:00:00
#SBATCH --array=0-7
#SBATCH -n 20
obs=("spherocity" "thrust" "transverse_spherocity" "transverse_thrust" "broadening" "isotropy" "nparticle:nparticle_fine" "mass")
#obs=("mass")

python evaluate_iterations.py --config config/dataset_v7_MCA3P_unfolddata_v3.json   --method multifold --migiter 0 1 --obs nparticle,${obs[${SLURM_ARRAY_TASK_ID}]} --eff-acc
