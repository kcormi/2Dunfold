#!/bin/env python
from optparse import OptionParser
import os, sys, re
import time
import glob
parser = OptionParser()

run_template = """#!/bin/bash
#SBATCH -p standard
#SBATCH --account=t3
#SBATCH --job-name=job-name
#SBATCH --mem=12G                     # memory 3GB (per job)
#SBATCH --time 10:00:00

for obs in nparticle:nparticle_fine mass spherocity thrust broadening transverse_spherocity transverse_thrust isotropy; do
#    python evaluate_iterations.py --config {DATASET}   --method multifold --migiter 0 1 --obs nparticle:nparticle_dim1_2D,${{obs}} --eff-acc --df-tag bs_{BSINDEX}
    python evaluate_iterations.py --config {DATASET}   --method multifold --migiter 0 1 --obs nparticle,${{obs}} --eff-acc --df-tag bs_{BSINDEX}
done
python merge_bs.py --input {MERGE} --output results_finebin_v7_MCA3P_unfoldCUEP8M1_maxweight10_minweightp1_coverage_v3_1/merge_bs_{BSINDEX}.csv
#python coverage.py --input config/merge_MCA3P_unfoldCUEP8M1_maxweight10_minweightp1_v3_stat_sys_bs.txt --output results_finebin_v7_MCA3P_unfoldCUEP8M1_maxweight10_minweightp1_coverage_v3_1/merge_bs.csv
python plot_bias.py
"""

bsindex_list = range(21,50)

for bsindex in bsindex_list:
  with open("config/dataset_v7_MCA3P_unfoldCUEP8M1_v3_bs_TEMPLATE.json",'r') as f_dataset_template:
    dataset = f_dataset_template.read()
    print(dataset)
    dataset = dataset.format(BSINDEX = bsindex)
  with open("config/merge_MCA3P_unfoldCUEP8M1_v3_bs_TEMPLATE.txt",'r') as f_merge_template:
    merge = f_merge_template.read().format(BSINDEX = bsindex)
  with open("config/plot_1d_v7_MCA3P_unfoldCUEP8M1_bs_TEMPLATE.json",'r') as f_plot_template:
    plot = f_plot_template.read().format(BSINDEX = bsindex)

  f_dataset_name = "config/dataset_v7_MCA3P_unfoldCUEP8M1_v3_sys_stat_bs_{}.json".format(bsindex)
  f_merge_name = "config/merge_MCA3P_unfoldCUEP8M1_v3_stat_sys_bs_{}.txt".format(bsindex)
  f_plot_name = "config/plot_1d_v7_MCA3P_unfoldCUEP8M1_v3_stat_sys_bs_{}.json".format(bsindex)
  f_run_name = "run_1D_MCA3P_unfoldCUEP8M1_stat_sys_bs_{}.sh".format(bsindex)
  script_dataset = open(f_dataset_name,'w')
  script_dataset.write(dataset)
  script_dataset.close()

  script_merge = open(f_merge_name,'w')
  script_merge.write(merge)
  script_merge.close()

  script_plot = open(f_plot_name,'w')
  script_plot.write(plot)
  script_plot.close()

  script_run = open(f_run_name,'w')
  script_run.write(run_template.format(DATASET=f_dataset_name,BSINDEX=bsindex,MERGE=f_merge_name,PLOT=f_plot_name))
  script_run.close()
  os.system("sbatch "+f_run_name)

