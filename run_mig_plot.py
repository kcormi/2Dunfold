#!/bin/env python
from optparse import OptionParser
import os, sys, re
import time
import glob
parser = OptionParser()

run_template_diff = """
python plot_migration.py --inputref results_finebin_v7_MCEPOS_unfoldCP1_1d_1p3M_eff_acc_ensemble4_coverage_stresstest_v2/unfold_{OBS1}_{OBS2}_nominal_optimize_multifold_bs_{BSREF}.root --inputcompare results_finebin_v7_MCEPOS_unfoldCP1_1d_1p3M_eff_acc_ensemble4_coverage_stresstest_v2/unfold_{OBS1}_{OBS2}_nominal_optimize_multifold_bs_{BSCOMPARE}.root --histref {HISTREF} --histcompare {HISTCOMPARE} --output results_finebin_v7_MCEPOS_unfoldCP1_1d_1p3M_eff_acc_ensemble4_coverage_stresstest_v2/plots_mig/{OUTPUT} --title "{TITLE}" --row {ROW} --column {COLUMN}
"""

run_template = """
python plot_migration.py --inputref results_finebin_v7_MCEPOS_unfoldCP1_1d_1p3M_eff_acc_ensemble4_coverage_stresstest_v2/unfold_{OBS1}_{OBS2}_nominal_optimize_multifold_bs_{BSREF}.root --histref {HISTREF}  --output results_finebin_v7_MCEPOS_unfoldCP1_1d_1p3M_eff_acc_ensemble4_coverage_stresstest_v2/plots_mig/{OUTPUT} --title "{TITLE}" --row {ROW} --column {COLUMN}
""" 
compare_list = ["CP5_migration_nui1p5","CP5_migration_nui-1p5",
                "CP1_migration_nui1p5","CP1_migration_nui-1p5",
                "CP1_gen_nui1p5","CP1_gen_nui-1p5",
                "CP5_gen_nui1p5","CP5_gen_nui-1p5"]

obs1 = "nparticle"
obs2 = "spherocity"
bsref = 1
bscompare = range(2,10)
#row_coarse_bin = ["reco_ntrk90.0-160.0","reco_ntrk70.0-90.0","reco_ntrk50.0-70.0","reco_ntrk30.0-50.0","reco_ntrk20.0-30.0","reco_ntrk10.0-20.0","reco_ntrk3.0-10.0"]
#column_coarse_bin = ["gen_nch3.0-10.0","gen_nch10.0-30.0","gen_nch30.0-50.0","gen_nch50.0-80.0","gen_nch80.0-140.0"]
row_coarse_bin = ["reco_ntrk3.0-160.0"]
column_coarse_bin = ["gen_nch3.0-140.0"]
histref = ""
title = ""
title_compare = ""
for rowbin in row_coarse_bin:
  for columnbin in column_coarse_bin:
    histref += "HistMig_MC_multifold_{COLUMN}_{ROW}_iter0,".format(COLUMN=columnbin,ROW=rowbin)
    title += "Migration matrix {{SAMPLE}} (spherocity)\n{COLUMN} {ROW},".format(COLUMN=columnbin,ROW=rowbin)
    title_compare += "Migration matrix difference {{SAMPLE}} vs nominal (spherocity)\n{COLUMN} {ROW},".format(COLUMN=columnbin,ROW=rowbin)
histref = histref[:-1]
title = title[:-1]
title_compare = title_compare[:-1]
histcompare = histref
output = "mig_nch_sph_{SAMPLE}"
output_compare = "migdiff_nch_sph_MC_{SAMPLE}"
for i,compare in enumerate(compare_list):
  command_diff = run_template_diff.format(OBS1=obs1,OBS2=obs2,BSREF=1,BSCOMPARE=bscompare[i],
                                HISTREF = histref.format(SAMPLE=compare),
                                HISTCOMPARE = histcompare.format(SAMPLE=compare),
                                OUTPUT = output_compare.format(SAMPLE=compare),
                                TITLE = title_compare.format(SAMPLE=compare),
                                ROW = len(row_coarse_bin),COLUMN = len(column_coarse_bin))
  print(command_diff)
  os.system(command_diff)
  command = run_template.format(OBS1=obs1,OBS2=obs2,BSREF=bscompare[i],
                                HISTREF = histref.format(SAMPLE=compare),
                                OUTPUT = output.format(SAMPLE=compare),
                                TITLE = title.format(SAMPLE=compare),
                                ROW = len(row_coarse_bin),COLUMN = len(column_coarse_bin))
  print(command)
  os.system(command)

command = run_template.format(OBS1=obs1,OBS2=obs2,BSREF=1,
                                HISTREF = histref.format(SAMPLE="MC"),
                                OUTPUT = output.format(SAMPLE="MC"),
                                TITLE = title.format(SAMPLE="MC"),
                                ROW = len(row_coarse_bin),COLUMN = len(column_coarse_bin))
print(command)
os.system(command)
