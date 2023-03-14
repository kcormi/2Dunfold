# Evaluate the unfolding results
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_102rc1/x86_64-centos7-gcc11-opt/setup.sh
python evaluate_iterations.py --config Config_checkunfold/Config_sph_1d_v7_MCEPOS_unfoldCP1_optimize_1p3M.json  --method multifold --migiter 0 1 2 10 20 --eff-acc
python evaluate_chi2_iterations.py --input results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc_ensemble4/unfold_nparticle_eta2p4pt05_pur_1d_spherocity_eta2p4pt05_pur_1d_nominal_optimize_multifold.root --output results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc_ensemble4/iter_nparticle_eta2p4pt05_pur_1d_spherocity_eta2p4pt05_pur_1d_nominal_optimize_multifold.root --config Config_plot/Config_sph_1d_v7_MCEPOS_unfoldCP1.json --plot --plotdir results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc/plots_multifold/
```
