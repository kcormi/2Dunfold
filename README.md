# Evaluate the unfolding results

To setup the local environment you need python3, which can be taken from a recent LCG release. We use:
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_102rc1/x86_64-centos7-gcc11-opt/setup.sh
```

To produce histograms and plots, you can run (showing spherocity as an example observable:
```
python evaluate_iterations.py --config config/dataset_v7_MCEPOS_unfoldCP1_optimize_1p3M.json  --method multifold --migiter 0 1 2 10 20 --eff-acc --obs nparticle,spherocity
python evaluate_chi2_iterations.py --input results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc_ensemble4/unfold_nparticle_spherocity_nominal_optimize_multifold.root --output results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc_ensemble4/iter_nparticle_spherocity_nominal_optimize_multifold.root --config config/plot_1d_v7_MCEPOS_unfoldCP1.json --plot --plotdir results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc_ensemble4/plots_multifold/ --obs nparticle,spherocity
```

to run for all observables you can do:

```
for obs in spherocity thrust transverse_spherocity transverse_thrust broadening isotropy; do
    python evaluate_iterations.py --config config/dataset_v7_MCEPOS_unfoldCP1_optimize_1p3M.json    --method multifold --migiter 0 1 2 10 20 --eff-acc --obs nparticle,${obs}
    python evaluate_chi2_iterations.py --input results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc_ensemble4/unfold_nparticle_${obs}_nominal_optimize_multifold.root --output results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc_ensemble4/iter    _nparticle_${obs}_nominal_optimize_multifold.root --config config/plot_1d_v7_MCEPOS_unfoldCP1.json  --plot --plotdir results_finebin_v7_MCEPOS_sysCP1_1d_1p3M_eff_acc_ensemble4/plots_multifold/ --obs nparticle,${obs}
done
```

