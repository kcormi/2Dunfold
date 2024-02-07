# Evaluate the unfolding results

To setup the local environment you need python3, which can be taken from a recent LCG release. We use:
```
source /cvmfs/sft.cern.ch/lcg/views/LCG_102rc1/x86_64-centos7-gcc11-opt/setup.sh
```
For the first time of run, install the `catpy` for calculation the ks-distance.
```
python -m venv catpy_env
git clone ssh://git@gitlab.cern.ch:7999/kcormier/catpy.git
cd catpy
python -m pip install .
cd ..
```
Next time do
```
source catpy_env/bin/activate
```
To produce histograms, you can run (showing spherocity as an example observable):
Reminder: each time you run this command, it will add new entries to the output dataframe, so in the end the histograms for different observables will be in the file. Issues may come out if you repeat the processing for the same observable because of duplications. To deleate the existing csv file and overwrite, use the `--df-overwrite` flag for the script.
```
python evaluate_iterations.py --config config/dataset_v7_MCA3P_unfolddata_v3.json   --method multifold --migiter 0 1 --obs nparticle,spherocity --eff-acc
```
To merge the results of systematic variations and get full uncertainty, or merge the results of toy experiments and calculate the covariance matrices:
```
python merge_bs.py --input config/merge_MCA3P_unfolddata_v3.txt --output results_finebin_v7_MCA3P_unfolddata_v3_1/merge.csv
```
To plot from the DataFrames:
```
python plot_iter.py --obs nparticle,spherocity --config config/plot_1d_v7_MCA3P_unfolddata_v3.json
```
To process and plot all the observables (remember to clean the output directory before):
```
bash run_1D_MCA3P_unfolddata_v3_sys.sh 
```
The first step of processing might take large memory when processing the bootstraps with multiprocessing. Batch jobs can be submitted  instead:
```
sbatch run_1D_MCA3P_unfolddata_v3_sys_batcharray.sh
```
to specify a different binning configuration for a variable simply add the name of the binning after the variable separated by a colon `:`. e.g. using `--obs nparticle,spherocity:coarse_sph_binning`, will use the default binning for `nparticle` as defined in the config file, and the `coarse_sph_binning` for `spherocity`. In the case `coarse_sph_binning` should be defined in the configuration file under the `binnings` heading.

Example:
```
python evaluate_iterations.py --config config/dataset_v7_MCA3P_unfolddata_v3.json   --method multifold --migiter 0 1 --obs nparticle,spherocity:spherocity_sym --eff-acc
```

To run the 2D unfolding of event shapes in slices of Nch:
```
bash run_1D_MCA3P_unfolddata_v3_2d.sh
```
The first step of processing in 2D can take a long time, better to run with batch:
```
sbatch run_1D_MCA3P_unfolddata_v3_2d_batcharray.sh 
```
To run the 2D unfolding of Nch in slices of event shapes:
```
bash run_1D_MCA3P_unfolddata_v3_2d_inverse.sh 
```
To submit the first step to batch:
```
sbatch run_1D_MCA3P_unfolddata_v3_2d_inverse_batcharray.sh 
```
# Validations

## Plots of reweighting for systematic variations
To run test on the gen-level reweighting
```
bash run_1D_MCA3P_genweightCH3.sh
bash run_1D_MCA3P_genweightCP1.sh
bash run_1D_MCA3P_genweightEPOS.sh
```
To run test on the gen- and reco-level reweighting from nomial MC to systematic variaitons.
```
bash run_1D_MCA3P_sysweightA3Ptrackdrop.sh
bash run_1D_MCA3P_sysweightCH3genweightA3P.sh
bash run_1D_MCA3P_sysweightCP1genweightA3P.sh
bash run_1D_MCA3P_sysweightEPOSgenweightA3P.sh
```

## Bias and coverage test:
To process the toy experiments:
```
python run_coverage_process.py
```
Then extract bias and coverage from the toys:
```
python coverage.py --input config/merge_MCA3P_unfoldCUEP8M1_maxweight10_minweightp1_v3_stat_sys_bs.txt --output results_finebin_v7_MCA3P_unfoldCUEP8M1_maxweight10_minweightp1_coverage_v3_1/merge_bs_dataunc.csv --dataunc
```
For plotting
```
python plot_bias.py 
```

## Unfold pseudo-data
Test the 1D unfolding with pseudo-data from A14, CP4 and CUETP8M1 tunes (also output plots of uncertainty decompostions and bottomline tests):
```
bash run_1D_MCA3P_unfoldA14_v3.sh
bash run_1D_MCA3P_unfoldCP5_v3.sh
bash run_1D_MCA3P_unfoldCUEP8M1_v3.sh
```
2D unfold of pseudo-data:
```
bash run_1D_MCA3P_unfoldCUEP8M1_v3_2d.sh
```
Unfold CP5 pseudo-data with alternative systematic templates from A14, CUETP8M1 and EPOS instead of nominal choices of CP1, Herwig CH3 and EPOS
```
bash run_1D_MCA3P_unfoldCP5_othersys_v3.sh
```

