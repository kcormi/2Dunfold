
#!/bin/bash
for obs in spherocity thrust transverse_spherocity transverse_thrust broadening isotropy nparticle:nparticle_fine mass; do
    python evaluate_iterations.py --config config/dataset_v7_MCA3P_sysweightEPOSgenweightA3P_v3.json  --method multifold --migiter 0 1 --obs nparticle,${obs} --step1
done
python merge_bs.py --input config/merge_MCA3P_sysweightEPOSgenweightA3P_v3.txt --output results_finebin_v7_MCA3P_sysweightEPOSgenweightA3P_1d_ensemble4_v3/merge.csv
python plot_iter.py --obs nparticle,spherocity nparticle,thrust nparticle,transverse_spherocity nparticle,transverse_thrust nparticle,broadening nparticle,isotropy nparticle,mass nparticle,nparticle:nparticle_fine --config config/plot_1d_v7_MCA3P_sysweightEPOSgenweightA3P.json
