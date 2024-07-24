
#!/bin/bash
for obs in spherocity:spherocity_sym thrust:thrust_sym transverse_spherocity:transverse_spherocity_sym transverse_thrust:transverse_thrust_sym broadening:broadening_sym isotropy:isotropy_sym mass:mass_sym nparticle:nparticle_fine chjet_deltaphi; do
    python evaluate_iterations.py --config config/dataset_v7_MCA3P_unfoldCUEP8M1_v7.json   --method multifold --migiter 0 1 --obs nparticle,${obs} --eff-acc
done
python merge_bs.py --input config/merge_MCA3P_unfoldCUEP8M1_v7.txt --output results_finebin_v7_MCA3P_unfoldCUETP8M1_v7_sym/merge.csv
python plot_iter.py --obs nparticle,spherocity:spherocity_sym nparticle,thrust:thrust_sym nparticle,transverse_spherocity:transverse_spherocity_sym nparticle,transverse_thrust:transverse_thrust_sym nparticle,broadening:broadening_sym nparticle,isotropy:isotropy_sym nparticle,mass:mass_sym nparticle,nparticle:nparticle_fine nparticle,chjet_deltaphi --config config/plot_1d_v7_MCA3P_unfoldCUEP8M1_v7.json
