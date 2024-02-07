
#!/bin/bash
for obs in spherocity thrust transverse_spherocity transverse_thrust broadening isotropy; do
    python evaluate_iterations.py --config config/dataset_v7_MCA3P_unfolddata_v3_2d_inverse.json   --method multifold --migiter 0 1 --obs ${obs}:${obs}_coarse,nparticle:nparticle_fine --eff-acc
done
python merge_bs.py --input config/merge_MCA3P_unfolddata_v3_2d_inverse.txt --output results_finebin_v7_MCA3P_unfolddata_v3_2d_inverse_1/merge.csv
python plot_iter.py --obs spherocity:spherocity_coarse,nparticle:nparticle_fine thrust:thrust_coarse,nparticle:nparticle_fine transverse_spherocity:transverse_spherocity_coarse,nparticle:nparticle_fine transverse_thrust:transverse_thrust_coarse,nparticle:nparticle_fine broadening:broadening_coarse,nparticle:nparticle_fine isotropy:isotropy_coarse,nparticle:nparticle_fine --config config/plot_1d_v7_MCA3P_unfolddata_v3_2d_inverse.json
