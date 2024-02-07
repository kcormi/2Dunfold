
#!/bin/bash
for obs in thrust transverse_spherocity transverse_thrust broadening isotropy mass; do
    python evaluate_iterations.py --config config/dataset_v7_MCA3P_unfolddata_v3_2d.json   --method multifold --migiter 0 1 --obs nparticle:nparticle_dim1_2D,${obs} --eff-acc
done
python merge_bs.py --input config/merge_MCA3P_unfolddata_v3_2d.txt --output results_finebin_v7_MCA3P_unfolddata_v3_2d_1/merge.csv
python plot_iter.py --obs nparticle:nparticle_dim1_2D,spherocity nparticle:nparticle_dim1_2D,thrust nparticle:nparticle_dim1_2D,transverse_spherocity nparticle:nparticle_dim1_2D,transverse_thrust nparticle:nparticle_dim1_2D,broadening nparticle:nparticle_dim1_2D,isotropy nparticle:nparticle_dim1_2D,mass --config config/plot_1d_v7_MCA3P_unfolddata_v3_2d.json
