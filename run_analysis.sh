for i in $(seq 1 5) ; 
do 

python measure_horizontal_displacement.py --topology ../1/prot.psf    --trajectory  ../$i/sum_full_prot.xtc --o out_$i.dat --ref_sel 'name CA and resid 418–454, 469–498, 511–532, 536–554, 563–598, 630–642, 656–712' --analysis_selections 'name CA and segid PROA and resid 114-183' 'name CA and segid PROD and resid 114-183' 'name CA and segid PROC and resid 114-183' 'name CA and segid PROB and resid 114-183'  --stride 50   --ref_top  ../1/ionized.psf --ref ../1/memb0.gro & 

python measure_vertical_displacement.py --topology ../1/prot.psf    --trajectory  ../$i/sum_full_prot.xtc --o out_$i.dat --ref_sel 'name CA and resid 418–454, 469–498, 511–532, 536–554, 563–598, 630–642, 656–712' --analysis_selections 'name CA and segid PROA and resid 114-183' 'name CA and segid PROD and resid 114-183' 'name CA and segid PROC and resid 114-183' 'name CA and segid PROB and resid 114-183'  --stride 50   --ref_top  ../1/ionized.psf --ref ../1/memb0.gro & 
    
done 

python pca.py --topology /mnt/ceph/users/mastore/unbiased_trpv1_second_round/double_toxin_bound/clustering/prot.psf --trajectory combined_trajectory.xtc --symmetry-list "segid PROD" "segid PROB" "segid PROA" "segid PROC"  -s 'backbone and segid PROA PROC PROD PROB' --save-proj --ref /mnt/ceph/users/mastore/unbiased_trpv1_second_round/double_toxin_bound/clustering/pc_analysis/reference_protein_structure.pdb  --dir-root pca/ --vis-multiply 1.5 -fn 50 -o mult_1.5 
