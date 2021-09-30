plink --pca header --vcf region_only_fin.vcf
#had to change np.float to np.float64 in the py script from https://github.com/diazale/gt-dimred/blob/master/scripts/general_umap_script.py
python gt-dimred/scripts/general_umap_script.py -dset plink.eigenvec -pc 5 -nn 10 -md 0.001 -outdir . -head T -log .
awk '{print $1}' plink.eigenvec | tail -n +2 | paste - plink.eigenvec_UMAP_PC5_NC2_NN10_MD0.001_euclidean_20210823_145238.txt > umap_ids.txt

#FOR SNPS, from write.csv(res.pca$x, file = "res.pca.txt") taken from PCA make via prcomp
sed 's/,/ /g' res.pca.txt |awk 'FS=" "{print $1}' | paste -d " " - res.pca.txt | sed -e 's/,/ /g' -e 's/"//g' -e 's/  PC/FID IID PC/' | cut -d " " -f 1-22 > snpPCs.txt

