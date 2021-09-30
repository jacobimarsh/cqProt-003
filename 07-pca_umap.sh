grep '#' fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf > fin_b51_173kb_only.vcf
#vcf2matrix
grep -m 1 -A 10000000000000 '#CHR' fin_b51_173kb_only.vcf | cut -f2,10- | sed -e 's/0|1/1/g' -e 's/1|0/1/g' -e 's/1|1/2/g' -e 's/0|0/0/g' > fin_b51_173kb_only.mtx

#UMAP
plink --pca header --vcf region_only_fin.vcf
#had to change np.float to np.float64 in the py script from https://github.com/diazale/gt-dimred/blob/master/scripts/general_umap_script.py
python gt-dimred/scripts/general_umap_script.py -dset plink.eigenvec -pc 7 -nn 15 -md 0.001 -outdir . -head T -log .
echo -e 'ID\tX\tY\tHap' > umap_ids.txt
awk '{print $1}' plink.eigenvec | tail -n +2 | paste - plink.eigenvec_UMAP_PC7_NC2_NN15_MD0.001_euclidean_20210930_153346.txt | sed 's/ /\t/' >> umap_ids.txt
#then manually put haplotype from Table S6 into the hap column using excel or other...

#Visualization
R
library(tidyverse)
library(factoextra)

regmtx <- read.csv("fin_b51_173kb_only.mtx", sep = "\t")
regmtx_2 <- regmtx[,-1]
rownames(regmtx_2) <- regmtx[,1]
res.pca.mark <- prcomp(t(regmtx_2), scale = T)
pca_eig <- fviz_eig(res.pca.mark, addlabels = T)
##visualize PCAs with below
#fviz_pca_ind(res.pca.mark,
#             col.ind = "contrib", # Color by the quality of representation
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE     # Avoid text overlapping
#             )

ggsave('pca_eig.pdf',
       pca_eig, 
       device = 'pdf', 
       dpi = 1800,
       height = 150,
       width = 150,
       units = "mm")

umap_ids3 <- read.csv("umap_ids.txt", sep = "\t")
cqUMAP <- ggplot(umap_ids3, aes(x = X, y = Y, color = Hap)) +
  geom_point() +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_minimal()
cqUMAP

ggsave('cqUMAP.pdf',
       cqUMAP, 
       device = 'pdf', 
       dpi = 1800,
       height = 150,
       width = 150,
       units = "mm")

##FOR SNPS (not included in manuscript), from write.csv(res.pca$x, file = "res.pca.txt") taken from PCA make via prcomp
#sed 's/,/ /g' res.pca.txt |awk 'FS=" "{print $1}' | paste -d " " - res.pca.txt | sed -e 's/,/ /g' -e 's/"//g' -e 's/  PC/FID IID PC/' | cut -d " " -f 1-22 > snpPCs.txt
