
#plink --pca header --vcf region_only_fin.vcf
##had to change np.float to np.float64 in the py script from https://github.com/diazale/gt-dimred/blob/master/scripts/general_umap_script.py
#python gt-dimred/scripts/general_umap_script.py -dset plink.eigenvec -pc 5 -nn 10 -md 0.001 -outdir . -head T -log .
#awk '{print $1}' plink.eigenvec | tail -n +2 | paste - plink.eigenvec_UMAP_PC5_NC2_NN10_MD0.001_euclidean_20210823_145238.txt > umap_ids.txt
library(tidyverse)
library(factoextra)

res.pca.mark <- prcomp(t(regmtx), scale = FALSE)
fviz_eig(res.pca.mark, addlabels = T)
##visualize PCAs with below
#fviz_pca_ind(res.pca.mark,
#             col.ind = "cos2", # Color by the quality of representation
#             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#             repel = TRUE     # Avoid text overlapping
#             )

umap_ids <- read.csv("C:/Users/21485753/Desktop/cqProt003/PCA/umap_ids.txt", sep = "\t")
cqUMAP <- ggplot(umap_ids, aes(x = X, y = Y, color = Hap)) +
  geom_point() +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  theme_minimal()

ggsave('cqUMAP.pdf',
       cqUMAP, 
       device = 'pdf', 
       dpi = 1800,
       height = 150,
       width = 150,
       units = "mm")

##FOR SNPS, from write.csv(res.pca$x, file = "res.pca.txt") taken from PCA make via prcomp
#sed 's/,/ /g' res.pca.txt |awk 'FS=" "{print $1}' | paste -d " " - res.pca.txt | sed -e 's/,/ /g' -e 's/"//g' -e 's/  PC/FID IID PC/' | cut -d " " -f 1-22 > snpPCs.txt

