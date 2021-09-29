##blocks
plink --vcf fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf --blocks no-pheno-req --allow-extra-chr --out fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf --blocks-max-kb 10000

##R2 with GWAS-SNP
echo 20_31632556 > 20_31632556.tmp
plink --vcf ../protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf --show-tags 20_31632556.tmp --list-all --tag-r2 0.9 --out 20_31632556
plink --r2 dprime --ld-snp-list 20_31632556.tags --vcf ../protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf --out 20_31632556_tagsdp --with-freqs --ld-window 100000000 --ld-window-kb 1300
echo -e 'POS\tMAF\tR2\tDP' > 20_31632556_tagsum.txt
grep 31632556 20_31632556_tagsdp.ld | sed 's/ \+/\t/g' | sed 's/20\t31632556\t[^\t]*\t[^\t]*\t//' | awk '$6>0.99||$5>0.9 {print $3,$4,$5,$6}' | sed 's/ /\t/g' | sort -r -n -k3 >> 20_31632556_tagsum.txt

#added a column to tagsum with whether it in haplotype or not (1|0) under a column called GWAS
R
library(tidyverse)
library(viridis)

gwas_tagsum <- read.csv("20_31632556_tagsum.txt", sep = '\t', header = T)

gwas_r2 <- gwas_tagsum %>% ggplot(aes(x = POS, y = R2, colour = GWAS)) +
  geom_point() +
  scale_color_gradient(low = 'black',
                      high = 'blue') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_vline(xintercept = c(31604127 ,31777346),
             color = "red",
             linetype = "dashed") +
  ylab(expression(paste("R"^"2"~" Linkage with GWAS-SNP"))) +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Position on Chr20")

gwas_dp <- gwas_tagsum %>% ggplot(aes(x = POS, y = DP, colour = GWAS)) +
  geom_point() +
  scale_color_gradient(low = 'black',
                      high = 'blue') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  geom_vline(xintercept = c(31604127 ,31777346),
             color = "red",
             linetype = "dashed") +
  ylab("D' Linkage with GWAS-SNP") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlab("Position on Chr20")


ggsave('gwas_r2.pdf',
       gwas_r2, 
       device = 'pdf', 
       dpi = 1800,
       height = 90,
       width = 160,
       units = "mm")

ggsave('gwas_dp.pdf',
       gwas_dp, 
       device = 'pdf', 
       dpi = 1800,
       height = 90,
       width = 160,
       units = "mm")
