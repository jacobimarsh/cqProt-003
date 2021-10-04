# cqProt-003

**00-filtering_indel2snp.sh**

From raw vcf to filtered, imputed SNPs and filtered, recoded InDels

**01-GWAS.R**

(Figures S1-4; Tables S1,S3) GWAS analysis example conducted using rMVP used for each input vcf

**02-missingness.R**

(Figure S5) Per-site missingness information from vcftools, before ggplotting

**03-LDBlockShow.sh**

(Figure 1; Figure S6) Inputs and LDBlockShow code for whole region linkage

**04-plink_linkage.sh**

(Figures S7-8; Table S4) Linkage of sites with GWAS-SNP, before ggplotting

**05-haplotypeminer**

(Figure S9; Table S6) Code for generating kinship, structure, hapmap and other inputs for HaplotypeMiner (05a-e); analysis run with chosen parameters (05f) and scripts to export log info to workable format from command line (05g). 05X were used to rapidly test and analyse different parameters, to aid in optimization with clutree (06). 

**06-haps_clustree.R**

**07-pca_umap.sh**

**08-haplotype_viz**

(Figures 2-5; Table S5

**09-pixy_ntdiv.sh**

**10-SV_mosdepth.sh**

**0X-site_phenotests.sh**
