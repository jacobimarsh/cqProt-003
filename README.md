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

(Figure S9; Table S6) Code for generating kinship, structure, hapmap and other inputs for HaplotypeMiner (05a-e); analysis run with chosen parameters (05f) and scripts to export log info to workable format from command line (05g). 05X were used to rapidly test and analyse different parameters, to aid in optimization with clutree (06)

**06-haps_clustree.R**

(Figure S10) Cluster visualization across different CT parameters, using clustree

**07-pca_umap.sh**

(Figures S11-12) Principal component analysis on individuals using factoExtra and PLINK, before visualizing haplotype cluster segmentation using UMAP

**08-haplotype_viz**

(Figures 2-5; Table S5) Custom visualization of results from haplotyping (05). Involves converting HaplotypeMiner output into inputs for visualization (05a), integrating supporting SNP info with phenotype data (05b); and the code to visualize marker/haplotype combinations (08c), summaries of marker groups in reference to protein (08d), summaries of haplotype groups in reference to oil and protein (08e), and to represent the demographic breakdowns for the different trinucleotide insertions at site 31727019 (08f). 08x was used for preliminary analysis to identify trends for haplotype and marker groups in parallel, it combines the key features of Figures 2-4 however is too large for publication, and still contains groups that were later pruned. 

**09-pixy_ntdiv.sh**

(Figures S13-14; Table S7) Variant calling from g.vcf to include invariant sites, allowing for pixy analysis, before visualizing dXY between haplotype groups as a heatmap

**10-SV_mosdepth.sh**

Mosdepth to identify the 304bp SV in 20G085100 using per base coverage information

**0X-site_phenotests.sh**

Miscelanous script used to test for allele frequencies and mean phenotype scores for populations possessing each allele at a given locus
