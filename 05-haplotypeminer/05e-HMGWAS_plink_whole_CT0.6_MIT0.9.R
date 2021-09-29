#Final haplotypeminer parameters used
library(snpStats)
library(HaplotypeMiner)
library(ggplot2)

paramsHMGWAS_plink_whole <- haplo_params(
   input_file = '../input/nohet_b51_173kb_only.hmp.txt',
   gene_db_file = '../input/my_gene_db.txt',
   chr_db_file = '../input/chr20_size.txt',
   structure_file = '../input/fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf.4.meanQ_HMformatted.txt',
   kinship_file = '../input/fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf_kinship_HMformatted.txt',
   gene_name = 'GWAS_plink_whole',
   R2_measure = 'r2vs',
   cluster_R2 = 'r2',
   max_missing_threshold = 0.6,
   max_het_threshold = 0.3,
   min_alt_threshold = 0.01,
   min_allele_count = 4,
   cluster_threshold = 0.6,
   max_marker_to_gene_distance = 2000000,
   max_flanking_pair_distance = 2*2000000,
   marker_independence_threshold = 0.9)

HMGWAS_plink_whole_haplotypes <- haplo_selection(paramsHMGWAS_plink_whole, verbose = TRUE)
str(HMGWAS_plink_whole_haplotypes, max.level = 1)
HMGWAS_plink_whole_haplotypes

graph_list_HMGWAS_plink_whole <- list('All_markers' = 'density',
                   'Filtered_markers' = c('matrix', 'distance', 'genotypes'),
                   'Clustered_markers' = c('matrix', 'genotypes'),
                   'Selected_clusters' = c('matrix', 'genotypes'),
                   'Selected_markers' = c('matrix', 'genotypes'),
                    'Haplotypes' = c('genotypes'))
haplo_output(HMGWAS_plink_whole_haplotypes, output_dir = '../analysis_output/HMGWAS_plink_whole_CT0.6_MIT0.9', graphs = graph_list_HMGWAS_plink_whole)

haplo_logfile(HMGWAS_plink_whole_haplotypes, to_file = TRUE)

genotype_plot(snp_data = HMGWAS_plink_whole_haplotypes$Filtered_markers,
              gene_pos = HMGWAS_plink_whole_haplotypes$Parameters$Gene_center,
              kept_markers = HMGWAS_plink_whole_haplotypes$Haplotypes$Markers,
              assignment = HMGWAS_plink_whole_haplotypes$Haplotypes$Assignment,
              name_order = FALSE) +
  theme(axis.text = element_text(size = 0.1))
                                                                                                                                                                                                                                             genotype_plot(snp_data = HMGWAS_plink_whole_haplotypes$Haplotypes,                                                                                                                                                                                         gene_pos = HMGWAS_plink_whole_haplotypes$Parameters$Gene_center,
              kept_markers = HMGWAS_plink_whole_haplotypes$Haplotypes$Markers,                                                                                                                                                                             assignment = NULL,
              name_order = TRUE)                                                                                                                                                                                                                                                                                                                                                                                                                                                          density_plot(snp_data = HMGWAS_plink_whole_haplotypes$All_markers,
             center_pos = HMGWAS_plink_whole_haplotypes$Parameters$Gene_center,
             chr_length = HMGWAS_plink_whole_haplotypes$Parameters$Chromosome_length)
                                                                                                                                                                                                                                             ld_plot(snp_data = HMGWAS_plink_whole_haplotypes$Filtered_markers,                                                                                                                                                                                   center_pos = HMGWAS_plink_whole_haplotypes$Parameters$Gene_center,
        kept_markers = HMGWAS_plink_whole_haplotypes$Haplotypes$Markers) +
  theme(axis.text = element_text(size = 5))
  
 ld_plot(snp_data = HMGWAS_plink_whole_haplotypes$Clustered_markers,                                                                                                                                                                                  center_pos = HMGWAS_plink_whole_haplotypes$Parameters$Gene_center,
        kept_markers = HMGWAS_plink_whole_haplotypes$Haplotypes$Markers)

distance_plot(snp_data = HMGWAS_plink_whole_haplotypes$Filtered_markers,
              center_pos = HMGWAS_plink_whole_haplotypes$Parameters$Gene_center,
              r2_threshold3 = HMGWAS_plink_whole_haplotypes$Parameters$Marker_independence_threshold)
