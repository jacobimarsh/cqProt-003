##PREP
echo '"CHROM" "POS" "p-value"' >> ingwas.txt
awk -F "," '31500000<$3 && $3<31950000 {print $2,$1,$8}' protein.FarmCPU.csv | sed 's/"//g' >> ingwas.txt

##RUN
LDBlockShow/bin/LDBlockShow -InVCF fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf.gz -OutPut GWAS_imp_R2 -OutPdf -Region 20:31500000:31950000 -SeleVar 2 -InGWAS ingwas.txt -InGFF in.gff3 -TopSite
LDBlockShow/bin/LDBlockShow -InVCF fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf.gz -OutPut GWAS_imp_DP_CI -OutPdf -Region 20:31550000:31850000 -SeleVar 1 -BlockType 1 -InGWAS ingwas.txt -InGFF in.gff3 -TopSite
##clean up in illustrator
