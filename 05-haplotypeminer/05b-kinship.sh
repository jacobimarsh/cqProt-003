 module load java
 /group/pawsey0149/jmarsh1/packages/tassel-5-standalone/run_pipeline.pl -Xmx100g -h fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.hmp.txt -KinshipPlugin -method Centered_IBS -endPlugin -export fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_kinship.txt -exportType SqrMatrix
 touch temp_empty.txt
 tail -n +4 fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_kinship.txt | awk '{print $1}' | tr "\n" "\t" >> temp_head.txt
 paste temp_empty.txt temp_head.txt >> fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf_kinship_HMformatted.txt
 tail -n +4 fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_kinship.txt >> fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf_kinship_HMformatted.txt
 rm temp_empty.txt
 rm temp_head.txt
