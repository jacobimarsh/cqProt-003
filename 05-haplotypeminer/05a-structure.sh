module load python/2.7.14
module load scipy
module load numpy/1.13.3
python /group/pawsey0149/jmarsh1/packages/fastStructure/structure.py -K 4 --input=/scratch/pawsey0149/jmarsh1/cqProt003/current_data/fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf --output=fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf_simple --prior=simple
grep '#CHROM' ../current_data/fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf >> temppopfile.txt
sed -i 's/^.*FORMAT //' temppopfile.txt
sed -i 's/  /\n/g' temppopfile.txt
echo '  p1  p4  p2  p3' >> fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf.4.meanQ_HMformatted.txt
paste temppopfile.txt fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf_simple.4.meanQ >> temp_fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf.4.meanQ_HMformatted.txt
cat temp_fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf.4.meanQ_HMformatted.txt >> fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf.4.meanQ_HMformatted.txt
rm a
rm temppopfile.txt
