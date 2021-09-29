module load java
## fill in the third column of the original VCF file
head -n 30000 gwas_Gm20.vcf|grep "#" >header

grep -v "#" gwas_Gm20.vcf|awk '$3=$1"_"$2' |sed 's/ /\t/g' >tem.vcf

cat header tem.vcf >gwas_Gm20.id.vcf

awk '$5 !~ /([[:alpha:]])+,[[:alpha:]]/{print}' gwas_Gm20.id.vcf |grep -v "*" > gwas_Gm20.biallic.vcf

~/biosoft/plink --vcf gwas_Gm20.biallic.vcf --maf 0.01 --geno 0.1 --recode  vcf-iid --out protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1 --allow-extra-chr

~/biosoft/plink --vcf gwas_Gm20.biallic.vcf --maf 0.05 --geno 0.1 --recode  vcf-iid --out protein_Gm20.SNPs.id.biallic_maf_0.05_geno_0.1 --allow-extra-chr

~/biosoft/plink --vcf gwas_Gm20.biallic.vcf --maf 0.01 --recode  vcf-iid --out protein_Gm20.SNPs.id.biallic_maf_0.01 --allow-extra-chr
