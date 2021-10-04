module load java
## fill in the third column of the original VCF file
head -n 30000 gwas_Gm20.vcf|grep "#" >header

grep -v "#" gwas_Gm20.vcf|awk '$3=$1"_"$2' |sed 's/ /\t/g' >tem.vcf

cat header tem.vcf >gwas_Gm20.id.vcf

awk '$5 !~ /([[:alpha:]])+,[[:alpha:]]/{print}' gwas_Gm20.id.vcf |grep -v "*" > gwas_Gm20.biallic.vcf

~/biosoft/plink --vcf gwas_Gm20.biallic.vcf --maf 0.01 --geno 0.1 --recode  vcf-iid --out protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1 --allow-extra-chr

#tested with maf at 0.05, same GWAS results
~/biosoft/plink --vcf gwas_Gm20.biallic.vcf --maf 0.05 --geno 0.1 --recode  vcf-iid --out protein_Gm20.SNPs.id.biallic_maf_0.05_geno_0.1 --allow-extra-chr

~/biosoft/plink --vcf gwas_Gm20.biallic.vcf --maf 0.01 --recode  vcf-iid --out protein_Gm20.SNPs.id.biallic_maf_0.01 --allow-extra-chr

#Imputation
java -jar beagle.18May20.d20.jar gt=protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf nthreads=22 out=fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01.vcf

##heterozygotes to missing + imputation
sed 's/0\/1/.\/./g' protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf >> nohet_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf
sed -i '/1\/0/.\/./g' nohet_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf
java -jar beagle.18May20.d20.jar gt=nohet_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf nthreads=22 out=nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01.vcf

mkdir bin
vcftools --vcf Gm20.INDELs.vcf --minQ 30 --out Gm20.qual.INDELS.vcf --recode
mv Gm20.qual.INDELS.vcf.log bin/
##FILTERED LINES (>30phred): 512907

vcftools --max-missing 0.9 --non-ref-af-any 0.01 --vcf Gm20.qual.INDELS.vcf.recode.vcf --recode --out qual_miss0.9_nra0.01
mv Gm20.qual.INDELS.recode.vcf bin/
mv qual_miss0.9_nra0.01.log bin/
##FILTERED LINES (<0.9 miss, nra>0.01): 430605 (429719 pheno_only)

##indels2snp reformatting

bcftools annotate -x ^FORMAT/GT qual_miss0.9_nra0.01.recode.vcf >> qual_miss_nra_gt_Gm20.INDELs.vcf
bcftools norm -m  - qual_miss_nra_gt_Gm20.INDELs.vcf | bcftools filter --include 'strlen(REF)<strlen(ALT)' >> norm_qual_miss_nra_gt_Gm20.INs.vcf
bcftools norm -m  - qual_miss_nra_gt_Gm20.INDELs.vcf | bcftools filter --include 'strlen(REF)>strlen(ALT)' >> norm_qual_miss_nra_gt_Gm20.DELs.vcf
mv qual_miss0.9_nra0.01.recode.vcf bin/
mv qual_miss_nra_gt_Gm20.INDELs.vcf bin/
##INs: 143934 (143017 pheno_only) ; OUTs: 428909 (427769 pheno_only)

bcftools norm -m + norm_qual_miss_nra_gt_Gm20.INs.vcf >> qual_miss_nra_gt_Gm20.INs.vcf
bcftools norm -m + norm_qual_miss_nra_gt_Gm20.DELs.vcf >> qual_miss_nra_gt_Gm20.DELs.vcf
mv norm_qual_miss_nra_gt_Gm20.INs.vcf bin/
mv norm_qual_miss_nra_gt_Gm20.DELs.vcf bin/
##INs: 104909 (104559 pheno_only) ; OUTs: 348347 (347696 pheno_only)

sed 's/\.\//m\//g' qual_miss_nra_gt_Gm20.DELs.vcf | sed 's/\/\./\/m/g' qual_miss_nra_gt_Gm20.DELs.vcf > mm_qual_miss_nra_gt_Gm20.DELs.vcf
sed 's/0\//n\//g' mm_qual_miss_nra_gt_Gm20.DELs.vcf | sed 's/\/0/\/n/g' > nn_mm_qual_miss_nra_gt_Gm20.DELs.vcf
sed -r 's/[0-9]+\//1\//g' nn_mm_qual_miss_nra_gt_Gm20.DELs.vcf | sed -r 's/\/[0-9]+/\/1/g' > 11_nn_mm_qual_miss_nra_gt_Gm20.DELs.vcf
sed 's/m\//\.\//g' 11_nn_mm_qual_miss_nra_gt_Gm20.DELs.vcf | sed 's/\/m/\/\./g' | sed 's/n\//0\//g' | sed 's/\/n/\/0/g' | sed 's/10\//1\//g' | sed 's/\/10/\/1/g' | sed 's/01\//1\//g' | sed 's/\/01/\/1/g' > 1alt_qual_miss_nra_gt_Gm20.DELs.vcf
sed 's/\.\t[^\t]*\t[^\t]*/.\tG\tC/' 1alt_qual_miss_nra_gt_Gm20.DELs.vcf > Aref_Talt_qual_miss_nra_gt_Gm20.DELs.vcf
bcftools annotate --set-id +'DEL_%POS' Aref_Talt_qual_miss_nra_gt_Gm20.DELs.vcf > id_Aref_Talt_qual_miss_nra_gt_Gm20.DELs.vcf
vcftools --vcf id_Aref_Talt_qual_miss_nra_gt_Gm20.DELs.vcf --maf 0.01 --recode --out maf_id_Aref_Talt_qual_miss_nra_gt_Gm20.DELs
mv maf_id_Aref_Talt_qual_miss_nra_gt_Gm20.DELs.recode.vcf filtered_biallelic_Gm20.DELs.vcf
mv qual_miss_nra_gt_Gm20.DELs.vcf bin/
mv maf_id_Aref_Talt_qual_miss_nra_gt_Gm20.DELs.log bin/
mv mm_qual_miss_nra_gt_Gm20.DELs.vcf bin/
mv nn_mm_qual_miss_nra_gt_Gm20.DELs.vcf bin/
mv 11_nn_mm_qual_miss_nra_gt_Gm20.DELs.vcf bin/
mv 1alt_qual_miss_nra_gt_Gm20.DELs.vcf bin/
mv Aref_Talt_qual_miss_nra_gt_Gm20.DELs.vcf bin/
mv id_Aref_Talt_qual_miss_nra_gt_Gm20.DELs.vcf bin/
##OUTs: 56065 (55512 pheno_only)

sed 's/\.\//m\//g' qual_miss_nra_gt_Gm20.INs.vcf | sed 's/\/\./\/m/g' qual_miss_nra_gt_Gm20.INs.vcf > mm_qual_miss_nra_gt_Gm20.INs.vcf
sed 's/0\//n\//g' mm_qual_miss_nra_gt_Gm20.INs.vcf | sed 's/\/0/\/n/g' > nn_mm_qual_miss_nra_gt_Gm20.INs.vcf
sed -r 's/[0-9]+\//1\//g' nn_mm_qual_miss_nra_gt_Gm20.INs.vcf | sed -r 's/\/[0-9]+/\/1/g' > 11_nn_mm_qual_miss_nra_gt_Gm20.INs.vcf
sed 's/m\//\.\//g' 11_nn_mm_qual_miss_nra_gt_Gm20.INs.vcf | sed 's/\/m/\/\./g' | sed 's/n\//0\//g' | sed 's/\/n/\/0/g' | sed 's/01\//1\//g' | sed 's/\/01/\/1/g' | sed -e 's/\/..$/\/1/' | sed 's/\t[^\t].\//\t1\//g' > 1alt_qual_miss_nra_gt_Gm20.INs.vcf
sed 's/\.\t[^\t]*\t[^\t]*/.\tA\tT/' 1alt_qual_miss_nra_gt_Gm20.INs.vcf > Aref_Talt_qual_miss_nra_gt_Gm20.INs.vcf
bcftools annotate --set-id +'INS_%POS' Aref_Talt_qual_miss_nra_gt_Gm20.INs.vcf > id_Aref_Talt_qual_miss_nra_gt_Gm20.INs.vcf
vcftools --recode --vcf id_Aref_Talt_qual_miss_nra_gt_Gm20.INs.vcf --out id_Aref_Talt_qual_miss_nra_gt_Gm20.INs
vcftools --vcf id_Aref_Talt_qual_miss_nra_gt_Gm20.INs.recode.vcf --maf 0.01 --recode --out maf_id_Aref_Talt_qual_miss_nra_gt_Gm20.INs
mv maf_id_Aref_Talt_qual_miss_nra_gt_Gm20.INs.recode.vcf filtered_biallelic_Gm20.INs.vcf
mv qual_miss_nra_gt_Gm20.INs.vcf bin/
mv maf_id_Aref_Talt_qual_miss_nra_gt_Gm20.INs.log bin/
mv mm_qual_miss_nra_gt_Gm20.INs.vcf bin/
mv nn_mm_qual_miss_nra_gt_Gm20.INs.vcf bin/
mv 11_nn_mm_qual_miss_nra_gt_Gm20.INs.vcf bin/
mv 1alt_qual_miss_nra_gt_Gm20.INs.vcf bin/
mv Aref_Talt_qual_miss_nra_gt_Gm20.INs.vcf bin/
mv id_Aref_Talt_qual_miss_nra_gt_Gm20.INs.vcf bin/
mv id_Aref_Talt_qual_miss_nra_gt_Gm20.INs.recode.vcf bin/
##INs: 35826 (35413 pheno_only)

grep ^20 filtered_biallelic_Gm20.INs.vcf | awk '{print $2}' > filtered_biallelic_Gm20.INs.sites.txt
grep ^20 filtered_biallelic_Gm20.DELs.vcf | awk '{print $2}' > filtered_biallelic_Gm20.DELs.sites.txt
cat filtered_biallelic_Gm20.INs.sites.txt filtered_biallelic_Gm20.DELs.sites.txt | sort -n | uniq -c | grep '2 ' | sed 's/.* //' > overlap_filtered_biallelic_Gm20.INs_DELs.txt
mv filtered_biallelic_Gm20.INs.sites.txt bin/
mv filtered_biallelic_Gm20.DELs.sites.txt bin/
wc -l overlap_filtered_biallelic_Gm20.INs_DELs.txt | awk '{print $1}' | echo Overlapping sites: $(</dev/stdin)
##OVERLAP: 4706 (13 in region) (4540 pheno_only)

bcftools query -l filtered_biallelic_Gm20.INs.vcf | sort > INs_samples.txt
bcftools view -S samples.txt protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf > sort_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf
mv protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf bin
grep '#' sort_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf > protein_INs.vcf
grep -v '#' sort_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf > protein_INs_nohead.txt
grep -v '#' filtered_biallelic_Gm20.INs.vcf >> protein_INs_nohead.txt
sort -n -k2 prot_INs_nohead.txt >> protein_INs.vcf
mv protein_INs_nohead.txt bin
grep '#' sort_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf > protein_DELs.vcf
grep -v '#' sort_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf > prot_DELs_nohead.txt
grep -v '#' filtered_biallelic_Gm20.DELs.vcf >> prot_DELs_nohead.txt
sort -n -k2 prot_DELs_nohead.txt >> protein_DELs.vcf
mv protein_DELs_nohead.txt bin
