##vcf2hapmap
sed 's/0\/1/.\/./g' protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf >> nohet_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf
sed -i 's/1\/0/.\/./g' nohet_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf
module load java
java -jar beagle.18May20.d20.jar gt=nohet_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf nthreads=22 out=nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1
bgzip -d nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf.gz
tassel-5-standalone/run_pipeline.pl -Xmx20G -vcf nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf -export nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1 -exportType HapmapDiploid
mv nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.hmp.txt nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_hmp.txt
sed -i 's/#//g' nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_hmp.txt
sed -i "s/^/Gm/" nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_hmp.txt
sed -i "s/Gmrs/rs/" nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_hmp.txt
grep -v $'\tAC\t\|\tAG\t\|\tAT\t\|\tCA\t\|\tCG\t\|\tCT\t\|\tGA\t\|\tGC\t\|\tGT\t\|\tTA\t\|\tTC\t\|\tTG\t\|\tA\t20\|\tC\t20\|\tG\t20\|\tT\t20' nohet_prefin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_hmp.txt >> nohet_fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_hmp.txt
##convert back
#tassel-5-standalone/run_pipeline.pl -Xmx20G -h nohet_fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_hmp.txt -export -exportType VCF
##get region only
head -n 1 nohet_fin_b51_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1_hmp.txt > nohet_b51_173kb_only.hmp.txt
awk '($4>31604127  && $4<31777346){print}' >> nohet_b51_173kb_only.hmp.txt
