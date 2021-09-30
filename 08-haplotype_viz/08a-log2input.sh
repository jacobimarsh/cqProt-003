##site

grep 'assigned to haplotype' Log.txt | sed 's/.*type //g' | sed 's/ : /,/' | sed 's/ *//' >Upsethapfile.txt

for i in {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O}; do grep -m 1 -A 4 ${i} HMGWAS_plink_whole_CT0.6_MIT0.9_breakdown.txt | grep -v Hap | sed -z 's/\n/,/g' | sed -z 's/,$/\n/'; done | paste -d , Upsethapfile.txt - | sed 's/ //g' > U_S_hapfilepre.txt

sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$3) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^B/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$4) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^C/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$5) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^D/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$6) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^E/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$7) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^F/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$8) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^G/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$9) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^H/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$10) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^I/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$11) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^J/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$12) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^K/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$13) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^L/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$14) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^M/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$15) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^N/ s/$/,{}/' U_S_hapfilepre.txt
sed 's/.\/.*NA\t//' Haplotypes.hmp.txt | awk '{if ($2!=$16) {print $1}}' | sed -z 's/\n/,/g' | sed 's/rs,//' | sed 's/,$//' | xargs -I {} sed -ie '/^O/ s/$/,{}/' U_S_hapfilepre.txt
#NOTE: Change ,'s to ;'s for Gm20_  and remove Gm20

##allele

grep -A 50 '$Haplotypes$Markers' HMGWAS_plink_whole_CT0.6MIT0.9.o | grep  | sed 's/._//' | sed 's/ .//' > U_S_allele.txt

cat U_S_allele.txt | while read line; do head -n 2 ${line}_phenosum.txt | tail -n 1 | awk '{print $5}' | xargs -I [] python -c "print([] - $(head -n 3 ${line}_phenosum.txt | tail -n 1 | awk '{print $5}'))";  done | paste -d "," U_S_allele.txt -  >> U_S_allele_temp.txt

cat U_S_allele.txt | while read line; do head -n 3 ${line}_phenosum.txt| tail -n 2 | awk '{sum+=$4} END {print sum}'; done | paste -d "," U_S_allele_temp.txt - > U_S_allele_temp2.txt
grep -m 1 '#CHR' protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf | sed 's/.*FORMAT\t//' | wc -w | xargs -I [] awk -F "," '{MB = $3; print MB/[]}' U_S_allele_temp2.txt | paste -d "," U_S_allele_temp2.txt - > U_S_allele_temp3.txt

cat U_S_site.txt | while read line ; do awk -v var="$line" -F '[[:space:]]+' '$1==var {print $4}' U_S_site_tags.tags.list; done | paste -d "," U_S_allele_temp3.txt - > U_S_allele_fin.txt

#NOTE: Add cleaned up phenotest results for set_size, filled out table with 0's/NAs, to commas, remove final comma, remove middle headers, add Type column

echo "Hap,ID,Prot" > allhappheno.txt
for i in {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O}; do cat hap${i}_pheno.txt | sed -e "s/^/${i}\t/" | sed 's/\t/,/g' ; done >> allhappheno.txt

#need to add lr/mc/oc etc.
