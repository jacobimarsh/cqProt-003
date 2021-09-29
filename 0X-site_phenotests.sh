#!/bin/bash
NAME=${1?Error: no input site}
INPUT=${2?Error: no input vcf}
 mkdir "${1}indel_pheno"
 grep -m 1 '#CHROM' ${2} > ${1}indel_pheno/${1}.vcf
 grep -m 1 ${1} ${2} >> ${1}indel_pheno/${1}.vcf
 awk '
 {
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
 }
 NF>p { p = NF }
 END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
 }' ${1}indel_pheno/${1}.vcf > ${1}indel_pheno/${1}_xpose.txt
 grep '0/0' ${1}indel_pheno/${1}_xpose.txt | sed 's/ .*//' > ${1}indel_pheno/${1}_ref_id.txt
 for b in {1..15}
 do
 grep "${b}/${b}" ${1}indel_pheno/${1}_xpose.txt | sed 's/ .*//' > ${1}indel_pheno/${1}_${b}_id.txt
 grep "\./${b}" ${1}indel_pheno/${1}_xpose.txt | sed 's/ .*//' > ${1}indel_pheno/${1}_hetmiss_id.txt
 grep "${b}/\." ${1}indel_pheno/${1}_xpose.txt | sed 's/ .*//' >> ${1}indel_pheno/${1}_hetmiss_id.txt
 done
 grep '\./\.' ${1}indel_pheno/${1}_xpose.txt | sed 's/ .*//' > ${1}indel_pheno/${1}_miss_id.txt
 grep '\./0' ${1}indel_pheno/${1}_xpose.txt | sed 's/ .*//' >> ${1}indel_pheno/${1}_hetmiss_id.txt
 grep '0/\.' ${1}indel_pheno/${1}_xpose.txt | sed 's/ .*//' >> ${1}indel_pheno/${1}_hetmiss_id.txt
 cat ${1}indel_pheno/*id.txt | grep -v -f - ${1}indel_pheno/${1}_xpose.txt | grep '/' | sed 's/ .*//' >> ${1}indel_pheno/${1}_het_id.txt
 echo -e 'Site\tAlleles\tnInd\tnIndPheno\tAvPheno' > ${1}indel_pheno/${1}_phensum.txt
 echo ${1} > ${1}indel_pheno/${1}_reftemp.txt
 grep 'REF ' ${1}indel_pheno/${1}_xpose.txt | sed 's/REF /REF_/' | sed -r 's/(.*)/\1\1/' >> ${1}indel_pheno/${1}_reftemp.txt
 wc -l ${1}indel_pheno/${1}_ref_id.txt | sed 's/ .*//g' >> ${1}indel_pheno/${1}_reftemp.txt
 grep -f ${1}indel_pheno/${1}_ref_id.txt allpheno.txt | awk '{ total += $2; count++ } END { print count, total/count }' | sed 's/ /\n/' >> ${1}indel_pheno/${1}_reftemp.txt
 cat ${1}indel_pheno/${1}_reftemp.txt | tr '\n' '\t' >> ${1}indel_pheno/${1}_phensum.txt
 echo '' >> ${1}indel_pheno/${1}_phensum.txt
 for c in {1..15}
 do
 echo ${1} > ${1}indel_pheno/${1}_${c}.txt
 grep 'ALT ' ${1}indel_pheno/${1}_xpose.txt | sed 's/ALT /ALT_/' | awk -F"[,]" "{print \$$c}" >> ${1}indel_pheno/${1}_${c}.txt
 wc -l ${1}indel_pheno/${1}_${c}_id.txt | sed 's/ .*//g' >> ${1}indel_pheno/${1}_${c}.txt
 grep -f ${1}indel_pheno/${1}_${c}_id.txt allpheno.txt | awk '{ total += $2; count++ } END { print count, total/count }' | sed 's/ /\n/' >> ${1}indel_pheno/${1}_${c}.txt
 cat ${1}indel_pheno/${1}_${c}.txt | tr '\n' '\t' >> ${1}indel_pheno/${1}_phensum.txt
 echo '' >> ${1}indel_pheno/${1}_phensum.txt
 done
 echo ${1} > ${1}indel_pheno/${1}_hettemp.txt
 echo 'HET' >> ${1}indel_pheno/${1}_hettemp.txt
 wc -l ${1}indel_pheno/${1}_het_id.txt | sed 's/ .*//g' >> ${1}indel_pheno/${1}_hettemp.txt
 grep -f ${1}indel_pheno/${1}_het_id.txt allpheno.txt | awk '{ total += $2; count++ } END { print count, total/count }' | sed 's/ /\n/' >>  ${1}indel_pheno/${1}_hettemp.txt
 cat ${1}indel_pheno/${1}_hettemp.txt | tr '\n' '\t' >> ${1}indel_pheno/${1}_phensum.txt
 echo '' >> ${1}indel_pheno/${1}_phensum.txt
 echo ${1} > ${1}indel_pheno/${1}_hetmisstemp.txt
 echo './N' >> ${1}indel_pheno/${1}_hetmisstemp.txt
 wc -l ${1}indel_pheno/${1}_hetmiss_id.txt | sed 's/ .*//g' >> ${1}indel_pheno/${1}_hetmisstemp.txt
 grep -f ${1}indel_pheno/${1}_hetmiss_id.txt allpheno.txt | awk '{ total += $2; count++ } END { print count, total/count }' | sed 's/ /\n/' >> ${1}indel_pheno/${1}_hetmisstemp.txt
 cat ${1}indel_pheno/${1}_hetmisstemp.txt | tr '\n' '\t' >> ${1}indel_pheno/${1}_phensum.txt
 echo '' >> ${1}indel_pheno/${1}_phensum.txt
 echo ${1} > ${1}indel_pheno/${1}_misstemp.txt
  echo './.' >> ${1}indel_pheno/${1}_misstemp.txt
 wc -l ${1}indel_pheno/${1}_miss_id.txt | sed 's/ .*//g' >> ${1}indel_pheno/${1}_misstemp.txt
 grep -f ${1}indel_pheno/${1}_miss_id.txt allpheno.txt | awk '{ total += $2; count++ } END { print count, total/count }' | sed 's/ /\n/' >>  ${1}indel_pheno/${1}_misstemp.txt
 cat ${1}indel_pheno/${1}_misstemp.txt | tr '\n' '\t' >> ${1}indel_pheno/${1}_phensum.txt
 echo '' >> ${1}indel_pheno/${1}_phensum.txt
 grep -v -P "\t\t0" ${1}indel_pheno/${1}_phensum.txt >> ${1}_phenosum.txt
 column -t ${1}_phenosum.txt
