#!/usr/bin/env bash
cd analysis_output/060421/
for FILE in *; do mkdir ${FILE}/groups; done
for FILE in *; do scp Rscripts/logfiles/${FILE}.o ${FILE}/; done
for FILE in *; do grep -A 10000 '$Haplotypes$Assignment' ${FILE}/${FILE}.o > ${FILE}/${FILE}_temp_assign.txt; awk '{print $2,$3}' ${FILE}/${FILE}_temp_assign.txt > ${FILE}/${FILE}_assignment.txt; done
for FILE in *; do
for i in {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BB,CC,DD,EE,FF,GG,HH,II,JJ,KK,LL,MM,NN,OO,PP,QQ,RR,SS,TT,UU,VV,WW,XX,YY,ZZ}; do grep -E ' '${i}$ ${FILE}/${FILE}_assignment.txt | awk '{print $1}' > ${FILE}/groups/hap${i}.txt;
if [ -s ${FILE}/groups/hap${i}.txt ]
then
grep -f ${FILE}/groups/hap${i}.txt ../../allpheno_nodash.txt > ${FILE}/groups/hap${i}_pheno.txt
wc -l ${FILE}/groups/hap${i}_pheno.txt > ${FILE}/groups/hap${i}_tempphenonum.txt
awk '{total += $2; count++} END {print total/count}' ${FILE}/groups/hap${i}_pheno.txt > ${FILE}/groups/hap${i}_tempphenoav.txt
echo ${i} | paste - ${FILE}/groups/hap${i}_tempphenonum.txt ${FILE}/groups/hap${i}_tempphenoav.txt | sed 's/ HM/\t/g' | sed 's/\/groups.*txt//g' > ${FILE}/groups/hap${i}phenohaptemp.txt
awk '{print $3,$1,$2,$4}' ${FILE}/groups/hap${i}phenohaptemp.txt | sed 's/ /\t/g' >> phenohap.txt
fi
done
done

##can use below for getting some country of origin info
#for FILE in *
#do
#for i in {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BB,CC,DD,EE,FF,GG,HH,II,JJ,KK,LL,MM,NN,OO,PP,QQ,RR,SS,TT,UU,VV,WW,XX,YY,ZZ};
#do
#grep -f ${FILE}/groups/hap${i}.txt ../../Wmv4_meta_info.txt | sort -k8 > ${FILE}/groups/info_hap${i}.txt;
#done
#grep -A 500 '$Haplotype$Assignment' ${FILE}/${FILE}.o > ${FILE}_summary.txt
##can improve
#for i in {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,AA,BB,CC,DD,EE,FF,GG,HH,II,JJ,KK,LL,MM,NN,OO,PP,QQ,RR,SS,TT,UU,VV,WW,XX,YY,ZZ};
#do
#echo "${i}Hap nSamples Group Country" >> ${FILE}_summary.txt;
#cat ${FILE}/groups/info_hap${i}.txt | awk '{print $14,$7}' | sort | uniq -c | sort -r >> ${FILE}_summary.txt;
#echo " nSamples Year" >> ${FILE}_summary.txt;
#cat ${FILE}/groups/info_hap${i}.txt | awk '{print $10}' | sort | uniq -c | sort -r >> ${FILE}_summary.txt;
#done
#done
