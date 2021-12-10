#conda create --name blast
conda activate blast
#conda install -c bioconda blast
module load bedtools

cat species.txt | while read line
do makeblastdb -dbtype nucl -in ${line}.fna
blastn -db ${line}.fna -query 20g085100.fa -outfmt 6 -word_size 11 | grep 'e-' | grep -v 'e-0' > ${line}_20g085100_blast.txt
sed -i "s/^/${line}\t/" ${line}_20g085100_blast.txt
done
##Results indicate the 304bp SV has a tonne of hits throughout many genomes, definitely an abundant TE

cat *20g085100_blast.txt | 
  awk 'FS="\t" {if ($8 > 4000 && $9 < 4400 && $NF > 149) print $0, OFS = "\t"}' > 304SV_150bs_hits.txt
  
cat *20g085100_blast.txt | 
  awk 'FS="\t" {if (!($8 > 4000 && $9 < 4400) && $NF > 149) print $0, OFS = "\t"}' > not_304SV_150bs_hits.txt  
