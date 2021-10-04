cat pop.txt | while read line
do mosdepth -b 1 -Q 10 ${line} ${line}.chr20.bam
zcat ${line}.regions.bed.gz | grep glyma.Wm82.gnm4.Gm20 > ${line}.chr20_regions.txt 
done
#probably better to do below by extracting a range with awk rather than this manual approach
cat pop.txt | while read line; do grep -A 5034 31724592 ${line}.chr20_regions.txt > ${line}.reg.tt; done
cat pop.txt | while read line; do echo ${line} >> 0sum.txt
awk '{print $4}' ${line}.reg.tt | uniq -c | sort -k1 | tail -n 1 >> 0sum.txt
done
grep -B 1 ' 0' 0sum.txt | grep -v ' 0' | grep -v -- -- > just0s.txt
grep ' 0' 0sum.txt | awk '{print $1}' | paste just0s.txt - > 0counts.txt
