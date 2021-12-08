cat pop.txt | while read line
do mosdepth -b 1 -Q 10 ${line} ${line}.chr20.bam
zcat ${line}.regions.bed.gz | grep glyma.Wm82.gnm4.Gm20 > ${line}.chr20_regions.txt 
done
#probably better to do below by extracting a range with awk rather than this manual approach
cat pop.txt | while read line; 
do grep -A 5034 31724592 ${line}.chr20_regions.txt > ${line}.reg.tt; 
done
#
cat pop.txt | while read line; 
do grep -A 302 31728620 ${line}.reg.tt | awk '{sum+=$NF; n++} END { if (n>0) print sum / n; }' | xargs -I [] echo -e "${line}\t[]"; 
done > 304bp_avs.txt
sed -i 's/-/./g' 304bp_avs.txt
sed -i 's/_/./g' 304bp_avs.txt
# Below gives a quick summary of individuals with '0' runs of coverage, 
# However some individuals have ~0-1 coverage in small sections of the region likely due to misassembly which significantly truncates the '0' counts
# For this reason we manually inspected individuals who had significant 0counts that were lower than 310bp
cat pop.txt | while read line; do echo ${line} >> 0sum.txt
awk '{print $4}' ${line}.reg.tt | uniq -c | sort -k1 | tail -n 1 >> 0sum.txt
done
grep -B 1 ' 0' 0sum.txt | grep -v ' 0' | grep -v -- -- > just0s.txt
grep ' 0' 0sum.txt | awk '{print $1}' | paste just0s.txt - > 0counts.txt
