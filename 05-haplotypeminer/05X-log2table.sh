##used to easily assess runs at different parameters in a directory
#!/usr/bin/env bash
echo -e 'Input\tCT\tPairD\tSNP2GeneD\tMIT\tnBiallelic_SNPs\tKept_for_clustering\tKept_from_clustering\tKept_markers\tUnassigned_inds\tHap_size\tnHaplotype\tnInds_per_hap' >> 060421_excel_Log.txt
for FILE in analysis_output/060421/*
do
echo ${FILE} | sed 's/analysis_output\///g' >> tempfile.txt
grep 'Gene_name : ' ${FILE}/Log.txt | sed 's/.*Gene_name : //g' >> tempfile.txt
echo >> tempfile.txt
grep 'Input_file : ' ${FILE}/Log.txt | sed 's/.*\///g' >> tempfile.txt
grep 'Marker_cluster_threshold : ' ${FILE}/Log.txt | sed 's/.*Marker_cluster_threshold : //g' >> tempfile.txt
grep 'Maximum_flanking_pair_distance : ' ${FILE}/Log.txt | sed 's/.*Maximum_flanking_pair_distance : //g' >> tempfile.txt
grep 'Maximum_marker_to_gene_distance : ' ${FILE}/Log.txt | sed 's/.*Maximum_marker_to_gene_distance : //g' >> tempfile.txt
grep 'Marker_independence_threshold : ' ${FILE}/Log.txt | sed 's/.*Marker_independence_threshold : //g' >> tempfile.txt
grep 'Number of biallelic markers : ' ${FILE}/Log.txt | sed 's/.*Number of biallelic markers : //g' >> tempfile.txt
grep 'Number of markers passing MAC filter (final number kept for clustering): ' ${FILE}/Log.txt | sed 's/.*Number of markers passing MAC filter (final number kept for clustering): //g' >> tempfile.txt
grep 'Total markers kept following clustering: ' ${FILE}/Log.txt | sed 's/.*Total markers kept following clustering: //g' >> tempfile.txt
grep 'Total markers kept: ' ${FILE}/Log.txt | sed 's/.*Total markers kept: //g' >> tempfile.txt
grep 'Number of individuals not unambiguously assigned a haplotype: ' ${FILE}/Log.txt | sed 's/.*Number of individuals not unambiguously assigned a haplotype: //g' >> tempfile.txt
grep 'Haplotype size (distance between two farthest markers): ' ${FILE}/Log.txt | sed 's/.*Haplotype size (distance between two farthest markers): //g' >> tempfile.txt
grep 'Number of distinct haplotypes : ' ${FILE}/Log.txt | sed 's/.*Number of distinct haplotypes : //g' >> tempfile.txt
grep 'Number of individuals assigned to haplotype' ${FILE}/Log.txt | awk '$NF > 4 {print $NF}' | sort -rn >> tempfile.txt
echo >> 060421_excel_Log.txt
cat tempfile.txt | tr "\n " "\t " >> 060421_excel_Log.txt
rm tempfile.txt
done
