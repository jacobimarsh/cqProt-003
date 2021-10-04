## call invariant sites
## have a list.txt with all the ID's, may need to fix list2.txt in later step

tac list.txt | while read line; do gatk --java-options "-Xmx10g -XX:ParallelGCThreads=27" GenotypeGVCFs -R glyma.Wm82.gnm4.4PTR.genome_main.fa -V ${line}.g.vcf -O justGm20_${line}.vcf -all-sites -L glyma.Wm82.gnm4.Gm20 --allow-old-rms-mapping-quality-annotation-data; done
tac list2.txt | sed 's/.gz//' | while read line ; do bcftools annotate -x ^FORMAT/GT ../${line}.gz -O v -o rdy_${line}; bgzip rdy_${line}; tabix -p vcf rdy_${line}.gz; done
bcftools merge -l list3.txt -m all -o my -O v --threads 25
bgzip -c my.vcf
tabix -p vcf my.vcf.gz

vcftools --gzvcf ../my.vcf.gz \
--max-maf 0 \
--remove-indels \
--max-missing 0.9 \
--recode --stdout | bgzip -c > test_invariant.vcf.gz
vcftools --gzvcf ../my.vcf.gz \
--mac 1 \ 
--max-missing 0.9 \
--minQ 30 \
--recode --stdout | bgzip -c > test_variant.vcf.gz

tabix -p vcf test_invariant.vcf.gz
tabix -p vcf test_variant.vcf.gz
bcftools concat \
--allow-overlaps \
test_variant.vcf.gz test_invariant.vcf.gz \
-O z -o test_filtered.vcf.gz

#set up allHaps.txt with id's in first col and hap in second

conda activate pixy
pixy --stats pi fst dxy \
--vcf test_filtered.vcf.gz \
--populations allHaps.txt \
--window_size 10000 \
--n_cores 20 \
--chromosomes 'glyma.Wm82.gnm4.Gm20' \
--output_prefix allHaps

conda activate pixy
pixy --stats pi fst dxy \
--vcf test_filtered.vcf.gz \
--populations allHaps.txt \
--window_size 10000 \
--n_cores 20 \
--chromosomes 'glyma.Wm82.gnm4.Gm20' \
--output_prefix allHaps_reg \
--interval_start 31604127 \
--interval_end 31777346

conda activate r3.6
R
library(tidyverse)
library(tidyverse)
library(pheatmap)
library(corrplot)
library(ComplexHeatmap)
library(circlize)
library(lazyeval)
#pi
allHaps_pi <- read.csv("allHaps_pi.txt", sep = "\t")
group_by(allHaps_pi, pop) %>% summarise(tot = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))
allHaps_reg_pi <- read.csv("allHaps_reg_pi.txt", sep = "\t")
group_by(allHaps_reg_pi, pop) %>% summarise(tot = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T)
#dXY
filt_Haps_DXY <- read.csv("Haps_dxy.txt", sep = "\t")
allDXY_sum <- group_by(filt_Haps_DXY, pop1, pop2) %>% summarise(tot = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))
DXY_premat <- spread(allDXY_sum, pop1, tot) %>% add_row(pop2 = "A",.before = T)%>% add_column(G = 0) %>% remove_rownames %>% column_to_rownames(var="pop2")%>% 
  mutate_each( funs_( interp( ~replace(., is.na(.),0) ) ) )
DXY_mat <- DXY_premat %>% select(order(colnames(DXY_premat))) %>% as.matrix.data.frame()
full_DXY_mat <- (t(DXY_mat) + DXY_mat)*1000
Heatmap(full_DXY_mat, cluster_rows = F, cluster_columns = F, row_names_side = "left", column_names_side = "top", column_names_rot = 0, 
        cell_fun = function(j, i, x, y, width, height, fill) 
  {grid.text(sprintf("%.1f", full_DXY_mat[i, j]), x, y, gp = gpar(fontsize = 10))},
  col = colorRamp2(c(0.5, 3, 5.5), c("blue", "white", "red")),
  name = "dXYE+03")

filt_Haps_regDXY <- read.csv("Haps_reg_dxy.txt", sep = "\t")
regDXY_sum <- group_by(filt_Haps_regDXY, pop1, pop2) %>% summarise(tot = sum(count_diffs, na.rm = T)/sum(count_comparisons, na.rm = T))
regDXY_premat <- spread(regDXY_sum, pop1, tot) %>% add_row(pop2 = "A",.before = T)%>% add_column(G = 0) %>% remove_rownames %>% column_to_rownames(var="pop2")%>% 
  mutate_each( funs_( interp( ~replace(., is.na(.),0) ) ) )
regDXY_mat <- regDXY_premat %>% select(order(colnames(regDXY_premat))) %>% as.matrix.data.frame()
full_regDXY_mat <- (t(regDXY_mat) + regDXY_mat)*1000
Heatmap(full_regDXY_mat, cluster_rows = F, cluster_columns = F, row_names_side = "left", column_names_side = "top", column_names_rot = 0, 
        cell_fun = function(j, i, x, y, width, height, fill) 
  {grid.text(sprintf("%.1f", full_regDXY_mat[i, j]), x, y, gp = gpar(fontsize = 10))},
  col = colorRamp2(c(0.5, 3, 5.5), c("blue", "white", "red")),
  name = "dXY*E3")
