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

R
library(tidyverse)
library(pals)
library(ggsci)

avcov <- read.csv("304bp_avs.txt", sep = "\t", header = F) %>% as.tibble()
haps <- read.csv("haps.txt", sep = "\t", header = F) %>% as.tibble()
hapavcov <- merge(x = avcov, 
               y = haps, 
               by.x = "V1", 
               by.y = "V2") %>% rename("ID" = "V1", "bpcov" = "V2", "hap" = "V1.y", "grp" = "V3")

npg_col = pal_npg("nrc")(9)
col_list2 <- c('wt'=npg_col[8],
              'lr' = npg_col[3],
              'ocult' =npg_col[2],
              'mcult' =npg_col[4])

bpcovplot <- ggplot(data = hapavcov %>% 
  mutate(grp=factor(grp,
                           levels = c('wt','lr','ocult','mcult')))
               ) +
  geom_jitter(aes(x = hap,
                  y = bpcov,
                  fill = grp),
                  alpha=0.4, 
              pch=21, 
              width=0.2) +
  scale_fill_manual("Grp", values = col_list2, labels = rev(c("Modern Cultivar", "Old Cultivar","Landrace" , "Wild-Type"))) +
  xlab("Haplotype Group") +
  ylab("Mean Depth Across SV") +
  theme_minimal() +
  theme(plot.title = element_text(hjust= 0.5),
        legend.title = element_blank())

ggsave('bpcovplot.pdf',
       bpcovplot, 
       device = 'pdf', 
       dpi = 1800,
       height = 90,
       width = 160,
       units = "mm")
