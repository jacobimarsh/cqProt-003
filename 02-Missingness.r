#vcftools --gzvcf Gm20.SNPs.id.biallic.vcf.gz --missing-site --out Gm20.SNPs.id.biallic.vcf
#reformat into two column file, c1 = position, c2 = missing proportion
library(tidyverse)
library(viridis)
missplot <- miss %>% ggplot(aes(x = POS, y = F_MISS, colour = F_MISS)) +
  geom_point() +
  ylab("Missing Frequency") +
  xlab("Position") + 
  scale_color_viridis_c('Missing Freq',
                        option = 3,
                        limits = c(0, 0.25), 
                        oob = scales::squish,
                        begin = 0,
                        end = 0.8) +
  geom_vline(xintercept = c(31604127 ,31777346),
             color = "red",
             linetype = "dashed") +
  theme_minimal()
  

ggsave('missplot.pdf',
       missplot, 
       device = 'pdf', 
       dpi = 1800,
       height = 90,
       width = 160,
       units = "mm")
