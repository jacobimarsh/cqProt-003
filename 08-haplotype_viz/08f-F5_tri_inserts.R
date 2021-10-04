library(tidyverse)
library(ggplot2)
library(forcats)
library(patchwork)
library(pals)

#input4 from TableS6: Hap     ID      Prot    Grp     Origin_country  GWAS-SNP        31727019_Allele How (How is '31727019_CopyNum')

alldat <- read.csv("input4.txt", sep = "\t")
col_list <- c('wt'=npg_col[8],
              'lr' = npg_col[3],
              'ocult' =npg_col[2],
              'mcult' =npg_col[4])
CNVnum_by_grp <- alldat %>% group_by(Grp) %>% filter(How != 'ref', How != 'mis') %>% count(How)

grpINS <- CNVnum_by_grp %>% 
  mutate(Grp=factor(Grp,
                           levels = rev(c('mcult','ocult','lr','wt')))) %>% 
  ggplot(aes(x = How, 
                     y = n, 
                     fill = Grp))  + 
  geom_bar(width =0.8, 
           stat = "identity",
           col = "black") + 
  labs(title = "Frequency of Insertions at 31727019", 
       y = "Individuals", 
       fill = "Domestication status", 
       x = 'Trinucleotide Repeat Count') +
   scale_fill_manual(values = col_list,
                     labels = rev(c("Modern Cultivar", "Old Cultivar","Landrace" , "Wild-Type"))) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.85,0.8)) 

grpINS

ggsave('grpINS2.pdf',
       grpINS, 
       device = 'pdf', 
       dpi = 1800,
       height = 90,
       width = 160,
       units = "mm")
