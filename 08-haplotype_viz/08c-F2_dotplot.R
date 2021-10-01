library(tidyverse)
library(dplyr)
library(splitstackshape)
library(viridis)
library(patchwork)
library(ggplot2)
library(vioplot)
library(scales)
library(pals)
library(ggsci)
library(ggalt)

#Set up inputs
#US_haps_fin3.txt has J-O trimmed already

HapsFile <- read.csv("U_S_haps_fin3.txt") %>% mutate(altalleles=ifelse(hap=='B','',altalleles),
         hap=gsub('B','A',hap),
         hap=gsub('C','B',hap),
         hap=gsub('D','C',hap),
         hap=gsub('E','D',hap),
         hap=gsub('F','E',hap),
         hap=gsub('G','F',hap),
         altalleles=ifelse(hap=='H','',altalleles),
         hap=gsub('H','F',hap),
         hap=gsub('I','G',hap)) %>% 
  group_by(hap) %>% 
  summarise(hap=hap,nInd=sum(nInd),wt=sum(wt),lr=sum(lr),oc=sum(oc),mc=sum(mc),altalleles=paste(altalleles,collapse = ';')) %>% 
  distinct() %>% 
  mutate(altalleles = gsub('^\\;|\\;$', '', altalleles), altalleles = ifelse(altalleles=='',NA,altalleles))
HapsFile %>% separate_rows(altalleles,sep=';') %>% 
  mutate(value=1)
HapAlMatrix <- HapsFile %>% separate_rows(altalleles,sep=';') %>% 
     mutate(value=1) %>% spread(altalleles,value,fill = 0) %>% select(-`<NA>`)
intersect <- HapAlMatrix %>% 
  as_tibble() %>% 
  gather(allele,
         present,
         7:ncol(.)) %>% 
  mutate(present=as.factor(present)) %>%
  mutate(allele = str_remove(allele, "altalleles_"))
intersect_lines <- intersect %>% 
  filter(present == 1) %>% 
  group_by(hap) %>% 
  mutate(allele = as.numeric(allele)) %>% 
  summarise(max = max(allele), 
            min = min(allele)) %>% 
  mutate(min = as.character(min),
         max = as.character(max))
         
#Run graph

E <- ggplot() +
  geom_segment(data = intersect_lines, col = 'grey', 
               aes(x = hap, 
                   xend = hap,
                   y = min, 
                   yend = max),
               size = 1.5)+
  geom_point(data = intersect,
             aes(hap,
                 as.character(allele),
                 fill =as.factor(present),
                 size = 2),
             col ='black', 
             pch = 21) +
  scale_fill_manual(values = c('white','black', 'white'))+
  theme_minimal()+
  theme(legend.position = 'none',
        plot.margin = unit(c(0,0,0,0), 
                           "cm"),
        plot.title = element_blank(),
        axis.text.y = element_text(size=10, face='bold', color = c("black", "black", "black", "black", "black", "red")),
        axis.text.x = element_text(size=10, face = 'bold', color = "black")) +
  ylab("Marker group") +
  xlab("Haplotype combination") +
  scale_y_discrete(position = "left", labels = c(paste0("M0",as.character(6:1))))
# filter(intersect_lines, !(hap %in% c('M01','M08')))
E
