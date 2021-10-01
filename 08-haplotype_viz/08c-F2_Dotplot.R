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
