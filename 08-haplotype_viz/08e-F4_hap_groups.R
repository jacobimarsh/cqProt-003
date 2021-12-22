library(tidyverse)
library(pals)
library(ggsci)
library(gridExtra)
library(ggpp)
library(grid)
library(gtable)

tabdat <- read.csv("U_S_haps_fin3.txt", header = T) %>%  #from log2input.sh
  mutate(altalleles=ifelse(hap=='B','',altalleles),
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
  mutate(altalleles = gsub('^\\;|\\;$', '', altalleles), altalleles = ifelse(altalleles=='',NA,altalleles)) %>% select(-altalleles) %>% 
  rename(c(Wild = wt, Landrace = lr, 'Old Ctvr' = oc, 'Mod Ctvr' = mc, 'Hap' = hap, Total = nInd))

TableS8 <- read.csv("TableS8.txt", sep = "\t") %>% filter(as.double(Oil) > 0) %>%  #from supp tables
  mutate(Oil = as.numeric(Oil))

TableS8_shad <- read.csv("TableS8.txt", sep = "\t") %>% filter(as.double(Oil) > 0) %>% 
  mutate(Oil = as.numeric(Oil)) %>% 
  distinct(Oil, Prot) %>% 
  mutate(Hap = paste("A","B","C","D","E","F","G", sep = "_")) %>% 
  separate(Hap,c(paste0(rep('Hap',7),c("A","B","C","D","E","F","G"))),sep='_') %>% 
  gather(dummy,Hap,3:9) %>% 
  select(-dummy)

npg_col = pal_npg("nrc")(9)
col_list <- c('wt'=npg_col[8],
              'lr' = npg_col[3],
              'ocult' =npg_col[2],
              'mcult' =npg_col[4])
          
lint <- tableGrob(tabdat, rows = NULL, theme = ttheme_gtstripes(
                       colhead = list(bg_params = list(fill = "white"),
                         fg_params=list(
                           col=c("black", "black", "#DC0000FF","#00A087FF","#4DBBD5FF","#3C5488FF"))),
        plot.margin = unit(c(0,0,0,0), "cm")))

lint1 <- gtable_add_grob(lint,
                         grobs = segmentsGrob(
                           x0 = unit(0,"npc"),
                           y0 = unit(1,"npc"),
                           x1 = unit(1,"npc"),
                           y1 = unit(1,"npc"),
                           gp = gpar(lwd = 1)),
                         t = 1, l = 1, r = ncol(lint))
lint2 <- gtable_add_grob(lint1,
                         grobs = segmentsGrob(
                           x0 = unit(0,"npc"),
                           y0 = unit(1,"npc"),
                           x1 = unit(1,"npc"),
                           y1 = unit(1,"npc"),
                           gp = gpar(lwd = 1)),
                         t = 2, l = 1, r = ncol(lint1))
lint3 <- gtable_add_grob(lint2,
                         grobs = segmentsGrob(
                           x0 = unit(0,"npc"),
                           y0 = unit(0,"npc"),
                           x1 = unit(1,"npc"),
                           y1 = unit(0, "npc"),
                           gp = gpar(lwd = 1)),
                         t = 8, l = 1, r = ncol(lint2))


lint3$widths <- unit(rep(0.85/ncol(lint3), ncol(lint3)), "npc")

duntab <- grid.arrange(lint3)

duntab$widths <- unit(rep(1.1/ncol(duntab), ncol(duntab)), "npc")

Aplot <- ggplot() +
  geom_point(data=TableS8_shad,aes(x = Oil, y = Prot),col='light grey',size = 1, alpha = 1) +
  geom_point(data=filter(TableS8,Grp=='lr', Hap == 'A'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.6) +
  geom_point(data=filter(TableS8,Grp=='wt', Hap == 'A'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.7) +
  geom_point(data=filter(TableS8,Grp=='ocult', Hap == 'A'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  geom_point(data=filter(TableS8,Grp=='mcult', Hap == 'A'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ylab("%Seed Protein")+
    theme_minimal() +
  ggtitle("A") +
  theme(strip.text = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        strip.background = element_rect(fill = "light grey"),
        legend.position = 'none',
        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("%Seed Oil")  +
  scale_color_manual(values = col_list)

Bplot <- ggplot() +
  geom_point(data=TableS8_shad,aes(x = Oil, y = Prot),col='light grey',size = 1, alpha = 1) +
  geom_point(data=filter(TableS8,Grp=='lr', Hap == 'B'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.6) +
  geom_point(data=filter(TableS8,Grp=='wt', Hap == 'B'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.7) +
  geom_point(data=filter(TableS8,Grp=='ocult', Hap == 'B'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  geom_point(data=filter(TableS8,Grp=='mcult', Hap == 'B'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ylab("Seed Protein [%]")+
    theme_minimal() +
  ggtitle("B") +
  theme(strip.text = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "light grey"),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        legend.position = 'none',
        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("%Seed Oil")  +
  scale_color_manual(values = col_list)

Cplot <- ggplot() +
  geom_point(data=TableS8_shad,aes(x = Oil, y = Prot),col='light grey',size = 1, alpha = 1) +
  geom_point(data=filter(TableS8,Grp=='lr', Hap == 'C'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.6) +
  geom_point(data=filter(TableS8,Grp=='wt', Hap == 'C'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.7) +
  geom_point(data=filter(TableS8,Grp=='ocult', Hap == 'C'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  geom_point(data=filter(TableS8,Grp=='mcult', Hap == 'C'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ylab("%Seed Protein")+
    theme_minimal() +
  ggtitle("C") +
  theme(strip.text = element_text(size = 12, face = 'bold'),
        axis.text = element_blank(),
        axis.title = element_blank(),
        strip.background = element_rect(fill = "light grey"),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("%Seed Oil")  +
  scale_color_manual(values = col_list)

Dplot <- ggplot() +
  geom_point(data=TableS8_shad,aes(x = Oil, y = Prot),col='light grey',size = 1, alpha = 1) +
  geom_point(data=filter(TableS8,Grp=='lr', Hap == 'D'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.6) +
  geom_point(data=filter(TableS8,Grp=='wt', Hap == 'D'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.7) +
  geom_point(data=filter(TableS8,Grp=='ocult', Hap == 'D'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  geom_point(data=filter(TableS8,Grp=='mcult', Hap == 'D'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ylab("%Seed Protein")+
    theme_minimal() +
  ggtitle("D") +
  theme(strip.text = element_text(size = 12, face = 'bold'),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "light grey"),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("%Seed Oil")  +
  scale_color_manual(values = col_list)

Eplot <- ggplot() +
  geom_point(data=TableS8_shad,aes(x = Oil, y = Prot),col='light grey',size = 1, alpha = 1) +
  geom_point(data=filter(TableS8,Grp=='lr', Hap == 'E'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.6) +
  geom_point(data=filter(TableS8,Grp=='wt', Hap == 'E'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.7) +
  geom_point(data=filter(TableS8,Grp=='ocult', Hap == 'E'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  geom_point(data=filter(TableS8,Grp=='mcult', Hap == 'E'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ylab("%Seed Protein")+
  ggtitle("E") +
    theme_minimal() +
  theme(strip.text = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),,
        strip.background = element_rect(fill = "light grey"),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("%Seed Oil")  +
  scale_color_manual(values = col_list)

Fplot <- ggplot() +
  geom_point(data=TableS8_shad,aes(x = Oil, y = Prot),col='light grey',size = 1, alpha = 1) +
  geom_point(data=filter(TableS8,Grp=='lr', Hap == 'F'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.6) +
  geom_point(data=filter(TableS8,Grp=='wt', Hap == 'F'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.7) +
  geom_point(data=filter(TableS8,Grp=='ocult', Hap == 'F'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  geom_point(data=filter(TableS8,Grp=='mcult', Hap == 'F'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ylab("%Seed Protein")+
    theme_minimal() +
  ggtitle("F") +
  theme(strip.text = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = "light grey"),
        legend.position = 'none',
        axis.title = element_text(size = 12),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  xlab("Seed Oil [%]")  +
  scale_color_manual(values = col_list)

Gplot <- ggplot() +
  geom_point(data=TableS8_shad,aes(x = Oil, y = Prot),col='light grey',size = 1, alpha = 1) +
  geom_point(data=filter(TableS8,Grp=='lr', Hap == 'G'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.6) +
  geom_point(data=filter(TableS8,Grp=='wt', Hap == 'G'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.7) +
  geom_point(data=filter(TableS8,Grp=='ocult', Hap == 'G'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  geom_point(data=filter(TableS8,Grp=='mcult', Hap == 'G'),aes(x = Oil, y = Prot, col = Grp, size = 1),size = 1.5, alpha = 0.5) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  ylab("%Seed Protein")+
    theme_minimal()+
  ggtitle("G") +
  theme(strip.text = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 10),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_rect(fill = "light grey"),
        legend.position = 'none',
        plot.margin = unit(c(0,0,0,0), "cm"),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) +
  xlab("%Seed Oil") +
  scale_color_manual(values = col_list)

layout <- "AHH
BCD
EFG"

allplot <- Aplot + Bplot + Cplot + Dplot + Eplot + Fplot + Gplot + duntab + plot_layout(design = layout)

ggsave("ProtXOil8.pdf",
       allplot,
       device = "pdf",
       dpi = 1800,
       height = 150,
       width = 180,
       units = 'mm')
