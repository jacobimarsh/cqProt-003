library(tidyverse)
library(data.table)
library(pals)
library(ggsci)
library(scales)

ts5 <- fread("tableS5.txt", sep = "\t") %>% as_tibble() %>% 
  mutate(MG=as.numeric(gsub('M0','',MG)))


p <- ggplot() + 
  geom_segment(data=filter(ts5,Type=='nrep'),aes(x=Pos,xend=Pos,y=MG-0.2,yend=MG+0.2),size=0.2)+
  geom_segment(data=filter(ts5,Type=='rep'),aes(x=Pos,xend=Pos,y=MG-0.2,yend=MG+0.2),col='red')+
  geom_point(data=filter(ts5,Type=='rep'),aes(Pos,MG), pch=23,fill='red',size=2)+
  geom_vline(xintercept = c(31604127 ,31777346),
             color = c("#ff9933"),
             linetype = "dashed")+
  geom_vline(xintercept = c(31724592 ,31729626),
             color = c("#66ccff"),
             linetype = "dashed")+
  scale_x_continuous(breaks= pretty_breaks(), labels = comma)+
  scale_y_reverse(breaks=1:6, labels=paste('M0',1:6,sep=''))+
  labs(x='Position', y='Marker group')+
  theme_minimal()+
  theme(axis.text.y = element_text(face='bold',color= c("red", "black", "black", "black", "black", "black")))+
  expand_limits(x=c(31450000,31900000))+
  NULL
  
ggsave('sup_plot.pdf',p,device = 'pdf',units = 'cm',height = 9,width = 16)  
