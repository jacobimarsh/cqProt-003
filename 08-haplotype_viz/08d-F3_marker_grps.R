PhenoSum <- read.csv("allpheno_resum3.txt") #Site,Type,Alleles,nInd,nIndPheno,AvPheno
AlleleFile <- read.csv("U_S_allele_fin3.txt") %>% filter(allele != 31628273 & allele != 31637882) #from log2input.sh with pruned sites

AlPhen <- merge(x = AlleleFile, 
                y = PhenoSum, 
                by.x = "allele", 
                by.y = "Site")

TagPercDiffs <- read.csv("percdifftags3.csv") %>% subset(select = -X) %>% filter(TAGGING >= 31604127 & TAGGING <= 31777346 & SNP != 31628273 & SNP != 31637882) #from tags_phenos.sh
AlleleCounts <- read.csv("ACAN_tagSNPs.txt") %>% separate(COUNTINFO, c("AC", "AF"), sep = ";" )%>% filter(SITE >= 31604127 & SITE <= 31777346) #from tags_phenos.sh

TagProp2 <- AlleleCounts %>% 
  mutate(AltAF = as.numeric(AC) / as.numeric(AF)) %>% 
  subset(select = -c(AC,AF))
TagAlleles2 <- right_join(TagPercDiffs, 
                          AlleleFile, 
                          by =c("SNP" = "allele"))
TagAlProp2 <- right_join(TagAlleles2, 
                         TagProp2, 
                         by = c("TAGGING" = "SITE")) %>% 
  mutate(percdiff = 100*percdiff) %>% filter(!is.na(SNP))

#Left plot (C)
C <- AlPhen %>% 
  mutate(Type=factor(Type,levels = c("REF","MISS","HET","HETMISS","ALT"))) %>%
  ggplot(aes(x = nInd, 
             y = as.character(allele),
             fill = Type,
             color = Type),
         ylim()) +
  geom_bar(aes(), 
           position = "stack", 
           stat = "identity",
           colour = "black",
           width = 0.8) +
  scale_x_reverse(breaks = scales::pretty_breaks(n = 5), 
                     expand = c(0,0)) +
  theme_void() + 
  theme(axis.text.y = element_blank(), 
        axis.title.x = element_text(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.key.size = unit(5, 
                               "mm"),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0.1), 
                           "cm"),
        plot.title = element_blank()) +
 scale_color_manual(values = c('black', 'black', 'black', 'black', 'black'), 
                    guide=F) +
 scale_fill_manual(values = rev(c('#440154FF', "#238A8DFF", "#73D055FF", "#FDE725FF", "#FFFFFF"))) +
  xlab("Allele count") +
  scale_y_discrete(position = "right", labels = c(paste0("M",as.character(20:10)), paste0("M0",as.character(7:1))))
  
C


# Right plot (D)
D <- ggplot() +
  geom_jitter(data = TagAlProp2,
              aes(x = abs(percdiff),
                  y = as.character(SNP),
                  fill = AltAF),
              alpha = 0.25,
              pch = 21,
              height = 0.25) +
  scale_fill_gradient('Minor allele frequency',
                      low = 'white',
                      high = '#440154FF') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "bold",
                                   size = 10, color = c("black", "black", "black", "black", "black", "red")),
        plot.margin = unit(c(0,0.1,0,0), 
                           "cm"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.position = "bottom",
        legend.key.size = unit(5, 
                               "mm"),
        plot.title = element_blank(),
        axis.text.x = element_text(face = "bold",
                                   size = 10),
        axis.title.x = element_text()) +
  xlab("Seed protein association [%]") +
  scale_y_discrete(position = "left", labels = c(paste0(" M0",as.character(6:1))))

D


CD <- C + D + plot_layout(design = "CD
                          CD")
CD

ggsave('CD3.pdf',
       CD, 
       device = 'pdf', 
       dpi = 1800,
       height = 97.875,
       width = 174,
       units = "mm")
