#cat ids.txt | while read line
#do for a in {1.0,0.9,0.8,0.7,0.6,0.5} 
#do for i in {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z} 
#do if grep -q ${line} HMGWAS_plink_whole_CT${a}_MIT0.9/groups/hap${i}.txt
#then echo ${line},CT${a},${i}
#fi
#done
#done
#done >> Id_haps.txt

##then just make a simple csv table with col 1=ID; col 2=PHEN

library(clustree)
library(tidyverse)

id_haps <- read.csv("Id_haps.txt")
wide_IH <- spread(id_haps, "CTCT", "HAP") %>% as_tibble()
nona_WIH <- wide_IH %>% replace(is.na(.), 'Z')
pheno <- read.csv("allphenono_.txt")
ph_NWIH <- left_join(nona_WIH,
                     pheno,
                     by="ID")
                     
phencol <- clustree(ph_NWIH, 
                    prefix = "CT", 
                    node_colour = "PHEN", 
                    node_colour_aggr = "mean", 
                    edge_width = 1, 
                    node_alpha = 0.9) + 
        scale_color_viridis_c(limits=c(max(top_frac(ph_NWIH,
                                              -0.1,
                                              PHEN)$PHEN),
                                 min(top_frac(ph_NWIH,
                                              0.1,
                                              PHEN)$PHEN)),
                                 oob = scales::squish, 
                                 begin = 0.25, 
                                 direction = -1, 
                                 name = "Pheno") +
        scale_edge_color_continuous(high = "black", 
                                    low = "grey80") +
        labs(edge_colour = "nIndividuals", 
             size = "nIndividuals") +
        guides(edge_alpha = FALSE)

ggsave('phencol_clustree.pdf',
       phencol, 
       device = 'pdf', 
       dpi = 300,
       height = 9,
       width = 16,
       units = "in")
