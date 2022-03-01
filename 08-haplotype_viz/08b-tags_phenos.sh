#These aren't the exact files used for inputs, but the methodology was the same

plink --vcf protein_INs.vcf --show-tags U_S_sites.txt --tag-r2 0.9 --list-all
awk '{print $1,$NF}' plink.tags.list > tags.list
sed 's/ /,/' tags.list | sed 's/|/;/g' >> semifile.txt
R
dat<- read.csv("semifile.txt")
dat_long <- dat %>% as_tibble %>% separate(TAGS,as.character(2:300),sep=";") %>% gather(col,TAGGING,2:ncol(.)) %>% select(-col) %>% na.omit()
write_csv(dat_long, 'snp_data_long.csv')
exit()
sed 's/,20_/,SNP,/' snp_data_long.csv | sed 's/,INS_/,INS,/' | sed 's/20_//' >> Type_snp_data_long.csv

plink --vcf protein_DELs.vcf --show-tags U_S_sites.txt --tag-r2 0.9 --list-all
awk '{print $1,$NF}' plink.tags.list > DELtags.list
sed 's/ /,/' DELtags.list | sed 's/|/;/g' >> DELsemifile.txt
R
library(tidyverse)
DELdat<- read.csv("DELsemifile.txt")
DELdat_long <- DELdat %>% as_tibble %>% separate(TAGS,as.character(2:300),sep=";") %>% gather(col,TAGGING,2:ncol(.)) %>% select(-col) %>% na.omit()
write_csv(DELdat_long, 'DELsnp_data_long.csv')
quit()
n
grep DEL DELsnp_data_long.csv | sed 's/,20_/,SNP,/'  | sed 's/,DEL_/,DEL,/' | sed 's/20_//' >> Type_snp_data_long.csv

sort -n -k1 -t ',' Type_snp_data_long.csv > sorted_Type_snp_data_long.csv

awk 'FS = "," {print $3}' sorted_Type_snp_data_long.csv | while read line; do indel_pheno.sh ${line} protein_INs.vcf; done
grep DEL sorted_Type_snp_data_long.csv | awk 'FS = "," {print $3}' | while read line; do indel_pheno.sh ${line} protein_DELs.vcf; done

awk 'FS = "," {print $3}' sorted_Type_snp_data_long.csv | while read line; do uniq ${line}_phenosum.txt >> sort_${line}_phenosum.txt  ; done

#not used awk 'FS = "," {print $3}' sorted_Type_snp_data_long.csv | while read line; do echo $line >> a; head -n 2 ${line}_phenosum.txt | tail -n 1 | awk '{print $5}' | xargs -I [] python -c "print([] - $(head -n 3 ${line}_phenosum.txt | tail -n 1 | awk '{print $5}'))";  done | paste -d "," Type_snp_data_long.csv - >> pdiff_Type_snp_data_long.csv


paste ${line}_phenosum.txt -; done > allphen_adType.txt
#maybe this grep SNP pdiff_Type_snp_data_long.csv | awk 'FS="," {print $3}' | while read line; do echo -e 'Type\nREF\nALT\nHET\nHETMISS\nMISS' | paste ${line}_phenosum.txt -; done > allphen_adType.txt
R
library(tidyverse)
doto <- read.csv("allphen_adType.txt", sep = "\t", row.names = NULL, na.strings="") %>% rename(tormv = AvPheno, AvPheno = nIndPheno, nIndPheno = nInd, nInd = Alleles, Alleles = Site, Site = row.names) %>% subset(select = -(tormv)) %>% filter(nIndPheno != "") %>% filter(Site != "Site") %>%distinct() %>% mutate(Alleles = str_replace(Alleles, "./.", "MISS")) %>% mutate(Alleles=ifelse(nchar(Alleles) == 2,'REF',Alleles)) %>% mutate(Alleles=ifelse(nchar(Alleles) == 1,'ALT',Alleles)) %>% subset(select = -(Alleles))
tormv <- doto %>% subset(select = -c(nInd,nIndPheno)) %>% spread(Alleles,AvPheno) %>% mutate(Type=ifelse(is.na(Type) == F,"RMV",Type))
justpdiff <- tormv[!grepl("RMV",tormv$Type),] %>% subset(select = -Type) %>% mutate(pdiff = as.numeric(as.character(REF)) - as.numeric(as.character(ALT))) %>% subset(select = c(Site,pdiff))
rest <- read.csv("sorted_Type_snp_data_long.csv")
all <- right_join(rest,justpdiff,by = c("TAGGING" = "Site")) %>% as_tibble()
write.csv(all, "pdifftags.csv")

###PERC DIFF####

R
library(tidyverse)
dote <- read.csv("allphen_adType.txt", sep = "\t", row.names = NULL, na.strings="") %>% rename(tormv = AvPheno, AvPheno = nIndPheno, nIndPheno = nInd, nInd = Alleles, Alleles = Site, Site = row.names) %>% subset(select = -(tormv)) %>% filter(nIndPheno != "") %>% filter(Site != "Site") %>%distinct() %>% mutate(Alleles = str_replace(Alleles, "./.", "MISS")) %>% mutate(Alleles=ifelse(nchar(Alleles) == 2,'REF',Alleles)) %>% mutate(Alleles=ifelse(nchar(Alleles) == 1,'ALT',Alleles))
tormv <- dote %>% subset(select = -c(nInd,nIndPheno)) %>% spread(Alleles,AvPheno) %>% mutate(Type=ifelse(is.na(Type) == F,"RMV",Type))
justpdiff <- tormv[!grepl("RMV",tormv$Type),] %>% subset(select = -Type) %>% mutate(percdiff = (as.numeric(as.character(ALT)) - as.numeric(as.character(REF)))/as.numeric(as.character(REF))) %>% subset(select = c(Site,percdiff)) %>% as_tibble()
rest <- read.csv("sorted_Type_snp_data_long.csv") %>% as_tibble()
rest$TAGGING <- as.character(rest$TAGGING)
all <- right_join(rest,justpdiff,by = c("TAGGING" = "Site")) %>% as_tibble()
write.csv(all, "percdifftags.csv")

##Stopped including INS/DEL
##GETTING alt AF 
echo "SITE,AC,AF" >> ACAN_tagSNPs.txt
awk 'FS = "," {print $4}' pdifftags.csv | tail -n +2 | grep -f - sort_protein_Gm20.SNPs.id.biallic_maf_0.01_geno_0.1.vcf | sed -e 's/^20\t//' |sed -e 's/\t20_.*AC/,AC/' | sed -e 's/\tGT.*//' | sed 's/..=//g' | sed 's/;/,/' | awk 'FS="," {if ($1>=31604127 && $1<=31777346) print}' >> ACAN_tagSNPs.txt
R
library(tidyverse)
AlleleCounts <- read.csv("ACAN_tagSNPs.txt")
