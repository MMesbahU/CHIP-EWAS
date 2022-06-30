
#### Samples for LD reference panel used in GSMR
#################### UKB 20k for LD Reference ####
## 20k random EUR UKB samples
ukb500k_eid <- data.table::fread("/Volumes/medpop_esp2/pradeep/UKBiobank/v3data/ukb7089_imp_chr3_v3_s487395.sample")
ukb200k <- read.csv("/Users/muddin/Documents/Project/UKB200k/CHIP_Calls_200k/April7_Shared/ukb200k_Age_Sex_Ethnicity.eid7089_31063.csv", header = T, stringsAsFactors = F)
ukb200k.EUR <- subset(ukb200k, ukb200k$ETHNIC=="1001")

ukb200k.EUR <- merge(ukb200k.EUR, ukb500k_eid, by.x = "eid_7089", by.y="ID_1")
set.seed(18072021)
ukb20k.EUR <- sample(x = ukb200k.EUR$eid_7089, size = 20000, replace = F)
write.table(ukb20k.EUR, "/Volumes/mesbah/ukb_chip/ukb_sample_random_20kEUR.tsv", row.names = F, col.names = F, quote = F)
###################################### 

###################################### Get cis-mQTL from goDMC ########
library(data.table)
library(httr)
library(dplyr)

# All CpGs
load("~/Documents/Project/Methylation/DNAm/Final_results/Discovery_EWAS/Replicated.CpGs.meta_ewas.rda")
all_fdr.cpgs <- unique(c(chip.meta_chs_aric.replicated$CPG.Labels, 
                         tet2.meta_chs_aric.replicated$CPG.Labels, 
                         dnmt3a.meta_chs_aric.replicated$CPG.Labels)) 
### Annotaation
annot <- fread("humanMethylation450k/ProductInfo/HumanMethylation450_15017482_v1-2.csv", 
               stringsAsFactors = F, header = T, 
               skip = 7, fill = TRUE, sep=",")

annot <- annot[, c(1,12,13)]
names(annot) <- c('CPG.Labels',"chr","Pos")
annot <- subset(annot, (annot$chr %in% c(1:22)) & !is.na(annot$Pos) )

all_fdr.cpgs.chr <- merge(all_fdr.cpgs, annot, by=1 )
names(all_fdr.cpgs.chr) <- c("CPG.Labels", "chr", "Pos")

########
res.mQTL.all <- vector("list", length(unique(all_fdr.cpgs.chr$chr)))

for (i in c(1:22)){
  res_list <- POST("http://api.godmc.org.uk/v0.1/query", body = list(
    cpgs = all_fdr.cpgs.chr$CPG.Labels[all_fdr.cpgs.chr$chr==i],
    cistrans = 'cis',
    clumped = '',
    pval = 5e-8), 
    encode = "json")
  res.mQTL.all[[i]] <- content(res_list) %>% lapply(., as_data_frame) %>% bind_rows; rm(res_list)
  
  cat("chr=",i,",mQTL=", nrow(res.mQTL.all[[i]]))
}

# save by chromosome
for(j in 1:22){
  d <- res.mQTL.all[[j]]
  d <- d[,c(10, 22, 1, 2, 12, 5, 24, 19, 23)]
  names(d) <- c("cpg", "SNP","A1", "A2", "freq", "b", "se", "p", "N") 
  fwrite(d, paste0("~/Documents/Project/Methylation/DNAm/Final_results/Discovery_EWAS/mqtl.chr",j,".tsv.gz"), row.names = F, quote = F, sep="\t")
}
######################## END 


