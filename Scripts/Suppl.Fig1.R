###################################### Suppl. Figure 1 ##################
######## plot Gene, sample, prevalence
############## Bar Plot CHIP Genes and Per Sample

library(ggplot2)
library(cowplot)
library(ggpubr)
theme_set(theme_cowplot())
library(readxl)
library(data.table)
library(dplyr)
load("~/Documents/Project/Methylation/CHS_CHIP_VAF_DNMT3A_TET2.rda")

vaf2020 <- fread("~/Documents/Project/Methylation/TOPMed CHIP calls 8-31-2020/TOPMed_variant_level_CHIPcalls_with_covariates_2020_08_31.tsv")

vaf2020CHS <- subset(vaf2020, vaf2020$Sample %in% CHS_CHIP_VAF_DNMT3A_TET2$NWD_ID); rm(vaf2020)

Dat_baseline_visit1 <- CHS_CHIP_VAF_DNMT3A_TET2[order(CHS_CHIP_VAF_DNMT3A_TET2$NWD_ID,CHS_CHIP_VAF_DNMT3A_TET2$agemeth, decreasing = FALSE), ]

Dat_baseline_v1 <- Dat_baseline_visit1[!duplicated(Dat_baseline_visit1$NWD_ID),]

CHIP_Gene_Table <- read_excel("~/Documents/Project/Methylation/DNAm/Replication_ARIC_April28/Aric_baseline_info.xlsx", 
                              sheet = "Gene_count")
CHIP_Gene_Table.ARIC <- subset(CHIP_Gene_Table, CHIP_Gene_Table$Cohort=="ARIC")

CHIP_Gene_Table.ARIC <- CHIP_Gene_Table.ARIC[,c(1,4,5)]
  # CHS
chip.count <- as.data.frame(table(vaf2020CHS$Sample))

dnmt3a_chip.chs <- subset(vaf2020CHS, vaf2020CHS$Sample %in% Dat_baseline_v1$NWD_ID[Dat_baseline_v1$hasDNMT3A==1])

dnmt3a_chip.chs <- subset(dnmt3a_chip.chs, !duplicated(dnmt3a_chip.chs$Sample))

tet2_chip.chs <- subset(vaf2020CHS, vaf2020CHS$Sample %in% Dat_baseline_v1$NWD_ID[Dat_baseline_v1$hasTET2==1])

tet2_chip.chs <- subset(tet2_chip.chs, !duplicated(tet2_chip.chs$Sample))

other_chip.chs <- subset(vaf2020CHS, vaf2020CHS$Sample %in% Dat_baseline_v1$NWD_ID[is.na(Dat_baseline_v1$haschip==1 & Dat_baseline_v1$hasDNMT3A==0 & Dat_baseline_v1$hasTET2==0)])

other_chip_noDnmtTet2.chs <- subset(other_chip.chs, !(other_chip.chs$Gene %in% c("DNMT3A", "TET2")) )

asxl1_chip.chs <- subset(other_chip.chs, other_chip.chs$Gene=="ASXL1")

dnmt3a_tet_asxl1.chs <- data.frame(rbind(dnmt3a_chip.chs, tet2_chip.chs, asxl1_chip.chs))

dnmt3a_tet_asxl1.chs <- dnmt3a_tet_asxl1.chs[!duplicated(dnmt3a_tet_asxl1.chs$Sample), c(6,8)]

dnmt3a_tet_asxl1.chs$Cohort <- "CHS"  

# 3 samples wtih both DNMT3A and TET2 CHIP
samids.to_include <- unique(other_chip.chs$Sample[!(other_chip.chs$Sample %in% other_chip_noDnmtTet2.chs$Sample)])
  # keeping the larger clone
other_chip_3samWtih_DNMT_TET2.chs <- subset(other_chip.chs, other_chip.chs$Sample %in% samids.to_include)
other_chip_3samWtih_DNMT_TET2.chs <- other_chip_3samWtih_DNMT_TET2.chs[order(other_chip_3samWtih_DNMT_TET2.chs$Sample, other_chip_3samWtih_DNMT_TET2.chs$VAF, decreasing = TRUE), ]
other_chip_3samWtih_DNMT_TET2.chs <- subset(other_chip_3samWtih_DNMT_TET2.chs, !duplicated(other_chip_3samWtih_DNMT_TET2.chs$Sample))

genes.CHS <- as.data.frame(rbind(dnmt3a_chip.chs[,c(1,6)], 
                           tet2_chip.chs[, c(1,6)], 
                           other_chip_noDnmtTet2.chs[,c(1,6)], 
                           other_chip_3samWtih_DNMT_TET2.chs[,c(1,6)]))
genes.CHS$Cohort <- "CHS"

CHIP_Gene_Table.CHS <- as.data.frame(table(genes.CHS$Gene))

names(CHIP_Gene_Table.CHS) <- c("Gene", "Freq")

CHIP_Gene_Table.CHS$Cohort <- "CHS"

CHIP_Gene_Table.CHS_ARIC <- as.data.frame(rbind(CHIP_Gene_Table.CHS, CHIP_Gene_Table.ARIC))   

CHIP_Gene_Table.CHS_ARIC$Cohort <- factor(CHIP_Gene_Table.CHS_ARIC$Cohort, levels = c("CHS", "ARIC"))

CHIP_Gene_Table.CHS_ARIC$Prop <- NA

CHIP_Gene_Table.CHS_ARIC$Prop[1:15] <- round(prop.table(CHIP_Gene_Table.CHS_ARIC$Freq[CHIP_Gene_Table.CHS_ARIC$Cohort=="CHS"])*100,1)

CHIP_Gene_Table.CHS_ARIC$Prop[16:nrow(CHIP_Gene_Table.CHS_ARIC)] <- round(prop.table(CHIP_Gene_Table.CHS_ARIC$Freq[CHIP_Gene_Table.CHS_ARIC$Cohort=="ARIC"])*100,1)

##

# Get prop
p1 <- ggplot(data=CHIP_Gene_Table.CHS_ARIC, aes(x=reorder(Gene, -Prop), 
                                       y=Prop, fill=Cohort)) + 
  xlab("") + ylab("Proportion of Individuals (%)") +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=Prop), vjust=-0.6, color="black",
            position = position_dodge(1.2), size=2.5) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,  hjust=1), 
        legend.position = "None", 
        plot.title = element_text(size = 20, face = "bold")) +
  scale_fill_manual(values=c("black", "gray45")) +  ggtitle("a)")

# VAF distribution
ARIC_gene_vaf_info_2022_01_26 <- readxl::read_excel("Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/review/ARIC_gene_vaf_info_2022-01-26.xlsx") 
ARIC_gene_vaf_info_2022_01_26$Cohort <- "ARIC"

vaf2020_chs_aric <- as.data.frame(rbind(dnmt3a_tet_asxl1.chs, ARIC_gene_vaf_info_2022_01_26))

vaf2020_chs_aric$Gene <- factor(vaf2020_chs_aric$Gene, levels = c("DNMT3A", "TET2", "ASXL1")) 
vaf2020_chs_aric$Cohort <- factor(vaf2020_chs_aric$Cohort, levels = c("CHS", "ARIC")) 


vaf2020_chs_aric  %>% group_by(Gene, Cohort) %>% summarise(median=median(VAF), mean=mean(VAF), Min=min(VAF), Max=max(VAF))

p2 <- ggplot(data=vaf2020_chs_aric, aes(x=Gene, y=VAF, fill=Cohort)) + 
  xlab("") + ylab("Variant Allele Fraction") + 
  geom_violin(trim = FALSE) + scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1,  hjust=1), 
        legend.title = element_blank(), 
        plot.title = element_text(size = 20, face = "bold")) + 
  scale_fill_manual(values=c("black", "gray45")) +  ggtitle("b)") +
  stat_summary(fun = "median", geom = "point",
               color = "white", 
               position = position_dodge(0.9))

##
# PerSampleCHIP_Table <- as.data.frame(table(as.data.frame(table(vaf_info$Sample))$Freq))
PerSampleCHIP_Table <- read_excel("~/Documents/Project/Methylation/DNAm/Replication_ARIC_April28/Aric_baseline_info.xlsx", 
                                  sheet = "Persample")
PerSampleCHIP_Table$Prop <- NA
PerSampleCHIP_Table$Prop[1:3] <- round(prop.table(PerSampleCHIP_Table$Freq[PerSampleCHIP_Table$Cohort=="ARIC"])*100,1) 
PerSampleCHIP_Table$Prop[4:6] <- round(prop.table(PerSampleCHIP_Table$Freq[PerSampleCHIP_Table$Cohort=="CHS"])*100,1) 
PerSampleCHIP_Table$Cohort <- factor(PerSampleCHIP_Table$Cohort, levels = c("CHS", "ARIC"))

# http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization
p3 <- ggplot(data=PerSampleCHIP_Table, aes(x=Var1, y=Prop, fill=Cohort)) + 
  xlab("Number of Variants") + ylab("Proportion of Individuals (%)") +
  geom_bar(stat="identity", width=0.8, 
           position="dodge",) +
  geom_text(aes(label=Prop), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=3.5) + 
  scale_fill_manual(values=c("black", "gray45")) + 
  theme(legend.position = "None", 
        plot.title = element_text(size = 20, face = "bold")) +  
  ggtitle("c)")

# Prevalence
Preval <- read_excel("~/Documents/Project/Methylation/DNAm/Replication_ARIC_April28/Aric_baseline_info.xlsx", 
                     sheet = "Preval")
	# Confidence interval
prop_CI <- function(x,n){
  my_prep <- x/n
  return(list( round( (my_prep - 1.96 * sqrt( (my_prep * (1-my_prep))/n ))*100,2), 
               round( (my_prep + 1.96 * sqrt( (my_prep * (1-my_prep))/n ) ) *100,2) ) )
}

Preval$Prev <- round(Preval$nCHIP/Preval$totN *100,2)
Preval$lowCI <- prop_CI(x = Preval$nCHIP, n = Preval$totN)[[1]]
Preval$uppCI <- prop_CI(x = Preval$nCHIP, n = Preval$totN)[[2]]
# Order
Preval$AgeBin <- factor(Preval$AgeBin, levels = c("41-50", "51-60", "61-70", "71-80", ">80"))
Preval$Cohort <- factor(Preval$Cohort, levels = c("CHS", "ARIC"))
#
p4 <- ggplot(data=Preval, aes(x=AgeBin, y=Prev, 
                              group=Cohort)) + 
  xlab("Age") + ylab("CHIP Prevalence (%)") + 
  geom_line(aes(linetype = Cohort, color= Cohort)) + 
  geom_text(aes(label=Prev), vjust=-1, color="black", 
            position = position_dodge(0.1), size=3.5) +
  geom_ribbon(aes(ymin = lowCI, ymax = uppCI, x = AgeBin), 
              alpha = 0.1)  +  ggtitle("d)") + 
  scale_color_manual(values=c("black", "gray45")) +
  theme(plot.title = element_text(size = 20, face = "bold"), 
        legend.title = element_blank())

# https://cran.r-project.org/web/packages/egg/vignettes/Ecosystem.html
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SupplementaryFigure1.png", 
    width=16, height=10, units= "in", res=300, pointsize = 4)

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2 )

dev.off()
######################### Supplementary Figure 1 END #################################

