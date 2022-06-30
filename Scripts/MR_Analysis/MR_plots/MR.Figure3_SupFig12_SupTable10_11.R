
 
################################ Supplementary Table 10 #####################
########  Summary statistics of Mendelian randomization between      ######## 
########  cis-mQTL and CAD. GCTA software GSMR package was used for  ######## 
########  the MR analysis. mQTL=methylation quantitative trait loci  ######## 
########  CAD=coronary artery disease.                               ########  
#############################################################################
#### Load libraries
library(data.table)
######## Read GSMR results 
gsmr <- fread("~/Documents/Project/Methylation/DNAm/Final_results/Figures/MR_updated/CHIP_CAD_MR.csv") 
# Format P 
gsmr$P <- formatC(gsmr$p, format = "E", digits = 2)
# Convert Beta to "OR" 
gsmr$OR <- round(exp(gsmr$b),2)
# 95% Confidence Interval around "OR" 
gsmr$CI_95percent <- paste(round(exp(gsmr$b-1.96*gsmr$se),2),
                           round(exp(gsmr$b + 1.96*gsmr$se),2),
                           sep="-")
# FDR 
gsmr$FDR <- p.adjust(gsmr$p, "BH")

# Annotate CpGs with CHIP/DNMT3A/TET2 EWAS association
#### Load Replicated CHIP/DNMT3A/TET2 CpGs
load("~/Documents/Project/Methylation/DNAm/Final_results/Discovery_EWAS/Replicated.CpGs.meta_ewas.chip_vaf10_dnmt3a_tet2.rda")
gsmr$DNMT3A <- ifelse(gsmr$Exposure %in% dnmt3a.meta_chs_aric.replicated$CPG.Labels, 1,0)
gsmr$TET2 <- ifelse(gsmr$Exposure %in% tet2.meta_chs_aric.replicated$CPG.Labels, 1,0)
gsmr$CHIP <- ifelse(gsmr$Exposure %in% chip.meta_chs_aric.replicated$CPG.Labels, 1,0)

gsmr$CHIP_category[gsmr$DNMT3A==1 & gsmr$TET2==0 & gsmr$CHIP==0] <- "DNMT3A"
gsmr$CHIP_category[gsmr$DNMT3A==1 & gsmr$TET2==1 & gsmr$CHIP==0] <- "DNMT3A/TET2"
gsmr$CHIP_category[gsmr$DNMT3A==1 & gsmr$TET2==1 & gsmr$CHIP==1] <- "CHIP/DNMT3A/TET2"
gsmr$CHIP_category[gsmr$DNMT3A==1 & gsmr$TET2==0 & gsmr$CHIP==1] <- "CHIP/DNMT3A"
gsmr$CHIP_category[gsmr$DNMT3A==0 & gsmr$TET2==1 & gsmr$CHIP==0] <- "TET2"
gsmr$CHIP_category[gsmr$DNMT3A==0 & gsmr$TET2==1 & gsmr$CHIP==1] <- "CHIP/TET2"
gsmr$CHIP_category[gsmr$DNMT3A==0 & gsmr$TET2==0 & gsmr$CHIP==1] <- "CHIP"


# DNMT3A
dnmt.gsmr <- merge(gsmr, dnmt3a.meta_chs_aric.replicated, 
                   by.x = "Exposure", by.y="CPG.Labels")
  # change in DNAm is associted with increased risk for CAD
# dnmt.gsmr <- subset(dnmt.gsmr, dnmt.gsmr$b*dnmt.gsmr$meta_chs_aric.Effect>0)
dnmt.gsmr$EWAS_Cat <- "DNMT3A" 

# TET2
tet.gsmr <- merge(gsmr, tet2.meta_chs_aric.replicated, 
                  by.x = "Exposure", by.y="CPG.Labels")
# tet.gsmr <- subset(tet.gsmr, tet.gsmr$b*tet.gsmr$meta_chs_aric.Effect>0)
tet.gsmr$EWAS_Cat <- "TET2" 

# CHIP
chip.gsmr <- merge(gsmr, chip.meta_chs_aric.replicated, 
                   by.x = "Exposure", by.y="CPG.Labels")
# chip.gsmr <- subset(chip.gsmr, chip.gsmr$b*chip.gsmr$meta_chs_aric.Effect>0)
chip.gsmr$EWAS_Cat <- "CHIP" 
  # combine EWAS and MR results
comb.gsmr <- as.data.frame(rbind(dnmt.gsmr, tet.gsmr, chip.gsmr)) 
comb.gsmr <- comb.gsmr[order(comb.gsmr$Exposure, comb.gsmr$meta_chs_aric.P.value, decreasing = FALSE), ]
comb.gsmr <- subset(comb.gsmr, !duplicated(comb.gsmr$Exposure))

  # Suppl. Table 10
write.csv(comb.gsmr, "~/Documents/Project/Methylation/DNAm/Final_results/Figures/MR_updated/Suppl.Table10.CHIP_CAD_MR.2580CpGs.csv", row.names = F)

  # Table for Fig 3: fdr<5% and increaseing CAD risk in MR
# comb.gsmr.fdr5 <- subset(comb.gsmr, comb.gsmr$FRD<0.05 & 
#                                comb.gsmr$b*comb.gsmr$meta_chs_aric.Effect>0)
#  write.csv(comb.gsmr.fdr5, "~/Documents/Project/Methylation/DNAm/Final_results/updated_MR.FDR05_51CpGs.csv", row.names = F)
# save.image(file = "~/Documents/Project/Methylation/DNAm/Final_results/updated_MR.CHIP_DNMT3A_TET2.Dec31.rda")
# load("/Volumes/medpop_esp2/mesbah/projects/DNAm_CHIP/May2022_Desktop/Final_results/updated_MR.CHIP_DNMT3A_TET2.Dec31.rda")
#############################################################################

###################### Figure 3: Forest plot #########################
############## Mendelian randomization analysis of CHIP-associated 
############## CpGs and CAD risk 
######################################################################
  # Read MR results for 51 FDR significant CpGs with increased risk for CAD
gsmr.fdr05 <- data.table::fread("~/Documents/Project/Methylation/DNAm/Final_results/updated_MR.FDR05_51CpGs.csv")
gsmr.fdr05$Exposure_Gene <- paste(gsmr.fdr05$Exposure, gsmr.fdr05$Gene, sep = " ")
gsmr.fdr05$Group <- ifelse(gsmr.fdr05$DNMT3A==1, "DNMT3A",
                           ifelse(gsmr.fdr05$TET2==1, "TET2",
                                  ifelse(gsmr.fdr05$TET2==0 & 
                                           gsmr.fdr05$DNMT3A==0 , "CHIP")))
table(gsmr.fdr05$Group, exclude = NULL)
gsmr.fdr05$lower <- exp(gsmr.fdr05$b - 1.96*gsmr.fdr05$se)
gsmr.fdr05$upper <- exp(gsmr.fdr05$b + 1.96*gsmr.fdr05$se)
names(gsmr.fdr05)
  # Sort by effect beta
gsmr.fdr05 <- gsmr.fdr05[order(gsmr.fdr05$b, decreasing = T),]

gsmr.fdr05$FDR_f <- formatC(x = gsmr.fdr05$FDR, digits = 2,format = "E")

  # Forestplot
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figure3.CHIP_CAD_MR_Forestplot.png",
    width=12, height=12, units= "in", res=300, pointsize = 4)

par(mfrow=c(1,1), mar= c(5, 4.5, 4, 2))

library(TwoSampleMR)

forest_plot_1_to_many(mr_res = gsmr.fdr05, b="b", se="se",
                      exponentiate=T, ao_slc=F, 
                      TraitM="Exposure_Gene", 
                      by="Group",
                      trans="identity",
                      xlab="OR for CAD per SD increase in DNAm (95% CI)",    
                      subheading_size=11,
                      col1_title="Risk factor",
                      col1_width=1.6,
                      col_text_size=4, colour_scheme="black",
                      addcols=c("meta_chs_aric.Direction", 
                                "nsnp","OR","FDR_f"),
                      addcol_widths=c(0.4, 0.3, 0.3,0.6),
                      addcol_titles=c("","","","")
)

dev.off()
# # text edited in inkscape and PowerPoint
####################################################################################


####################### Supplementary Figure 12:######################## 
##### Scatterplots depicting SNP effects on 
##### exposures (CpGs) vs. SNP effects on outcome (CAD). 
#### GSMR plots
# https://yanglab.westlake.edu.cn/software/gcta/res/gsmr_plot.r
source("gsmr_plot.r")

# Plot Bonferroni significant MR
comb.gsmr.bonferroni <- subset(gsmr.fdr05, gsmr.fdr05$p<0.05/2580)
pdf("~/Documents/Project/Methylation/DNAm/Final_results/Figures/Suppl.Fig12.GSMR_plots.pdf")
for(i in 1:nrow(comb.gsmr.bonferroni)){
  chr <- comb.gsmr.bonferroni$Chr[i]
  cpg <- comb.gsmr.bonferroni$Exposure[i]
  d <- read_gsmr_data(paste0("~/Documents/Project/Methylation/DNAm/Final_results/MR_CAD.Jul28/gsmr_mQTL.chr",chr,".cad.results.eff_plot.gz"))
  plot_gsmr_effect(d, cpg, "CAD", colors()[75])
}
dev.off()

####################################################################################

### Suppl. Table 11
################################ Supplementary Table 11 #####################
######## Partially independent (LD r2<0.05) cis-mQTL and cis-eQTL    ######## 
######## summary statistics from Min et al. 2021 Nat Genet           ######## 
########  (http://mqtldb.godmc.org.uk/)                              ######## 
######## and Vosa et al 2021 Nat Genet (https://www.eqtlgen.org/),   ########  
######## respectively. These mQTL were used as instrumental variable ########  
######## in Mendelian randomization in GSMR analysis.                ######## 
######## mQTL=methylation quantitative trait loci;                   ######## 
######## eQTL=expression quantitative trait loci.                    ######## 
#############################################################################
################ 

#### Extract mQTL SNPs
source("gsmr_plot.r")
comb.gsmr.fdr5 <- fread("~/Documents/Project/Methylation/DNAm/Final_results/final_mr51.comb.gsmr.fdr5.csv")
dat <- comb.gsmr.fdr5[,c(1,7)]

for(chr in unique(dat$Chr)){
  print(chr)

  d <- read_gsmr_data(paste0("~/Documents/Project/Methylation/DNAm/Final_results/MR_CAD.Jul28/gsmr_mQTL.chr",chr,".cad.results.eff_plot.gz"))

  my_cpgs <- subset(dat, dat$Chr==chr)[,1]
  for (i in 1:length(my_cpgs)) {
    cpg <- my_cpgs[i]
    print(cpg)
    cat(cpg, gsmr_snp_effect(d, expo_str = cpg, outcome_str = "CAD")$snp,
        sep = "\n", 
        file = paste0("~/Documents/Project/Methylation/DNAm/Final_results/MR_CAD.Jul28/",cpg,".mQTL.snps.chr",chr,".txt") )
  }
}

###### Extract eQTL
eQTLs <- data.table::fread("/Volumes/mesbah/DNAm/eQTLGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz", header=T)
# file_list <- system("ls -l ~/Documents/Project/Methylation/DNAm/Final_results/MR_CAD.Jul28/*.mQTL.snps.chr*.txt | awk '{print $NF}' ", intern = T)
for(chr in unique(dat$Chr)){
  print(chr)
  d <- fread(paste0("~/Documents/Project/Methylation/DNAm/Final_results/Discovery_EWAS/mqtl.chr",chr,".tsv.gz"))
  my_cpgs <- subset(dat, dat$Chr==chr)[,1]
  d <- subset(d, d$cpg %in% my_cpgs$Exposure)
  for (i in 1:nrow(my_cpgs)) {
    cpg <- my_cpgs[i]
    
    snp_list <- fread(paste0("~/Documents/Project/Methylation/DNAm/Final_results/MR_CAD.Jul28/",cpg,".mQTL.snps.chr",chr,".txt"), header=T)
    
    my_mQTL <- subset(d, (d$cpg %in% cpg) & (d$SNP %in% snp_list[[1]]))
    
    my_mQTL_eQTL <- merge(my_mQTL, eQTLs, by="SNP", all.x = T)
    write.csv(my_mQTL_eQTL, paste0("~/Documents/Project/Methylation/DNAm/Final_results/MR_CAD.Jul28/",cpg,".my_mQTL_eQTL.chr",chr,".csv"),row.names = F)
    
  }
}

#############################################################################


