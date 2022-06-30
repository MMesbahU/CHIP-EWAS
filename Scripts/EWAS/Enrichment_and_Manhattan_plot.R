################ Baseline characteristics: CHS vs ARIC #########3
  ## CHS
# load("/Volumes/medpop_esp2/mesbah/projects/DNAm_CHIP/May2022_Desktop/etcs/CHS_CHIP_VAF_DNMT3A_TET2.rda")
  # CHS vs ARIC baseine info
load("~/Documents/Project/Methylation/Revision/R001/table1_info_CHS_vs_ARIC.rda")

  # Keeping the 1st DNAm observation for repeat samples 
Dat_baseline_visit1 <- CHS_CHIP_VAF_DNMT3A_TET2[order(CHS_CHIP_VAF_DNMT3A_TET2$NWD_ID,CHS_CHIP_VAF_DNMT3A_TET2$agemeth, decreasing = FALSE), ]
Dat_baseline_v1 <- Dat_baseline_visit1[!duplicated(Dat_baseline_visit1$NWD_ID),]; rm(Dat_baseline_visit1)

Chars_baseline <- fread("/Volumes/medpop_esp2/mesbah/projects/DNAm_CHIP/May2022_Desktop/etcs/20210329_CHS.ARIC_baseline.tsv")
Dat_baseline_v1.pheno <- merge(Dat_baseline_v1, Chars_baseline, by.x="NWD_ID", by.y="NWDID")
# Table 1 CHS
# Race
table(Dat_baseline_v1.pheno$RACE)
# AA  EA 
# 280 302 
# Female
table(Dat_baseline_v1.pheno$Gender, Dat_baseline_v1.pheno$RACE, exclude = NULL)  
#   AA  EA
# 0 177 179
# 1 103 123
# Age
round(sd(Dat_baseline_v1.pheno$ageblood),1)
# 5.2
round(mean(Dat_baseline_v1.pheno$ageblood),1)
# 73.6
# EA Male==1
round(mean(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1]),1)
round(range(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1]),1)
# EA FeMale
round(mean(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0]),1)
round(range(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0]),1)

# AA Male=1
round(mean(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1]),1)
round(range(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1]),1)
# AA FeMale=0
round(mean(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0]),1)
round(range(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0]),1)

# Smoking
table(Dat_baseline_v1.pheno$Ever_smoker_baseline, Dat_baseline_v1.pheno$RACE, exclude = NULL)
#   AA  EA
# 0    126 135
# 1    154 166
# <NA>   0   1
round(166/301 * 100)
round(154/280 * 100)
round(prop.table(table(Dat_baseline_v1.pheno$Ever_smoker_baseline[ Dat_baseline_v1.pheno$RACE=="AA"]))*100)
round(prop.table(table(Dat_baseline_v1.pheno$Ever_smoker_baseline[ Dat_baseline_v1.pheno$RACE=="EA"]))*100)

## Smoking by sex and race
table(Dat_baseline_v1.pheno$Ever_smoker_baseline[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$Ever_smoker_baseline[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], exclude = NULL)

table(Dat_baseline_v1.pheno$Ever_smoker_baseline[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$Ever_smoker_baseline[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], exclude = NULL)


# BMI
round(summary(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1]),1)
round(sd(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], na.rm = T),1)
round(summary(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0]),1)
round(sd(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], na.rm = T),1)

round(summary(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1]),1)
round(sd(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], na.rm = T),1)
round(summary(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0]),1)
round(sd(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], na.rm = T),1)

# prev CAD all
Dat_baseline_v1.pheno$prevCAD <- ifelse(Dat_baseline_v1.pheno$ageblood>=Dat_baseline_v1.pheno$t_CAD_all_age & Dat_baseline_v1.pheno$t_CAD_all==1, 1,0)

table(Dat_baseline_v1.pheno$prevCAD [Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$prevCAD [Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], exclude = NULL)

table(Dat_baseline_v1.pheno$prevCAD [Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$prevCAD [Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], exclude = NULL)


table(Dat_baseline_v1.pheno[(Dat_baseline_v1.pheno$ageblood>=Dat_baseline_v1.pheno$t_CAD_all_age), c("NWD_ID","ageblood", "t_CAD_all_age", "RACE")]$RACE, exclude = NULL)
# AA EA 
# 3 11 
table(Dat_baseline_v1.pheno$ageblood>=Dat_baseline_v1.pheno$t_CAD_all_age & Dat_baseline_v1.pheno$t_CAD_all==1, Dat_baseline_v1.pheno$RACE, exclude = NULL)
#   AA  EA
# FALSE 277 291
# TRUE    3  11
round(3/280 * 100,1)
round(11/302 * 100,1)

# T2D
table(Dat_baseline_v1.pheno$T2D[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$T2D[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], exclude = NULL)
table(Dat_baseline_v1.pheno$T2D[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$T2D[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], exclude = NULL)

# table(Dat_baseline_v1.pheno$T2D, Dat_baseline_v1.pheno$RACE, exclude = NULL)
#   AA  EA
# 0 224 254
# 1  56  48
round(56/280 * 100,1)
round(48/302 * 100,1)

# CHIP
table(Dat_baseline_v1.pheno$haschip[Dat_baseline_v1.pheno$Gender==1], Dat_baseline_v1.pheno$RACE[Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$haschip[Dat_baseline_v1.pheno$Gender==0], Dat_baseline_v1.pheno$RACE[Dat_baseline_v1.pheno$Gender==0], exclude = NULL)

# table(Dat_baseline_v1.pheno$haschip, Dat_baseline_v1.pheno$RACE, exclude = NULL)
#   AA  EA
# 0 239 257
# 1  41  45
round(41/(41+239) *100,1)
round(45/(45+257) *100,1)
# VAF10
table(Dat_baseline_v1.pheno$hasVAF10[Dat_baseline_v1.pheno$Gender==1], Dat_baseline_v1.pheno$RACE[Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$hasVAF10[Dat_baseline_v1.pheno$Gender==0], Dat_baseline_v1.pheno$RACE[Dat_baseline_v1.pheno$Gender==0], exclude = NULL)

# table(Dat_baseline_v1.pheno$hasVAF10, Dat_baseline_v1.pheno$RACE, exclude = NULL)
#       AA  EA
# 0    239 257
# 1     36  40
# <NA>   5   5
round(36/(36+5+239) *100,1)
round(40/(40+5+257) *100,1)
# DNMT3A
# table(Dat_baseline_v1.pheno$hasDNMT3A, Dat_baseline_v1.pheno$RACE, exclude = NULL)
#       AA  EA
# 0    239 257
# 1     21  14
# <NA>  20  31
round(21/(21+20+239) *100,1)
round(12/(14+31+257) *100,1)
table(Dat_baseline_v1.pheno$hasDNMT3A[Dat_baseline_v1.pheno$Gender==1], Dat_baseline_v1.pheno$RACE[Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$hasDNMT3A[Dat_baseline_v1.pheno$Gender==0], Dat_baseline_v1.pheno$RACE[Dat_baseline_v1.pheno$Gender==0], exclude = NULL)

# TET2
# table(Dat_baseline_v1.pheno$hasTET2, Dat_baseline_v1.pheno$RACE, exclude = NULL)
#       AA  EA
# 0    239 257
# 1      6  12
# <NA>  35  33
round(6/(6+35+239) *100,1)
round(12/(12+33+257) *100,1)
table(Dat_baseline_v1.pheno$hasTET2[Dat_baseline_v1.pheno$Gender==1], Dat_baseline_v1.pheno$RACE[Dat_baseline_v1.pheno$Gender==1], exclude = NULL)
table(Dat_baseline_v1.pheno$hasTET2[Dat_baseline_v1.pheno$Gender==0], Dat_baseline_v1.pheno$RACE[Dat_baseline_v1.pheno$Gender==0], exclude = NULL)

######################


  ## ARIC 
# aric_baseline <- fread("~/Documents/Project/Methylation/Revision/R001/ARIC_dat_table1.csv")
table(aric_baseline$prvchd, exclude = NULL)

### Compare CHS vs ARIC 
# Age: Two sided Welch Two Sample t-test
  # Male
t.test(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$age[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$age[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
  # Female
t.test(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$age[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$ageblood[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$age[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# Ever Smoked: Two sided Welch Two Sample t-test
  # Male
t.test(Dat_baseline_v1.pheno$Ever_smoker_baseline[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$eversmk[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$Ever_smoker_baseline[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$eversmk[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
  # Female
t.test(Dat_baseline_v1.pheno$Ever_smoker_baseline[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$eversmk[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$Ever_smoker_baseline[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$eversmk[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# BMI: Two sided Welch Two Sample t-test
  # Male
t.test(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$bmi[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$bmi[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
  # Female
t.test(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$bmi[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$Bmi[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$bmi[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# prevCAD: Two sided Welch Two Sample t-test
  # Male
t.test(Dat_baseline_v1.pheno$prevCAD[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$prvchd[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$prevCAD[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$prvchd[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
  # Female
t.test(Dat_baseline_v1.pheno$prevCAD[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$prvchd[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$prevCAD[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$prvchd[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# prevT2D: Two sided Welch Two Sample t-test
  # Male
t.test(Dat_baseline_v1.pheno$T2D[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$diabts[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$T2D[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$diabts[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
  # Female
t.test(Dat_baseline_v1.pheno$T2D[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$diabts[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$T2D[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$diabts[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# Any CHIP: Two sided Welch Two Sample t-test
  # Male
t.test(Dat_baseline_v1.pheno$haschip[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$chip[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$haschip[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$chip[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
  # Female
t.test(Dat_baseline_v1.pheno$haschip[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$chip[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$haschip[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$chip[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# Expanded CHIP: Two sided Welch Two Sample t-test
# Male
t.test(Dat_baseline_v1.pheno$hasVAF10[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$afgt10[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$hasVAF10[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$afgt10[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
# Female
t.test(Dat_baseline_v1.pheno$hasVAF10[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$afgt10[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$hasVAF10[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$afgt10[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# DNMT3A: Two sided Welch Two Sample t-test
  # Male
t.test(Dat_baseline_v1.pheno$hasDNMT3A[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$dnmt3a[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$hasDNMT3A[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$dnmt3a[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
  # Female
t.test(Dat_baseline_v1.pheno$hasDNMT3A[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$dnmt3a[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$hasDNMT3A[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$dnmt3a[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# TET2: Two sided Welch Two Sample t-test
  # Male
t.test(Dat_baseline_v1.pheno$hasTET2[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$tet2[aric_baseline$race=="B" & aric_baseline$gender=="M"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$hasTET2[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==1], 
       aric_baseline$tet2[aric_baseline$race=="W" & aric_baseline$gender=="M"], alternative = "t")$p.value
# Female
t.test(Dat_baseline_v1.pheno$hasTET2[Dat_baseline_v1.pheno$RACE=="AA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$tet2[aric_baseline$race=="B" & aric_baseline$gender=="F"], alternative = "t")$p.value
t.test(Dat_baseline_v1.pheno$hasTET2[Dat_baseline_v1.pheno$RACE=="EA" & Dat_baseline_v1.pheno$Gender==0], 
       aric_baseline$tet2[aric_baseline$race=="W" & aric_baseline$gender=="F"], alternative = "t")$p.value

# save(Dat_baseline_v1.pheno, aric_baseline, file = "~/Documents/Project/Methylation/Revision/R001/table1_info_CHS_vs_ARIC.rda")

############################################

################ Replication in ARIC ####################################
## Load CHS meta-EWAS
# load("/Volumes/mesbah/projects/DNAm_CHIP/May2022_Desktop/Final_results/Discovery_EWAS/chs.meta_ewas.chip_vaf10_dnmt3a_tet2.rda")

load("/Volumes/medpop_esp2/mesbah/projects/DNAm_CHIP/May2022_Desktop/Final_results/Discovery_EWAS/chs.meta_ewas.chip_vaf10_dnmt3a_tet2.rda")

#########
# TET2 ##
#########
# Suppl. Table 3 (old)
suppl_tabl3 <- fread("~/Documents/Project/Methylation/Revision/R001/SupplT3.TET2.csv")

  ## Meta-ARIC and meta-CHS
load("/Volumes/medpop_esp2/mesbah/projects/DNAm_CHIP/May2022_Desktop/Final_results/Discovery_EWAS/meta_tet2.chs_aric_rep_combined.rda")
tet2.meta_chs_aric_rep$meta4_FDR <- p.adjust(tet2.meta_chs_aric_rep$meta_chs_aric.P.value, "BH")

   ## replicated in meta-ARIC
tet2.meta_chs_aric_rep$repl_metaARIC <- ifelse( (tet2.meta_chs_aric_rep$aric_meta.FDR<0.05 & (tet2.meta_chs_aric_rep$chs_meta.Effect * tet2.meta_chs_aric_rep$aric_meta.Effect) >0 ) ,1,0 )
## replicated in ARIC-AA
tet2.meta_chs_aric_rep$repl_ARIC_AA <- ifelse( tet2.meta_chs_aric_rep$CPG.Labels %in% meta.tet2_CHS.aricAA$CPG.Labels[meta.tet2_CHS.aricAA$repl_ARIC_AA_fdr==1] ,1,0 )
table(tet2.meta_chs_aric_rep$repl_ARIC_AA)
# 1423
## replicated in ARIC-EA
tet2.meta_chs_aric_rep$repl_ARIC_EA <- ifelse( tet2.meta_chs_aric_rep$CPG.Labels %in% meta.tet2_CHS.aricEA$CPG.Labels[meta.tet2_CHS.aricEA$repl_ARIC_EA_fdr==1] ,1,0 )
table(tet2.meta_chs_aric_rep$repl_ARIC_EA)
# 50
## merge AA, EA
tet2.meta_chs_aric_rep <- merge(tet2.meta_chs_aric_rep, tet2.result.aa, by.x="CPG.Labels", by.y = "aric_AA.results.CPG.Labels", all.x = T)
tet2.meta_chs_aric_rep <- merge(tet2.meta_chs_aric_rep, tet2.result.ea, by.x="CPG.Labels", by.y = "aric.EA.results.CPG.Labels", all.x = T)
table(tet2.meta_chs_aric_rep$repl_metaARIC)

  # annotation
tet2.meta_chs_aric_rep <- merge(tet2.meta_chs_aric_rep, annot[,c(1,12,13,17,22:33)], by.x="CPG.Labels", by.y="IlmnID")
tet2.meta_chs_aric_rep$Strand[tet2.meta_chs_aric_rep$Strand=="F"] = "+"
tet2.meta_chs_aric_rep$Strand[tet2.meta_chs_aric_rep$Strand=="R"] = "-"

write.csv(tet2.meta_chs_aric_rep, "~/Documents/Project/Methylation/Revision/R001/Jun22_2022.tet2.meta_chs_aric_rep.csv", row.names = F)

#########################
### ARIC: AA only
load("/Volumes/GoogleDrive-110230452486495871691/My Drive/CHIP and DNA Methylation/DNAm/Replication_ARIC_April28/aric.AA.chip.vaf.tet2.dnmt3a.results.04232021.rda")

names(tet2.result.aa) <- paste0("aric_AA.",names(tet2.result.aa))
meta.tet2_CHS.aricAA <- merge(meta.tet2, tet2.result.aa, by.x="CPG.Labels", by.y = "aric_AA.results.CPG.Labels")
meta.tet2_CHS.aricAA <- subset(meta.tet2_CHS.aricAA, meta.tet2_CHS.aricAA$chs_meta.FDR<0.05)
meta.tet2_CHS.aricAA$aric_AA.results.fdr <- p.adjust(p = meta.tet2_CHS.aricAA$aric_AA.results.P.value, 
                                                  method = "BH")

meta.tet2_CHS.aricAA$repl_ARIC_AA_fdr <- ifelse(meta.tet2_CHS.aricAA$chs_meta.Effect * 
                                                  meta.tet2_CHS.aricAA$aric_AA.coefficients.effect.size>0 & 
                                                  meta.tet2_CHS.aricAA$aric_AA.results.fdr<0.05,
                                                1,0) 
# fdr<0.05
(table(meta.tet2_CHS.aricAA$repl_ARIC_AA_fdr, 
                       exclude = NULL))
# 1423/11798 == 12.06%
meta.tet2_CHS.aricAA$repl_metaARIC <- ifelse(meta.tet2_CHS.aricAA$CPG.Labels %in% suppl_tabl3$CpG, 1, 0)
table(meta.tet2_CHS.aricAA$repl_metaARIC, meta.tet2_CHS.aricAA$repl_ARIC_AA_fdr, exclude = NULL)
#     0     1
# 0 10204   115
# 1   171  1308
meta.tet2_CHS.aricAA_nMetaARIC <- subset(meta.tet2_CHS.aricAA, meta.tet2_CHS.aricAA$repl_ARIC_AA_fdr==1 | meta.tet2_CHS.aricAA$repl_metaARIC==1) 
meta.tet2_CHS.aricAA_nMetaARIC <- merge(meta.tet2_CHS.aricAA_nMetaARIC, annot, by.x = "CPG.Labels", by.y="IlmnID")
####

## ARIC EA only
load("/Volumes/GoogleDrive-110230452486495871691/My Drive/CHIP and DNA Methylation/DNAm/Replication_ARIC_April28/aric.EA.chip.vaf.tet2.dnmt3a.results.04232021.rda")
names(tet2.result.ea) <- paste0("aric.EA.",names(tet2.result.ea))
meta.tet2_CHS.aricEA <- merge(meta.tet2, tet2.result.ea, by.x="CPG.Labels", by.y = "aric.EA.results.CPG.Labels")
meta.tet2_CHS.aricEA <- subset(meta.tet2_CHS.aricEA, 
                               meta.tet2_CHS.aricEA$chs_meta.FDR<0.05)
meta.tet2_CHS.aricEA$aric.EA.results.fdr <- p.adjust(p = meta.tet2_CHS.aricEA$aric.EA.results.P.value, method = "BH")
meta.tet2_CHS.aricEA$repl_ARIC_EA_fdr <- ifelse(meta.tet2_CHS.aricEA$chs_meta.Effect * meta.tet2_CHS.aricEA$aric.EA.coefficients.effect.size>0 & 
                                                        meta.tet2_CHS.aricEA$aric.EA.results.fdr<0.05, 1,0)
table(meta.tet2_CHS.aricEA$repl_ARIC_EA_fdr)
# 50/11785 == 0.42%

## Overlap with meta-ARIC repl
  # AA vs meta-ARIC
table(meta.tet2_CHS.aricAA$CPG.Labels[meta.tet2_CHS.aricAA$repl_ARIC_AA_fdr==1] %in% tet2.meta_chs_aric.replicated$CPG.Labels)
# FALSE  TRUE 
# 115  1308
  # EA vs meta-ARIC
table(meta.tet2_CHS.aricEA$CPG.Labels[meta.tet2_CHS.aricEA$repl_ARIC_EA_fdr==1] %in% tet2.meta_chs_aric.replicated$CPG.Labels)
# FALSE  TRUE 
# 34    16
  # AA vs EA
table(meta.tet2_CHS.aricEA$CPG.Labels[meta.tet2_CHS.aricEA$repl_ARIC_EA_fdr==1] %in% meta.tet2_CHS.aricAA$CPG.Labels[meta.tet2_CHS.aricAA$repl_ARIC_AA_fdr==1])
# FALSE  TRUE 
# 41     9
#############################

  # AA vs Tulstrup
table(meta.tet2_CHS.aricAA$CPG.Labels[meta.tet2_CHS.aricAA$repl_ARIC_AA_fdr==1] %in% 
        replicated_meta.tet2.DAN_ewas.fdr05$CPG.Labels)
# FALSE  TRUE 
# 228  1195 
  # EA vs Tulstrup
table(meta.tet2_CHS.aricEA$CPG.Labels[meta.tet2_CHS.aricEA$repl_ARIC_EA_fdr==1] %in% 
        replicated_meta.tet2.DAN_ewas.fdr05$CPG.Labels)
# FALSE  TRUE 
# 24    26
  # meta-ARIC vs Tulstrup
table(tet2.meta_chs_aric.replicated$CPG.Labels %in% 
        replicated_meta.tet2.DAN_ewas.fdr05$CPG.Labels)
# FALSE  TRUE 
# 221  1258 
########################################################################

################ TET2 Replication in Tulstrup et al. 2021 ####################################

# Load meta-EWAS CHS
load("/Volumes/mesbah/projects/DNAm_CHIP/May2022_Desktop/Final_results/Discovery_EWAS/chs.meta_ewas.chip_vaf10_dnmt3a_tet2.rda")

# Load replicated CpGs in ARIC 
load("/Volumes/mesbah/projects/DNAm_CHIP/May2022_Desktop/Final_results/Discovery_EWAS/Replicated.CpGs.meta_ewas.chip_vaf10_dnmt3a_tet2.rda")


## TET2 EWAS from Tulstrup et al. 2021
library(readxl)
Tulstrup.tet2.ewas <- read_excel("/Volumes/GoogleDrive-110230452486495871691/My Drive/Manuscript/Methylation_CHIP/Resoures_published/TET2_DNAm/TET2 mutation status.41467_2021_26093_MOESM5_ESM.xlsx")
table(meta.tet2$CPG.Labels %in% Tulstrup.tet2.ewas$probe) 
# FALSE   TRUE 
# 55768 422893 
meta.tet2.Tulstrup_ewas <- merge(meta.tet2, Tulstrup.tet2.ewas, by.x="CPG.Labels", 
                                 by.y="probe")

meta.tet2.Tulstrup_ewas.fdr05 <- subset(meta.tet2.Tulstrup_ewas, 
                                   meta.tet2.Tulstrup_ewas$chs_meta.FDR<0.05)

meta.tet2.Tulstrup_ewas.fdr05$FDR_tulstrup <- p.adjust(p = meta.tet2.Tulstrup_ewas.fdr05$pval, 
                                                  method = "BH")
  # FDR<0.05 and concordant effect direction
table(meta.tet2.Tulstrup_ewas.fdr05$FDR_tulstrup<0.05 & 
        (meta.tet2.Tulstrup_ewas.fdr05$chs_meta.Effect * 
           meta.tet2.Tulstrup_ewas.fdr05$estimate>0), 
      exclude = NULL)
# FALSE  TRUE 
# 4067  6943 
round(prop.table(table(meta.tet2.Tulstrup_ewas.fdr05$FDR_tulstrup<0.05 & 
                         (meta.tet2.Tulstrup_ewas.fdr05$chs_meta.Effect * 
                            meta.tet2.Tulstrup_ewas.fdr05$estimate>0), 
                       exclude = NULL))*100,1)
# FALSE  TRUE 
# 36.9  63.1 
meta.tet2.Tulstrup_ewas.fdr05$replicated_FDR05 <- ifelse(meta.tet2.Tulstrup_ewas.fdr05$FDR_tulstrup<0.05 & 
                                                           (meta.tet2.Tulstrup_ewas.fdr05$chs_meta.Effect*meta.tet2.Tulstrup_ewas.fdr05$estimate>0),1,0)

## Replicated in Danish study
replicated_meta.tet2.Tulstrup_ewas <- subset(meta.tet2.Tulstrup_ewas.fdr05, 
                                             meta.tet2.Tulstrup_ewas.fdr05$replicated_FDR05==1)
# Add Illumina 450k annotation
annot <- fread("/Volumes/GoogleDrive-110230452486495871691/My Drive/CHIP and DNA Methylation/humanMethylation450k/ProductInfo/HumanMethylation450_15017482_v1-2.csv",
               stringsAsFactors = F, header = T,
               skip = 7, fill = TRUE, sep=",")

replicated_meta.tet2.Tulstrup_ewas <- merge(replicated_meta.tet2.Tulstrup_ewas, 
                                            annot, by.x="CPG.Labels", 
                                            by.y="IlmnID")

write.csv(replicated_meta.tet2.Tulstrup_ewas, "Documents/Project/Methylation/Revision/replicated_meta.tet2.Tulstrup_ewas.csv", row.names = F)

## plot Avg. DNAm 
cor(meta.tet2.Tulstrup_ewas$chs_meta.Effect, 
    meta.tet2.Tulstrup_ewas$estimate)
## all matched 422,893 cpgs 0.2688976
cor(meta.tet2.Tulstrup_ewas$chs_meta.Freq1, 
    meta.tet2.Tulstrup_ewas$mean_beta_controls)
# 0.9981009
cor(meta.tet2.Tulstrup_ewas.fdr05$chs_meta.Effect, 
    meta.tet2.Tulstrup_ewas.fdr05$estimate)
# FDR significant CHS cpGs 0.4981827
cor(meta.tet2.Tulstrup_ewas.fdr05$chs_meta.Freq1, 
    meta.tet2.Tulstrup_ewas.fdr05$mean_beta_controls)
#  0.9927734
cor(replicated_meta.tet2.Tulstrup_ewas$chs_meta.Freq1, 
    replicated_meta.tet2.Tulstrup_ewas$mean_beta_controls)
#  0.9910702
cor(replicated_meta.tet2.Tulstrup_ewas$chs_meta.Effect, 
    replicated_meta.tet2.Tulstrup_ewas$estimate)
# replicated CpGs in Dans 0.7321157
plot(replicated_meta.tet2.Tulstrup_ewas$chs_meta.Freq1 ~ 
       replicated_meta.tet2.Tulstrup_ewas$mean_beta_controls)
plot(replicated_meta.tet2.Tulstrup_ewas$chs_meta.Effect ~ 
       replicated_meta.tet2.Tulstrup_ewas$estimate)
########################################################################


############################ GO Enrichment ##############################
  # ARIC replicated CpGs
load(file = "/Volumes/mesbah/projects/DNAm_CHIP/May2022_Desktop/Final_results/Discovery_EWAS/go_kegg_enrichment.July23_2021.rda")
library(missMethyl) 
  # CHIP
goCHIP.ARIC_rep <- gometh(sig.cpg = chip.meta_chs_aric.replicated$CPG.Labels,
                         all.cpg = meta.chip$CPG.Labels,
                         array.type = "450K", collection = "GO",
                         plot.bias = F, fract.counts = T, prior.prob = T)

table(goCHIP.ARIC_rep$FDR<0.05)
topGSA(goCHIP.ARIC_rep, number = 10)
table(goCHIP.ARIC_rep$ONTOLOGY, exclude = NULL)
    # DNMT3A
goDNMT3A.ARIC_rep <- gometh(sig.cpg = dnmt3a.meta_chs_aric.replicated$CPG.Labels,
                          all.cpg = meta.dnmt3a$CPG.Labels,
                          array.type = "450K", collection = "GO",
                          plot.bias = F, fract.counts = T, prior.prob = T)
goDNMT3A.ARIC_rep$GO_Identifier <- row.names(goDNMT3A.ARIC_rep) 
goDNMT3A.ARIC_rep_neg <- gometh(sig.cpg = dnmt3a.meta_chs_aric.replicated$CPG.Labels[dnmt3a.meta_chs_aric.replicated$chs_meta.Effect<0],
                            all.cpg = meta.dnmt3a$CPG.Labels,
                            array.type = "450K", collection = "GO",
                            plot.bias = F, fract.counts = T, prior.prob = T)

table(goDNMT3A.ARIC_rep$FDR<0.05)
table(goDNMT3A.ARIC_rep$P.DE<2.22e-6)

topGSA(goDNMT3A.ARIC_rep, number = 10)
table(goDNMT3A.ARIC_rep$ONTOLOGY, exclude = NULL)
write.csv(goDNMT3A.ARIC_rep, "Documents/Project/Methylation/Revision/goDNMT3A.ARIC_rep.csv", row.names = F)

# table(goDNMT3A.ARIC_rep_neg$TERM[goDNMT3A.ARIC_rep_neg$FDR<=0.05] %in% goTET2.DAN_rep$TERM[goTET2.DAN_rep$FDR<0.05])


  # TET2
goTET2.ARIC_rep <- gometh(sig.cpg = tet2.meta_chs_aric.replicated$CPG.Labels,
                          all.cpg = meta.tet2$CPG.Labels,
                          array.type = "450K", collection = "GO",
                          plot.bias = F, fract.counts = T, prior.prob = T)
goTET2.ARIC_rep$GO_Identifier <- row.names(goTET2.ARIC_rep) 

table(goTET2.ARIC_rep$TERM[goTET2.ARIC_rep$FDR<0.05] %in% goTET2.DAN_rep$TERM[goTET2.DAN_rep$FDR<0.05])

goTET2.ARIC_rep_pos <- gometh(sig.cpg = tet2.meta_chs_aric.replicated$CPG.Labels[tet2.meta_chs_aric.replicated$chs_meta.Effect>0],
                          all.cpg = meta.tet2$CPG.Labels,
                          array.type = "450K", collection = "GO",
                          plot.bias = F, fract.counts = T, prior.prob = T)

table(goTET2.ARIC_rep_pos$TERM[goTET2.ARIC_rep_pos$FDR<0.05] %in% goTET2.DAN_rep$TERM[goTET2.DAN_rep$FDR<0.05])

table(goTET2.ARIC_rep$FDR<0.05)
topGSA(goTET2.ARIC_rep, number = 10)
table(goTET2.ARIC_rep$ONTOLOGY, exclude = NULL)

# write.csv(goTET2.ARIC_rep, "Documents/Project/Methylation/Revision/goTET2.ARIC_rep.csv", row.names = F)
  # Exp
# goExpCHIP.ARIC_rep <- gometh(sig.cpg = vaf10.meta_chs_aric.replicated$CPG.Labels,
#                           all.cpg = meta.vaf10$CPG.Labels,
#                           array.type = "450K", collection = "GO",
#                           plot.bias = F, fract.counts = T, prior.prob = T)
# 
# table(goExpCHIP.ARIC_rep$FDR<0.05)
# topGSA(goExpCHIP.ARIC_rep, number = 10)
# table(goExpCHIP.ARIC_rep$ONTOLOGY, exclude = NULL)
#########################################################################


############################## Fig 1a EWAS Manhattan ######################
## Load Discovery EWAS 
load("~/Documents/Project/Methylation/DNAm/Final_results/Discovery_EWAS/final.CHS_meta_EWAS.annot.rda")

# # Manhattan plot
source("~/Documents/Project/Methylation/ewas_manhattan.R")
my_manhattan.reflect(x = m.chip, cpgname = annot$IlmnID,
                     chr = annot$CHR, pos = annot$MAPINFO,
                     main.title = "(a) Manhattan plot for association between DNAm and Any CHIP",
                     file.type = "eps",
                     save.plot = "~/Documents/Project/Manuscript/Methylation_CHIP/Pics/Figure2a.chip_metaEWAS_chs")
# 
my_manhattan.reflect(x = m.vaf, cpgname = annot$IlmnID,
                     chr = annot$CHR, pos = annot$MAPINFO,
                     main.title = "(b) Manhattan plot for association between DNAm and Expanded CHIP",
                     file.type = "eps",
                     save.plot = "~/Documents/Project/Manuscript/Methylation_CHIP/Pics/Figure2b.Expandedchip_metaEWAS_chs")

# 
my_manhattan.reflect(x = m.dnmt3a, cpgname = annot$IlmnID,
                     chr = annot$CHR, pos = annot$MAPINFO,
                     main.title = "(c) Manhattan plot for association between DNAm and DNMT3A CHIP",
                     file.type = "eps",
                     save.plot = "~/Documents/Project/Manuscript/Methylation_CHIP/Pics/Figure2c.DNMT3Achip_metaEWAS_chs")

my_manhattan.reflect(x = m.tet2, cpgname = annot$IlmnID,
                     chr = annot$CHR, pos = annot$MAPINFO,
                     main.title = "(d) Manhattan plot for association between DNAm and TET2 CHIP",
                     file.type = "eps",
                     save.plot = "~/Documents/Project/Manuscript/Methylation_CHIP/Pics/Figure2d.TET2chip_metaEWAS_chs")
# 
########
### Inflation
round(median(qchisq(1 - m.chip$P.value, 1))/qchisq(0.5, 1),2)
# 1.11
round(median(qchisq(1 - m.vaf$P.value, 1))/qchisq(0.5, 1),2)
# 1.4
round(median(qchisq(1 - m.dnmt3a$P.value, 1))/qchisq(0.5, 1),2)
# 0.92
round(median(qchisq(1 - m.tet2$P.value, 1))/qchisq(0.5, 1), 2)
# 1.45


############################## Fig 1b-d Volcano plots######################
# Volcano plot
# png("~/Documents/Project/Methylation/DNAm/Final_results/Figures/S1ad.metaEWAS.volcano_plot_chs.highlight_replicated.png", 
#     width=6, height=6, units= "in", res=300, 
#     pointsize = 5)
# par(mfrow=c(2,2), mar= c(5, 4.5, 4, 2))
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/Fig2b.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
# http://rfunction.com/archives/1302
par(mfrow=c(1,1), mar= c(5, 5, 4, 2))
plot(m.chip$Effect,
     -log10(m.chip$P.value), xlab = expression("Effect Size" ~ (beta) ), 
     ylab=expression(-log10 ~ italic("(P)")), 
     main="", font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
     cex.axis=1.5, font.axis=2, font.main=2)
abline(h=-log10(max(m.chip$P.value[m.chip$FDR<0.05])), lty="dashed")
points(m.chip$Effect[m.chip$CPG.Labels%in%chip.meta_chs_aric.replicated$CPG.Labels], 
       -log10(m.chip$P.value)[m.chip$CPG.Labels%in%chip.meta_chs_aric.replicated$CPG.Labels], 
       col="red")
legend("topright",c("CHIP EWAS"), text.col="red", 
       cex=1.8, text.font = 2, bty="n", bg="n")
dev.off()

# "gray45"
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/Supplementary_Figure_4b.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5); par(mfrow=c(1,1), mar= c(5, 5, 4, 2))
plot(m.vaf$Effect,
     -log10(m.vaf$P.value), xlab = expression("Effect Size" ~ (beta) ), 
     ylab=expression(-log10 ~ italic("(P)") ), 
     main="", font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
     cex.axis=1.5, font.axis=2, font.main=2)
abline(h=-log10(max(m.vaf$P.value[m.vaf$FDR<0.05])), lty="dashed")
points(m.vaf$Effect[m.vaf$CPG.Labels%in%vaf10.meta_chs_aric.replicated$CPG.Labels], 
       -log10(m.vaf$P.value)[m.vaf$CPG.Labels%in%vaf10.meta_chs_aric.replicated$CPG.Labels], 
       col="gray45")
# legend("topright",c("Expanded CHIP EWAS"), text.col="gray45", 
#        cex=1.8, text.font = 2, bty="n", bg="n")
dev.off()
# DNMT3A
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/Fig2c.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1,1), mar= c(5, 5, 4, 2))
plot(m.dnmt3a$Effect,
     -log10(m.dnmt3a$P.value), xlab = expression("Effect Size" ~ (beta) ), 
     ylab=expression(-log10 ~ italic("(P)")), 
     main="", font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
     cex.axis=1.5, font.axis=2, font.main=2)
abline(h=-log10(max(m.dnmt3a$P.value[m.dnmt3a$FDR<0.05])), lty="dashed")
points(m.dnmt3a$Effect[m.dnmt3a$CPG.Labels%in%dnmt3a.meta_chs_aric.replicated$CPG.Labels], 
       -log10(m.dnmt3a$P.value)[m.dnmt3a$CPG.Labels%in%dnmt3a.meta_chs_aric.replicated$CPG.Labels], 
       col="cornflowerblue")
legend("topright",legend = expression(italic("DNMT3A") ~ EWAS), text.col="cornflowerblue", 
       cex=2, text.font = 2, bty="n", bg="n")
dev.off()

# TET2
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/Fig2d.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1,1), mar= c(5, 5, 4, 2))
plot(m.tet2$Effect,
     -log10(m.tet2$P.value), xlab = expression("Effect Size" ~ (beta) ), 
     ylab=expression(-log10 ~ italic("(P)")), 
     main="", font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
     cex.axis=1.5, font.axis=2, font.main=2) 
abline(h=-log10(max(m.tet2$P.value[m.tet2$FDR<0.05])), lty="dashed")
points(m.tet2$Effect[m.tet2$CPG.Labels%in%tet2.meta_chs_aric.replicated$CPG.Labels], 
       -log10(m.tet2$P.value)[m.tet2$CPG.Labels%in%tet2.meta_chs_aric.replicated$CPG.Labels], 
       col="seagreen")
legend("topleft",legend = expression(italic("TET2") ~ EWAS), 
       text.col="seagreen", 
       cex=2, text.font = 2, bty="n", bg="n")

# par(mfrow=c(1,1), mar= c(5, 4.5, 4, 2))
dev.off()

###########################################################################

##################### QQ plot
# QQ plot
library(qqman)

png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SuppFig3a.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1, 1), mar= c(5, 5, 4, 2))
qq(m.chip$P.value, main="" , font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
   cex.axis=1.5, font.axis=2, font.main=2)
dev.off()


png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SuppFig4c.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1, 1), mar= c(5, 5, 4, 2))
qq(m.vaf$P.value, main="" , font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
   cex.axis=1.5, font.axis=2, font.main=2)
dev.off()

png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SuppFig3b.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1, 1), mar= c(5, 5, 4, 2))
qq(m.dnmt3a$P.value, main="" , font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
   cex.axis=1.5, font.axis=2, font.main=2)
dev.off()

png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SuppFig3c.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1, 1), mar= c(5, 5, 4, 2))
qq(m.tet2$P.value, main="" , font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
   cex.axis=1.5, font.axis=2, font.main=2)
dev.off()

# Zoomed QQ
# Supp Figs 3d-f
# 0-1 range
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SuppFig3d_chip.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1, 1), mar= c(5, 5, 4, 2))
qq(m.chip$P.value, main="" , font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
   cex.axis=1.5, font.axis=2, font.main=2, xlim=c(0, 1), ylim=c(0,1))
# m=0.30103
# m <- -log10(qchisq(0.5, 1))
m = -log10(.5)
lines(c(m,m), c(-0.1,m), lty=2, col="blue")
lines(c(-0.1,m), c(m,m), lty=2, col="blue")
dev.off()


# Dnmt3a 3e
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SuppFig3e_dnmt3a.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1, 1), mar= c(5, 5, 4, 2))
qq(m.dnmt3a$P.value, main="" , font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
   cex.axis=1.5, font.axis=2, font.main=2, xlim=c(0, 1), ylim=c(0,1))
m=-log10(.5)
# m <- -log10(qchisq(0.5, 1))
lines(c(m,m), c(-0.1,m), lty=2, col="blue")
lines(c(-0.1,m), c(m,m), lty=2, col="blue")
dev.off()


# Tet2 Supp Fig3f
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SuppFig3f_tet2.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1, 1), mar= c(5, 5, 4, 2))
qq(m.tet2$P.value, main="" , font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
   cex.axis=1.5, font.axis=2, font.main=2, xlim=c(0, 1), ylim=c(0,1))
m = -log10(.5)
# m <- -log10(qchisq(0.5, 1))
lines(c(m,m), c(-0.1,m), lty=2, col="blue")
lines(c(-0.1,m), c(m,m), lty=2, col="blue")
dev.off()


# Expanded CHIP;  Supp Figs 4d
png("~/Documents/Project/Methylation/DNAm/Manuscript_for_review/Updated/Figures/SuppFig4d_expchip.png", 
    width=3, height=3, units= "in", res=300, 
    pointsize = 5)
par(mfrow=c(1, 1), mar= c(5, 5, 4, 2))
qq(m.vaf$P.value, main="" , font.lab=2, font=2, cex.lab=1.8, cex.main=2, 
   cex.axis=1.5, font.axis=2, font.main=2, xlim=c(0, 1), ylim=c(0,1))
m=-log10(.5)
# m <- -log10(qchisq(0.5, 1))
lines(c(m,m), c(-0.1,m), lty=2, col="blue")
lines(c(-0.1,m), c(m,m), lty=2, col="blue")
dev.off()

###########################################################################


############################## Fig 1e Venn diagram ######################
library(VennDiagram)
replicated.venn.plot <- venn.diagram(
  x = list(
    "Any CHIP" = chip.meta_chs_aric.replicated$CPG.Labels,
    TET2 = tet2.meta_chs_aric.replicated$CPG.Labels,
    "Expanded CHIP" = vaf10.meta_chs_aric.replicated$CPG.Labels,
    DNMT3A = dnmt3a.meta_chs_aric.replicated$CPG.Labels
  ),
  filename = "~/Documents/Project/Methylation/DNAm/Final_results/Figures/Figure4.Overlap_Replicated_CpGs.tiff",
  col = "black",
  fill = c("red", "seagreen", "gray45", "cornflowerblue"),
  alpha = 0.50,
  cat.col = c("red", "seagreen", "gray45", "cornflowerblue"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05
)
#########################################################################

############################### UpSet Plot ##############################
library(UpSetR)
load(file = "/Volumes/mesbah/projects/DNAm_CHIP/May2022_Desktop/Final_results/Discovery_EWAS/go_kegg_enrichment.July23_2021.rda")

## Format data for upset plot
all_replicated_cpg <- unique(c(chip.meta_chs_aric.replicated$CPG.Labels, 
                               vaf10.meta_chs_aric.replicated$CPG.Labels, 
                               dnmt3a.meta_chs_aric.replicated$CPG.Labels,
                               tet2.meta_chs_aric.replicated$CPG.Labels))

dat_cpg <- as.data.frame(matrix(NA, 
                                nrow = length(all_replicated_cpg), 
                                ncol = 5))

names(dat_cpg) <- c("CpG", "Any CHIP", "Expanded CHIP", "DNMT3A", "TET2")

dat_cpg$CpG <- all_replicated_cpg
dat_cpg$`Any CHIP` <- ifelse(dat_cpg$CpG %in% chip.meta_chs_aric.replicated$CPG.Labels, 1,0)
dat_cpg$`Expanded CHIP` <- ifelse(dat_cpg$CpG %in% vaf10.meta_chs_aric.replicated$CPG.Labels, 1,0)
dat_cpg$DNMT3A <- ifelse(dat_cpg$CpG %in% dnmt3a.meta_chs_aric.replicated$CPG.Labels, 1,0)
dat_cpg$TET2 <- ifelse(dat_cpg$CpG %in% tet2.meta_chs_aric.replicated$CPG.Labels, 1,0)

table(dat_cpg$`Any CHIP`, dat_cpg$`Expanded CHIP`,
      dat_cpg$DNMT3A, dat_cpg$TET2)

pdf("~/Documents/Project/Methylation/Revision/upset_fig1e.pdf",
    width=7, height=5, bg="transparent")

upset(dat_cpg, sets = c("Any CHIP", "Expanded CHIP", "DNMT3A", "TET2"),
      order.by ="freq", decreasing = F, 
      mb.ratio = c(0.7, 0.3),
      sets.bar.color=c("red", "gray45", 
                       "cornflowerblue", "seagreen"))

dev.off()
#########################################################################


###########################################################################





############################### CpG Island Enrichment #####################
###### CPG Island
setwd("~/Documents/Project/Methylation/DNAm/Final_results/")
load(file = "~/Documents/Project/Methylation/DNAm/Final_results/Discovery_EWAS/go_kegg_enrichment.July23_2021.rda")
rm(list = ls()[3:42])
load("/Volumes/mesbah/projects/DNAm_CHIP/DNAm/Final_results/Discovery_EWAS/final.CHS_meta_EWAS.annot.rda")

chip.relation2cpgIsland <- cbind.data.frame(
  category=c("Open Sea", "Island", "Shore", "Shelf"),
  CpGs_in_mylist=c(table(m.chip$Relation_to_UCSC_CpG_Island[m.chip$CPG.Labels%in% chip.meta_chs_aric.replicated$CPG.Labels]=="")[["TRUE"]],table(m.chip$Relation_to_UCSC_CpG_Island[m.chip$CPG.Labels%in% chip.meta_chs_aric.replicated$CPG.Labels]=="Island")[["TRUE"]], table(grepl(pattern = "Shore",x = m.chip$Relation_to_UCSC_CpG_Island[m.chip$CPG.Labels%in% chip.meta_chs_aric.replicated$CPG.Labels]))[["TRUE"]], table(grepl(pattern = "Shelf",x = m.chip$Relation_to_UCSC_CpG_Island[m.chip$CPG.Labels%in% chip.meta_chs_aric.replicated$CPG.Labels]))[["TRUE"]]),
  CpGs_in_path = c(table(m.chip$Relation_to_UCSC_CpG_Island=="")[["TRUE"]], table(m.chip$Relation_to_UCSC_CpG_Island=="Island")[["TRUE"]], table(grepl(pattern = "Shore",x = m.chip$Relation_to_UCSC_CpG_Island))[["TRUE"]], table(grepl(pattern = "Shelf",x = m.chip$Relation_to_UCSC_CpG_Island))[["TRUE"]]))
chip.relation2cpgIsland$OR <- NA
chip.relation2cpgIsland$CI_95percent <- NA
chip.relation2cpgIsland$pval <- NA
for(i in 1:nrow(chip.relation2cpgIsland)) {
  chip.relation2cpgIsland$pval[i] <- formatC(fisher.test(matrix(c(chip.relation2cpgIsland[i,2], nrow(chip.meta_chs_aric.replicated) - chip.relation2cpgIsland[i,2], chip.relation2cpgIsland[i,3] - chip.relation2cpgIsland[i,2], nrow(m.chip) - chip.relation2cpgIsland[i,3] - nrow(chip.meta_chs_aric.replicated)), nrow = 2) )[[1]],digits = 1, format = "E")
  chip.relation2cpgIsland$CI_95percent[i] <- paste(round(fisher.test(matrix(c(chip.relation2cpgIsland[i,2], nrow(chip.meta_chs_aric.replicated) - chip.relation2cpgIsland[i,2], chip.relation2cpgIsland[i,3] - chip.relation2cpgIsland[i,2], nrow(m.chip) - chip.relation2cpgIsland[i,3] - nrow(chip.meta_chs_aric.replicated)), nrow = 2) )[[2]][1], 2), 
                                                   round(fisher.test(matrix(c(chip.relation2cpgIsland[i,2], nrow(chip.meta_chs_aric.replicated) - chip.relation2cpgIsland[i,2], chip.relation2cpgIsland[i,3] - chip.relation2cpgIsland[i,2], nrow(m.chip) - chip.relation2cpgIsland[i,3] - nrow(chip.meta_chs_aric.replicated)), nrow = 2) )[[2]][2], 2), sep="-")
  chip.relation2cpgIsland$OR[i] <- round(fisher.test(matrix(c(chip.relation2cpgIsland[i,2], nrow(chip.meta_chs_aric.replicated) - chip.relation2cpgIsland[i,2], chip.relation2cpgIsland[i,3] - chip.relation2cpgIsland[i,2], nrow(m.chip) - chip.relation2cpgIsland[i,3] - nrow(chip.meta_chs_aric.replicated)), nrow = 2) )[[3]],2)
}
chip.relation2cpgIsland$CHIP_cat <- "Any CHIP"
# VAF
vaf.relation2cpgIsland <- cbind.data.frame(
  category=c("Open Sea", "Island", "Shore", "Shelf"),
  CpGs_in_mylist=c(table(m.vaf$Relation_to_UCSC_CpG_Island[m.vaf$CPG.Labels %in% vaf10.meta_chs_aric.replicated$CPG.Labels]=="")[["TRUE"]],table(m.vaf$Relation_to_UCSC_CpG_Island[m.vaf$CPG.Labels %in% vaf10.meta_chs_aric.replicated$CPG.Labels]=="Island")[["TRUE"]], table(grepl(pattern = "Shore",x = m.vaf$Relation_to_UCSC_CpG_Island[m.vaf$CPG.Labels %in% vaf10.meta_chs_aric.replicated$CPG.Labels]))[["TRUE"]], table(grepl(pattern = "Shelf",x = m.vaf$Relation_to_UCSC_CpG_Island[m.vaf$CPG.Labels %in% vaf10.meta_chs_aric.replicated$CPG.Labels]))[["TRUE"]]),
  CpGs_in_path = c(table(m.vaf$Relation_to_UCSC_CpG_Island=="")[["TRUE"]], table(m.vaf$Relation_to_UCSC_CpG_Island=="Island")[["TRUE"]], table(grepl(pattern = "Shore",x = m.vaf$Relation_to_UCSC_CpG_Island))[["TRUE"]], table(grepl(pattern = "Shelf",x = m.vaf$Relation_to_UCSC_CpG_Island))[["TRUE"]]))
vaf.relation2cpgIsland$OR <- NA
vaf.relation2cpgIsland$CI_95percent <- NA
vaf.relation2cpgIsland$pval <- NA
for(i in 1:nrow(vaf.relation2cpgIsland)) {
  vaf.relation2cpgIsland$pval[i] <- formatC(fisher.test(matrix(c(vaf.relation2cpgIsland[i,2], nrow(vaf10.meta_chs_aric.replicated) - vaf.relation2cpgIsland[i,2], vaf.relation2cpgIsland[i,3] - vaf.relation2cpgIsland[i,2], nrow(m.vaf) - vaf.relation2cpgIsland[i,3] - nrow(vaf10.meta_chs_aric.replicated)), nrow = 2) )[[1]], digits = 1, format = "E")
  vaf.relation2cpgIsland$CI_95percent[i] <- paste(round(fisher.test(matrix(c(vaf.relation2cpgIsland[i,2], nrow(vaf10.meta_chs_aric.replicated) - vaf.relation2cpgIsland[i,2], vaf.relation2cpgIsland[i,3] - vaf.relation2cpgIsland[i,2], nrow(m.vaf) - vaf.relation2cpgIsland[i,3] - nrow(vaf10.meta_chs_aric.replicated)), nrow = 2) )[[2]][1],2),
                                                  round(fisher.test(matrix(c(vaf.relation2cpgIsland[i,2], nrow(vaf10.meta_chs_aric.replicated) - vaf.relation2cpgIsland[i,2], vaf.relation2cpgIsland[i,3] - vaf.relation2cpgIsland[i,2], nrow(m.vaf) - vaf.relation2cpgIsland[i,3] - nrow(vaf10.meta_chs_aric.replicated)), nrow = 2) )[[2]][2],2),sep="-")
  vaf.relation2cpgIsland$OR[i] <- round(fisher.test(matrix(c(vaf.relation2cpgIsland[i,2], nrow(vaf10.meta_chs_aric.replicated) - vaf.relation2cpgIsland[i,2], vaf.relation2cpgIsland[i,3] - vaf.relation2cpgIsland[i,2], nrow(m.vaf) - vaf.relation2cpgIsland[i,3] - nrow(vaf10.meta_chs_aric.replicated)), nrow = 2) )[[3]],2)
}
vaf.relation2cpgIsland$CHIP_cat <- "Expanded CHIP"
# DNMT3A
dnmt3a.relation2cpgIsland <- cbind.data.frame(
  category=c("Open Sea", "Island", "Shore", "Shelf"),
  CpGs_in_mylist=c(table(m.dnmt3a$Relation_to_UCSC_CpG_Island[m.dnmt3a$CPG.Labels %in% dnmt3a.meta_chs_aric.replicated$CPG.Labels]=="")[["TRUE"]],table(m.dnmt3a$Relation_to_UCSC_CpG_Island[m.dnmt3a$CPG.Labels %in% dnmt3a.meta_chs_aric.replicated$CPG.Labels]=="Island")[["TRUE"]], table(grepl(pattern = "Shore",x = m.dnmt3a$Relation_to_UCSC_CpG_Island[m.dnmt3a$CPG.Labels %in% dnmt3a.meta_chs_aric.replicated$CPG.Labels]))[["TRUE"]], table(grepl(pattern = "Shelf",x = m.dnmt3a$Relation_to_UCSC_CpG_Island[m.dnmt3a$CPG.Labels %in% dnmt3a.meta_chs_aric.replicated$CPG.Labels]))[["TRUE"]]),
  CpGs_in_path = c(table(m.dnmt3a$Relation_to_UCSC_CpG_Island=="")[["TRUE"]], table(m.dnmt3a$Relation_to_UCSC_CpG_Island=="Island")[["TRUE"]], table(grepl(pattern = "Shore",x = m.dnmt3a$Relation_to_UCSC_CpG_Island))[["TRUE"]], table(grepl(pattern = "Shelf",x = m.dnmt3a$Relation_to_UCSC_CpG_Island))[["TRUE"]]))
dnmt3a.relation2cpgIsland$OR <- NA
dnmt3a.relation2cpgIsland$CI_95percent <- NA
dnmt3a.relation2cpgIsland$pval <- NA
for(i in 1:nrow(dnmt3a.relation2cpgIsland)) {
  dnmt3a.relation2cpgIsland$pval[i] <- formatC(fisher.test(matrix(c(dnmt3a.relation2cpgIsland[i,2], nrow(dnmt3a.meta_chs_aric.replicated) - dnmt3a.relation2cpgIsland[i,2], dnmt3a.relation2cpgIsland[i,3] - dnmt3a.relation2cpgIsland[i,2], nrow(m.dnmt3a) - dnmt3a.relation2cpgIsland[i,3] - nrow(dnmt3a.meta_chs_aric.replicated)), nrow = 2) )[[1]], digits = 1,format = "E")
  dnmt3a.relation2cpgIsland$CI_95percent[i] <- paste(round(fisher.test(matrix(c(dnmt3a.relation2cpgIsland[i,2], nrow(dnmt3a.meta_chs_aric.replicated) - dnmt3a.relation2cpgIsland[i,2], dnmt3a.relation2cpgIsland[i,3] - dnmt3a.relation2cpgIsland[i,2], nrow(m.dnmt3a) - dnmt3a.relation2cpgIsland[i,3] - nrow(dnmt3a.meta_chs_aric.replicated)), nrow = 2) )[[2]][1],2),
                                                     round(fisher.test(matrix(c(dnmt3a.relation2cpgIsland[i,2], nrow(dnmt3a.meta_chs_aric.replicated) - dnmt3a.relation2cpgIsland[i,2], dnmt3a.relation2cpgIsland[i,3] - dnmt3a.relation2cpgIsland[i,2], nrow(m.dnmt3a) - dnmt3a.relation2cpgIsland[i,3] - nrow(dnmt3a.meta_chs_aric.replicated)), nrow = 2) )[[2]][2],2),sep = "-")
  dnmt3a.relation2cpgIsland$OR[i] <- round(fisher.test(matrix(c(dnmt3a.relation2cpgIsland[i,2], nrow(dnmt3a.meta_chs_aric.replicated) - dnmt3a.relation2cpgIsland[i,2], dnmt3a.relation2cpgIsland[i,3] - dnmt3a.relation2cpgIsland[i,2], nrow(m.dnmt3a) - dnmt3a.relation2cpgIsland[i,3] - nrow(dnmt3a.meta_chs_aric.replicated)), nrow = 2) )[[3]],2)
}
dnmt3a.relation2cpgIsland$CHIP_cat <- "DNMT3A"
# Tet2
tet2.relation2cpgIsland <- cbind.data.frame(
  category=c("Open Sea", "Island", "Shore", "Shelf"),
  CpGs_in_mylist=c(table(m.tet2$Relation_to_UCSC_CpG_Island[m.tet2$CPG.Labels %in% tet2.meta_chs_aric.replicated$CPG.Labels]=="")[["TRUE"]],table(m.tet2$Relation_to_UCSC_CpG_Island[m.tet2$CPG.Labels %in% tet2.meta_chs_aric.replicated$CPG.Labels]=="Island")[["TRUE"]], table(grepl(pattern = "Shore",x = m.tet2$Relation_to_UCSC_CpG_Island[m.tet2$CPG.Labels %in% tet2.meta_chs_aric.replicated$CPG.Labels]))[["TRUE"]], table(grepl(pattern = "Shelf",x = m.tet2$Relation_to_UCSC_CpG_Island[m.tet2$CPG.Labels %in% tet2.meta_chs_aric.replicated$CPG.Labels]))[["TRUE"]]),
  CpGs_in_path = c(table(m.tet2$Relation_to_UCSC_CpG_Island=="")[["TRUE"]], table(m.tet2$Relation_to_UCSC_CpG_Island=="Island")[["TRUE"]], table(grepl(pattern = "Shore",x = m.tet2$Relation_to_UCSC_CpG_Island))[["TRUE"]], table(grepl(pattern = "Shelf",x = m.tet2$Relation_to_UCSC_CpG_Island))[["TRUE"]]))
tet2.relation2cpgIsland$OR <- NA
tet2.relation2cpgIsland$CI_95percent <- NA
# tet2.relation2cpgIsland$uCI <- NA
tet2.relation2cpgIsland$pval <- NA
for(i in 1:nrow(tet2.relation2cpgIsland)) {
  tet2.relation2cpgIsland$pval[i] <- formatC(fisher.test(matrix(c(tet2.relation2cpgIsland[i,2], nrow(tet2.meta_chs_aric.replicated) - tet2.relation2cpgIsland[i,2], tet2.relation2cpgIsland[i,3] - tet2.relation2cpgIsland[i,2], nrow(m.tet2) - tet2.relation2cpgIsland[i,3] - nrow(tet2.meta_chs_aric.replicated)), nrow = 2) )[[1]], format = "E", digits = 1)
  tet2.relation2cpgIsland$OR[i] <- round(fisher.test(matrix(c(tet2.relation2cpgIsland[i,2], nrow(tet2.meta_chs_aric.replicated) - tet2.relation2cpgIsland[i,2], tet2.relation2cpgIsland[i,3] - tet2.relation2cpgIsland[i,2], nrow(m.tet2) - tet2.relation2cpgIsland[i,3] - nrow(tet2.meta_chs_aric.replicated)), nrow = 2) )[[3]],2)
  tet2.relation2cpgIsland$CI_95percent[i] <- paste(round(fisher.test(matrix(c(tet2.relation2cpgIsland[i,2], nrow(tet2.meta_chs_aric.replicated) - tet2.relation2cpgIsland[i,2], tet2.relation2cpgIsland[i,3] - tet2.relation2cpgIsland[i,2], nrow(m.tet2) - tet2.relation2cpgIsland[i,3] - nrow(tet2.meta_chs_aric.replicated)), nrow = 2) )[[2]][1],2), round(fisher.test(matrix(c(tet2.relation2cpgIsland[i,2], nrow(tet2.meta_chs_aric.replicated) - tet2.relation2cpgIsland[i,2], tet2.relation2cpgIsland[i,3] - tet2.relation2cpgIsland[i,2], nrow(m.tet2) - tet2.relation2cpgIsland[i,3] - nrow(tet2.meta_chs_aric.replicated)), nrow = 2) )[[2]][2],2),sep = "-")
  
}
tet2.relation2cpgIsland$CHIP_cat <- "TET2"

# combined Enrichment of cpg types
relation2cpgIsland <- rbind.data.frame(chip.relation2cpgIsland, vaf.relation2cpgIsland, dnmt3a.relation2cpgIsland, tet2.relation2cpgIsland)


write.csv(relation2cpgIsland, "Final.chs.enrichment_of_CpGs_in_Island_shore_shelf_openSea.csv", row.names = F)
#########



### DMR enrichment
## DMR
dnmt3a.dmr <- merge(as.data.frame(table(dnmt3a.meta_chs_aric.replicated$DMR)), as.data.frame(table(m.dnmt3a$DMR)), by="Var1" )
names(dnmt3a.dmr) <- c("DMR", "CpGs_in_my_lits", "CpGs_in_path")
dnmt3a.dmr$OR <- NA
dnmt3a.dmr$CI_95percent <- NA
dnmt3a.dmr$pval <- NA
for(i in 1:nrow(dnmt3a.dmr)) {
  dnmt3a.dmr$pval[i] <- formatC(fisher.test(matrix(c(dnmt3a.dmr[i,2], nrow(dnmt3a.meta_chs_aric.replicated) - dnmt3a.dmr[i,2], dnmt3a.dmr[i,3] - dnmt3a.dmr[i,2], nrow(m.dnmt3a) - dnmt3a.dmr[i,3] - nrow(dnmt3a.meta_chs_aric.replicated)), nrow = 2) )[[1]],format = "E", digits = 1)
  dnmt3a.dmr$CI_95percent[i] <- paste(round(fisher.test(matrix(c(dnmt3a.dmr[i,2], nrow(dnmt3a.meta_chs_aric.replicated) - dnmt3a.dmr[i,2], dnmt3a.dmr[i,3] - dnmt3a.dmr[i,2], nrow(m.dnmt3a) - dnmt3a.dmr[i,3] - nrow(dnmt3a.meta_chs_aric.replicated)), nrow = 2) )[[2]][1], 2),
                                      round(fisher.test(matrix(c(dnmt3a.dmr[i,2], nrow(dnmt3a.meta_chs_aric.replicated) - dnmt3a.dmr[i,2], dnmt3a.dmr[i,3] - dnmt3a.dmr[i,2], nrow(m.dnmt3a) - dnmt3a.dmr[i,3] - nrow(dnmt3a.meta_chs_aric.replicated)), nrow = 2) )[[2]][2],2), sep="-")
  dnmt3a.dmr$OR[i] <- round(fisher.test(matrix(c(dnmt3a.dmr[i,2], nrow(dnmt3a.meta_chs_aric.replicated) - dnmt3a.dmr[i,2], dnmt3a.dmr[i,3] - dnmt3a.dmr[i,2], nrow(m.dnmt3a) - dnmt3a.dmr[i,3] - nrow(dnmt3a.meta_chs_aric.replicated)), nrow = 2) )[[3]],2)
}
dnmt3a.dmr$CHIP_cat <- "DNMT3A"
# Benjamini & Hochberg (1995) ("BH" or its alias "fdr")
#dnmt3a.dmr$FDR <- p.adjust(dnmt3a.dmr$pval, method = "BH")
#TET
tet2.dmr <- merge(as.data.frame(table(tet2.meta_chs_aric.replicated$DMR)), as.data.frame(table(m.tet2$DMR)), by="Var1" )
names(tet2.dmr) <- c("DMR", "CpGs_in_my_lits", "CpGs_in_path")
tet2.dmr$OR <- NA
tet2.dmr$CI_95percent <- NA
tet2.dmr$pval <- NA
for(i in 1:nrow(tet2.dmr)) {
  tet2.dmr$pval[i] <- formatC(fisher.test(matrix(c(tet2.dmr[i,2], nrow(tet2.meta_chs_aric.replicated) - tet2.dmr[i,2], tet2.dmr[i,3] - tet2.dmr[i,2], nrow(m.tet2) - tet2.dmr[i,3] - nrow(tet2.meta_chs_aric.replicated)), nrow = 2) )[[1]],format = "E", digits = 1)
  tet2.dmr$CI_95percent[i] <- paste(round(fisher.test(matrix(c(tet2.dmr[i,2], nrow(tet2.meta_chs_aric.replicated) - tet2.dmr[i,2], tet2.dmr[i,3] - tet2.dmr[i,2], nrow(m.tet2) - tet2.dmr[i,3] - nrow(tet2.meta_chs_aric.replicated)), nrow = 2) )[[2]][1],2), 
                                    round(fisher.test(matrix(c(tet2.dmr[i,2], nrow(tet2.meta_chs_aric.replicated) - tet2.dmr[i,2], tet2.dmr[i,3] - tet2.dmr[i,2], nrow(m.tet2) - tet2.dmr[i,3] - nrow(tet2.meta_chs_aric.replicated)), nrow = 2) )[[2]][2],2),sep="-")
  tet2.dmr$OR[i] <- round(fisher.test(matrix(c(tet2.dmr[i,2], nrow(tet2.meta_chs_aric.replicated) - tet2.dmr[i,2], tet2.dmr[i,3] - tet2.dmr[i,2], nrow(m.tet2) - tet2.dmr[i,3] - nrow(tet2.meta_chs_aric.replicated)), nrow = 2) )[[3]],2)
}
tet2.dmr$CHIP_cat <- "TET2"

# Any CHIP
chip.dmr <- merge(as.data.frame(table(chip.meta_chs_aric.replicated$DMR)), as.data.frame(table(m.chip$DMR)), by="Var1" )
names(chip.dmr) <- c("DMR", "CpGs_in_my_lits", "CpGs_in_path")
chip.dmr$OR <- NA
chip.dmr$CI_95percent <- NA
chip.dmr$pval <- NA
for(i in 1:nrow(chip.dmr)) {
  chip.dmr$pval[i] <- formatC(fisher.test(matrix(c(chip.dmr[i,2], nrow(chip.meta_chs_aric.replicated) - chip.dmr[i,2], chip.dmr[i,3] - chip.dmr[i,2], nrow(m.chip) - chip.dmr[i,3] - nrow(chip.meta_chs_aric.replicated)), nrow = 2) )[[1]],format = "E", digits = 1)
  chip.dmr$CI_95percent[i] <- paste(round(fisher.test(matrix(c(chip.dmr[i,2], nrow(chip.meta_chs_aric.replicated) - chip.dmr[i,2], chip.dmr[i,3] - chip.dmr[i,2], nrow(m.chip) - chip.dmr[i,3] - nrow(chip.meta_chs_aric.replicated)), nrow = 2) )[[2]][1],2), 
                                    round(fisher.test(matrix(c(chip.dmr[i,2], nrow(chip.meta_chs_aric.replicated) - chip.dmr[i,2], chip.dmr[i,3] - chip.dmr[i,2], nrow(m.chip) - chip.dmr[i,3] - nrow(chip.meta_chs_aric.replicated)), nrow = 2) )[[2]][2],2),sep="-")
  chip.dmr$OR[i] <- round(fisher.test(matrix(c(chip.dmr[i,2], nrow(chip.meta_chs_aric.replicated) - chip.dmr[i,2], chip.dmr[i,3] - chip.dmr[i,2], nrow(m.chip) - chip.dmr[i,3] - nrow(chip.meta_chs_aric.replicated)), nrow = 2) )[[3]],2)
}
chip.dmr$CHIP_cat <- "Any CHIP"

# Expanded CHIP
vaf10.dmr <- merge(as.data.frame(table(vaf10.meta_chs_aric.replicated$DMR)), as.data.frame(table(m.vaf$DMR)), by="Var1" )
names(vaf10.dmr) <- c("DMR", "CpGs_in_my_lits", "CpGs_in_path")
vaf10.dmr$OR <- NA
vaf10.dmr$CI_95percent <- NA
vaf10.dmr$pval <- NA
for(i in 1:nrow(vaf10.dmr)) {
  vaf10.dmr$pval[i] <- formatC(fisher.test(matrix(c(vaf10.dmr[i,2], nrow(vaf10.meta_chs_aric.replicated) - vaf10.dmr[i,2], vaf10.dmr[i,3] - vaf10.dmr[i,2], nrow(m.vaf) - vaf10.dmr[i,3] - nrow(vaf10.meta_chs_aric.replicated)), nrow = 2) )[[1]],format = "E", digits = 1)
  vaf10.dmr$CI_95percent[i] <- paste(round(fisher.test(matrix(c(vaf10.dmr[i,2], nrow(vaf10.meta_chs_aric.replicated) - vaf10.dmr[i,2], vaf10.dmr[i,3] - vaf10.dmr[i,2], nrow(m.vaf) - vaf10.dmr[i,3] - nrow(vaf10.meta_chs_aric.replicated)), nrow = 2) )[[2]][1],2), 
                                     round(fisher.test(matrix(c(vaf10.dmr[i,2], nrow(vaf10.meta_chs_aric.replicated) - vaf10.dmr[i,2], vaf10.dmr[i,3] - vaf10.dmr[i,2], nrow(m.vaf) - vaf10.dmr[i,3] - nrow(vaf10.meta_chs_aric.replicated)), nrow = 2) )[[2]][2],2),sep="-")
  vaf10.dmr$OR[i] <- round(fisher.test(matrix(c(vaf10.dmr[i,2], nrow(vaf10.meta_chs_aric.replicated) - vaf10.dmr[i,2], vaf10.dmr[i,3] - vaf10.dmr[i,2], nrow(m.vaf) - vaf10.dmr[i,3] - nrow(vaf10.meta_chs_aric.replicated)), nrow = 2) )[[3]],2)
}
vaf10.dmr$CHIP_cat <- "Expanded CHIP"

# combined DMR Enrichment of cpg types
cpgs_in_dmr <- rbind.data.frame(chip.dmr, vaf10.dmr, dnmt3a.dmr, tet2.dmr)

write.csv(cpgs_in_dmr, "Final.cpgs_in_dmr.csv", row.names = F)

write.csv(tulstrup.tet2.dmr, "Documents//Project/Methylation/Revision/tulstrup.tet2.dmr.csv")
write.csv(tulstrup.tet2.relation2cpgIsland, "Documents//Project/Methylation/Revision/tulstrup.tet2.relation2cpgIsland.csv")

save.image(file = "Documents/Project/Methylation/Revision/Latest_DNAm_enrichemnt_data.rda")
###########################################################################
