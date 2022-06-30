
############

ARGs <- commandArgs(TRUE)

methylation_file <- ARGs[1]

output_filename <- ARGs[2]

sample_info_n_age <- ARGs[3]

Ancestry_2_test <- ARGs[4]

#################
# Load Libraries  
library(data.table)

library(CpGassoc)

# Load Methylation sample info files with Meth_Age 
load(sample_info_n_age)

#### 
CHIP <- subset(CHS_CHIP_VAF_DNMT3A_TET2, (CHS_CHIP_VAF_DNMT3A_TET2$Ancestry== Ancestry_2_test) & !is.na(CHS_CHIP_VAF_DNMT3A_TET2$hasVAF10) )

cat("VAF10 status:", table(CHIP$hasVAF10),"\n")

cat ("Ancestry :", Ancestry_2_test, "\t")

cat("Unique VAF10 carriers:",length(unique(CHIP$NWD_ID[CHIP$hasVAF10==0])), length(unique(CHIP$NWD_ID[CHIP$hasVAF10==1])), "\n") 

# Load SWAN (or Raw) data
methylation_dat <- fread(methylation_file, header=T, stringsAsFactors=F, sep=',')

  # CHIP
CHIP_sample_list <- names(methylation_dat) %in% CHIP$sample_id.y
CHIP_final_methylation_dat <- methylation_dat [ , ..CHIP_sample_list ]

  # Reorder
CHIP_ordered_methylation_dat0 <- setcolorder( CHIP_final_methylation_dat, CHIP$sample_id.y )
  # Convert DT to DF
CHIP_ordered_methylation_dat0 <- setDF(CHIP_ordered_methylation_dat0)
  # Add row names from Column 1
.rowNamesDF(CHIP_ordered_methylation_dat0, make.names=F) <- methylation_dat$V1
CHIP_ordered_methylation_dat <- CHIP_ordered_methylation_dat0
#
CHIP_sanity_check = as.data.frame(matrix(NA, nrow=length(names(CHIP_ordered_methylation_dat)), ncol=1))
for(i in 1:nrow(CHIP_sanity_check)) { 
  CHIP_sanity_check$V1[i] <- names(CHIP_ordered_methylation_dat)[i] == CHIP$sample_id.y[i]
}

table(CHIP_sanity_check)


############ table should show only TRUE values
# clean unnecessary files
rm(CHIP_final_methylation_dat) 
rm(methylation_dat)
gc()
############################# Run Association ###########################
# 1. EWAS1: CHIP as numeric 0/1 by Ancestry separetely
# DNAm ~ VAF10 + Age_Meth + Age_Meth2 + Gender + Batch + 5_cell_types
VAF10_ewas_age2_EA_AA <- cpg.assoc(CHIP_ordered_methylation_dat, CHIP$hasVAF10, ~CHIP$agemeth + CHIP$agemeth_sqr + CHIP$Gender + CHIP$Batch_factor + CHIP$CD8T + CHIP$CD4T + CHIP$NK + CHIP$Bcell + CHIP$Mono, chip.id=CHIP$NWD_ID, random=TRUE, large.data = TRUE)

  # Save CHIP ewas results
save(VAF10_ewas_age2_EA_AA, file = output_filename)

############################################
  # Extract CpGs for Selected Samples

# load(all_samples)

# CHIP_cpgs_swan <- methylation_dat[methylation_dat$V1 %in% CHIP_ewas_age2$FDR.sig$CPG.Labels, ]

# CHIP_cpgs_colmn_list <- names(CHIP_cpgs_swan) %in% c("V1", CHS_CHIP_DNMT3A_TET2$sample_id.y)

# CHIP_cpgs_swan <- CHIP_cpgs_swan[, ..CHIP_cpgs_colmn_list]

# save(CHIP_cpgs_swan, file = cpg_beta_out)

############################################


