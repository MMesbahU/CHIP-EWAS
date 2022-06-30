
############
ARGs <- commandArgs(TRUE)

methylation_file <- ARGs[1]

output_filename <- ARGs[2]

sample_info_n_age <- ARGs[3]

Ancestry_2_test <- as.character(ARGs[4])
#################
# Load Libraries  
library(data.table)

library(CpGassoc)

# Load Methylation sample info files with Meth_Age 
load(sample_info_n_age)

#### 
tet2 <- subset(CHS_TET2, CHS_TET2$Ancestry == Ancestry_2_test)
cat("TET2 status:",table(tet2$hasTET2),"\n")
cat ("Ancestry :", Ancestry_2_test, "\t")
cat("Unique TET2 carriers:",length(unique(tet2$NWD_ID[tet2$hasTET2==0])), length(unique(tet2$NWD_ID[tet2$hasTET2==1])), "\n") 

# Load SWAN (or Raw) data
methylation_dat <- fread(methylation_file, header=T, stringsAsFactors=F, sep=',' , fill = TRUE)

  # TET2
tet2_sample_list <- names(methylation_dat) %in% tet2$sample_id.y
tet2_final_methylation_dat <- methylation_dat [ , ..tet2_sample_list ]

  # Reorder
tet2_ordered_methylation_dat0 <- setcolorder( tet2_final_methylation_dat, tet2$sample_id.y )
  # Convert DT to DF
tet2_ordered_methylation_dat0 <- setDF(tet2_ordered_methylation_dat0)
  # Add row names from Column 1
.rowNamesDF(tet2_ordered_methylation_dat0, make.names=F) <- methylation_dat$V1
tet2_ordered_methylation_dat <- tet2_ordered_methylation_dat0
#
tet2_sanity_check = as.data.frame(matrix(NA, nrow=length(names(tet2_ordered_methylation_dat)), ncol=1))
for(i in 1:nrow(tet2_sanity_check)) { 
  tet2_sanity_check$V1[i] <- names(tet2_ordered_methylation_dat)[i] == tet2$sample_id.y[i]
}

table(tet2_sanity_check)


############ table should show only TRUE values
# clean unnecessary files
rm(tet2_final_methylation_dat) 
# rm(methylation_dat)

############################# Run Association ###########################
# 1. EWAS1: CHIP as numeric 0/1; Stratified by Ancestry
# DNAm ~ CHIP_Gene + Age_Meth + Age_Meth2 + Gender + Batch + 5_cell_types
tet2_ewas_age2_EA_AA <- cpg.assoc(tet2_ordered_methylation_dat, 
                                  tet2$hasTET2, ~tet2$agemeth + 
                                    tet2$agemeth_sqr + 
                                    tet2$Gender + tet2$Batch_factor + 
                                    tet2$CD8T + tet2$CD4T + tet2$NK + 
                                    tet2$Bcell + tet2$Mono, 
                                  chip.id=tet2$NWD_ID, 
                                  random=TRUE, large.data = TRUE)

  # Save CHIP ewas results
save(tet2_ewas_age2_EA_AA, file = output_filename)
############################################

