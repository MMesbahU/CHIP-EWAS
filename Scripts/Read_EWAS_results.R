############################################################################
### Read_EWAS_results.R - prepares EWAS output for functional annotation 
### Read in meta-EWAS results for CHS (discovery) and ARIC (replication)
### Output: intermediate_OurResults.rda
### Author: Karen Conneely
############################################################################

############################################################################
### CHS results (as output tables from METAL)
### Throughout, d prefix indicates results from the DNMT3A CHIP EWAS
### while t prefix indicates results from the TET2 CHIP EWAS
############################################################################

d_ewas = read.delim('/Users/kconnee/Desktop/CHIP/SuppTab1/chs_AA_EA.meta_ewas.DNMT3A.TBL',as.is=T)
t_ewas = read.delim('/Users/kconnee/Desktop/CHIP/SuppTab1/chs_AA_EA.meta_ewas.TET2.TBL',as.is=T)

### Add column for Benjamini-Hochberg FDR estimate
d_ewas[,16] = p.adjust(d_ewas$P.value,method="BH")
t_ewas[,16] = p.adjust(t_ewas$P.value,method="BH")
colnames(d_ewas)[16] = colnames(t_ewas)[16] = "FDR"

### Extract summary statistics useful for functional annotation
dnmt3a_ewas = d_ewas[,c("MarkerName","Effect","StdErr","P.value","FDR")]
tet2_ewas = t_ewas[,c("MarkerName","Effect","StdErr","P.value","FDR")]

############################################################################
### Load sets of CpGs that replicated in ARIC EWAS
### (as R dataset containing summary statistics ONLY for replicated CpGs)
############################################################################

load("/Users/kconnee/Desktop/CHIP/SuppTab1/Replicated.CpGs.meta_ewas.chip_vaf10_dnmt3a_tet2.rda")
dnmt3a_rep = dnmt3a.meta_chs_aric.replicated
tet2_rep = tet2.meta_chs_aric.replicated

### Add column to CHS summary statistics indicating TRUE/FALSE for replication in ARIC
dnmt3a_ewas[,6] = (dnmt3a_ewas[,1]%in%dnmt3a_rep[,1])
tet2_ewas[,6] = (tet2_ewas[,1]%in%tet2_rep[,1])

### combine DNMT3A and TET2 results into a single data.frame of 478,687 CpG sites, and save 

our_res = cbind(dnmt3a_ewas,tet2_ewas)
colnames(our_res)[c(6,12)] = c("RepD","RepT")
save(our_res,file="intermediate_OurResults.rda")

