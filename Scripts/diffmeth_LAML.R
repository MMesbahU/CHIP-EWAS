########################################################################################################################
### diffmeth_LAML.R: Perform EWAS of AML by driver mutation, to facilitate comparison to our results
### AML EWAS is performed below using data from TCGA et al., N Engl J Med. 2013 368(22):2059-74
### An intermediate file is produced that is then used by the program Figsef.R
### Code is commented out because it requires large files and doesn't need to be repeated - just use the intermediate file.
### Output: diffmeth_LAML_intermediate.rda
### Author: Karen Conneely
########################################################################################################################

########################################################################################################################
### Read in phenotype data from published supplement (TCGA et al., NEJM 2013)
########################################################################################################################

#samp = read.csv('SuppTable01.csv',as.is=T)
#patients = samp[,c("TCGA.Patient.ID","Sex","Race","Age","FAB","Cytogenetic.Classification","DNMT3A","TET2")]

########################################################################################################################
### Create indicator variable for mutation type
########################################################################################################################

#patients[,9] = NA
#colnames(patients)[9] = "mutation"
#attach(patients)

#patients[DNMT3A==""&TET2=="",9]="neither"
#patients[DNMT3A!=""&TET2=="",9]="DNMT3A"
#patients[DNMT3A==""&TET2!="",9]="TET2"
#patients[DNMT3A!=""&TET2!="",9]="both"

#k=order(patients[,1])
#patients = patients[k,]

########################################################################################################################
### Read in DNAm data from published supplement (TCGA et al., NEJM 2013)
########################################################################################################################

#load('LAML_betas.rda')
#k=order(colnames(laml_beta))
#laml_beta=laml_beta[,k]

########################################################################################################################
### Extract subset of phenotype data with DNAm (sub) and verify that the patients and sort order are the same
########################################################################################################################

#sub = patients[patients[,1]%in%colnames(laml_beta),]
#print(table(colnames(laml_beta)==sub[,1]))

########################################################################################################################
### Perform EWAS by mutation type, comparing LAML patients with DNMT3A or TET2 mutation to those with other mutations
########################################################################################################################

#library(CpGassoc)

#keep = sub$mutation%in%c("DNMT3A","neither")
#ewas_d = cpg.assoc(laml_beta[,keep],(sub$mutation[keep]=="DNMT3A"),covariates=sub[keep,c("Sex","Age")])

#keep2 = sub$mutation%in%c("TET2","neither")
#ewas_t = cpg.assoc(laml_beta[,keep2],(sub$mutation[keep2]=="TET2"),covariates=sub[keep2,c("Sex","Age")])

#res = cbind(ewas_d$results,ewas_t$results)

#save(res,file="diffmeth_LAML_intermediate.rda")