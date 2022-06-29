################################################################################################
### Figs3ef.R: Compare DNMT3A and TET2 CHIP EWAS results to an EWAS of those mutations in LAML
### Output: Figs3ef.jpg
### Author: Karen Conneely
################################################################################################

#####################################################################################################
### Reads in results from our LAML EWAS of TCGA data (intermediate file produced by diffmeth_LAML.R)
### Reads in intermediate results from our CHIP EWAS
#####################################################################################################
load(file="diffmeth_LAML_intermediate.rda")
load(file="intermediate_OurResults.rda")

###Merge results files
allres=merge(our_res,res,by=1)

###Examine overlap in contingency tables.  

print(table(allres[,6],allres[,16]<.05))
print(fisher.test(table(allres[,6],allres[,16]<.05)))
#DNMT3A CHIP (replicated) vs. LAML DNMT3A (FDR<.05):   OR=19.61717 
#  FALSE 376751  11633
#  TRUE    2067   1252

print(table(allres[,12],allres[,22]<.05))
print(fisher.test(table(allres[,12],allres[,22]<.05)))
#TET2 CHIP (replicated) vs. LAML TET2 (FDR<.05): OR=0
#  FALSE 390482     12
#  TRUE    1209      0

print(table(allres[,6],allres[,14]<.05))
print(fisher.test(table(allres[,6],allres[,14]<.05)))
##DNMT3A CHIP (replicated) vs. LAML DNMT3A (p<.05): OR=10.52225
#  FALSE 322303  66081
#  TRUE    1051   2268

print(table(allres[,12],allres[,20]<.05))
print(fisher.test(table(allres[,12],allres[,20]<.05)))
#TET2 CHIP (replicated) vs. LAML TET2 (p<.05): OR=6.695799 
#  FALSE 365204  25290
#  TRUE     826    383


### Correlation across all four EWAS (DNMT3A and TET2 CHIP, LAML with DNMT3A and TET2)
cormat = cor(cbind(cbind(allres[,c(2,8)]/allres[,c(3,9)]),allres[,c(13,19)]),use="complete.obs")

### Create Figures 3e & f
par(mfcol=c(1,2),xpd=FALSE,mgp=c(2.3,1,0),font.lab=2)
temp=allres[allres[,5]<.05,]
txt <- format(c(cormat[1,3], 0.123456789), digits = 2)[1]
    txt <- paste0("r = ", txt)
    
colmat = matrix(1,dim(temp)[1],1)
colmat[temp[,6]] = "cornflower blue"
plot(temp[,2]/temp[,3], temp[,13], col=colmat, cex=.3, xlim=c(-10,10), ylim=c(-10,10), font=2, font.lab=2, xlab=expression(italic("DNMT3A")~"CHIP EWAS"), ylab = expression(italic("DNMT3A")~"mutation in AML"))
abline(0,0)
lines(c(0,0),c(-12,11))
text(-6.5,7.5,txt,font=2)
temp=allres[allres[,11]<.05,]
txt <- format(c(cormat[2,4], 0.123456789), digits = 2)[1]
    txt <- paste0("r = ", txt)

colmat = matrix(1,dim(temp)[1],1)
colmat[temp[,12]] = "seagreen"
plot(temp[,8]/temp[,9], temp[,19], col=colmat, cex=.3, xlim=c(-10,10), ylim=c(-10,10), font=2, font.lab=2, xlab=expression(italic("TET2")~"CHIP EWAS"), ylab = expression(italic("TET2")~"mutation in AML"))

abline(0,0)
lines(c(0,0),c(-12,12))
text(-6.5,7.5,txt,font=2)


dev.copy2pdf(file="Figs3ef.pdf")

