#####################################################################################################
### Fig3d.R - create dataset to upload to eFORGE tool, and plot results
### Makes CpG lists that are then provided to eFORGE, with ENCODE option selected
### Compare methylation profiles across cell types for CpG sites replicated in DNMT3a and TET2 EWAS 
### Output: Figd.jpg
### Author: Karen Conneely
#####################################################################################################

#########################################################################################################
### Step 1: Read in lists of successfully replicated CpGs, and format into table for upload to eFORGE
#########################################################################################################

load("/Users/kconnee/Desktop/CHIP/SuppTab1/Replicated.CpGs.meta_ewas.chip_vaf10_dnmt3a_tet2.rda")

dnmt3a_rep = dnmt3a.meta_chs_aric.replicated
tet2_rep = tet2.meta_chs_aric.replicated

### Sort by meta-analysis p-value so we can easily grab the 1000 most significant

k=order(dnmt3a_rep$meta_chs_aric.P.value)
dnmt3a_rep=dnmt3a_rep[k,]

k=order(tet2_rep$meta_chs_aric.P.value)
tet2_rep=tet2_rep[k,]

write.table(dnmt3a_rep[1:1000,1],file='d3.csv',quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(tet2_rep[1:1000,1],file='t2.csv',quote=FALSE,row.names=FALSE,col.names=FALSE)

#########################################################################################################
###Step 2: upload CpG list to eFORGE, and run with ENCODE option selected
#########################################################################################################

#########################################################################################################
### Step 3: Read in results from eFORGE site (named DNMT3A.450K and TET2.450K), and create a barplot
#########################################################################################################

    d_enrich = read.delim('/Users/kconnee/Desktop/CHIP/SuppTab1/DNMT3A.450k.encode.chart.tsv', as.is=T)
    t_enrich = read.delim('/Users/kconnee/Desktop/CHIP/SuppTab1/TET2.450k.encode.chart.tsv', as.is=T)

    ### Double-check that both output tables have the same cell types in the same rows, prior to cbinding
    print(table(d_enrich[,c(3:6,8)]==t_enrich[,c(3:6,8)]))
  ###  TRUE 
  ###   625 

    ### Concatenate results for plotting to compare DNMT3A to TET2 results
    comp_final = cbind(d_enrich[,c(3:6,8,1,2,9)],t_enrich[,c(1,2,9)])
    
    print(table(comp_final$Tissue))
 
keep=(comp_final$Tissue=="Blood"&substr(comp_final$Cell,1,2)!="GM")
temp = as.matrix(-log10(comp_final[keep,c(7,10)]))
row.names(temp)=comp_final$Cell[keep]
colnames(temp) = c("DNMT3A","TET2")
temp2=t(temp)

library(yarrr)
blue50 = transparent(orig.col = "cornflowerblue", trans.val = 0.5, maxColorValue = 255)
green50 = transparent(orig.col = "seagreen", trans.val = 0.5, maxColorValue = 255)
blue75 = transparent(orig.col = "cornflowerblue", trans.val = 0.25, maxColorValue = 255)
green75 = transparent(orig.col = "seagreen", trans.val = 0.25, maxColorValue = 255)

### Restrict to just blood cell types for barplot of -log(enrichment P-values)
keep=(comp_final$Tissue=="Blood"&substr(comp_final$Cell,1,2)!="GM")
temp = as.matrix(-log10(comp_final[keep,c(7,10)]))
row.names(temp)=comp_final$Cell[keep]
colnames(temp) = c("DNMT3A","TET2")
temp2=t(temp)
temp3=temp2[,c(4,2,3,1,11,12)]

coln = colnames(temp3)
coln2=coln
coln2[coln=="CD34+"] = "HSCs"
coln2[coln=="CD14+"] = "Monocytes"
coln2[coln=="CD20+"] = "B cells"
coln2[coln=="Adult_CD4+"] = "Naive T"
colnames(temp3) = coln2

jpeg(file="Fig3d.jpg")
par(xpd=NA,font.axis=2,font.lab=2,cex.axis=1.2,cex.lab=1.2)
barplot(temp3,beside=T,las=2,col=c("cornflower blue","seagreen"),ylab=NULL,mar=c(10, 4, 2, 2),ylim=c(0,100),xaxt="n")
legend("topright",c(expression(italic("DNMT3A")-"associated CpGs"),expression(italic("TET2")-"associated CpGs")),fill=c("cornflower blue","seagreen"),cex=1.2)
text(x = 3*(1:6)-1, y = -5, labels = coln2, xpd = NA, srt = 35, adj = 0.965, cex = 1.2,font=2)

abline(-log10(.05/12),0,lty=2)
lines(c(0,18),c(0,0))

dev.off()
