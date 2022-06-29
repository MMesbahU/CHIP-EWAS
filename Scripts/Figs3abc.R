#####################################################################################################
### Figs3abc.R - analyses to create Figures 3a-c
### Compare methylation profiles across cell types for CpG sites replicated in DNMT3a and TET2 EWAS 
### Output: Figs3abc.jpg
### Author: Karen Conneely
#####################################################################################################

### Read in WGS summary data (processed from raw Blueprint data downloaded from GEO using read_WGS.R)
load("/Users/kconnee/Desktop/CHIP/Blueprint/wgs_hsc.rda")
load("/Users/kconnee/Desktop/CHIP/Blueprint/wgs_blood.rda")
load("/Users/kconnee/Desktop/CHIP/Blueprint/wgs_myeloid.rda")

###########################################################################################################
### Create small data.frame ("small") containing relevant Illumina annotation (chromosome, position, etc.)
###########################################################################################################

library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ga=getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
small=cbind(ga[,c("chr","pos","Relation_to_Island","DMR","Enhancer","DHS","Regulatory_Feature_Group","UCSC_RefGene_Group","UCSC_RefGene_Name")])

###########################################################################################################
### Convert our Illumina annotation file (data.frame "small") into GRanges object
###########################################################################################################

library(GenomicRanges)
annot = as.data.frame(small[,c("chr","pos")])
gr_annot = GRanges(seqnames=annot$chr,ranges=IRanges(start=annot$pos,end=annot$pos),names=row.names(annot))

###########################################################################################################
### Perform liftover to convert to hg38 coordinates (using annotation file downloaded from UCSC browser)
###########################################################################################################

library(rtracklayer)
chain <- import.chain("~/Desktop/Annotations/hg19ToHg38.over.chain")
annot38 = liftOver(gr_annot,chain)

############################################################################
### Convert WGS summmary data into GRanges objects
############################################################################

lymphoid_annot = GRanges(seqnames=sum_wgs[,1],ranges=IRanges(start=sum_wgs[,2],end=sum_wgs[,2]))
myeloid_annot = GRanges(seqnames=sum_myeloid[,1],ranges=IRanges(start=sum_myeloid[,2],end=sum_myeloid[,2]))
hsc_annot = GRanges(seqnames=sum_hsc[,1],ranges=IRanges(start=sum_hsc[,2],end=sum_hsc[,2]))

###########################################################################################################
### Create matrices containing positions of CpGs that appear both in Illumina 450K and in each WGS dataset
###########################################################################################################

overlap_lymphoid=as.matrix(findOverlaps(annot38,lymphoid_annot))
overlap_myeloid=as.matrix(findOverlaps(annot38,myeloid_annot))
overlap_hsc=as.matrix(findOverlaps(annot38,hsc_annot))

#####################################################################################
### Use this positional data to create "functional annotation" matrices 
### Matrices contain WGS methylation proportions and other annotation for each CpG
#####################################################################################
annot_lymphoid = cbind(annot[overlap_lymphoid[,1],],sum_wgs[overlap_lymphoid[,2],])
annot_myeloid = cbind(annot[overlap_myeloid[,1],],sum_myeloid[overlap_myeloid[,2],])
annot_hsc = cbind(annot[overlap_hsc[,1],],sum_hsc[overlap_hsc[,2],])

###############################################################################
### Load the EWAS results file and merge with functional annotation mattrices
#############################################################################
load(file="intermediate_OurResults.rda")

hsc = merge(our_res,annot_hsc[,c(5:7)],by.x=1,by.y=0,all.x=TRUE)
lymphoid = merge(our_res,annot_lymphoid[,c(7:12)],by.x=1,by.y=0,all.x=TRUE)
myeloid = merge(our_res,annot_myeloid[,c(7:12)],by.x=1,by.y=0,all.x=TRUE)

#######################################################################################################
### Within lymphoid and myeloid cell types, add columns to summarize across more granular cell types
### Column 19: sum of methylated reads across all cell types
### Column 20: sum of total (U + M) reads across all cell types
### Column 21: estimated methylation proportion in lymphoid cells (averaged over B, CD4, and CD8 cells)
### or in myeloid cells (averaged over monocytes and neutrophils).  This strategy was taken due to sparse coverage.
#######################################################################################################

lymphoid[,19] = apply(cbind(lymphoid[,c("MB","M4","M8")],myeloid$Mnk),1,sum,na.rm=T)
lymphoid[,20] = apply(cbind(lymphoid[,c("NB","N4","N8")],myeloid$Nnk),1,sum,na.rm=T)
lymphoid[,21] = lymphoid[,19]/lymphoid[,20]

myeloid[,19] = apply(myeloid[,c("Mmono","Mneu")],1,sum,na.rm=T)
myeloid[,20] = apply(myeloid[,c("Nmono","Nneu")],1,sum,na.rm=T)
myeloid[,21] = myeloid[,19]/myeloid[,20]

#######################################################################################################
### Create matrices with EWAS results (columns 2-3) from DNMT3a (d3) or TET2 (t2) EWAS, 
### and DNA methylation proportions for HSC, Lymphoid, and Myeloid cells (columns 4-6)
#######################################################################################################

d3 = cbind(hsc[,c(1,5,6,15)],lymphoid[,21],myeloid[,21])
t2 = cbind(hsc[,c(1,11,12,15)],lymphoid[,21],myeloid[,21])

colnames(d3)[4:6] = c("HSC","Lymphoid","Myeloid")
colnames(t2)[4:6] = c("HSC","Lymphoid","Myeloid")

############################################################################
###Get medians to include in text of manuscript
############################################################################
apply(d3[!(d3$RepD|t2$RepT),4:6],2,median,na.rm=T)
###     HSC  Lymphoid   Myeloid 
###.  0.7500000 0.7500000 0.6666667 

apply(d3[d3$RepD,4:6],2,median,na.rm=T)
###      HSC  Lymphoid   Myeloid 
### 0.3333333 0.6428571 0.4545455 

apply(t2[t2$RepT,4:6],2,median,na.rm=T)
###1.0000000 0.8461538 0.4000000 

#####################################################################################################
### Plot distributions of mean methylation for all CpGs vs. those associated with DNMT3A or TET2 CHIP
#####################################################################################################
library(vioplot)

colpal = c("grey","cornflowerblue","seagreen")

#####################################################################################################
###Create functions to facilitate plotting this for each of three cell types
#####################################################################################################

### "prep" function called below to enable plotting of only replicated CpGs in each category by setting others to missing
prep <- function(x) {
x[d3$RepD|t2$RepT,1]=NA
x[!d3$RepD,2]=NA
x[!t2$RepT,3]=NA
x
}

### "makepanel" function is where all the plotting happens
makepanel <- function(celltype,colnum) {
mye_toplot = d3[,c(colnum,colnum,colnum)]
names(mye_toplot) = c("\nAll CpGs \non array",expression("\n"~italic("DNMT3A")~"-\nassociated"),expression("\n"~italic("TET2")~"-\nassociated"))
mye2=prep(mye_toplot)

vioplot(mye_toplot[1:10,],plotCentre="point",col=0,ylim=c(0,1),ylab=NULL,main=celltype,rectCol="white", side="left",lineCol=0,border=0,xaxt="n")
for (i in 1:3) {
	temp=na.omit(mye2[,i])
    vioplot(temp,at=i-0.1,plotCentre="line",col=colpal[i],ylim=c(0,1),rectCol="white", side="left",add=T,xaxt="n")
text(x = 1:3, y = par("usr")[3]-.032, labels = c("All CpGs", expression(italic("DNMT3A-")), expression(italic("TET2-"))), xpd = NA, cex = 2.5)
text(x = 1:3, y = par("usr")[3]-.09, labels = c("on array", "associated", "associated"), xpd = NA, cex = 2.5)

    if (i==1) {temp=sample(temp,50000)}
    stripchart(temp, method="jitter", col=colpal[i], vertical=TRUE,pch=19, cex=.2, add=TRUE,at=i+.1)
}
}

jpeg('Figs3abc.jpg',width=1440,height=480)
par(mfcol=c(1,3),cex.main=3,cex.axis=2.5,cex.lab=2.5,xpd=NA,font.axis=2,font.lab=2,mgp=c(3,2,0),lheight=0.7)
makepanel("Myeloid cells",6)
makepanel("Lymphoid cells",5)
makepanel("HSCs",4)
dev.off()


