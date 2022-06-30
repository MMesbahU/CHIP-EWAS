### Rscript /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/format_4Metal.R /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/AA.CHIP_EWAS.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/AA_CHIP_EWAS.SWAN.tsv

## Rscript /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/format_4Metal.R /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/EA.CHIP_EWAS.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/EA_CHIP_EWAS.SWAN.tsv

############
ARGs <- commandArgs(TRUE)

ewas_input <- ARGs[1]

output_filename <- ARGs[2]

#################
# Load Libraries  
# library(data.table)

# library(CpGassoc)
loadRData <- function(fileName){
#loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
}

##
d <- loadRData(fileName=ewas_input)

res <- d$results

res$adj.intercept <- d$coefficients$adj.intercept

res$effect.size <- d$coefficients$effect.size

res$std.error <- d$coefficients$std.error

res$Ref <- "A"

res$Alt <- "T"

res <- subset(res, !is.na(res$T.statistic) )
##
write.table(res, output_filename, row.name=F, quote = F, sep="\t")


