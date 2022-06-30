#!/bin/bash

# qsub -R y -l h_rt=20:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasTET2 -N EA_ewas_tet2_swan /broad/hptmp/mesbah/DNAm/April5_2021/ewasTET2/Ancestry_stratified.TET2.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasTET2/Ancestry_stratified.TET2.EWAS.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasTET2/TET2_EWAS.EA.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_TET2.rda EA
# qsub -R y -l h_rt=20:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasTET2 -N AA_ewas_tet2_swan /broad/hptmp/mesbah/DNAm/April5_2021/ewasTET2/Ancestry_stratified.TET2.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasTET2/Ancestry_stratified.TET2.EWAS.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasTET2/TET2_EWAS.AA.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_TET2.rda AA

source /broad/software/scripts/useuse

use R-4.0
###
MyEWASscript=${1}

methylation=${2}

outFile=${3}

sample_info=${4}

Ancestry2Test=${5}
#### Clock Time: Start
echo -e "Job started at: $(date)"

Job_START=$(date +%s)
####
Rscript ${MyEWASscript} ${methylation} ${outFile} ${sample_info} ${Ancestry2Test}
############################# Clock Time: END
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
##############################

