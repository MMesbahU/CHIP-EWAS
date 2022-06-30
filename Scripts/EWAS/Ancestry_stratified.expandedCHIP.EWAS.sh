#!/bin/bash

# AA: qsub -R y -l h_rt=20:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasVAF10 -N AA_ewas_VAF10_swan /broad/hptmp/mesbah/DNAm/April5_2021/ewasVAF10/Ancestry_stratified.expandedCHIP.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasVAF10/Ancestry_stratified.expandedCHIP.EWAS.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasVAF10/AA.VAF10_EWAS.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_VAF_DNMT3A_TET2.rda "AA"

# EA: qsub -R y -l h_rt=20:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasVAF10 -N EA_ewas_VAF10_swan /broad/hptmp/mesbah/DNAm/April5_2021/ewasVAF10/Ancestry_stratified.expandedCHIP.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/Ancestry_stratified.expandedCHIP.EWAS.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasVAF10/EA.VAF10_EWAS.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_VAF_DNMT3A_TET2.rda "EA"

## 
source /broad/software/scripts/useuse


use R-4.0

MyEWASscript=${1}

methylation=${2}

outFile=${3}

sample_info=${4}

ancestry_2_test=${5}


#### Clock Time: Start
echo -e "Job started at: $(date)"

Job_START=$(date +%s)
####
Rscript ${MyEWASscript} ${methylation} ${outFile} ${sample_info} ${ancestry_2_test}


############################# Clock Time: END
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
##############################

