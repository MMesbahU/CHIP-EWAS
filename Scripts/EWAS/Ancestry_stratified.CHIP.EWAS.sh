#!/bin/bash

### Modified:
# AA: qsub -R y -l h_rt=20:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP -N AA_ewas_CHIP_swan_v2 /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/Ancestry_stratified.CHIP.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/Ancestry_stratified.CHIP.EWAS.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/AA.CHIP_EWAS.swan.v2.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_DNMT3A_TET2.rda "AA"
# EA: qsub -R y -l h_rt=20:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP -N EA_ewas_CHIP_swan_v2 /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/Ancestry_stratified.CHIP.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/Ancestry_stratified.CHIP.EWAS.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasCHIP/EA.CHIP_EWAS.swan.v2.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_DNMT3A_TET2.rda "EA"


## 
source /broad/software/scripts/useuse


use R-4.0

MyEWASscript=${1}

methylation=${2}

outFile=${3}

sample_info=${4}

ancestry_2_test=${5}

# all_CHIP=${6}

# out_CpGs_beta=${7}

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

