#!/bin/bash

# swan: qsub -R y -l h_rt=40:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge -N comb_AAEA_ewas_AGE_swan_noCHIP /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge/Age.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge/Age.EWAS_noCHIPadj.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge/comb_AAEA.AGE_EWAS.noCHIP.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_DNMT3A_TET2.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_DNMT3A_TET2.rda /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge/AGE_EWAS_noCHIP_CpG_Betas_swan.rda

# swan: qsub -R y -l h_rt=40:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge -N comb_AAEA_ewas_AGE_swan_withCHIP /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge/Age.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge/Age.EWAS_withCHIPadj.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge/comb_AAEA.AGE_EWAS.withCHIP.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_DNMT3A_TET2.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_DNMT3A_TET2.rda /broad/hptmp/mesbah/DNAm/April5_2021/ewasAge/AGE_EWAS_withCHIP_CpG_Betas_swan.rda



source /broad/software/scripts/useuse


use R-4.0

MyEWASscript=${1}

methylation=${2}

outFile=${3}

sample_info=${4}

all_CHIP=${5}

out_CpGs_beta=${6}

#### Clock Time: Start
echo -e "Job started at: $(date)"

Job_START=$(date +%s)
####

Rscript ${MyEWASscript} ${methylation} ${outFile} ${sample_info} ${all_CHIP} ${out_CpGs_beta}


############################# Clock Time: END
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
##############################

