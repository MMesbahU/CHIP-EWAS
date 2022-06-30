#!/bin/bash

# qsub -R y -l h_rt=20:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A -N EA_ewas_dnmt3a_swan /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A/Ancestry_stratified.DNMT3A.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A/Ancestry_stratified.DNMT3A.EWAS.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A/DNMT3A_EWAS.EA.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_DNMT3A.rda "EA" /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_DNMT3A_TET2.rda /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A/EA.DNMT3A_EWAS_CpG_Betas_SWAN.rda

# qsub -R y -l h_rt=20:00:00 -l h_vmem=20G -wd /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A -N AA_ewas_dnmt3a_swan /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A/Ancestry_stratified.DNMT3A.EWAS.sh /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A/Ancestry_stratified.DNMT3A.EWAS.R /medpop/esp2/mesbah/methCHIP/CHS_methylation_SWAN_2019May02.csv.gz /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A/DNMT3A_EWAS.AA.swan.rda /broad/hptmp/mesbah/DNAm/April5_2021/CHS_DNMT3A.rda "AA" /broad/hptmp/mesbah/DNAm/April5_2021/CHS_CHIP_DNMT3A_TET2.rda /broad/hptmp/mesbah/DNAm/April5_2021/ewasDNMT3A/AA.DNMT3A_EWAS_CpG_Betas_SWAN.rda


source /broad/software/scripts/useuse


use R-4.0

MyEWASscript=${1}

methylation=${2}

outFile=${3}

sample_info=${4}

Ancestry2Test=${5}

all_CHIP=${6}

out_CpGs_beta=${7}

#### Clock Time: Start
echo -e "Job started at: $(date)"

Job_START=$(date +%s)
####

Rscript ${MyEWASscript} ${methylation} ${outFile} ${sample_info} ${Ancestry2Test} ${all_CHIP} ${out_CpGs_beta}


############################# Clock Time: END
echo "Job ended at: $(date)" 

Job_END=$(date +%s)

echo $(( Job_END - Job_START)) | awk '{print "Total run time " int($1/3600)"H:"int(($1%3600)/60)"M:"int($1%60)"S"}'
##############################

