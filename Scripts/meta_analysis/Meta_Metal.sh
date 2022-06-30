#!/bin/bash

source /broad/software/scripts/useuse

use Tabix


# qsub -R y -wd /medpop/esp2/mesbah/CHIP_GWAS/metaAnalysis -pe smp 1 -binding linear:1 -l h_vmem=30G -l h_rt=20:00:00 -N MetaMetal_Oct20_2020 /medpop/esp2/mesbah/CHIP_GWAS/metaAnalysis/Meta_Metal.sh /medpop/esp2/mesbah/CHIP_GWAS/metaAnalysis/directive.METAL.txt /medpop/esp2/mesbah/CHIP_GWAS/metaAnalysis/MetaMetal_logs

metal_directive=${1}
temp_log_file=${2}
metal_out=${3} 

Metal=/medpop/esp2/mesbah/tools/METAL/build/bin/metal

${Metal} ${metal_directive} 1>>${temp_log_file}.log 2>>${temp_log_file}.err

bgzip ${metal_out}


