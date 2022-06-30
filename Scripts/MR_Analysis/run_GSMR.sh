#!/bin/bash

# Run: qsub -wd /broad/hptmp/mesbah/ukb_chip/ukb20k_plink/tmpdir -R y -t 1-22 -l h_vmem=20G -l h_rt=10:00:00 -pe smp 4 -binding linear:4 -N prep_ukb20k /broad/hptmp/mesbah/ukb_chip/ukb20k_plink/convertBgen.sh 4


nNodes=${1}
chr=${SGE_TASK_ID}
### 
/medpop/esp2/mesbah/tools/qcTool/bin/qctool_v2.0.7 \
	-g ukb_imp_chr${chr}_v3.bgen \
	-s ukb7089_imp_chr3_v3_s487395.sample \
	-threads ${nNodes} \
	-incl-samples ukb_sample_random_20kEUR.tsv \
	-og ukb20k_plink/ukb_imp_chr${chr}_v3 \
	-ofiletype binary_ped


## Prep input for GSMR
# ls -lhv ukb20k_plink/ukb_imp_chr*_v3.bim | awk '{print $NF}' | sed 's:.bim::g' > ukb20k_plink_file.txt

# echo -e "CAD\tMR_CAD_Jul17_2021/Outcome_CAD2017.tsv.gz" > gsmr_cad_outcome.txt


# Run GSMR



## Exposure files
# for chrom in {1..22}; do paste -d ' ' <(ls -lht /broad/hptmp/mesbah/DNAm/goDMC/mqtl_chr${chrom}/gsmr.mQTL.*.txt | awk '{print $NF}' | sed 's:gsmr.mQTL.:\t:g' | cut -f2 | sed 's:.txt::g') <(ls -lht /broad/hptmp/mesbah/DNAm/goDMC/mqtl_chr${chrom}/gsmr.mQTL.*.txt | awk '{print $NF}') > /broad/hptmp/mesbah/DNAm/goDMC/chr${chrom}_mQTL.exposures.txt; done 



# vanderHarst  2017
# for chr in {1..22}; do echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" >chr${chr}.CAD_UKBIOBANK.txt && zcat CAD_UKBIOBANK.gz | awk -v chrom=${chr} '$2==chrom{ print $11"\t"$4"\t"$5"\t"$10"\t"$6"\t"$7"\t"$8"\t"$9}' >> chr${chr}.CAD_UKBIOBANK.txt; done

# for chrom in {1..22}; do paste -d ' ' <(echo "CAD") <(ls -lht /broad/hptmp/mesbah/DNAm/CAD_GWAS/chr${chrom}.CAD_UKBIOBANK.txt | awk '{print $NF}') > /broad/hptmp/mesbah/DNAm/CAD_GWAS/chr${chrom}.CAD_vanderHarst2017.outcome.txt; done 


### Run GSMR
for chr in {1..22}
do 
	qsub -wd /broad/hptmp/mesbah/DNAm/vanderHarst2017.MR/tmpdirs -R y -b y -N chr${chr}.GSMR_CAD_allCpGs.26July2021 -l h_rt=10:00:00 -l h_vmem=20G -pe smp 2 -binding linear:2 /medpop/esp2/mesbah/tools/gcta_1.93.2beta/gcta64 --bfile /broad/hptmp/mesbah/ukb_chip/ukb20k_plink/ukb_imp_chr${chr}_v3 --gsmr-file /broad/hptmp/mesbah/DNAm/goDMC/chr${chr}_mQTL.exposures.txt /broad/hptmp/mesbah/DNAm/CAD_GWAS/chr${chr}.CAD_vanderHarst2017.outcome.txt --gsmr-direction 2 --clump-r2 0.05 --gwas-thresh 5e-8 --gsmr-snp-min 5 --diff-freq 0.2 --effect-plot --out /broad/hptmp/mesbah/DNAm/vanderHarst2017.MR/chr${chr}.cad_vanderHarst2017_mQTL.gsmr_results --threads 2

done


