# Clonal hematopoiesis of indeterminate potential, DNA methylation, and risk for coronary artery disease
Code used to generate figures and tables in Uddin, M. M. et al. 2022 [Clonal hematopoiesis of indeterminate potential, DNA methylation, and risk for coronary artery disease.](https://www.nature.com/articles/s41467-022-33093-3)

## CHIP calling
CHIP calling was done using the pipeline described in [Bick,  et  al. 2021](https://www.nature.com/articles/s41586-020-2819-2)  (Terra workflow:https://app.terra.bio/#workspaces/terra-outreach/CHIP-Detection-Mutect2)

## EWAS
EWAS was performed using **CpGassoc** R package [Scripts](https://github.com/MMesbahU/CHIP-EWAS/tree/main/Scripts/EWAS)

## Meta-Analysis
Meta-analysis was performed using **METAL** software [Scripts](https://github.com/MMesbahU/CHIP-EWAS/tree/main/Scripts/meta_analysis)

## Mendelian randomization (MR)
Mendelian randomization was performed using GSMR from GCTA software [Scripts](https://github.com/MMesbahU/CHIP-EWAS/tree/main/Scripts/MR_Analysis)

## Enrichment analyses
Enrichment and functional annotation analyses were performed using custom R scripts (coming soon)

## Citation
* Uddin, M. M., N. Q. Nguyen, B. Yu, J. Brody, A. Pampana, T. Nakao, M. Fornage, J. Bressler, N. Sotoodehnia, J. Weinstock, M. Honigberg, D. Nachun, R. Bhattacharya, G. Griffin, V. Chander, R. Gibbs, J. Rotter, C. Liu, A. Baccarelli, D. Chasman, E. Whitsel, D. Kiel, J. Murabito, E. Boerwinkle, B. Ebert, S. Jaiswal, J. Floyd, A. Bick, C. Ballantyne, B. Psaty, P. Natarajan and K. Conneely (2022). ["Clonal hematopoiesis of indeterminate potential, DNA methylation, and risk for coronary artery disease."](https://www.nature.com/articles/s41467-022-33093-3) *Nat Commun.*
