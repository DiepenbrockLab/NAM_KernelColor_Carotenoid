# NAM_KernelColor_Carotenoid
Scripts for LaPorte et al. (2021)

A short description of the included scripts:
- 1b: Calculates adjusted P-values with the multtest library using the imputed markers (GBS SNPs), transformed BLUES, and output from the joint linkage analysis.
- 2a: Consolidates GWAS results into a single file.
- 3a: Creates tabular summary using 1b results, imputed markers, and BLUEs
- 3b: Finalizes tabular summary, including support interval information.
- 3c: Updates tabular summary of GWAS results
- RandomForest10Fam: this script includes the code for the Random Forest, which predicts kernel color based on carotenoid values. Included also in this script is the Linear Model run on the same data, as well as the parameter optimization (grid search). 
