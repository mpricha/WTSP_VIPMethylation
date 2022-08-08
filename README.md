# WTSP_VIPMethylation

This respository contains the scripts and data to reproduce the analyses from Prichard et al., 2022 (under review).

This data analysis can be reproduced using the R script and files included here. The script and analyses were conducted in R Studio (RStudio 2022.07.1+554 "Spotted Wakerobin"). The script skips the step of reading in the individual read counts that are the output from bismark. The code for that step is included but not intended to be run on the files in this respository. An R Markdown draft of the script will be available soon.

Raw bisulfite-converted DNA sequencing data will be available on Dryad at this link after peer review: https://doi.org/10.5061/dryad.547d7wm9s  
Raw RNA sequencing data is available at: https://www-ncbi-nlm-nih-gov.proxy.library.emory.edu/bioproject/657006

The FASTA DNA file and the QC-Mapping-CallingMethylation.txt files are available to help with reproducing the quality control, read mapping, allele-assignment, and methylation count steps leading up to the .cov.gz outputs from Bismark (i.e., L48-61 of the R script). 
