# WTSP_VIPMethylation

This respository contains the scripts and data to reproduce the analyses from Prichard et al., 2022 (under review).

This data analysis can be reproduced using the MeanMethyaltion_byallele.csv, adult-dfcg.csv, and chick-dfcg.csv files with the R scripts included here.

The WTSP_VIPMethylation.html file includes the entire analyses for both age groups starting from the raw sequencing data. It requires the following files:
  * **cpg_manifest.txt** A list of each CpG in the sequenced region. This includes all of the polymorphic CpGs and a factor distinguishing the CpGs that are shared between both alleles, or which allele it is present on.
  * **files_manifest.txt** A list of all of the bilsulfite converted DNA sequencing file names (matches the file names in the link below).
  * **files_rna.txt** A list of all the RNA sequencing files (does not match the file names in the link below).
  * **HypNormalizationallele.txt** The normalization factors for the RNA sequencing data - by individual / sequencing run. 

Raw bisulfite-converted DNA sequencing data is available at: 
Raw RNA sequencing data is available at: https://www-ncbi-nlm-nih-gov.proxy.library.emory.edu/bioproject/657006

