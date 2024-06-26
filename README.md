# Quantile regression for genome-wide association studies (QR GWAS)



## Install dependent packages in R
quantreg [https://cran.r-project.org/package=quantreg](https://cran.r-project.org/package=quantreg)

QRank [https://CRAN.R-project.org/package=QRank](https://CRAN.R-project.org/package=QRank)

data.table [https://CRAN.R-project.org/package=data.table](https://CRAN.R-project.org/package=data.table)

dplyr [https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)



## Example
We provide a toy dataset in [example](/example). 
The input data are 1) the [genotype matrix](example/example.genotype.tsv) that includes N sample = 503 individuals (column "IID") 
and N variant = 10 SNPs (columns starting with "1-555") from the [1000 Genome Phase 3](https://www.internationalgenome.org/category/phase-3/)
and 2) the [phenotype table](example/example.phenotype.tsv) that includes one simulated quantitative trait (column "PHENOTYPE") and two simulated covariates (columns "COVAR1" and "COVAR2") for the N sample = 503 individuals (column "IID").

Run the R script [example.qr.R](example.qr.R) to perform QR GWAS (single variant tests) for the phenotype and genotype data in [example](/example). 

The [expected output](example/example.sumstat.tsv) is a tab-delimited text file with the QR summary statistics. 
- Column "ID": variant ID from the genotype data.
- Column "P_QR": the integrated QR p-value across multiple quantile levels.
- Columns from "P_QR0.1" to "P_QR0.9": quantile-specific QR p-value for the quantile levels 0.1, 0.2, ..., 0.9.
- Columns from "BETA_QR0.1" to "BETA_QR0.9": quantile-specific QR effect size for the quantile levels 0.1, 0.2, ..., 0.9.

