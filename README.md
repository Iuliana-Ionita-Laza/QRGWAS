# Quantile regression for genome-wide association studies (QR GWAS)

## System requirements

The code has been tested with R versions 3.5.3 and 4.2.2 on the Unix-like operating systems including macOS 14 and CentOS Linux 7. 

## Install dependent packages in R

quantreg [https://cran.r-project.org/package=quantreg](https://cran.r-project.org/package=quantreg)

QRank [https://CRAN.R-project.org/package=QRank](https://CRAN.R-project.org/package=QRank)

data.table [https://CRAN.R-project.org/package=data.table](https://CRAN.R-project.org/package=data.table)

dplyr [https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)

## Example

We provide a toy dataset in [example](/example). 

Input files: 
- A [genotype matrix](example/example.genotype.tsv) that includes N sample = 503 individuals (column "IID") 
and N variant = 10 SNPs (columns starting with "1-555") from the [1000 Genome Phase 3](https://www.internationalgenome.org/category/phase-3/)
- A [phenotype table](example/example.phenotype.tsv) that includes one simulated quantitative trait (column "PHENOTYPE") and one simulated covariate (column "COVAR1") for the N sample = 503 individuals (column "IID").

Run the R script [example.qr.R](example.qr.R) to perform QR GWAS (single variant tests) for the phenotype and genotype data in [example](/example). 

The [expected output](example/example.sumstat.tsv) is a tab-delimited text file with the QR summary statistics. 
- Column "ID": variant ID from the genotype data.
- Column "P_QR": the integrated QR p-value across multiple quantile levels.
- Columns from "P_QR0.1" to "P_QR0.9": quantile-specific QR p-value for the quantile levels 0.1, 0.2, ..., 0.9.
- Columns from "BETA_QR0.1" to "BETA_QR0.9": quantile-specific QR effect size for the quantile levels 0.1, 0.2, ..., 0.9.

## Performance for large data

QR GWAS analysis across nine quantile levels for a dataset with 325K samples and 8.5M variants typically requires approximately 1,946 CPU hours on a single CPU core. To speed up the process, it is highly recommended to divide the whole genome genotype data into scattered segments of variants. For instance, by splitting the genome-wide data into 1,000 scattered chunks, the analysis of 325K samples and 8.5K variants can be completed in about 1.9 hours, on avearge approximately 0.2 hours per quantile level.
