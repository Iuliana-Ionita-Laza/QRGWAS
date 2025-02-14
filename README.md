# Quantile regression for genome-wide association studies (QR GWAS)

This repository provides code and instructions for conducting QR GWAS for quantitative traits as described in our paper ["Genome-wide discovery for biomarkers using quantile regression at biobank scale"](https://www.nature.com/articles/s41467-024-50726-x).



## Install dependent packages in R

quantreg [https://cran.r-project.org/package=quantreg](https://cran.r-project.org/package=quantreg)

QRank [https://CRAN.R-project.org/package=QRank](https://CRAN.R-project.org/package=QRank)

data.table [https://CRAN.R-project.org/package=data.table](https://CRAN.R-project.org/package=data.table)

dplyr [https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)



## Example

We provide a toy dataset in [example](/example). 
The input data includes 
1) the [genotype matrix](example/example.genotype.tsv) that contains N sample = 503 individuals (column "IID") 
and N variant = 10 SNPs (columns starting with "1-555") from the [1000 Genome Phase 3](https://www.internationalgenome.org/category/phase-3/).
2) the [phenotype table](example/example.phenotype.tsv) that contains one simulated quantitative trait (column "PHENOTYPE") and two simulated covariates (columns "COVAR1" and "COVAR2") for the N sample = 503 individuals (column "IID").

Run the R script [example.qr.R](example.qr.R) to perform QR GWAS (single variant tests) for the phenotype and genotype data in [example](/example). Set ```is.effect.estimated = T``` to enable the estimation of quantile-specific effect size (disabled by default).

The [expected output](example/example.sumstat.tsv) is a tab-delimited text file with the QR summary statistics. 
- Column "ID": variant ID from the genotype data.
- Column "P_QR": the integrated QR p-value across multiple quantile levels.
- Columns from "P_QR0.1" to "P_QR0.9": quantile-specific QR p-value for the quantile levels 0.1, 0.2, ..., 0.9.
- Columns from "BETA_QR0.1" to "BETA_QR0.9": quantile-specific QR effect size for the quantile levels 0.1, 0.2, ..., 0.9.



## Suggestions for genome-wide analysis

We recommend splitting genome-wide genotype data into smaller variant subsets or small genomic regions to improve computational efficiency. Data scattering reduces the computational burden on a single machine and enables parallel execution across multiple machines, significantly enhancing scalability. 

For example, in a dataset with 325,306 individuals and 8,572,925 variants, performing QR GWAS for a single quantile level took 216 CPU hours using R. To optimize this process, we split the genome into 1,000 segments, each containing ~8.6K variants. Each segment required ~1.9 CPU hours to complete QR testing across 9 quantile levels when executed on a single CPU core. By leveraging data scatter and parallel computing, we distributed the workload across 1,000 CPU cores. This reduced the wall-clock time to ~1.9 hours for 9 quantile levels, making the QR GWAS feasible for biobank-scale data.


Users may consider this data partitioning strategry and determine the number of genomic chunks based on data volume and available computing resources. The R code for QR GWAS requires loading a genotype matrix in memory. Loading whole-genome or chromosome-wide genotype data into R consumes significant memory and processing time. It would be much easier to load a smaller chunk of genotype data into memory and then perform QR on the small subsets. 

Users can directly import the commonly used PLINK/VCF format data into R using packages such as [genio](https://github.com/OchoaLab/genio) and [seqminer](https://github.com/zhanxw/seqminer). Below, we demonstrate how to convert a genotype data chunk into an R-friendly format using [PLINK2](https://www.cog-genomics.org/plink/2.0/).


Assume the genotype data for chromosome 1 is stored in a PLINK binary fileset: genotype_chr1.fam, genotype_chr1.bim, and genotype_chr1.bed. If we aim to perform QR GWAS on a 3 Mb genomic region, we can extract the data as follows.

```bash
plink2 \
--bfile genotype_chr1 \
--chr 1 \
--from-bp 1 \
--to-bp 3000000
--export A \
--out genotype_scatter_chr1_1_3000000
```

This generates a plain text file, "genotype_scatter_chr1_1_3000000.raw", which follows the same format as the [example genotype matrix](example/example.genotype.tsv) we provided. Such file can be loaded into R using the code from [example.qr.R](example.qr.R) (replacing example.genotype.tsv with genotype_scatter_chr1_1_3000000.raw). 

By default, PLINK2 exports the additive coding of REF allele to the resulting ".raw" file. If needed, the allele coding can be modified using the [--export-allele](https://www.cog-genomics.org/plink/2.0/data#export) option.



## Citation

Chen Wang, Tianying Wang, Krzysztof Kiryluk, Ying Wei, Hugues Aschard, and Iuliana Ionita-Laza. [Genome-wide discovery for biomarkers using quantile regression at biobank scale](https://www.nature.com/articles/s41467-024-50726-x). Nature Communications, 15(1), 6460.


