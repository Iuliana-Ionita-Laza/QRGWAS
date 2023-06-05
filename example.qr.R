# example of QR GWAS
library(data.table)
library(tidyverse)
library(QRank)
library(quantreg)
options(scipen = 999)

## QRANK
## https://pubmed.ncbi.nlm.nih.gov/28334222/
## https://CRAN.R-project.org/package=QRank
## Quantile regression
## https://doi.org/10.1201/9781315120256
## https://cran.r-project.org/package=quantreg
QRank <- function(gene, snp, cov = NULL, tau = c(0.25, 0.5, 0.75)) {
  ltau = length(tau)
  x = as.matrix(snp)
  y = as.matrix(gene)
  cov = as.matrix(cov)
  zz = cbind(rep(1, nrow(y)), cov)
  
  VN = matrix(0, nrow = ltau, ncol = ltau) 
  for(i in 1:ltau) {
    for(j in 1:ltau) {
      VN[i,j] = min(tau[i], tau[j]) - tau[i] * tau[j]     
    }
  }
  xstar = lm(x~zz-1)$residual
  SN = NULL
  for(i in 1:ltau) {
    qreg = rq.fit.pfn(zz, y, tau = tau[i])
    qreg$residual = y - zz %*% as.matrix(qreg$coefficients)
    qreg$dual = 1 * (qreg$residual > 0)
    ranks = qreg$dual - (1 - tau[i])
    Sn = as.matrix(t(xstar) %*% (ranks))
    SN = c(SN,Sn)
  } 
  
  VN2 = matrix(outer(VN,t(xstar)%*% xstar,"*"), nrow = ltau)
  pvalue1 = pchisq(SN^2/diag(VN2), 1, lower.tail = F)
  names(pvalue1) = tau
  
  e = solve(chol(VN2))
  SN2 = t(e) %*% SN
  pvalue = pchisq(sum(SN2^2), ltau, lower.tail = F)
  
  result = list(composite.pvalue = pvalue, quantile.specific.pvalue = pvalue1, tau = tau)
  return(result)
}

## Cauchy combination
## https://pubmed.ncbi.nlm.nih.gov/33012899/
## https://github.com/yaowuliu/ACAT
cauchy.meta <- function(pvals) {
  # Check input
  pvals = pvals[!is.na(pvals)]
  if (length(pvals) == 0) {
    return(NA)
  }
  # pvals[pvals == 0] = 2.2E-308
  # Convert to Cauchy
  cauchy = 1 / (pvals * pi)
  cauchy[pvals >= 1e-15] = tanpi(0.5 - pvals[pvals >= 1e-15])
  stats = mean(cauchy)
  p = pcauchy(q = stats, lower.tail = F)
  return(p)
}



####################################################################################################



## quantile levels
## 0 < values < 1
qntl = (1:9) / 10


## phenotype table
## row: individual
## column: phenotypic variable
pheno = fread(file = paste0("example/example.phenotype.tsv"),
              header = T, sep = "\t", data.table = F, stringsAsFactors = FALSE)
#### one phenotype
# pheno %>% select(PHENOTYPE)
#### one covariate
# pheno %>% select(COVAR1)



## genotype matrix
## row: individual
## column: variant
geno = fread(file = paste0("example/example.genotype.tsv"),
             header = T, sep = "\t", data.table = F, stringsAsFactors = FALSE)
geno.mat = as.matrix(geno %>% select(-(FID:PHENOTYPE)))



## QR single variant test
set.seed(256)
df.p = data.frame()
for (i in 1:ncol(geno.mat)) {
  p.qr = QRank(gene = pheno %>% select(PHENOTYPE), 
               cov = pheno %>% select(COVAR1), 
               snp = geno.mat[, i], 
               tau = qntl)
  df.p = rbind(df.p, data.frame(p.qr$quantile.specific.pvalue %>% t()))
}
colnames(df.p) = paste0("P QR ", qntl)
df.p = df.p %>% mutate(`P QR` = apply(X = df.p, MARGIN = 1, FUN = cauchy.meta), .before = 1)
df.p = df.p %>% mutate(ID = colnames(geno.mat), .before = 1)



## QR summary statistics
## row: variant
## column: p-value
fwrite(x = df.p, file = paste0("example/example.sumstat.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)


