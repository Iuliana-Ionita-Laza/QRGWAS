# example of QR GWAS
library(data.table)
library(dplyr)
library(QRank)
library(quantreg)
options(scipen = 999)

## ACAT (Cauchy combination)
## reference => https://pubmed.ncbi.nlm.nih.gov/33012899/
## source code => https://github.com/yaowuliu/ACAT
acat <- function(pvals) {
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

## QRank
## reference => https://pubmed.ncbi.nlm.nih.gov/28334222/
## source code => https://CRAN.R-project.org/package=QRank
## Quantile regression
## reference => https://doi.org/10.1201/9781315120256
## source code => https://cran.r-project.org/package=quantreg
qrank.test <- function(phe, cov, snp, fit.residual, tau = c(0.25, 0.5, 0.75)) {
  ltau = length(tau)
  x = as.matrix(snp)
  y = as.matrix(phe)
  zz = cbind(rep(1, nrow(y)), as.matrix(cov))
  
  VN = matrix(0, nrow = ltau, ncol = ltau) 
  for(i in 1:ltau) {
    for(j in 1:ltau) {
      VN[i,j] = min(tau[i], tau[j]) - tau[i] * tau[j]     
    }
  }
  
  xstar = lm(x ~ zz - 1)$residual
  
  SN = NULL
  for(i in 1:ltau) {
    dual = 1 * ((fit.residual %>% filter(Quantile == tau[i]) %>% pull(Residual)) > 0)
    ranks = dual - (1 - tau[i])
    Sn = as.matrix(t(xstar) %*% (ranks))
    SN = c(SN,Sn)
  } 
  
  VN2 = matrix(outer(VN,t(xstar)%*% xstar,"*"), nrow = ltau)
  pvalue = pchisq(SN^2/diag(VN2), 1, lower.tail = F)
    
  # e = solve(chol(VN2))
  # SN2 = t(e) %*% SN
  # pvalue.composite = pchisq(sum(SN2^2), ltau, lower.tail = F)
  pvalue.cauchy = acat(pvalue)
  names(pvalue) = tau
  
  result = list(tau = tau, quantile.pvalue = pvalue, cauchy.pvalue = pvalue.cauchy)
  return(result)
}



####################################################################################################



## quantile levels
## 0 < values < 1
qntl = (1:9) / 10


## phenotype table
## row: individual
## column: phenotypic variable
pheno = fread(file = paste0("example/example.phenotype.tsv"),
              header = T, sep = "\t", data.table = F, stringsAsFactors = F)
#### one phenotype
# pheno %>% select(PHENOTYPE)
#### one covariate
# pheno %>% select(COVAR1)



## genotype matrix
## row: individual
## column: variant
geno = fread(file = paste0("example/example.genotype.tsv"),
             header = T, sep = "\t", data.table = F, stringsAsFactors = F)
geno.mat = as.matrix(geno %>% select(-(FID:PHENOTYPE)))



## QR single variant test
#### null model fitting
set.seed(seed = 256)
df.res = data.frame()
fml.null = paste0("PHENOTYPE ~ COVAR1")
for (idx.qntl in 1:length(qntl)) {
  fit.null = rq(data = pheno %>% select(PHENOTYPE, COVAR1), 
                formula = as.formula(fml.null), 
                tau = qntl[idx.qntl], 
                method = "pfn")
  df.res = rbind(df.res,
                 data.frame(Quantile = qntl[idx.qntl],
                            pheno %>% select(IID),
                            Residual = unname(as.matrix(pheno %>% select(PHENOTYPE)) - as.matrix(pheno %>% select(COVAR1) %>% mutate(1, .before = 1)) %*% as.matrix(fit.null$coefficients))))
}

#### score test
is.effect.estimated = F # whether estimation of quantile-specific effect size is required
set.seed(seed = 256)
df.qr = data.frame()
for (idx.geno in 1:ncol(geno.mat)) {
  df.qreg = data.frame(pheno %>% select(PHENOTYPE, COVAR1),
                       GENOTYPE = geno.mat[, idx.geno],
                       check.names = F, stringsAsFactors = F)
  
  beta = double()
  if (is.effect.estimated) {
    fml.qreg = paste0("PHENOTYPE ~ COVAR1 + GENOTYPE")
    for (idx.qntl in 1:length(qntl)) {
      qreg = rq(data = df.qreg, formula = as.formula(fml.qreg), tau = qntl[idx.qntl], method = "pfn")
      beta = c(beta, qreg$coefficients["GENOTYPE"])
    }
  }

  p = qrank.test(phe = df.qreg %>% select(PHENOTYPE), 
                 cov = df.qreg %>% select(COVAR1), 
                 snp = df.qreg %>% select(GENOTYPE), 
                 fit.residual = df.res,
                 tau = qntl)
  df.qr = rbind(df.qr,
                data.frame(p$cauchy.pvalue,
                           t(p$quantile.pvalue),
                           t(beta)))
}

colnames(df.qr) = c("P_QR",
                    paste0("P_QR", qntl),
                    ifelse(is.effect.estimated, paste0("BETA_QR", qntl), NA)) %>% 
  head(n = 1 + length(qntl) + length(beta))
df.qr = df.qr %>% mutate(ID = colnames(geno.mat), .before = 1)



## QR summary statistics
## row: variant
## column: p-value # and beta
fwrite(x = df.qr, file = paste0("example/example.sumstat.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)


