library(tidyverse)
library(qvalue)
library(ggman)

p <- scan("scripts/Pc35_pvalues_only.txt")
qobj <- qvalue_truncp(p)
max(qobj$pvalues[qobj$qvalues <= 0.05])
qobj
qobj$pvalues[qobj$qvalues <= 0.05]

p <- scan("/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/Pc63_pvalues_only.txt")
qobj <- qvalue_truncp(p)
max(qobj$pvalues[qobj$qvalues <= 0.05])
min(qobj$qvalues)

0.05/860632

gwas <- read_tsv('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/Pc63.txt', col_names = FALSE)
gwas
ggman(gwas, snp='X2', bp='X4', chrom='X3',pvalue='X7')

chroms <- gwas %>% filter(X7 < 7.8256e-05) %>%
  group_by(X3) %>%
  summarise(n()) %>%
  print(n = Inf)
