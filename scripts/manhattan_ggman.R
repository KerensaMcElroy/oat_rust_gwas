library(tidyverse)
library(qvalue)
library(ggman)

pc_stats <- read_tsv(file='analysis/tassel/12NC29/interim/binary/diff_0/mlm_nocomp2.txt')

p_vals <- pc_stats %>% 
  filter(!is.na(p) & !is.nan(p)) %>%
  select(p)

oat_diff <- pc_stats$Trait[1]

#p <- scan("scripts/Pc35_pvalues_only.txt")
qobj <- qvalue_truncp(p)
max(qobj$pvalues[qobj$qvalues <= 0.05])
qobj
qobj$pvalues[qobj$qvalues <= 0.05]

p <- scan("scripts/Pc63_pvalues_only.txt")
qobj <- qvalue_truncp(p)
max(qobj$pvalues[qobj$qvalues <= 0.05])
min(qobj$qvalues)

0.05/860632

gwas <- read_tsv('scripts/Pc63.txt', col_names = FALSE)
gwas
ggman(gwas, snp='X2', bp='X4', chrom='X3',pvalue='X7')
