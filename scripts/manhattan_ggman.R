library(tidyverse)
library(qvalue)

p <- scan("scripts/Pc35_pvalues_only.txt")
qobj <- qvalue_truncp(p)
max(qobj$pvalues[qobj$qvalues <= 0.05])
qobj
qobj$pvalues[qobj$qvalues <= 0.05]

p <- scan("scripts/Pc63_pvalues_only.txt")
qobj <- qvalue_truncp(p)
max(qobj$pvalues[qobj$qvalues <= 0.05])
min(qobj$qvalues)
