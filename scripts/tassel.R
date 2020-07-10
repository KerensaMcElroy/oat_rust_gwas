options(java.parameters = c("-Xmx10g", "-Xms10g"))

library(rTASSEL)
library(tidyverse)

out <- '/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/'
startLogger(fullPath = paste0(out, 'logs'), fileName = 'rtassel.log')
tasGenoVCF <- readGenotypeTableFromPath(
  path = paste0(out, 'data/all_isolates_12SD80.biallelic.maf_missing_filter.vcf'),
  sortPositions = TRUE
)

tasPheno <- readPhenotypeFromPath(
  path = paste0(out, 'data/pheno_complete_avg_selected.txt')
)

tasGenoPheno <- readGenotypePhenotype(
  genoPathOrObj = tasGenoVCF,
  phenoPathDFOrObj = tasPheno
)

tasGenoPheno

tasKin <- kinshipMatrix(tasObj = tasGenoPheno)

tasMLM <- assocModelFitter(tasObj = tasGenoPheno,             # <- our prior TASSEL object
                           formula = Pc35 ~ .,               # <- run only EarHT
                           fitMarkers = TRUE,                 # <- set this to TRUE for GLM
                           kinship = tasKin,                  # <- our prior kinship object
                           fastAssociation = FALSE
)

plot <- manhattanPlot(tasMLM$MLM_Stats, 'Pc35', 0.05, showSigMarkers = TRUE)
?manhattanPlot
manhattanPlot
