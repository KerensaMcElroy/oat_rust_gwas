module load tassel/5.2.63


#run_pipeline.pl -Xms10g -Xmx10g -SortGenotypeFilePlugin -inputFile /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD80.biallelic.maf_missing_filter.vcf -outputFile /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD80.biallelic.maf_missing_filter.sorted.vcf -fileType VCF > tassel.log

#run_pipeline.pl -Xms10g -Xmx10g -importGuess /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD80.biallelic.maf_missing_filter.sorted.vcf -KinshipPlugin -method Centered_IBS -endPlugin -export /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD80.biallelic.maf_missing_filter.sorted.txt -exportType SqrMatrix

#run_pipeline.pl -Xms10g -Xmx10g -importGuess /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD80.biallelic.maf_missing_filter.sorted.vcf -PrincipalComponentsPlugin -covariance true -endPlugin -export /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD80.biallelic.maf_missing_filter.sorted.pca

#run_pipeline.pl -Xmx10g -Xms10g -fork1 -importGuess /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD80.biallelic.maf_missing_filter.sorted.vcf -filterAlign -filterAlignMinFreq 0.05 -fork2 -importGuess /scratch1/mce03b/2020-02-28_GWAS/data/pheno_complete_avg_selected.txt -fork3 -importGuess /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD801.txt -combine4 -input1 -input2 -input3 -intersect -fork5 -importGuess /scratch1/mce03b/2020-02-28_GWAS/data/all_isolates_12SD80.biallelic.maf_missing_filter.sorted.txt -combine6 -input5 -input4 -mlm -mlmVarCompEst P3D -mlmCompressionLevel Optimum -mlmMaxP 0.01 -export /scratch1/mce03b/2020-02-28_GWAS/data/mlm

grep 'Pc35' /scratch1/mce03b/2020-02-28_GWAS/data/mlm7.txt > Pc35.txt
awk '$1 == "Pc35" { print $7 }' /scratch1/mce03b/2020-02-28_GWAS/data/mlm7.txt > Pc35_pvalues.txt
grep -vwE "(NaN)" Pc35_pvalues.txt > Pc35_pvalues_only.txt







