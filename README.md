# Genome wide association study for oat crown rust

## workplan

### Project background

* Review literature to understand pathosystem (done)
* Meet with E. Heningsen and E. Nazareno to understand the steps that were taken to identify Avr candidate genes in the oat crown rust fungus (done)

### Data collation

* reference genomes
  - 12SD80 and 12NC29 from NCBI, check integrity with md5sum (done)
  - corresponding annotations from https://raw.githubusercontent.com/figueroalab/Pca-genome/master/Pca_assemblies/ (checksums not available) (done)
  - download script in github (done)

* sequence reads
  - 60 isolates 2015/1990 raw
  - 60 isolates 2015/1990 aligned
  - 80 isolates buckthorn nursery 2016/17/18 raw
  - 50 isolates 2017 raw
 
* metadata
  - 60 isolates 2015/1990 (data/TableS1_final.xlsx) (done)
  - 80 isolates buckthorn nursery 2016/17/18 (data/Figueroa_Project_006_NovaSeq_Summary.xlsx) (done)
  - 50 isolates 2017 (data/Figueroa_Project_006_NovaSeq_Summary.xlsx) (done)
 
* phenotypes
  - 60 isolates 2015/1990 (data/1990_2015_phenotypes.txt) (done)
  - 80 isolates buckthorn nursery 2016/17/18 (data/2016_Pca_isolates_phenotype_completeset.xlsx, data/2017-2018_BT_isolates.xlsx) 
      * 2016 (done)
      * 2017/18 requires importing and averaging 9 point scale replicate results (due 08.05.2020 - done)
  - 50 isolates 2017 (data/Copy of OCR2017Survey.xlsx) (only subset sequenced? And only one replicate? Check with M) (due 08.05.2020 - done, two missing)
 
* snps 12SD80
  - 60 isolates 2015/1990 (subset of data/all_isolates_12SD80.filter.vcf, data/all_isolates_12SD80.biallelic.maf_missing_filter.vcf)
  - 80 isolates buckthorn nursery 2016/17/18 (data/Buckthorn_variants.vcf.tar.gz)
  - 50 isolates 2017 (subset of data/SA_SOU_NOR_variants.vcf.tar.gz)
  - request and check md5sums (due 08.05.2020 - done)

* review allele frequency plots to discard contaminated isolates (due 13.05.2020)

* tidy and merge all metadata and phenotype data using R script (80% complete - due 13.05.2020)

### Read processing

* establish snakemake pipeline
  - establish bioinformatics software and settings consistant with 12SD80 pipeline (trimmomatic > bwa > filter > dedup > freebayes) (done)
  - complete Daniel Collin's snakemake training (done)
  - work through WGS snakemake tutorial (due 15.05.2020 - done)
  - write and test snakemake pipeline to align small number of files to 12NC29 and call SNPs (due 22.05.2020 - done)
  - run full analysis on all read files against 12NC29 reference using snakemake pipeline (due 05.06.2020)

### Analysis

* establish GWAS analysis script / pipeline in R (due 05.06.2020)
  - binary classifier e.g. GGLM or LG
  - assess model fit
  - Manhattan plots
  
* test pipeline performance on 1990/2015 12SD80 data. Compare with tassal analysis (due 10.06.2020)

* run pipeline on all data (12SD80, 12NC29) for all 40 differential lines to find candidate genomic regions for new Avr genes (due 19.06.2020)

* closely examine identified regions for candidate new Avr genes (e.g. effectorP, comparative genomics with other rusts, etc) (due 10.07.2020)

* explore and compare alternative GWAS approaches (due 19.09.2020)
  - Machine learning
  - reference free e.g. kmers


### Write up

* draft manuscript
  - create google doc for manuscript (due 08.05.2020)
  - update methods (ongoing - due 10.07.2020)
  - prepare figures (ongoing - due 10.07.2020)
  - reference documents as read / add to introduction (ongoing - due 17.07.2020)
  - write analysis / discussion / conclusion (due 31.07.2020)
* feedback (due 07.08.2020)
* submission (due 14.08.2020)


