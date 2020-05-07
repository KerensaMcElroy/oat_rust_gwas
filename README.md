# Genome wide association study for oat crown rust

## workplan

### Project background

* Review literature to understand pathosystem - done
* Meet with E. Heningsen and E. Nazareno to understand the steps that were taken to identify Avr candidate genes in the oat crown rust fungus - done

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
 - 60 isolates 2015/1990
 - 80 isolates buckthorn nursery 2016/17/18
 - 50 isolates 2017
 
* phenotypes
 - 60 isolates 2015/1990
 - 80 isolates buckthorn nursery 2016/17/18
 - 50 isolates 2017
 
* snps 12SD80
 - 60 isolates 2015/1990 
 - 80 isolates buckthorn nursery 2016/17/18
 - 50 isolates 2017

* review allele frequency plots to discard contaminated isolates

* tidy and merge all metadata and phenotype data using R script

### Read processing

* establish snakemake pipeline
 - establish bioinformatics software and settings consistant with 12SD80 pipeline (trimmomatic > bwa > filter > dedup > freebayes) (done)
 - complete Daniel Collin's snakemake training (done)
 - work through WGS snakemake tutorial
 - write and test snakemake pipeline to align small number of files to 12NC29 and call SNPs
 - run full analysis on all read files against 12NC29 reference using snakemake pipeline

### Analysis

* establish GWAS analysis script / pipeline in R
  - binary classifier e.g. GGLM or LG
  - assess model fit
  - Manhattan plots
  
* test pipeline performance on 1990/2015 12SD80 data. Compare with tassal analysis

* run pipeline on all data (12SD80, 12NC29) for all 40 differential lines to find candidate genomic regions for new Avr genes

* closely examine identified regions for candidate new Avr genes (e.g. effectorP, comparative genomics with other rusts, etc)

* explore and compare alternative GWAS approaches
- Machine learning
- reference free e.g. kmers


### Write up

- draft manuscript
- feedback
- submission


