#!/bin/bash

export BIG=/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS


#genome files


#12NC29
#wget -nc -P ${BIG}/data ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/873/275/GCA_002873275.1_ASM287327v1/md5checksums.txt
#wget -nc -P ${BIG}/data ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/873/275/GCA_002873275.1_ASM287327v1/*genomic.fna.gz
#grep 'genomic.fna.gz' ${BIG}/data/md5checksums.txt > ${BIG}/data/genome_md5.txt
#
##12SD80
#wget -nc -P ${BIG}/data ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/873/125/GCA_002873125.1_ASM287312v1/md5checksums.txt
#wget -nc -P ${BIG}/data ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/873/125/GCA_002873125.1_ASM287312v1/*genomic.fna.gz
#grep 'genomic.fna.gz' ${BIG}/data/md5checksums.txt >> ${BIG}/data/genome_md5.txt
#
#cd ${BIG}/data/
#md5sum -c genome_md5.txt > integrity_check.txt 2>&1
#
#cd -
#
##annotations
#wget -nc -P ${BIG}/data https://raw.githubusercontent.com/figueroalab/Pca-genome/master/Pca_assemblies/Puccinia_coronata_avenae_12NC29.gff3
#wget -nc -P ${BIG}/data https://raw.githubusercontent.com/figueroalab/Pca-genome/master/Pca_assemblies/Puccinia_coronata_avenae_12NC29.proteins.fa
#wget -nc -P ${BIG}/data https://raw.githubusercontent.com/figueroalab/Pca-genome/master/Pca_assemblies/Puccinia_coronata_avenae_12SD80.proteins.fa
#wget -nc -P ${BIG}/data https://raw.githubusercontent.com/figueroalab/Pca-genome/master/Pca_assemblies/Puccinia_coronata_avenae_12SD80.gff3
#
#snps
wget -nc -P ${BIG}/data https://s3.msi.umn.edu/nifa_crown_rust/all_isolates_12SD80.filter.vcf
wget -nc -P ${BIG}/data https://s3.msi.umn.edu/nifa_crown_rust/all_isolates_12SD80.biallelic.maf_missing_filter.vcf
