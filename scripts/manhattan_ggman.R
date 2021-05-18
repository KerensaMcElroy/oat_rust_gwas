library(tidyverse)
library(ggman)
library(qqman)

#funtion for producing tables of significant peaks
sig_table <- function(genome, oat_diff, pval, nchrom, chroms, ext){
  pc_stats <- read_tsv(file=paste(oat_diff, '/mlm_nocomp2.txt', sep='')) 
  #print(pc_stats)
  for_table <- pc_stats %>%
      select(c(Marker, Pos, Chr, p)) %>%    
      left_join(chroms) %>%
      filter(!is.na(p) & !is.nan(p) & p < pval) %>%
      group_by(Chr) %>%
      filter(n()>nchrom) %>%
      ungroup()%>%
      select(c(Marker, Chromosome, Pos, p)) %>%
      print()
  diff <- pc_stats$Trait[1]
  write_tsv(for_table,path = paste('/datasets/work/af-crown-rust/work/2020-02-28_GWAS/results/tables/', genome, '_', diff, ext, sep=''), col_names = TRUE, append=FALSE)
}

#function for producing manhattan plots
my_manhattan <- function(in_file, chroms, genome){
  pc_stats <- read_tsv(file=in_file) %>%
    left_join(chroms) %>%
    filter(!is.na(p) & !is.nan(p)) %>%
    rename('snp' = 'Marker') %>%
    rename('bp' = 'Pos') %>%
    rename('chrom' = 'Chromosome') %>%
    rename('pvalue' = 'p') %>%
    mutate(snp = paste(str_extract(chrom, "[123456789]+.*F"),bp,sep='_'))
  
  to_plot <- pc_stats %>%
    select(c(snp, chrom, bp, pvalue))
  
  to_plot <- as.data.frame(to_plot)
  
  oat_diff <- pc_stats$Trait[1]
  
  plot <- ggman(to_plot, title=oat_diff, sigLine = -log10(6.5747e-08), ylim=c(0,30)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  snps_label <- to_plot %>%
    filter(pvalue < 6.5747e-08) %>%
    group_by(chrom) %>%
    filter(n()>2) %>%
    filter(pvalue==min(pvalue)) %>%
    filter(row_number()==1) %>%
    mutate(label = str_extract(chrom, "[123456789]+.*F")) %>%
    ungroup()
  
  snps_label <- as.data.frame(snps_label)
  
  if (nrow(snps_label) > 0) {
    snps_highlight <- to_plot %>%
      filter(chrom %in% snps_label$chrom)
    plot <- ggmanLabel(plot, labelDfm = snps_label, snp = "snp", label = "label", type = 'text')
    plot <- ggmanHighlight(plot, highlight = snps_highlight$snp, size=.1)
  } else {
    plot <- plot + scale_colour_grey(start = 0.5,end = 0.6)
  }
  ggsave(plot, width=6, height=4, device = 'jpeg', filename = paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/', genome, '_nine_man_', oat_diff,'.jpeg', sep=''))
  jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/', genome, '_nine_qq_', oat_diff,'.jpeg', sep=''), width = 3, height = 3, units = 'in', res = 300)
  qplt <- qq(to_plot$pvalue, main=oat_diff, cex=.2)
  dev.off()
}

#list the output and location of tassel results
folders_nine_12SD80 <- list.files(path='/datasets/work/af-crown-rust/work/2020-02-28_GWAS/analysis/freebayes/12SD80/nine', full.names=TRUE)
folders_nine_12NC29 <- list.files(path='/datasets/work/af-crown-rust/work/2020-02-28_GWAS/analysis/freebayes/12NC29/interim/nine', full.names=TRUE)
folders_bin_12NC29 <- list.files(path='/datasets/work/af-crown-rust/work/2020-02-28_GWAS/analysis/freebayes/12NC29/interim/binary', full.names=TRUE)

#contig names for genomes
chroms_12SD80 <- read_delim(file='/datasets/work/af-crown-rust/work/2020-02-28_GWAS/data/chroms_12SD80.txt', col_names = FALSE, delim=' ') %>%
  select(c(Chr = X1, Chromosome = X11)) %>%
  mutate(Chromosome = str_remove(Chromosome, ","))

chroms_12NC29 <- read_delim(file='/datasets/work/af-crown-rust/work/2020-02-28_GWAS/data/chroms.txt', col_names = FALSE, delim=' ') %>%
  select(c(Chr = X1, Chromosome = X11)) %>%
  mutate(Chromosome = str_remove(Chromosome, ","))

#bonferroni calculations
sites <- read_tsv(file='/datasets/work/af-crown-rust/work/2020-02-28_GWAS/analysis/freebayes/12SD80/nine/diff_1/mlm_nocomp2.txt') %>%
summarise(n())
cut_off_12SD80 <- 0.05/pull(sites,1)
#total number of sites for 12DS80 is 792835. Bonferroni cutoff:6.306482e-08
sites <- read_tsv(file='/datasets/work/af-crown-rust/work/2020-02-28_GWAS/analysis/freebayes/12NC29/interim/nine/diff_1/mlm_nocomp2.txt') %>%
  summarise(n())
cut_off_12NC29 <- 0.05/pull(sites,1)
#total number of sites for 12DS80 is 760492. Bonferroni cutoff:6.574691e-08

#generate tables
for(i in folders_nine_12SD80) {
    sig_table(genome='12SD80',oat_diff=i,pval=cut_off_12SD80,nchrom=2,chroms=chroms_12SD80, '_nine.tsv')
}

for(i in folders_nine_12NC29) {
  sig_table(genome='12NC29',oat_diff=i,pval=cut_off_12NC29,nchrom=2,chroms=chroms_12NC29, '_nine.tsv')
}

for(i in folders_bin_12NC29) {
  sig_table(genome='12NC29',oat_diff=i,pval=cut_off_12NC29,nchrom=2,chroms=chroms_12NC29, '_binary.tsv')
}
for (i in folders_nine) {
  my_manhattan(paste(i,'/mlm_nocomp2.txt', sep=''), chroms, '12SD80')
}

manhat_bw <- function(diff) {
  pc_stats <- read_tsv(file=paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/analysis/freebayes/12NC29/interim/nine/', diff, '/mlm_nocomp2.txt', sep='')) 

    
  to_plot <- pc_stats %>%
    select(c(Marker, Pos, Chr, p)) %>%    
    left_join(chroms) %>%
    filter(!is.na(p) & !is.nan(p))
  
  to_plot <- as.data.frame(to_plot)
  
  oat_diff <- pc_stats$Trait[1]
  
  plot <- ggman(to_plot,snp='Marker', bp='Pos',chrom='Chromosome', pvalue='p',title=oat_diff, sigLine = -log10(6.5747e-08), ylim=c(0,30)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    scale_colour_grey(start = 0.5,end = 0.6)
  return(plot)
}
792,835 
#total number of sites is 760491. Bonferroni cutoff: 6.5747e-08
0.05/760491

## code for zooming in on each peak and deciding whether to highlight

# belle - not required (no peaks)

# h548
h548 <- manhat_bw('diff_1')
peak='h548_7F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(h548, '000007F,',start=200000, end=206000, gene.tracks = FALSE, size=.2, title = peak)
dev.off()
# outcome: highlight 7F

# hi_fi
hifi <- manhat_bw('diff_2')
peak='hifi_179F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(hifi, '000179F,', start=206500, end=209500,gene.tracks = FALSE, size=.2, title = peak)
dev.off()

peak='hifi_314F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(hifi, '000314F,', gene.tracks = FALSE, size=.2, title = peak) #start=206500, end=209500
dev.off()

peak='hifi_310F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(hifi, '000310F,', gene.tracks = FALSE, size=.2, title = peak) #start=206500, end=209500
dev.off()

peak='hifi_566F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(hifi, '000566F,', gene.tracks = FALSE, size=.2, title = peak) #start=206500, end=209500
dev.off()

peak='hifi_589F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(hifi, '000589F,', gene.tracks = FALSE, size=.2, title = peak) #start=206500, end=209500
dev.off()

# outcome: 
#179F probable. 
#310 over LHS contig; possibly split contig. Quite spread out.
#314F over most of contig; possibly split contig. Quite spread out.
#566 possible. No SNPs in first half of contig
#589 unlikely.

## pc35
pc35 <- manhat_bw('diff_7')

peak='pc35_6F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc35, '000006F,', gene.tracks = FALSE, size=.2, title = peak, start=345000, end=415000)
dev.off()

# outcome
# 6F good

## pc 38
pc38 <- manhat_bw('diff_9')

peak='pc38_9F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000009F,', gene.tracks = FALSE, size=.2, title = peak, start=180000, end=205000)
dev.off()

peak='pc38_16F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000016F,', gene.tracks = FALSE, size=.2, title = peak, start=30000, end=80000)
dev.off()

peak='pc38_20F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000020F,', gene.tracks = FALSE, size=.2, title = peak, start=60000, end=160000)
dev.off()

peak='pc38_32F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000032F,', gene.tracks = FALSE, size=.2, title = peak, start=300000, end=400000)
dev.off()

peak='pc38_56F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000056F,', gene.tracks = FALSE, size=.2, title = peak, start=90000, end=120000)
dev.off()

peak='pc38_72F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000072F,', gene.tracks = FALSE, size=.2, title = peak, start=200000, end=230000)
dev.off()

peak='pc38_177F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000177F,', gene.tracks = FALSE, size=.2, title = peak, start=120000, end=180000)
dev.off()

peak='pc38_197F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000197F,', gene.tracks = FALSE, size=.2, title = peak, start=150000, end=180000)
dev.off()

peak='pc38_199F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc35, '000199F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)
dev.off()

peak='pc38_260F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc35, '000260F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)
dev.off()

peak='pc38_268F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000268F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)
dev.off()

peak='pc38_277F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000277F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)
dev.off()

peak='pc38_338F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000338F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)
dev.off()

peak='pc38_359F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000359F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)
dev.off()

peak='pc38_403F'
jpeg(paste('/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/results/images/zoom/12NC29_nine_man_', peak,'.jpeg', sep=''), width = 4, height = 3, units = 'in', res = 300)
ggmanZoom(pc38, '000403F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)
dev.off()

ggmanZoom(pc38, '000314F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)

# outcome
# 9F good
# 13F possible
# 16F probably, but gaps nearby
# 32F great
# 20F possible, gaps nearby and dispersed
# 56F probable - highlight
# 72F probable
# 173 remove
# 177 great
# 197 possible
# 199 possible
# 260 great
# 268 good, RHS chromosome
# 277 good
# 338 good
# 359 okay, entire chromosome. Large gaps in coverage.
# 403F good
# 310 possibly needs adding. Inspected due to comparison with PC39

## pc40
pc40 <- manhat_bw('diff_11')
ggmanZoom(pc40, '000436F,', gene.tracks = FALSE, size=.2, title = peak)#, start=150000, end=180000)

# 18F false
# 60F false
# 191 probably false
# 436F false
