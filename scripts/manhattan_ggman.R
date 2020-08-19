library(tidyverse)
library(ggman)
library(qqman)

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

folders_nine <- list.files(path='/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/analysis/freebayes/12NC29/interim/nine', full.names=TRUE)

chroms <- read_delim(file='/datasets/work/AF_CROWN_RUST_WORK/2020-02-28_GWAS/data/chroms.txt', col_names = FALSE, delim=' ') %>%
  select(c(Chr = X1, Chromosome = X11)) 

for (i in folders_nine) {
  my_manhattan(paste(i,'/mlm_nocomp2.txt', sep=''), chroms, '12NC29')
}

#total number of sites is 760491. Bonferroni cutoff: 6.5747e-08
0.05/760491


