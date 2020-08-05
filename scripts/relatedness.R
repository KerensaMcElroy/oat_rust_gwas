library(tidyverse)
library(ggrepel)
clones <- c('15FL1-2','15Fl1-4',
            '15MN10-4','15MN10-5',
            '15ND19-2', '15ND19-5',
            '15ND20-3', '15ND20-4',
            '15NE8-4', '15NE8-5',
            '15NE9-1', '15NE9-3',
            '15SD11-1','15SD11-2',
            '15SD30-1','15SD30-3',
            '16MN102-1','16MN102-3',
            '16MN107-1','16MN107-3',
            '16MN110-1','16MN110-2',
            '16MN93-1', '16MN93-2',
            '16MN94-1','16MN94-2',
            '16MN98-2','16MN98-3',
            '16MN99-1','16MN99-2','16MN99-3',
            '17FL20-1','17FL20-2','17FL20-3',
            '17FL23-1','17FL23-2',
            '17IL76-1','17IL76-3',
            '17LA50-2','17LA50-3',
            '17MN145-1','17MN145-2',
            '17MN146-1','17MN146-2',
            '17MN150-1','17MN150-2',
            '17MS65-1','17MS65-3',
            '17SD126-1','17SD126-2','17SD126-3')

pca <- read_delim(file='results/tables/12NC29_all1.txt', delim='\t',skip=2,col_names = TRUE)
meta <- read_delim(file='data/sample_phenos.tsv', delim='\t')

pca <- pca %>% left_join(meta, by=c('Taxa'='sample')) %>%
  mutate(lab = ifelse(Taxa %in% clones, 1, 0))

pca %>% ggplot(aes(x=PC1, y=PC2, colour=factor(year), shape=factor(buckthorn), label=Taxa))+
  geom_point() + 
 # xlim(c(-15,45))+
  #ylim(c(-60,60))+
  scale_colour_manual(values=c('#edf8fb', '#b3cde3', '#8c96c6', '#8856a7', '#810f7c'))+
  geom_text_repel(data = subset(pca, lab==1),size=3,colour='black',segment.size = .1)
subset(pca,lab==1)
