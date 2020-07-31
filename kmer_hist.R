library(tidyverse)

single <- read_delim(file='results/rust_abundance/15MN27-3_abundance.txt', delim=' ',col_names=c('count','kmer_no'), col_types='nn')
single <- single %>% 
  mutate(kmers = count*kmer_no)

iso15ND19_2 <- ggplot(single, aes(x=count, y=kmer_no))+
  xlim(0,100)+
  ylim(1, 5000000)+
  geom_point()

iso15MN27_3 <- ggplot(single, aes(x=count, y=kmer_no))+
  xlim(0,100)+
  ylim(1, 5000000)+
  geom_point()
iso15MN27_3
ggsave(file='results/images/kmer_15ND19_2.pdf',iso15ND19_2)
print(single, n=40)

(sum(single$kmers) - single$kmers[1] - single$kmers[2])/30

data_path <- "results/rust_abundance"   # path to the data
files <- dir(data_path, pattern = '*abundance.txt') # get file names

data <- files %>%
  setNames(nm = .) %>% 
  map_dfr(~ read_delim(file.path(data_path, .), delim=' ', col_names=c('count','kmer_no'), col_types='nn'), .id = "isolate") %>%
  separate(isolate, into = c('isolate', NA), sep = '_abund') %>%
  group_by(isolate)

plot_combined <- data %>% ggplot(aes(x=count, y=kmer_no, colour = isolate))+
  xlim(1,100)+
  ylim(0, 5000000)+
  geom_line(show.legend = FALSE)
ggsave(file='results/images/kmer_combined.pdf',plot_combined)

plot_facet <- data %>% ggplot(aes(x=count, y=kmer_no))+
  xlim(1,100)+
  ylim(0, 5000000)+
  geom_line(show.legend = FALSE)+
  facet_wrap(vars(isolate)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(file='results/images/kmer_facet.pdf', plot_facet)
