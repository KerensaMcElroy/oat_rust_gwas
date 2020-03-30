library(tidyverse)
library(readxl)
library(janitor)

#column names are split across two rows. Extract them separately first
headers <- read_excel(path = 'data/Figueroa_Project_006_NovaSeq_Summary.xlsx',
                    skip = 9, col_names = FALSE, n_max = 2, 
                    col_types = "text", na = "")

headers <- data.frame(t(headers)) 
colnames(headers) <- c('one','two')

headers <- headers %>% mutate(one = as.character(one), 
                              two = as.character(two)) %>%
                       mutate(two = replace_na(two, '_')) 

headers <- data.frame(t(headers))

column_labels <- headers %>% summarize_all(str_c, collapse = "_")

headers <-  unname(unlist(column_labels))

#now read the rest of the file ,'state','sample_no','pustule'
lanes <- read_excel(path = 'data/Figueroa_Project_006_NovaSeq_Summary.xlsx',
                      skip = 11, col_names = headers) %>%
  clean_names() %>%
  separate(sample, c('year', 'remainder'),
           remove = FALSE, sep = "[:alpha:]")

#read in South African data and clean names to match lane df

SA <- read_excel(path = 'data/BoshoffWHP_OatCrownRust_SA2018_1.xlsx', skip = 1) %>%
  clean_names() %>%
  mutate(sample_id = str_replace(sample_id, '/', '-')) %>%
  mutate(sample_id = str_replace(sample_id, ' ', ''))

setdiff(SA$sample_id, lanes$sample)

#Pca75 not in lanes df. PCA84-1-1 ambiguous, awaiting further info
