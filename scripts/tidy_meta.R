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

headers %>% unname(unlist(column_labels[1,]))

#now read the rest of the file
lanes <- read_excel(path = 'data/Figueroa_Project_006_NovaSeq_Summary.xlsx',
                      skip = 11, col_names = headers) %>%
  clean_names()
