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
           remove = FALSE, sep = "[:alpha:]+") %>%
  mutate(year = case_when(as.numeric(year) <= 20 ~ str_c('20', year),     #reformat for later coalesce with SA
                          as.numeric(year) > 20 & as.numeric(year) <= 99 ~ str_c('19', year))) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(sample = if_else(sample == "PCA84-1-1_B03", "Pca84-1_H07_S85", sample)) %>%
  mutate(sample = if_else(sample == "PCA84-1-1_H07", "PCA84-1-1_B03_S102", sample))

#read in South African data and clean names to match lane df

SA <- read_excel(path = 'data/BoshoffWHP_OatCrownRust_SA2018_1.xlsx', skip = 1) %>%
  clean_names() %>%
  mutate(sample_id = str_replace(sample_id, '/', '-')) %>%
  mutate(sample_id = str_replace(sample_id, ' ', '')) 

#manually correcting sample names. See email chain from Eva.
SA[11,2] <- "Pca84-1_H07_S85"
SA[32,2] <- "PCA84-1-1_B03_S102"
setdiff(SA$sample_id, lanes$sample)

#Pca75 not in lanes df. PCA84-1-1 ambiguous, awaiting further info
#Update 21.4.19: ignoring PCA84-1-1 for now as it is a SA isolate.

#Import and clean phenotype data
#note: scores are on scale from 0 to 9, and are the mean of 
#two observations. Individual values for each observation are lower in 
#the excel sheet.

phenos <- read_excel(path = 'data/2016_Pca_isolates_phenotype_completeset.xlsx',
                    sheet = 3, n_max = 41) %>%
  column_to_rownames(var = "Differential Line") %>%
  select(-1) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble()

#merging data so far

merged <- full_join(lanes, SA, by = c("sample" = "sample_id")) %>%
  mutate(year = coalesce(year.x, year.y)) %>%
  select(-year.x, -year.y) %>%
  mutate(pustule = case_when(str_detect(sample, "[:digit:]{2}[:alpha:]{2}[:digit:]+-[123]$") ~
                               str_extract(sample, '[123]$'))) %>%
  mutate(state = case_when(str_detect(sample,"[:digit:]{2}[:alpha:]{2,4}[:digit:]*-[:digit:]+$") ~
                             str_extract(sample,'[:alpha:]{2}'))) %>% 
  mutate(buckthorn = if_else(str_detect(sample,"[:digit:]{2}[:alpha:]{2}BT-[:digit:]+$"), TRUE, FALSE)) %>%
  mutate(buckthorn = if_else(year == 2016 & !is.na(state), TRUE, buckthorn)) %>%
  full_join(phenos, by = "sample")

#melania's existing data for the published 60 isolates from 1990 and 2015

pub <- read_excel(path = 'data/TableS1_final.xlsx', skip = 1) %>%
  mutate(year = str_extract(`Isolate ID`, '^[:digit:]{2}')) %>%
  mutate(year = case_when(as.numeric(year) <= 20 ~ str_c('20', year),     #reformat for later coalesce with SA
                          as.numeric(year) > 20 & as.numeric(year) <= 99 ~ str_c('19', year))) %>%
  mutate(year = as.numeric(year)) %>%
  mutate(state = str_extract(`Isolate ID`,'[:alpha:]{2}')) %>%
  mutate(buckthorn = if_else(str_detect(`Location`,"Buckthorn") | str_detect(Location, "BT"), TRUE, FALSE))
  
#merge this in again...
merged <- full_join(merged, pub, by =c("sample" = "Isolate ID", "year", "state", "buckthorn"))

