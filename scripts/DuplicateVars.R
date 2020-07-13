suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = T)
Infile <- args[1]

dupvar <- read_tsv(Infile, col_names = T) 
dupvar %>% 
  separate(IDS, c('Var1', 'Var2', 'Var3'), sep = " ") %>% 
  select(Var1) %>% 
  write_tsv(paste0(Infile, '.delete'))
