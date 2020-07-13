#!/usr/bin/env Rscript

library(readr)
suppressPackageStartupMessages(library(dplyr))

message("Reading phenos")
pheno <- snakemake@input[['pheno']] %>%
  read_tsv(guess_max = 100000)

message("Reading eigenvectors")
eigens <- snakemake@input[['eigenvec']] %>%
  read_delim(delim=" ", col_names = c("FID", "IID", paste0("jointPC", 1:10)),
             col_types = "ccdddddddddd")

message("Adding eigenvectors to phenos")
pheno %>%
  left_join(eigens, by = c("FID", "IID")) %>%
  filter(!is.na(jointPC1)) %>%
  write_tsv(snakemake@output[[1]])


