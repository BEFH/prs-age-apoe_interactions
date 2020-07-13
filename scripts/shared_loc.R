#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
library(readr)

lz <- read_tsv(snakemake@input[["lz"]],
               col_names = c("chrom", "lzsnp", "pos"),
               col_types = "ic-i--")

adgc <- read_tsv(snakemake@input[["adgc"]],
               col_names = c("chrom", "adgcsnp", "pos"),
               col_types = "ic-i--")

isect <- inner_join(lz, adgc, by = c("chrom", "pos"))

write_lines(isect[["lzsnp"]], snakemake@output[["lz"]])
write_lines(isect[["adgcsnp"]], snakemake@output[["adgc"]])
