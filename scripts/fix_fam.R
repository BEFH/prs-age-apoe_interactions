#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)

old <- snakemake@input[["allfam"]]
newstem <- snakemake@params[["ins"]]
outstem <- snakemake@params[["out_"]]

col.n <- c("FID", "IID", "PID", "MID")
col.nn <- c("none", "FIDIID", "PID_new", "MID_new")
col.no <- c("FID", "IID", "PID", "MID", "Sex", "Phe")
col.t <- "cccc--"

sexcol <- snakemake@params[["sex"]]
phecol <- snakemake@params[["status"]]
pcol <- as.list(rep("c", 4))
names(pcol) <- c("FID", "IID", sexcol, phecol)
pcol <- do.call(cols_only, pcol)

phenos <- snakemake@input[["phenos"]] %>%
  read_tsv(col_types = pcol) %>%
  rename(Sex = !!sexcol, Phe = !!phecol)

new_fam <- newstem %>%
  paste0(".fam") %>%
  read_table2( col_names = col.nn, col_types = col.t) %>%
  separate(FIDIID, c("FID_new", "IID_new"), sep = "_",
           remove = F, extra = "merge")

read_table2(old, col_names = col.n, col_types = col.t) %>%
  unite("FIDIID", FID, IID, sep = "_", remove = F) %>%
  left_join(phenos, by = c("FID", "IID")) %>%
  right_join(new_fam, by = "FIDIID") %>%
  mutate(IID = ifelse(is.na(IID), IID_new, IID)) %>%
  mutate(FID = ifelse(is.na(FID), FID_new, FID)) %>%
  mutate(PID = ifelse(is.na(PID), PID_new, PID)) %>%
  mutate(MID = ifelse(is.na(MID), MID_new, MID)) %>%
  mutate(Sex = ifelse(is.na(Sex), 0, Sex)) %>%
  mutate(Phe = ifelse(is.na(Phe), -9, Phe)) %>%
  select(!!col.no) %>%
  write_tsv(paste0(outstem, ".fam"), col_names = F)

file.copy(paste0(newstem, ".bim"), paste0(outstem, ".bim"))
file.copy(paste0(newstem, ".bed"), paste0(outstem, ".bed"))

