#!/usr/bin/env Rscript
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
library(stringr)
suppressPackageStartupMessages(library(tidyr))


#logging:
logf <- snakemake@log[[1]]
if (file.exists(logf)) file.remove(logf)
file.create(logf)
logs <- file(logf)
sink(logs, append = T)
sink(logs, append = T, type = "message")

mk_ikeep <- function(df, ikeep) {
  df %>% summarise(N = n(), .groups = "drop") %>% print
  df %>%
    ungroup %>%
    select(FID, IID) %>%
    write_tsv(snakemake@output[[ikeep]], col_names = F)
}

message("Reading phenotypes")
phenos <- read_tsv(snakemake@input[[1]], guess_max = 1000000)

p_withpheno <- phenos %>%
  filter(!ADGC_omit & !is.na(APOE)) %>%
  mutate(status_cc = ifelse(status == 1, "control", "case")) %>%
  group_by(status_cc)

message("Has phenotypes:")
p_withpheno %>% mk_ikeep("ikeep_withIGAP")

message("Has phenotypes and passes QC:")

noremrel <- c( #do not remove samples that are related to future studies
  "ADC4-NACC897662", "ADC4-NACC813784", "ADC5-NACC447248", "ADC6-NACC623541",
  "ADC7-NACC817085", "ADC7-NACC008802", "ADC7-NACC064202", "ADC7-NACC323196"
)
extraremrel <- c("ADC7-NACC407705", "WHICAP-RM0135", "WHICAP-RX0174")

p_QCpass <- p_withpheno %>% filter(
  ((!rel_omit & !QC_omit) | IID %in% noremrel) & !(IID %in% extraremrel)
)
p_QCpass %>% mk_ikeep("ikeep_withIGAP_qc")

## Not in Lambert
# ACT2
# ADC 4-7
# BIOCARD
# CHAP
# EAS
# MTV
# NBB
# ROSMAP2
# WASHU2
# WHICAP

message("Has phenotypes, isn't in IGAP and passes QC:")
p_QCpass_noIGAP <- p_QCpass %>% filter(
  cohort %in% c("ACT2", "ADC4", "ADC5", "ADC6", "ADC7", "BIOCARD", "CHAP",
                "EAS", "MTV", "NBB", "ROSMAP2", "WASHU2", "WHICAP"))
p_QCpass_noIGAP %>% mk_ikeep("ikeep_noIGAP_qc")

p_QCpass_noIGAP %>% ungroup %>% select(-status_cc) %>% write_tsv(snakemake@output[["pheno"]])

message("APOE genotypes of retained controls over 90")
p_QCpass_noIGAP %>%
  ungroup %>%
  filter(status == 1 & aaoaae2 > 90) %>%
  count(APOE) %>%
  mutate(APOE = as.character(APOE)) %>%
  add_row(APOE = "total", n = sum(.$n))
