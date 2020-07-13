#!/usr/bin/env Rscript
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
library(stringr)
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))

#logging:
logs <- file(snakemake@log[[1]])
sink(logs, append = T)
sink(logs, append = T, type = "message")

irem <- snakemake@input[["lz_cases"]] %>%
    read_tsv(col_names = c("FID", "IID"), col_types = "cc")

"Processing LZ phenos" %>% message
yun_phenos <- read_tsv(snakemake@input[["lz"]], col_types = "ccciciidd") %>%
  rename(yun_status = status) %>%
  mutate(status = ifelse(yun_status == "AD", 2,
    ifelse(yun_status %in% c("90+", "Control"), 1, NA)),
         consortium = "northwell") %>%
  anti_join(irem, by = c("FID", "IID"))

adgc_cols <- cols(
  .default = col_double(),
  FID = col_character(),
  IID = col_character(),
  cohort = col_character(),
  famstudy = col_skip(),
  rel_omit = col_logical(),
  ADSP_omit = col_logical(),
  ADSPibd_omit = col_logical(),
  ADGC_omit = col_logical(),
  QC_omit = col_logical(),
  QC_reason = col_logical(),
  keep_allQCibd = col_skip(),
  keep_noADSP = col_skip(),
  keep_ageAPOE = col_skip(),
  aaoaae_orig = col_skip(),
  InADSP = col_skip(),
  ADSP_Type = col_skip()
)

prs_status <- function(age, status) {
  ifelse(age >= 90 & status == 1,
         1, ifelse(status == 1, NA, status))
}

prs_status_under90 <- function(age, status) {
  ifelse(age < 90 & status == 1,
         1, ifelse(status == 1, NA, status))
}

"Reading ADGC phenos" %>% message
adgc_phenos <- read_tsv(snakemake@input[["adgc"]], col_types = adgc_cols) %>%
  select(-contains("PC"), -contains("pc"), -contains("age"), -contains("famfile")) %>%
  mutate(consortium = "ADGC")


#PCA <- read_tsv(snakemake@input[["PCA"]],
#  col_types = cols(.default = "d", FID = "c", IID = "c"))

"Merging and further processing phenotypes" %>% message
phenos <- bind_rows(adgc_phenos, yun_phenos) %>%
  #left_join(PCA, by = c("FID", "IID")) %>%
  #filter(!is.na(jointPC1)) %>%
  rename(all_status = status) %>%
  mutate(status = prs_status(aaoaae2, all_status),
         status_under = prs_status_under90(aaoaae2, all_status)) %>%
  mutate(status_nonorthwell = ifelse(consortium == "northwell", NA, all_status))

phenos %>% write_tsv(snakemake@output[["phenos"]])

message("all ages")
phenos %>% count(all_status, consortium)

message("with controls only 90+")
phenos %>% count(status, consortium)

message("with controls only 89-")
phenos %>% count(status_under, consortium)

message("APOE-33 only, with controls only 90+")
phenos %>% filter(APOE == 33) %>% count(status, consortium)


histo <- phenos %>%
  ggplot(aes(aaoaae2)) + geom_histogram() +
  facet_wrap(vars(consortium, status)) + theme_bw()

ggsave(snakemake@output[["hist"]], plot = histo)
