library(readr)
library(dplyr)
library(magrittr)

phenos <- read_tsv(snakemake@input[["phenos"]], guess_max = 100000)
rawscores <- read_table2(snakemake@input[["scores"]], guess_max = 50000)

names(rawscores)[!(names(rawscores) %in% c("FID", "IID"))] %<>%
  stringr::str_replace("-", "_minus") %>% paste0("PRS_", .)

# function to add pheno columns
modscores <- . %>%
  left_join(phenos, by=c("FID", "IID")) %>%
  mutate(APOE = forcats::fct_infreq(as.factor(APOE)),
         status10 = status -1,
         status = ifelse(status == 1, "control",
                  ifelse(status == 2, "case", NA))) %>%
  mutate(all_status10 = all_status -1,
         all_status = ifelse(all_status == 1, "control",
                      ifelse(all_status == 2, "case", NA))) %>%
  mutate(status_superage = ifelse(all_status10 == 1, "case",
                                  ifelse(all_status10 == 0,
                                         ifelse(aaoaae2 >= 90, "90+ control","89- control"),
                                                all_status)))


scores <- rawscores %>% modscores

scores_std <- rawscores %>%
  mutate_if(is.numeric, scale) %>%
  modscores

save(list = c("scores", "scores_std"), file = snakemake@output[[1]])
