#!/usr/bin/env Rscript

pacman::p_load(
  dplyr, readr, tidyr, magrittr, forcats, stringr,
  survival, survminer, broom, ggplot2, here)

if (exists("snakemake")) {
  load(here(snakemake@input[[1]]))
  setwd(here(snakemake@params[["outdir"]]))
  #logging:
  logf <- here(snakemake@log[[1]])
  if (file.exists(logf)) file.remove(logf)
  file.create(logf)
  logs <- file(logf)
  sink(logs, append = T)
  sink(logs, append = T, type = "message")
} else {
  load(here("output/scores+phenos_withAPOE.Rdata"))
  setwd(here("analysis/primary"))
}

unpack_kv <- function(df, col) { #Unpack kv pairs like from summary
  dplyr::mutate(df, KVpairs = str_split(!!sym(col), ", ")) %>%
  tidyr::unnest(KVpairs) %>%
  tidyr::separate(KVpairs, into = c("key", "value"), sep = "=") %>%
  tidyr::pivot_wider(names_from = key, values_from = value) %>%
  dplyr::select(-!!col)
}

summary_with_function <- function(model) { #puts full model in summary call
  smry <- summary(model)
  frml <- model$formula %>%
    deparse %>%
    paste(collapse = " ") %>%
    str_squish
  smry$call <- smry$call %>%
    deparse %>%
    str_replace("formula = .+(?=,)", frml) %>%
    str2lang
  smry
}

add_quanta <- function(dataf, tail = 0.1) {
  lab_percent <- Vectorize(function(percentile, tail) {
    if (percentile <= tail) {
      sprintf("Bottom %s of PRS", pct)
    } else if (percentile >= 1 - tail) {
      sprintf("Top %s of PRS", pct)
    } else {
      sprintf("Middle %s of PRS", pct2)
    }
  })
  pct <- round(100 * tail) %>% sprintf("%i%%", .)
  pct2 <- round(100 * (1 - 2 * tail)) %>% sprintf("%i%%", .)
  dataf %>%
  mutate(percentile = percent_rank(SCORE)) %>%
  mutate(prs_top = lab_percent(percentile, tail))
}

prs_quantified <- scores_std %>%
  rename(age = aaoaae2) %>%
  filter(APOE %in% c(22, 23, 33, 34, 44)) %>% # & cohort != "LZ"
  mutate(apoe4_carrier = ifelse(APOE %in% c(34, 44),
                                "APOE4 Carrier", "APOE4 Non-carrier")) %>%
  mutate(SCORE = PRS_1e_minus05) %>%
  add_quanta(0.1) %>%
  mutate(prs_top = as.factor(prs_top) %>% fct_infreq())

surv_form_strat <- formula("Surv(age, all_status10) ~ jointPC1 + jointPC2 + strata(prs_top, apoe4_carrier, shortlabel = F)")
surv_form_ustrat <- formula("Surv(age, all_status10) ~ jointPC1 + jointPC2 + prs_top + strata(cohort, apoe4_carrier)")

message("Stratified survival Cox proportional hazards model:")
surv_mod <- survival::coxph(surv_form_strat, data = prs_quantified)
surv_mod %>% summary_with_function %>% print

message("unstratified survival Cox proportional hazards model:")
surv_mod_unstrat <- survival::coxph(surv_form_ustrat, data = prs_quantified)
surv_mod_unstrat %>% summary_with_function %>% print

message("Test for proportionality of residuals:")
cox.zph(surv_mod_unstrat) %>% print

pdf("plots/coxzph.pdf")
ggcoxzph(cox.zph(surv_mod_unstrat))
dev.off()

message("unstratified survival ANOVA:")
anova(surv_mod_unstrat) %>% print

kmplot <- surv_mod %>%
  survfit %>%
  tidy %>%
  unpack_kv("strata") %>%
  ggplot() +
    aes(x = time, y = estimate, color = apoe4_carrier, linetype = prs_top) +
    geom_step() +
    geom_point(aes(x = ifelse(n.censor > 0, time, -1)), shape = 3) +
    scale_y_continuous(breaks = seq(0, 1, 0.2),
                       minor_breaks = seq(0.1, 0.9, 0.2)) +
    coord_cartesian(xlim = c(60, 110), ylim = c(0, 1.01), expand = F) +
    theme_bw() +
    labs(x = "Age", y = "Proportion Unaffected",
         linetype = "PRS stratum", color = "APOE genotype") +
    scale_color_manual(values = c("#B45ADC", "#E66100")) +
    scale_linetype_manual(values =
      c("Middle 80% of PRS" = "solid",
        "Bottom 10% of PRS" = "dotted",
        "Top 10% of PRS" = "longdash")) +
    theme(legend.position = c(0.15, 0.2),
          plot.margin = unit(c(5.5, 11, 5.5, 5.5), "points"))

ggsave("plots/kmplot.pdf", plot = kmplot, height = 7, width = 7)
ggsave("plots/kmplot.png", plot = kmplot, height = 7, width = 7)
