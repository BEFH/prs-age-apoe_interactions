#!/usr/bin/env Rscript

pacman::p_load(dplyr, readr, tidyr, magrittr, broom, forcats, purrr, car,
               ggplot2, viridis, ggsignif, scales, pROC, sjPlot, sjmisc,
               rstatix, ggstance, knitr)

if (exists("snakemake")) {
  scoresfile <- snakemake@input[["scores"]]
  wdir <- snakemake@params[["outdir"]]
  thresh_standard <- snakemake@input[["thresh_standard"]]
  thresh_superager <- snakemake@input[["thresh_superager"]]
  sens <- snakemake@wildcards[["sensitivity"]]
} else {
  scoresfile <- "output/scores+phenos_noAPOE.Rdata"
  wdir <- "testing"
  thresh_standard <- "output/PRS_noAPOE_nonorthwell.prsice"
  thresh_superager <- "output/PRS_noAPOE_superagers.prsice"
  sens <- "primary"
}

owd <- getwd()
load(scoresfile)
setwd(wdir)

capt <- . %>% capture.output(file = statsfile, append = T, split = T)
msgcapt <- function(msg) {
  padded <- sprintf("%s\n\n", msg)
  cat(padded, file = statsfile, append = T)
  cat(padded)
}
capt_only <- . %>% sprintf("%s\n\n", .) %>% cat(file = statsfile, append = T)
sumcapt <- . %>% summary %>% capt

opencon <- function(fname) {
  if (file.exists(fname)) file.remove(fname)
  a <- file.create(fname)
  file(fname, open = "wt")
}

if (sens == "agesensitivity") {
  deaths <- paste0("/sc/arion/projects/LOAD/Data/ADGC/ADGC_2018/",
    "Brian_pipeline_withfamily_and_adc8/ADGC.phenotypes.withexclusions.tsv") %>%
    read_tsv(col_types = cols_only(FID = "c", IID = "c", age_death = "d"))

  simple_ages <- . %>%
    left_join(deaths, by = c("FID", "IID")) %>%
    filter(!is.na(aaoaae) | consortium == "northwell") %>%
    select(-age_death)

  scores %<>% simple_ages
  scores_std %<>% simple_ages
}

# %% Figure 1
statsfile <- opencon("stats/fig1stats.txt")

residualize <- function(scores_std) {
  confound_model <- lm(PRS_1e_minus05 ~ sex + jointPC1 + jointPC2 + jointPC3 +
    jointPC4 + jointPC5 + jointPC6 + jointPC7 + jointPC8 + jointPC9 + jointPC10,
    data = scores_std)
  mutate(scores_std, residualized = predict(confound_model) %>% unlist)
}

#run ANOVA and pairwise t-tests
message("residualized PRS")
capt_only("Fig 1a:\n\n\nresidualized PRS")

residualized <- residualize(scores_std)

residualized %$%
  aov(residualized ~ status_superage) %>%
  sumcapt

msgcapt("\nPairwise comparisons using t-tests wiht pooled SD:")
ttests <- residualized %>%
  filter(!is.na(status_superage)) %>%
  pairwise_t_test(residualized ~ status_superage,
    pool.sd = T, p.adjust.method = "bonferroni")

ttests %>% kable %>% capt

ttests_plot <- ttests %>%
  unite("pair", group2, group1, sep = " vs. ") %>%
  mutate(p.adj.signif = ifelse(
    p.adj.signif == "****", "***", p.adj.signif))

get_signif <- function(x) {
  ttests_plot %>%
    filter(pair == x) %>%
    pull(p.adj.signif)
}

offset <- 0.006
distance <- 0.001

distplot <- scores %>%
  filter(!is.na(status_superage)) %>%
  mutate(Status = status_superage %>%
    factor(levels = c("case", "89- control", "90+ control"))) %>%
  ggplot(aes(fill = Status, y = PRS_1e_minus05, x = Status)) +
  geom_violin() + geom_boxplot(width = .1, fill = "white") + theme_bw() +
  geom_signif(comparisons = list(c("case", "90+ control")),
    y_position = offset + 3 * distance,
    annotation = get_signif("case vs. 90+ control")) +
  geom_signif(comparisons = list(c("case", "89- control")),
    y_position = offset + 2 * distance,
    annotation = get_signif("case vs. 89- control")) +
  geom_signif(comparisons = list(c("89- control", "90+ control")),
    y_position = offset + distance,
    annotation = get_signif("90+ control vs. 89- control")) +
  ylab("PRS") + ggtitle("PRS distributions") +
  scale_fill_manual(values = c(
    case = "#785EF0", "89- control" = "#DC267F", "90+ control" = "#FE6100"))

ggsave("plots/prsdists.png", plot = distplot, width = 7, height = 7)
ggsave("plots/prsdists.pdf", plot = distplot, width = 7, height = 7)


## OR plots

add_quanta <- function(df, tail = 0.1) {
  df %>%
  mutate(percentile = percent_rank(SCORE)) %>%
  mutate(prs_top = ifelse(percentile <= tail, 0,
                          ifelse(percentile >= 1 - tail, 1, -1))) %>%
  filter(prs_top >= 0) %>%
  mutate(cohort = as_factor(cohort))
}

make_oddsratios <- function(scores, model_extras = "") {
  ten_pcs <- paste(paste0("jointPC", 1:10), collapse = " + ")

  mod_extremes <- function(df, model, model_extras) {
    mod_e <- sprintf("status10 ~ %s + %s + sex%s",
                     "prs_top", ten_pcs, model_extras)
    mod <- glm(mod_e, family = binomial, data = df)
    mod %>%
      tidy %>%
      mutate(lower = exp(estimate - (1.96 * std.error)),
      upper = exp(estimate + (1.96 * std.error)),
      estimate = exp(estimate),
      std.error = exp(std.error)) %>%
      mutate(PRS_model = model)
  }

  test_prs_extremes <- . %>%
    purrr::imap(mod_extremes, model_extras) %>%
    bind_rows %>%
    filter(term == "prs_top") %>%
    mutate(OR_type = "Extremes") %>%
    select(-term)

  mod_continuous <- function(df, model, model_extras) {
    mod_c <- sprintf("status10 ~ %s + %s + sex%s", "SCORE",
      ten_pcs, model_extras)
    glm(mod_c, family = binomial, data = df) %>%
     tidy %>%
     mutate(lower = exp(estimate - (1.96 * std.error)),
     upper = exp(estimate + (1.96 * std.error)),
     estimate = exp(estimate),
     std.error = exp(std.error)) %>%
     mutate(PRS_model = model)
  }

  test_prs_continuous <- . %>%
    purrr::imap(mod_continuous, model_extras) %>%
    bind_rows %>%
    filter(term == "SCORE") %>%
    mutate(OR_type = "Continuous") %>%
    select(-term)

    #select(OR = estimate,
    #  "OR (-95% CI)" = lower, "OR (-95% CI)" = upper,
    #  SE = std.error, P = p.value, "T value" = statistic)

  test_prs_both <- . %>%
    {bind_rows(test_prs_continuous(.), test_prs_extremes(.))}

  oddsratios <- scores %>%
    split(scores$model) %>%
    lapply(add_quanta, 0.1) %>%
    test_prs_both

  oddsratio_labs <- oddsratios %>%
    select(PRS_model, OR_type, estimate, lower, upper) %>%
    pivot_wider(names_from = OR_type,
                values_from = c(estimate, lower, upper)) %>%
    mutate(lab = sprintf(
     "%s\nContOR = %.2f (%.2f - %.2f)\nExtrOR = %.2f (%.2f - %.2f)",
     PRS_model, estimate_Continuous, lower_Continuous, upper_Continuous,
     estimate_Extremes, lower_Extremes, upper_Extremes)) %>%
    select(PRS_model, lab)

  oddsratios %>%
    left_join(oddsratio_labs, by = "PRS_model")
}

plot_oddsratios <- function(oddsratios) {
  oddsratios %>% ggplot() +
  aes(y = PRS_model, x = estimate,
      xmin = lower, xmax = upper, color = OR_type) +
  geom_point(position = ggstance::position_dodgev(-0.3)) +
  geom_errorbarh(height = .1, position = ggstance::position_dodgev(-0.3)) +
  geom_vline(xintercept = 1, lty = 2) +  # add a dotted line at x=1 after flip
  ylab("PRS Model") + xlab("OR (95% CI) [scale is logarithmic]") +
  theme_bw() + # use a white background
  scale_color_manual(values = c(Continuous = "#648FFF", Extremes = "#FFB000")) +
  labs(color = "OR Type")
}

#suppress warnings from dodgev
ggsave_nododgewarn <- function(...) {
  h <- function(w) {
    if( any( grepl( "position_dodgev", w) ) ) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(ggsave(...), warning = h)
}

e4lab <- . %>%
  mutate(apoe4_carrier = ifelse(APOE %in% c(34, 44),
                                "APOE4 Carrier", "APOE4 Non-carrier") %>%
         factor(levels = c("APOE4 Non-carrier", "APOE4 Carrier")))

scores_std_with89 <- scores_std %>%
  mutate(status_89minus = ifelse(
    status_superage %in% c("case", "89- control"), all_status10, NA)) %>%
  e4lab

make_oddsratios2 <- . %>%
  mutate(SCORE = PRS_1e_minus05) %>%
  pivot_longer(cols = c("status_89minus", "status10"),
               names_to = "model", values_to = "s") %>%
  mutate(
    model = dplyr::recode(model,
      "status_89minus" = "89- ctrl", "status10" = "90+ ctrl"),
    status10 = s) %>%
  make_oddsratios() %>%
  mutate(PRS_model = factor(PRS_model, levels = c("90+ ctrl", "89- ctrl")))

oddsratio_noagecont <- scores_std_with89 %>% make_oddsratios2

ortab <- oddsratio_noagecont %>% select(-lab)
ortab %>% write_tsv("stats/OR_table.tsv")

msgcapt("Fig 1b:\n\nOdds Ratios")
ortab %>% kable %>% capt

or_breaks <- c(0:10)

oddsratio_fig <- oddsratio_noagecont %>%
  plot_oddsratios() +
  coord_cartesian(xlim = c(1, 10)) +
  scale_x_continuous(trans = log_trans(), breaks = or_breaks,
    labels = c(as.character(or_breaks)))

ggsave_nododgewarn("plots/OR_nocont.png",
  plot = oddsratio_fig, width = 7, height = 1.5)
ggsave_nododgewarn("plots/OR_nocont.pdf",
  plot = oddsratio_fig, width = 7, height = 1.5)

make_oddsratios_strat <- function(df, stratum) {
  make_oddsratios2(df) %>%
    mutate(Stratum = stratum)
}

oddsratio_noagecont_strat <- scores_std_with89 %>%
  filter(APOE != 24) %>%
  split(., .$apoe4_carrier) %>%
  imap(make_oddsratios_strat) %>%
  bind_rows

oddsratio_fig_strat <- oddsratio_noagecont_strat %>%
  plot_oddsratios() +
  coord_cartesian(xlim = c(1, 10)) +
  scale_x_continuous(trans = log_trans(), breaks = or_breaks,
    labels = c(as.character(or_breaks))) +
  facet_grid(Stratum ~ .)

ggsave_nododgewarn("plots/OR_strat.png",
  plot = oddsratio_fig_strat, width = 7, height = 3)
ggsave_nododgewarn("plots/OR_strat.pdf",
  plot = oddsratio_fig_strat, width = 7, height = 3)

fig1 <- gridExtra::arrangeGrob(
  distplot + ggtitle("(a) PRS distributions"),
  oddsratio_fig + ggtitle("(b) PRS effect sizes"),
  layout_matrix = rbind(
   c(1),
   c(1),
   c(1),
   c(2)))

ggsave_nododgewarn("plots/fig1.png", fig1, width = 7, height = 7)
ggsave_nododgewarn("plots/fig1.pdf", fig1, width = 7, height = 7)

close(statsfile)

# %% Figure 2

statsfile <- opencon("stats/fig2stats.txt")

msgcapt("Figure 2")
msgcapt("Figure 2a")
msgcapt("Age is significant for status:")

apcommon <- scores_std %>%
  e4lab %>%
  filter(APOE != 24) %>%
  mutate(Status = all_status) %>%
  filter(!is.na(all_status10 + aaoaae2))

ageplot <- apcommon %>%
  ggplot() + aes(x = aaoaae2, y = PRS_1e_minus05, color = Status) +
  geom_point(size = 0.01) + geom_smooth(method = "lm") +
  xlab("Age") + ylab("PRS") +
  ggtitle("(a) Effect of age and APOE on PRS") +
  theme_bw() + facet_grid(~apoe4_carrier) +
  scale_color_manual(values = c(case = "#DC3220", control = "#005AB5"))

ggsave("plots/ageplot.png", plot = ageplot, height = 7, width = 7)
ggsave("plots/ageplot.pdf", plot = ageplot, height = 7, width = 7)


ageplot2 <- apcommon %>%
  ggplot() + aes(y = aaoaae2, x = PRS_1e_minus05, color = Status) +
  geom_point(size = 0.01) + geom_smooth(method = "lm") +
  ylab("Age") + xlab("PRS") +
  ggtitle("(a) Effect of age and APOE on PRS") +
  theme_bw() + facet_grid(~apoe4_carrier)

lm(PRS_1e_minus05 ~ aaoaae2 + APOE + sex + jointPC1 + jointPC2 + jointPC3 +
  jointPC4 + jointPC5 + jointPC6 + jointPC7 + jointPC8 + jointPC9 + jointPC10,
  data = scores_std) %>%
    sumcapt

msgcapt("Even when accounting for case-control status")
age_lm <- lm(PRS_1e_minus05 ~ aaoaae2 + APOE + all_status + sex +
  jointPC1 + jointPC2 + jointPC3 + jointPC4 + jointPC5 + jointPC6 + jointPC7 +
  jointPC8 + jointPC9 + jointPC10, data = scores_std)
age_lm %>% sumcapt

msgcapt("No colinearity:")
age_lm %>%
  vif %>%
  as_tibble(rownames = "coef") %>%
  kable %>%
  capt

age_strat <- scores_std %>%
  mutate(E4strat = ifelse(APOE %in% c(34, 44), "APOE4 carrier",
                          ifelse(APOE %in% c(22, 23, 33),
                                 "APOE4 noncarrier", NA))) %>%
  mutate(E4strat = E4strat %>%
           as.factor %>%
           fct_relevel(c("APOE4 noncarrier", "APOE4 carrier")))

case_age_strat <- age_strat %>%
  filter(all_status10 == 1 & !is.na(E4strat))

ctrl_age_strat <- age_strat %>%
  filter(all_status10 == 0 & !is.na(E4strat))

case_age <- scores_std %>%
  filter(all_status10 == 1)

msgcapt("\n\nFig 2b: Case Associations:")
msgcapt("Age is associated with PRS:")
lm(PRS_1e_minus05 ~ aaoaae2 + E4strat + sex + jointPC1 + jointPC2 + jointPC3 +
  jointPC4 + jointPC5 + jointPC6 + jointPC7 + jointPC8 + jointPC9 + jointPC10,
  data = case_age_strat) %>%
    sumcapt

msgcapt("Interaction: Age is associated with PRS in E4 only:")
prs_age_lm_strat <- lm(PRS_1e_minus05 ~ aaoaae2 * E4strat + sex + jointPC1 +
  jointPC2 + jointPC3 + jointPC4 + jointPC5 + jointPC6 + jointPC7 + jointPC8 +
  jointPC9 + jointPC10, data = case_age_strat)
prs_age_lm_strat %>% sumcapt

## YUN PUT T-TESTS HERE ##
msgcapt("T-test with all samples")

cases333444 <- age_strat %>%
  mutate(SCORE = PRS_1e_minus05) %>%
  add_quanta %>%
  filter(all_status10 == 1 & !is.na(E4strat))

casese4 <- cases333444 %>% filter(E4strat == "APOE4 carrier")

bottom_ages_e4 <- casese4 %>% filter(prs_top == 1) %>% pull(aaoaae2)
top_ages_e4 <- casese4 %>% filter(prs_top == 0) %>% pull(aaoaae2)

ttestalle4case <- t.test(bottom_ages_e4, top_ages_e4)

ttestalle4case %>% capt

msgcapt("T-test with e3e4 samples")
cases34 <- casese4 %>% filter(APOE == 34)


bottom_ages_e3e4 <- cases34 %>% filter(prs_top == 1) %>% pull(aaoaae2)
top_ages_e3e4 <- cases34 %>% filter(prs_top == 0) %>% pull(aaoaae2)

tteste3e4case <- t.test(bottom_ages_e3e4, top_ages_e3e4)

tteste3e4case %>% capt

msgcapt("\n\nFig 2c: Control Associations:")
msgcapt("Age is associated with PRS:")
lm(PRS_1e_minus05 ~ aaoaae2 + E4strat + sex + jointPC1 + jointPC2 + jointPC3 +
  jointPC4 + jointPC5 + jointPC6 + jointPC7 + jointPC8 + jointPC9 + jointPC10,
  data = ctrl_age_strat) %>%
    sumcapt

msgcapt("\nInteraction: There is no interaction in controls:")
prs_age_lm_strat_ctrl <- lm(PRS_1e_minus05 ~ aaoaae2 * E4strat + sex +
  jointPC1 + jointPC2 + jointPC3 + jointPC4 + jointPC5 + jointPC6 + jointPC7 +
  jointPC8 + jointPC9 + jointPC10, data = ctrl_age_strat)
prs_age_lm_strat_ctrl %>% sumcapt

stratpal <- c("#E66100", "#B45ADC")

intplot_common <- function(plt) {
  plt +
    xlab("Age") + ylab("PRS") +
    theme_bw(base_size = 8) +
    scale_color_manual(values = stratpal) +
    scale_fill_manual(values = stratpal) +
    coord_cartesian(ylim = c(-0.5, 0.5))
}

cont <- plot_model(prs_age_lm_strat_ctrl, type = "int") %>% intplot_common +
  ggtitle("(c) Age x APOE interaction in controls") +
  theme(legend.position = "none")

case <- plot_model(prs_age_lm_strat, type = "int") %>% intplot_common +
  ggtitle("(b) Age x APOE interaction in cases") +
  labs(color = "Stratum") +
  theme(legend.position = "bottom")

fig2 <- gridExtra::arrangeGrob(
  ageplot + theme_bw(base_size = 8) + theme(legend.position = "bottom"),
  case,
  cont,
  layout_matrix = rbind(
   c(1, 1, 1, 2, 2),
   c(1, 1, 1, 2, 2),
   c(1, 1, 1, 2, 2),
   c(1, 1, 1, 2, 2),
   c(1, 1, 1, 2, 2),
   c(1, 1, 1, 3, 3),
   c(1, 1, 1, 3, 3),
   c(1, 1, 1, 3, 3),
   c(1, 1, 1, 3, 3)))


fig2 %>% ggsave("plots/fig2.png", ., width = 8, height = 5)
fig2 %>% ggsave("plots/fig2.pdf", ., width = 8, height = 5)

status_glm_sa <- glm(all_status10 ~ aaoaae2 * PRS_1e_minus05 * E4strat + sex +
  jointPC1 + jointPC2 + jointPC3 + jointPC4 + jointPC5 + jointPC6 + jointPC7 +
  jointPC8 + jointPC9 + jointPC10, data = age_strat)
status_glm <- glm(all_status10 ~ aaoaae2 + PRS_1e_minus05 + E4strat + sex +
  jointPC1 + jointPC2 + jointPC3 + jointPC4 + jointPC5 + jointPC6 + jointPC7 +
  jointPC8 + jointPC9 + jointPC10, data = age_strat)

msgcapt("\n\nOveral AD status associations:")
msgcapt("All associated with status:")
status_glm %>% sumcapt
msgcapt("Associations depend on interactions:")
status_glm_sa %>% sumcapt

msgcapt("No colinearity:")
status_glm %>% vif %>% as_tibble(rownames = "coef") %>% kable %>% capt

# %% Supplements

## Threshold plot

setwd(owd)

prs_superager <- read_tsv(
  thresh_superager,
  col_types = "-dddddi") %>%
  mutate(optimization = "superager")

prs_standard <- read_tsv(
  thresh_standard,
  col_types = "-dddddi") %>%
  mutate(optimization = "standard")

setwd(wdir)

thresholds <- bind_rows(prs_superager, prs_standard)

threshhold_plt <- thresholds %>% ggplot() + aes(x = log10(Threshold), y = R2) +
  geom_line(aes(linetype = optimization)) +
  geom_point(aes(color = -log10(P), shape = optimization)) +
  theme_bw() + scale_color_viridis(option = "B", end = 0.85, begin = 0.1)

ggsave("plots/lineplot.png", plot = threshhold_plt)
ggsave("plots/lineplot.pdf", plot = threshhold_plt)

scores %>%
  filter(!is.na(status_superage)) %>%
  ggplot() +
  aes(fill = factor(consortium), y = PRS_1e_minus05, x = status_superage) +
  geom_violin()
ggsave("PRSdists_byconsortium_10e-5_rerun.png")

## ROCs


ROCs <- list(standard = roc(predictor = scores$PRS_1e_minus05,
                            response = scores$all_status10),
             superager = roc(predictor = scores$PRS_1e_minus05,
                             response = scores$status10))

ROCsp <- list(standard = roc(predictor = scores$PRS_1e_minus05,
                             response = scores$all_status10,
                             partial.auc = c(1, 0.75)),
             superager = roc(predictor = scores$PRS_1e_minus05,
                             response = scores$status10,
                             partial.auc = c(1, 0.75)))

ROCsp

has.partial.auc(ROCsp$standard)

roc.test(ROCs$standard, ROCs$superager, reuse.auc = FALSE)

names(ROCs) <- sprintf("%s: AUC = %.2f", names(ROCs), lapply(ROCs, pROC::auc))

roc.test(ROCsp$standard, ROCsp$superager, reuse.auc = FALSE)

ggroc(ROCs) +
    geom_abline(intercept = 1, slope = 1, linetype = 2,
      size = 0.1, color = "black") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", legend.direction = "vertical") +
    ggtitle(label = "ROCs for AD Risk Models") +
    labs(color = "Risk Model")

ggsave("plots/ROCs.png")
ggsave("plots/ROCs.pdf")

#ORs

oddsratio_agecont <- scores_std %>%
  mutate(SCORE = PRS_1e_minus05) %>%
  pivot_longer(cols = c("status10", "all_status10"),
               names_to = "model", values_to = "s") %>%
  rename(status10 = s) %>%
  mutate(model = dplyr::recode(model,
    "all_status10" = "All ctrl (age cont)",
    "status10" = "90+ ctrl (age cont)")) %>%
  make_oddsratios(" + aaoaae2")

oddsratios <- bind_rows(oddsratio_noagecont, oddsratio_agecont)

plot_oddsratios(oddsratios)
ggsave("plots/OR.png", width = 7, height = 2.5)
ggsave("plots/OR.pdf", width = 7, height = 2.5)
