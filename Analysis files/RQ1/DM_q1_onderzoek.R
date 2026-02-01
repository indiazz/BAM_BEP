#analysis of sensitivity univariate DM on scoring rule

#libraries
library(gridExtra)
library(dplyr)
library(tidyr)
library(purrr)


#load data
#load("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/_settingpval2_ECC_Q_6_obsmodel_6_fcmodel_6.Rdata")
#GCA
load("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/_settingpval2GCA_6_obsmodel_6_fcmodel_6.Rdata")
#SSH
#load("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/_settingpval2_SSH_6_obsmodel_6_fcmodel_6.Rdata")
#clean data frame
dfmc_n <- dfmc %>%
  filter(!score %in% c("crps_list", "param_list", "chosenCopula_list", "timing_list"))

#set alpha
alpha <- 0.05

#ADD DM decision meaning
add_dm_decision <- function(df, alpha = 0.05) {
  df %>%
    mutate(
      decision = case_when(
        is.na(p_value) ~ "NA",
        p_value < alpha & value < 0 ~ "ModelBetter",
        p_value < alpha & value > 0 ~ "BenchmarkBetter",
        TRUE ~ "NoDifference"
      )
    )
}

# disagreement table with counts (score,model)
disagreement_counts_score_model <- function(df_decision) {
  df_decision %>%
    group_by(score, model, decision) %>%
    summarise(n = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = decision,
      values_from = n,
      values_fill = 0
    ) %>%
    arrange(score, model)
}

# disagreement table with proportions (score,model)

disagreement_props_score_model <- function(df_decision) {
  df_decision %>%
    group_by(score, model, decision) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(score, model) %>%
    mutate(total = sum(n)) %>%
    ungroup() %>%
    mutate(prop = n / total) %>%
    select(score, model, decision, prop) %>%
    tidyr::pivot_wider(
      names_from  = decision,
      values_from = prop,
      values_fill = 0
    ) %>%
    arrange(score, model)
}

# disagreement table with counts (scores)
disagreement_counts_score <- function(df_decision) {
  df_decision %>%
    group_by(score, decision) %>%
    summarise(n = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = decision,
      values_from = n,
      values_fill = 0
    ) %>%
    arrange(score)
}

# disagreement table with proportions (scores)
disagreement_props_score <- function(df_decision) {
  df_decision %>%
    group_by(score, decision) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(score) %>%
    mutate(total = sum(n)) %>%
    ungroup() %>%
    mutate(prop = n / total) %>%
    select(score, decision, prop) %>%
    tidyr::pivot_wider(
      names_from  = decision,
      values_from = prop,
      values_fill = 0
    ) %>%
    arrange(score)
}

# add decisions to df
dfmc_dec <- add_dm_decision(dfmc_n, alpha = alpha)

# create disagreement tables
tab_counts_score_model <- disagreement_counts_score_model(dfmc_dec)
tab_props_score_model  <- disagreement_props_score_model(dfmc_dec)
tab_counts_score <- disagreement_counts_score(dfmc_dec)
tab_props_score  <- disagreement_props_score(dfmc_dec)

#wide format of results
df_cases <- dfmc_dec %>%
  group_by(file, copula, rho0, eps, sigma, rho, d, repetition, model, score) %>%
  mutate(mc_id = row_number()) %>%
  ungroup()

df_decisions_wide <- df_cases %>%
  select(file, copula, rho0, eps, sigma, rho, d, repetition, model, mc_id, score, decision) %>%
  pivot_wider(
    names_from  = score,
    values_from = decision
  )

#pairwise disagreement overall
pairwise_disagreement_overall <- function(df_wide) {
  score_cols <- setdiff(names(df_wide), c("file", "copula", "rho0", "eps", "sigma", "rho", "d", "repetition", "model", "mc_id"))
  pairs_of_scores <- t(combn(score_cols, 2))
  
  map_dfr(seq_len(nrow(pairs_of_scores)), function(i) {
    score_1 <- pairs_of_scores[i, 1]
    score_2 <- pairs_of_scores[i, 2]
    
    sub <- df_wide %>%
      select(all_of(c("file", "copula", "rho0", "eps", "sigma", "rho",
                      "d", "repetition", "model", "mc_id", score_1, score_2))) %>%
      filter(!is.na(.data[[score_1]]), !is.na(.data[[score_2]]))
    
    agree <- mean(sub[[score_1]] == sub[[score_2]])
    
    rank_reversal <- mean(
      (sub[[score_1]] == "ModelBetter"     & sub[[score_2]] == "BenchmarkBetter") |
        (sub[[score_1]] == "BenchmarkBetter" & sub[[score_2]] == "ModelBetter")
    )
    
    tibble(
      score1         = score_1,
      score2         = score_2,
      n_cases        = nrow(sub),
      agreement           = agree,
      disagree = 1 - agree,
      rank_reversal  = rank_reversal
    )
  })
}

pairwise_overall <- pairwise_disagreement_overall(df_decisions_wide)

pairwise_overall
