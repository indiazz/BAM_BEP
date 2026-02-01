#libraries
library(dplyr)
library(tidyr)
library(stringr)

#load data
load("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/dfmc_bootstrapped_artificial.Rdata")
#long version
df_long <- dfmc_bootstrapped_artificial  %>%
  rename_with(~ str_replace(.x, "^p_bootstrap_", "pbootstrap_")) %>%
  pivot_longer(
    cols = matches("^(bootstrap|pbootstrap)_\\d+$"),
    names_to = c(".value", "bootstrap_num"),
    names_pattern = "^(bootstrap|pbootstrap)_(\\d+)$"
  ) %>%
  rename(p_bootstrap = pbootstrap) %>%
  mutate(bootstrap_num = as.integer(bootstrap_num)) %>%
  dplyr::select(model, score, bootstrap_num, p_bootstrap, bootstrap)

df_long_nice_names <- df_long %>%
  dplyr::rename(
    bootstrap_id = bootstrap_num,
    p_value = p_bootstrap,
    value    = bootstrap

  )

#add decision
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
df_dec <- add_dm_decision(df_long_nice_names)
# wide verison
df_decisions_wide <- df_dec %>%
  dplyr::select(model, score, bootstrap_id, decision) %>%
  pivot_wider(
    names_from  = score,
    values_from = decision
  )


compute_pairwise_disagreement_overall <- function(df_wide) {
    score_cols <- setdiff(names(df_wide),c("model", "bootstrap_id"))
    score_pairs <- t(combn(score_cols, 2))
    map_dfr(seq_len(nrow(score_pairs)), function(i) {
      score_1 <- score_pairs[i, 1]
      score_2 <- score_pairs[i, 2]
      
      sub <- df_wide %>%
        dplyr::select(all_of(c("model", "bootstrap_id", score_1, score_2))) %>%
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
        agreement         = agree,
        disagreement = 1 - agree,
        rank_reversal  = rank_reversal
      )
    })
  }

pairwise_overall <- compute_pairwise_disagreement_overall(df_decisions_wide)

pairwise_overall


#disagreement table proportions
disagreement_props_score <- function(df_dec) {
  df_dec %>%
    group_by(score, decision) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(score) %>%
    mutate(total = sum(n)) %>%
    ungroup() %>%
    mutate(prop = n / total) %>%
    dplyr::select(score, decision, prop) %>%
    tidyr::pivot_wider(names_from  = decision, values_from = prop, values_fill = 0) %>%
    arrange(score)
}
prop <- disagreement_props_score(df_dec)