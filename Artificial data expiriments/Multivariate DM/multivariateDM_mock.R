rm(list=ls())

#libraries
library(multDM)
#load data and sources
load("C:/Users/20202943/Documents/BAM BEP/CODE/artificial data experiments/data files/res_setting_mock_with_timewindow_1.Rdata")
dir <- "C:/Users/20202943/Documents/BAM BEP/Simulation code/Multivariate DM code/"
source(paste0(dir, "MDM_scoretest.R"))
source(paste0(dir, "nieuwe MDM selection.R"))

# when DM values from COBASE paper wanted to use:
# for (lst_name in c("crps_list", "es_list", "vs0w_list", "vs0_list","vs1w_list", "vs1_list")) {
#   
#   l <- get(lst_name)
#   
#   names(l)[names(l) == "COBASE-Gumbel-Q"] <- "Gumbel_shuffle"
#   names(l)[names(l) == "COBASE-Frank-Q"] <- "Frank_shuffle"
#   names(l)[names(l) == "COBASE-Clayton-Q"] <- "Clayton_shuffle"
#   names(l)[names(l) == "COBASE-CopGCA-Q"] <- "CopGCA_shuffle"
#   names(l)[names(l) == "ECC-Q"] <- "ecc.q"
#   names(l)[names(l) == "SSh-I14-Q"] <- "ssh"
#   names(l)[names(l) == "EMOS-Q"] <- "emos.q"
#   names(l)[names(l) == "ECC-R"] <- "ecc.r"
#   assign(lst_name, l)
# }


crps_list        <- res$crps_list
es_list          <- res$es_list
vs1_list         <- res$vs1_list
vs1w_list        <- res$vs1w_list
vs0_list         <- res$vs0_list
vs0w_list        <- res$vs0w_list


#choose input models for MDM computations

#input_models <- c("emos.q","Gumbel_shuffle","Gumbel", "Frank_shuffle","Frank","Clayton_shuffle","Clayton", "ecc.q","ecc.r", "ssh","CopGCA","CopGCA_shuffle")
#input_models <- c("ecc.q","ecc.r", "ssh")
#input_models <- c("Gumbel_shuffle","Gumbel", "Frank_shuffle","Frank","Clayton_shuffle","Clayton", "CopGCA","CopGCA_shuffle")
#input_models <- c("ecc.q", "ecc.s", "ssh", "decc.q")
input_models <- c("Frank","Frank_shuffle", "Gumbel","Gumbel_shuffle", "Clayton","Clayton_shuffle")
#input_models <- c("Frank","Frank_shuffle", "Gumbel","Gumbel_shuffle", "Clayton","Clayton_shuffle","ecc.q", "ecc.s", "ssh", "decc.q")


score_lists <- list(
  vs0  = vs0_list,
  vs0w = vs0w_list,
  vs1  = vs1_list,
  vs1w = vs1w_list,
  ES   = es_list
)

chunks <- vector("list", length(score_lists))
names(chunks) <- names(score_lists)

for (sc in names(score_lists)) {
#only take input models
  S_filter <- score_lists[[sc]][input_models]
  
  S_matrix <- do.call(rbind, S_filter)
  rownames(S_matrix) <- input_models

  q_val <- floor(ncol(S_matrix)^(1/3))
  
# run MDM using MDM_on_scores_chain2
  result <- MDM_on_scores_chain2(S_matrix = S_matrix,q= q_val,statistic = "S") 
  
  result_selection <- selection(alpha= 0.05,q= q_val,statistic = "Sc",S_matrix  = S_matrix)
  
  outcomes <- result_selection$outcomes
  best_id  <- which.min(outcomes[, "Mean loss"])
  best_model <- rownames(outcomes)[best_id]
  
  T_value <- if (!is.null(result$parameter["T"])) {
    as.numeric(result$parameter["T"])
  } else {
    ncol(S_matrix)
  }
  
  df_value <- if (!is.null(result$parameter["df"])) {
    as.numeric(result$parameter["df"])
  } else {
    NA_real_
  }
  
#output
  chunks[[sc]] <- data.frame(
    score       = sc,
    T           = T_value,
    q_value     = q_val,
    df          = df_value,
    Statistic   = as.numeric(result$statistic),
    p_value     = as.numeric(result$p.value),
    best_model  = best_model,
    stringsAsFactors = FALSE
  )
}
output_all <- do.call(rbind, chunks)
