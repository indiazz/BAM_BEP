#libraries
library(gridExtra)
library(dplyr)
library(multDM)
library(dplyr)
library(tidyr)
#MDM test:
MC_reps <- 100
fName <- paste0("res_setting6_")
flist <- list.files("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/")
existing <- as.numeric(sapply(flist, FUN = function(x) as.numeric(strsplit(strsplit(x, fName)[[1]][2], ".Rdata"))))
existing
load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/",fName, "1.Rdata"))

#load needed functions
dir <- "C:/Users/20202943/Documents/BAM BEP/Multivariate DM code/"
source(paste0(dir, "MDM_scoretest.R"))
#source(paste0(dir, "modelselectionMDM.R"))
source(paste0(dir, "nieuwe MDM selection.R"))

#selecting settings
input_scores <- names(res)
input_scores <- c("es_list","vs1_list","vs1w_list","vs0_list","vs0w_list")
input_models <- names(res$es_list)
#choice of all models, parametric and non-parametric
all_models <- setdiff(input_models, c("GOF"))
non_parametric <- c("ecc.q", "ecc.s", "ssh", "decc.q")
parametric <-c("Frank","Frank_shuffle", "Gumbel","Gumbel_shuffle", "Clayton","Clayton_shuffle", "gca")

input_models <- all_models #parametric or non_parametric

rho0 <- c(0.25, 0.5, 0.75)
rho <- c(0.25, 0.5, 0.75)
s <- setNames(vector("list", length(input_scores)), input_scores)
chunks <- list()
cf <- 0 
out <- list(); k <- 0
#MDM for all cases
for (r0 in rho0){
  for(r in rho){
    for(ID in existing[!is.na(existing)]){
      load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/", fName, ID,".Rdata"))
      print(ID)
      for (score in input_scores){
        print(score)
        s  <- res[[score]]
        for (rr in seq_len(MC_reps)){
          ES_mat <- do.call(rbind, lapply(input_models, function(m) s[[m]][rr, ]))
          rownames(ES_mat) <- input_models
          keep <- apply(!is.na(ES_mat), 2, all)
          S_matrix <- ES_mat[, keep, drop = FALSE]
          #oud <- try(MDM_on_scores_chain(S_matrix, q = floor(ncol(S_matrix)^(1/3)), small_sample = TRUE))
          result <- MDM_on_scores_chain2(S_matrix = S_matrix, q = floor(ncol(S_matrix)^(1/3)),statistic = "Sc")
          result_selection <- selection(alpha = 0.05, q = floor(ncol(S_matrix)^(1/3)), statistic = "Sc", S_matrix = S_matrix)
          outcomes <- result_selection$outcomes
          best_model_id <- which.min(outcomes[, "Mean loss"])
          best_model <-  rownames(outcomes)[best_model_id]
          df_value <- if (!is.null(result$parameter["df"])) as.numeric(result$parameter["df"]) else (nrow(S) - 1)
          q_value  <- if (!is.null(result$parameter["q"]))  as.numeric(result$parameter["q"])  else q_use
          T_value  <- if (!is.null(result$parameter["T"]))  as.numeric(result$parameter["T"])  else ncol(S)
          
          k <- k +1
          out[[k]] <- data.frame( file =ID, score = score, rr = rr, rho = r, rho0 = r0, T = T_value, q_value, df = df_value, Statistic = as.numeric(result$statistic), p_value = as.numeric(result$p.value), best_model = best_model, stringsAsFactors = FALSE )
        }
      }
 
  if (length(out)) {
    cf <- cf + 1
    chunks[[cf]] <- do.call(rbind, out)
    print(cf)
      }
    }
  }
}

output <- do.call(rbind, chunks)

save(output, file = paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/","_MDM_output_with_selection_allmodels_SC_","_fcmodel_",6,".Rdata"))

# RQ3 research
#load outputs of MDM test
load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/_MDM_output_with_selection_parametric_SC_fcmodel_6.Rdata"))
output_para <- output
load(paste0("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/_MDM_output_with_selection_nonparametric__fcmodel_6.Rdata"))
output_nonpara <- output
load("C:/Users/20202943/Documents/BAM BEP/setting 6/goed/data/_MDM_output_with_selection_nieuworgineel2__fcmodel_6.Rdata")
output_all <- output

alpha <- 0.05

prep_output <- function(df) {
  df %>%
    mutate(
      file = as.character(file),
      score = as.character(score),
      rr = as.integer(rr),
      rho = as.numeric(rho),
      rho0 = as.numeric(rho0),
      Statistic = as.numeric(Statistic),
      p_value = as.numeric(p_value),
      reject = p_value < alpha,
      misspec = abs(rho - rho0)
    )
}


#prepped data
out_par    <- prep_output(output_para)
out_nonpar <- prep_output(output_nonpara)
out_all <- prep_output(output_all)

# elimination process winners agreement

win_share <- function(df) {
  df %>%
    group_by(score, misspec, best_model) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(score, misspec) %>%
    mutate(share = n / sum(n)) %>%
    ungroup()
}

win_par    <- win_share(out_par)
win_nonpar <- win_share(out_nonpar)
win_all <- win_share(out_all)

#wide format of elimination proces winners
win_wide_par <- out_par %>%
  dplyr::select(file, rr, score, best_model) %>%
  distinct() %>%   # safety
  pivot_wider(
    names_from  = score,
    values_from = best_model
  )

win_wide_nonpar <- out_nonpar %>%
  dplyr::select(file, rr, score, best_model) %>%
  distinct() %>%   # safety
  pivot_wider(
    names_from  = score,
    values_from = best_model
  )
win_wide_all <- out_all %>%
  dplyr::select(file, rr, score, best_model) %>%
  distinct() %>%   # safety
  pivot_wider(
    names_from  = score,
    values_from = best_model
  )

#agreement on the winners of the elimination process
agreement <- function(df, pairs_of_scores) {
  purrr::map_dfr(pairs_of_scores, function(p) {
    score_1 <- p[1]
    score_2 <- p[2]
    
    tibble(
      comparison = paste(score_1, "vs", score_2),
      agreement  = mean(df[[score_1]] == df[[score_2]], na.rm = TRUE)
    )})
}

pairs_of_scores <- list(
  c("es_list", "vs0_list"), c("es_list", "vs1_list"), c("es_list", "vs0w_list"), c("es_list", "vs1w_list"), c("vs0_list", "vs0w_list"), c("vs0_list", "vs1w_list"), c("vs0_list", "vs1_list"), c("vs1_list", "vs0w_list"), c("vs1_list", "vs1w_list"), c("vs1_list", "vs0_list"), c("vs1w_list", "vs0w_list"))

# winner agreement tables presented in report
tab_agreement_par <- agreement(win_wide_par, pairs_of_scores)
tab_agreement_nonpar <- agreement(win_wide_nonpar, pairs_of_scores)
tab_agreement_all <- agreement(win_wide_all, pairs_of_scores)


# check: winners in wide form by misspecification
win_wide_misspec_par <- out_par %>%
  mutate(misspec = abs(rho - rho0)) %>%
  dplyr::select(file, rr, misspec, score, best_model) %>%
  distinct() %>%
  pivot_wider(
    names_from  = score,
    values_from = best_model
  )
win_wide_misspec_nonpar <- out_nonpar %>%
  mutate(misspec = abs(rho - rho0)) %>%
  dplyr::select(file, rr, misspec, score, best_model) %>%
  distinct() %>%
  pivot_wider(
    names_from  = score,
    values_from = best_model
  )
win_wide_misspec_nall <- out_all %>%
  mutate(misspec = abs(rho - rho0)) %>%
  dplyr::select(file, rr, misspec, score, best_model) %>%
  distinct() %>%
  pivot_wider(
    names_from  = score,
    values_from = best_model
  )

#agreement on the winners per misspecification of the elimination process
agreement_misspec_par <- win_wide_misspec_par %>%
  group_by(misspec) %>%
  summarise(
    agree_ES_VS05 = mean(es_list == vs0_list, na.rm = TRUE),
    agree_ES_VS1  = mean(es_list == vs1_list, na.rm = TRUE),
    agree_ES_VS0w  = mean(es_list == vs0w_list, na.rm = TRUE),
    agree_ES_VS1w  = mean(es_list == vs1w_list, na.rm = TRUE),
    agree_VS0_VS0w  = mean(vs0_list == vs0w_list, na.rm = TRUE),
    agree_VS0_VS1w  = mean(vs0_list == vs1w_list, na.rm = TRUE),
    agree_VS0_VS1  = mean(vs0_list == vs1_list, na.rm = TRUE),
    agree_VS1_VS1w  = mean(vs1_list == vs0w_list, na.rm = TRUE),
    agree_VS1_VS0  = mean(vs1_list == vs0_list, na.rm = TRUE),
    agree_VS1w_VS0w  = mean(vs1w_list == vs0w_list, na.rm = TRUE),
    .groups = "drop"
  )


agreement_misspec_nonpar <- winner_wide_misspec_nonpar %>%
  group_by(misspec) %>%
  summarise(
    agree_ES_VS05 = mean(es_list == vs0_list, na.rm = TRUE),
    agree_ES_VS1  = mean(es_list == vs1_list, na.rm = TRUE),
    agree_ES_VS0w  = mean(es_list == vs0w_list, na.rm = TRUE),
    agree_ES_VS1w  = mean(es_list == vs1w_list, na.rm = TRUE),
    agree_VS0_VS0w  = mean(vs0_list == vs0w_list, na.rm = TRUE),
    agree_VS0_VS1w  = mean(vs0_list == vs1w_list, na.rm = TRUE),
    agree_VS0_VS1  = mean(vs0_list == vs1_list, na.rm = TRUE),
    agree_VS1_VS1w  = mean(vs1_list == vs0w_list, na.rm = TRUE),
    agree_VS1_VS0  = mean(vs1_list == vs0_list, na.rm = TRUE),
    agree_VS1w_VS0w  = mean(vs1w_list == vs0w_list, na.rm = TRUE),
    .groups = "drop"
  )

agreement_misspec_all <- winner_wide_misspec_nall %>%
  group_by(misspec) %>%
  summarise(
    agree_ES_VS05 = mean(es_list == vs0_list, na.rm = TRUE),
    agree_ES_VS1  = mean(es_list == vs1_list, na.rm = TRUE),
    agree_ES_VS0w  = mean(es_list == vs0w_list, na.rm = TRUE),
    agree_ES_VS1w  = mean(es_list == vs1w_list, na.rm = TRUE),
    agree_VS0_VS0w  = mean(vs0_list == vs0w_list, na.rm = TRUE),
    agree_VS0_VS1w  = mean(vs0_list == vs1w_list, na.rm = TRUE),
    agree_VS0_VS1  = mean(vs0_list == vs1_list, na.rm = TRUE),
    agree_VS1_VS1w  = mean(vs1_list == vs0w_list, na.rm = TRUE),
    agree_VS1_VS0  = mean(vs1_list == vs0_list, na.rm = TRUE),
    agree_VS1w_VS0w  = mean(vs1w_list == vs0w_list, na.rm = TRUE),
    .groups = "drop"
  )

# agreement on rejections of MDM test

wide_rject_par <- out_par %>%
  dplyr::select(file, rr, score, reject) %>%
  distinct() %>%   
  pivot_wider(
    names_from  = score,
    values_from = reject
  )

wide_rject_nonpar <- out_nonpar %>%
  dplyr::select(file, rr, score, reject) %>%
  distinct() %>%   
  pivot_wider(
    names_from  = score,
    values_from = reject
  )

wide_rject_all <- out_all %>%
  dplyr::select(file, rr, score, reject) %>%
  distinct() %>%   
  pivot_wider(
    names_from  = score,
    values_from = reject
  )


#rejection agreement tables presented in report
rejection_agreement_par <- agreement(wide_rject_par, pairs_of_scores)
rejection_agreement_nonpar <- agreement(wide_rject_nonpar, pairs_of_scores)
rejection_agreement_all <- agreement(wide_rject_all, pairs_of_scores)
