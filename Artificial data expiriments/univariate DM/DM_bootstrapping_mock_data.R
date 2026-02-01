#bootstrapping method as used in the COBASE paper
#libraries
library(forecast)
library(future)
library(future.apply)
library(progressr)
#load data
load("C:/Users/20202943/Documents/BAM BEP/artificial data experiments/data files/res_setting_mock_with_timewindow_1.Rdata")
scores <- res
#benchmark
b <- c("ecc.q")


 compute_dfmc <- function(res, benchmark, count = 0, bootstrap = TRUE) {

  input_models <- names(res$es_list)
  input_models <- input_models[!input_models %in% c(benchmark, "EMOS")]
  input_models <- Filter(function(m) !is.logical(res$es_list[[m]]), input_models)

  input_scores <- setdiff(names(res), c("timing_list", "param_list", "chosenCopula_list"))

  nout <- 240
  set.seed(count)

  if (bootstrap) {
    sample_indices <- sample(1:nout, size = nout, replace = TRUE)
    stat_col <- paste0("bootstrap_", count)
    p_col    <- paste0("p_bootstrap_", count)
  } else {
    sample_indices <- 1:nout
    stat_col <- "full"
    p_col    <- "p_full"
  }

  # Preallocate
  dfmc <- expand.grid(model = input_models, score = input_scores, stringsAsFactors = FALSE)
  dfmc[[stat_col]] <- NA_real_
  dfmc[[p_col]]    <- NA_real_

  dm_power <- 1
  varestimator <- "bartlett"

  for (this_model in input_models) {
    for (this_score in input_scores) {

      ind <- which(dfmc$model == this_model & dfmc$score == this_score)

      tmp_stat <- 0
      tmp_p    <- 0


      if (grepl("[0-9]+$", this_score)) {

        stationNumber <- as.numeric(gsub("[^0-9]", "", gsub("[^-]+-", "", this_score)))
        n <- match(stationNumber, res$stations)

        e1 <- res[["crps_list"]][[this_model]][, n][sample_indices]
        e2 <- res[["crps_list"]][[benchmark]][, n][sample_indices]

        if (sum(abs(e1 - e2)) != 0) {
          tryDM <- try(dm.test(e1 = e1, e2 = e2, h = 1, power = dm_power, varestimator = varestimator),
                       silent = TRUE)
          if (!inherits(tryDM, "try-error")) {
            tmp_stat <- as.numeric(tryDM$statistic)
            tmp_p    <- as.numeric(tryDM$p.value)
          }
        }

      } else if (this_score == "crps_list") {

        e1 <- res[[this_score]][[this_model]][sample_indices]
        e2 <- res[[this_score]][[benchmark]][sample_indices]

        if (sum(abs(e1 - e2)) != 0) {
          tryDM <- try(dm.test(e1 = e1, e2 = e2, h = 1, power = dm_power, varestimator = varestimator),
                       silent = TRUE)
          if (!inherits(tryDM, "try-error")) {
            tmp_stat <- as.numeric(tryDM$statistic)
            tmp_p    <- as.numeric(tryDM$p.value)
          }
        }

      } else {

        e1 <- res[[this_score]][[this_model]][1, sample_indices]
        e2 <- res[[this_score]][[benchmark]][1, sample_indices]

        tryDM <- try(dm.test(e1 = e1, e2 = e2, h = 1, power = dm_power, varestimator = varestimator),
                     silent = TRUE)
        if (!inherits(tryDM, "try-error")) {
          tmp_stat <- as.numeric(tryDM$statistic)
          tmp_p    <- as.numeric(tryDM$p.value)
        }
      }

      dfmc[[stat_col]][ind] <- tmp_stat
      dfmc[[p_col]][ind]    <- tmp_p
    }
  }

  dfmc
}
compute_DM_scores <- function(res, benchmarks, parallelization = TRUE) {

  n_bootstrap_samples <- 100
  dfmc_bootstrapped_artificial <- data.frame()

  future::plan(future::sequential)

  # ---- Progress bar ----
  progressr::handlers(global = TRUE)
  progressr::handlers("cli")
  options(future.globals.maxSize = 2 * 1024^3)

  for (benchmark in benchmarks) {

    cat("In Postprocessing/Utilities/DM_util.R: Computing DM scores for",
        benchmark, "as the benchmark,\n")

    # Compute 1..n bootstrap runs + full sample
    data_frames <- progressr::with_progress({
      p <- progressr::progressor(along = 1:(n_bootstrap_samples + 1))

      future.apply::future_lapply(1:(n_bootstrap_samples + 1), function(i) {

        df <- if (i <= n_bootstrap_samples) {
          compute_dfmc(res = res, benchmark = benchmark, count = i, bootstrap = TRUE)
        } else {
          compute_dfmc(res = res, benchmark = benchmark, bootstrap = FALSE)
        }

        p()
        df
      }, future.seed = TRUE)
    })

    # Merge all runs by (model, score) only
    dfmc <- Reduce(
      f = function(x, y) merge(x, y, by = c("model", "score"), all = TRUE),
      x = data_frames
    )

    # Add benchmark label
    dfmc$benchmark <- benchmark

    # Append
    dfmc_bootstrapped_artificial <- rbind(dfmc_bootstrapped_artificial, dfmc)
  }

  # Save (adjust path if you prefer)
  fName <- "C:/Users/20202943/Documents/BAM BEP/artificial data experiments/data files/dfmc_bootstrapped_artificial.Rdata"
  message("In Postprocessing/Utilities/DM_util.R: Saving scores in: ", fName)

  if (!dir.exists(dirname(fName))) {
    dir.create(dirname(fName), recursive = TRUE)
  }

  save(dfmc_bootstrapped_artificial, file = fName)

  return(dfmc_bootstrapped_artificial)
}




dfmc_bootstrapped_artificial <- compute_DM_scores(res = res, benchmarks = b, parallelization = FALSE)
