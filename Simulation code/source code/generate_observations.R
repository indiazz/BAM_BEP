# code to generate multivariate observations 
# Note that parts of the notation may differ from the notation in the paper
#   ... in particular regarding the numbering of the settings which was changed during the revision (2 -> S1; 3 -> 2; 4 -> 3A, new: 3B)

# input:
#   model: numeric, indicating which model is used for the observations
#   nout: number of multivariate observations to be generated as evaluation period 
#   ninit: additional initial training period for model estimation purposes
#   d: dimension of the multivariate vectors
#   ... additional parameters, depending on the chosen model

# output:
#   list of two arrays "obs_init" and "obs" containing the observations
#   dimensions: ninit, d; and nout; d
#   first dimension: forecast instance; second dimension: dimension in multivariate setting

#obs <- generate_obs(model = 7, nout = 100, ninit = 50,d =2, rho0 = 0.1)

generate_obs <- function(model, nout, ninit, d, ...){
  require(MASS)
  
  # check input 
  if(any(!is.numeric(c(nout, ninit, d)))){
    stop("Input 'nout', 'ninit' and 'd' need to be numeric of length 1")
  }
  
  # initialize output arrays
  obs_init <- array(NA, dim = c(ninit, d))
  obs <- array(NA, dim = c(nout, d))
  
  # Setting 7 (multivariate Gaussian distribution as in Lerch et al.)
 
  if (model == 7) {
    # check if appropriate additional parameters are given
    input <- list(...)
    ind_check <- match(c("rho0"), names(input), nomatch = 0)
    if(ind_check == 0){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste("rho0", collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign additional parameters from input
    rho0 <- input$rho0
    
    # correlation matrix
    S <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:d){
        S[i,j] <- rho0^(abs(i-j))
      }
    }
    
    # mean vector
    mu <- rep(0, d)
    
    # observations
    obs_init <- mvrnorm(n = ninit, mu = mu, Sigma = S)
    obs <- mvrnorm(n = nout, mu = mu, Sigma = S)
  }
  return(list("obs_init" = obs_init, "obs" = obs))
}
# obs <- generate_obs2( data = data_mock,
#                       file_name               = "Mock_data",
#                      observation_columns     = c("obs"),
#                       ensemble_regex          = c("^M_"),
#                       store_uvpp = FALSE,
#                       output_dim_standard     = 10)

#generation of observations for artificial data Flos paper:
generate_obs2 <- function(data, file_name, observation_columns, ensemble_regex, store_uvpp = TRUE, ...){
      load_dimension_transform <- function(stations, observation_columns) {
      get_dimension <- function(stat, obs) {
        stat_idx <- match(stat, stations)
        obs_idx <- match(obs, observation_columns)
        
        return(stat_idx + length(stations) * (obs_idx - 1))
      }
      
      return(get_dimension)
    }
    print("running transform data...")
    # Stations and days from the data
    stations <- unique(data$station)
    days <- sort(unique(data$td))
    trainingDays <- 60
    nout <- length(days) - trainingDays
    d <- length(stations) * length(observation_columns)
    get_dimension <- load_dimension_transform(stations, observation_columns)
    
    # Ensemble members
    m <- sum(grepl(ensemble_regex[1], names(data))) # Each observation should have the same number of ens members
    
    ##############################################
    ## Ensemble and observation data structures ##
    ##############################################
    
    cat("Formatting ensemble and observations...\n")
    
   
    obs_init <- array(NA, dim = c(trainingDays, d))
    obs <- array(NA, dim = c(nout, d))
    
    # Dimensions consists of groups of (station, obs_type)
    for (var_idx in 1:length(observation_columns))
    {
      # Column names for the ensemble members
      obs_column <- observation_columns[var_idx]
      
      for (station_nr in stations)
      {
        # The dimension of the variable
        dd <- get_dimension(station_nr, obs_column)
        
        for (day in 1:length(days))
        {
          # Extract forecast for day and dim = index
          dat <- subset(data, td == days[day] & station == station_nr)
          if (day <= trainingDays) {
            obs_init[day, dd] <- unlist(dat[obs_column], use.names = FALSE)
          } else {
            obs[day - trainingDays, dd] <- unlist(dat[obs_column], use.names = FALSE)
          }
        }
      }
    }
    return(list(
      "obs"         = obs,
      "obs_init"    = obs_init))
}