# Load R packages
# library(MCMCpack)
# library(miscTools)
# library(mvnfast)
# library(sn)
library(Rcpp)
library(RcppArmadillo)

# This code must be in your current directory or you can change the path.
# Load cpp code
sourceCpp("EDPsource_cpp.cpp")
# Load R code
source("EDPsource_r.R")

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
Scenario = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

txt.title = paste0("Results/EDP_time_result.txt")
if (Scenario == 1) {
  df = data.frame(matrix(ncol = 5, nrow = 0))
  df_col_names = c("Scn","N250","N500","N1000","N2500")
  colnames(df) = df_col_names
  write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
}

# Specify the seed so can repeat simulation later if necessary
set.seed(1)

# ------------------------------------------------------------------------------
# Define of constants (adjust to fit the data generating scenario) -------------
# ------------------------------------------------------------------------------

esttype = "median"
save_cluster = FALSE

# Define number of observations per draw when using Monte Carlo integration to
# integrate out confounders
D = 1e3

# the number of MC integration per iteration for a MCMC chain
gibbs_thin = 1e2

# Define number of Monte Carlo draws when using Monte Carlo integration to
# integrate out confounders
num_MC = 1e2

# num_MC_prior
num_MC_prior = 1e5

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Define number of observations for each dataset
N = 250
{
  # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
  if (N>999){
    # for big N (e.g. N = 1000), can use less samples (check posterior)
    gibbs_iter = 2e4
    gibbs_burnin = 2e4
  } else if(N>499){
    # for N = 500
    gibbs_iter = 2e4
    gibbs_burnin = 5e4
  } else{
    # for small N (e.g. N = 250), need more samples
    gibbs_iter = 5e4
    gibbs_burnin = 5e4
  }

  # ------------------------------------------------------------------------------
  # End Definition of Constants --------------------------------------------------
  # ------------------------------------------------------------------------------

  # Load data for the specific dataset by run_ID
  temp_data = generate_data(N,Scenario)

  # Load outcome
  Y = temp_data$Y

  # Load mediation
  M = temp_data$M

  # Load covariates
  V = temp_data$V

  # Load treatment
  Z = temp_data$Z

  # Load covariates
  C = temp_data$C

  # Load true causal effect given other variable
  CE_true = temp_data$CE_true
}
N250time = system.time(N250EDPmediation_result <-
                         EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
                                      D, num_MC, num_MC_prior, esttype, save_cluster))
N250time

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Define number of observations for each dataset
N = 500
{
  # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
  if (N>999){
    # for big N (e.g. N = 1000), can use less samples (check posterior)
    gibbs_iter = 2e4
    gibbs_burnin = 2e4
  } else if(N>499){
    # for N = 500
    gibbs_iter = 2e4
    gibbs_burnin = 5e4
  } else{
    # for small N (e.g. N = 250), need more samples
    gibbs_iter = 5e4
    gibbs_burnin = 5e4
  }

  # ------------------------------------------------------------------------------
  # End Definition of Constants --------------------------------------------------
  # ------------------------------------------------------------------------------

  # Load data for the specific dataset by run_ID
  temp_data = generate_data(N,Scenario)

  # Load outcome
  Y = temp_data$Y

  # Load mediation
  M = temp_data$M

  # Load covariates
  V = temp_data$V

  # Load treatment
  Z = temp_data$Z

  # Load covariates
  C = temp_data$C

  # Load true causal effect given other variable
  CE_true = temp_data$CE_true
}
N500time = system.time(N500EDPmediation_result <-
                         EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
                                      D, num_MC, num_MC_prior, esttype, save_cluster))
N500time

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Define number of observations for each dataset
N = 1000
{
  # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
  if (N>999){
    # for big N (e.g. N = 1000), can use less samples (check posterior)
    gibbs_iter = 2e4
    gibbs_burnin = 2e4
  } else if(N>499){
    # for N = 500
    gibbs_iter = 2e4
    gibbs_burnin = 5e4
  } else{
    # for small N (e.g. N = 250), need more samples
    gibbs_iter = 5e4
    gibbs_burnin = 5e4
  }

  # ------------------------------------------------------------------------------
  # End Definition of Constants --------------------------------------------------
  # ------------------------------------------------------------------------------

  # Load data for the specific dataset by run_ID
  temp_data = generate_data(N,Scenario)

  # Load outcome
  Y = temp_data$Y

  # Load mediation
  M = temp_data$M

  # Load covariates
  V = temp_data$V

  # Load treatment
  Z = temp_data$Z

  # Load covariates
  C = temp_data$C

  # Load true causal effect given other variable
  CE_true = temp_data$CE_true
}
N1000time = system.time(N1000EDPmediation_result <-
                          EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
                                       D, num_MC, num_MC_prior, esttype, save_cluster))
N1000time

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Define number of observations for each dataset
N = 2500
{
  # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
  if (N>999){
    # for big N (e.g. N = 1000), can use less samples (check posterior)
    gibbs_iter = 2e4
    gibbs_burnin = 2e4
  } else if(N>499){
    # for N = 500
    gibbs_iter = 2e4
    gibbs_burnin = 5e4
  } else{
    # for small N (e.g. N = 250), need more samples
    gibbs_iter = 5e4
    gibbs_burnin = 5e4
  }

  # ------------------------------------------------------------------------------
  # End Definition of Constants --------------------------------------------------
  # ------------------------------------------------------------------------------

  # Load data for the specific dataset by run_ID
  temp_data = generate_data(N,Scenario)

  # Load outcome
  Y = temp_data$Y

  # Load mediation
  M = temp_data$M

  # Load covariates
  V = temp_data$V

  # Load treatment
  Z = temp_data$Z

  # Load covariates
  C = temp_data$C

  # Load true causal effect given other variable
  CE_true = temp_data$CE_true
}
N2500time = system.time(N2500EDPmediation_result <-
                          EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
                                       D, num_MC, num_MC_prior, esttype, save_cluster))
N2500time

allinfo = data.frame(t(c(Scenario,N250time[1],N500time[1],N1000time[1],N2500time[1])))
write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
