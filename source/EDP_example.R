# Load R packages
# library(MCMCpack)
# library(miscTools)
# library(mvnfast)
# library(sn)
library(Rcpp)
library(RcppArmadillo)

# library(devtools)
# devtools::install_github("RcppCore/RcppArmadillo")
setwd("C:/Users/WooJung/Documents/Rproject/EDPpostMediation/source")

# This code must be in your current directory or you can change the path.
# Load cpp code
sourceCpp("EDPsource_cpp.cpp")
# Load R code
source("EDPsource_r.R")

Scenario = 1

# {
#   # esttype = "mean"
#   esttype = "median"
#   save_cluster = TRUE
#   # save_cluster = FALSE
#   
#   # Extract ID for simulated dataset (specific to LSF computing cluster)
#   # Note: The LSB_JOBINDEX is specified in the bsub command using the -J
#   # option
#   run_ID = 1
#   
#   # Specify the seed so can repeat simulation later if necessary
#   set.seed(run_ID)
#   
#   # ------------------------------------------------------------------------------
#   # Define of constants (adjust to fit the data generating scenario) -------------
#   # ------------------------------------------------------------------------------
#   
#   # Define number of observations per draw when using Monte Carlo integration to
#   # integrate out confounders
#   D = 1e3
#   
#   # the number of MC integration per iteration for a MCMC chain
#   # gibbs_thin = 10
#   # gibbs_thin = 50
#   gibbs_thin = 1e2
#   # gibbs_num = 1e3
#   # gibbs_thin = floor(gibbs_iter/gibbs_num)
#   
#   # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
#   # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
#   if (N>999){
#     # for big N (e.g. N = 1000), can use less samples (check posterior)
#     gibbs_iter = 2e4
#     gibbs_burnin = 2e4
#   } else if(N>499){
#     # for N = 500
#     gibbs_iter = 2e4
#     gibbs_burnin = 5e4
#   } else{
#     # for small N (e.g. N = 250), need more samples
#     gibbs_iter = 5e4
#     gibbs_burnin = 5e4
#   }
#   
#   # Define number of Monte Carlo draws when using Monte Carlo integration to
#   # integrate out confounders
#   # num_MC = 1e4
#   # num_MC = 1e3
#   num_MC = 1e2
#   
#   # num_MC_prior
#   num_MC_prior = 1e5
#   # num_MC_prior = 1e4
#   
#   # ------------------------------------------------------------------------------
#   # End Definition of Constants --------------------------------------------------
#   # ------------------------------------------------------------------------------
#   
#   # Load data for the specific dataset by run_ID
#   temp_data = generate_data(N,Scenario)
#   
#   # Load outcome
#   Y = temp_data$Y
#   
#   # Load mediation
#   M = temp_data$M
#   
#   # Load covariates
#   V = temp_data$V
#   
#   # Load treatment
#   Z = temp_data$Z
#   
#   # Load covariates
#   C = temp_data$C
#   
#   # Load true causal effect given other variable
#   CE_true = temp_data$CE_true
# }

# EDPmediation_result = EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
#                                    D, num_MC, num_MC_prior, esttype, save_cluster)
# NIE_table = EDPmediation_result$NIE_result_mc
# NDE_table = EDPmediation_result$NDE_result_mc
# ATE_table = EDPmediation_result$ATE_result_mc
# 
# allinfo = rbind(c(CE_true[1]-CE_true[2],CE_true[2]-CE_true[3],CE_true[1]-CE_true[3]),
#                 matrix(c(NIE_table[-2],NDE_table[-2],ATE_table[-2]),nrow = 4))
# rownames(allinfo) = c("True","Est","q025","q975","CIlength")
# colnames(allinfo) = c("NIE","NDE","ATE")
# allinfo

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Define number of observations for each dataset
# N = 250
# {
#   # esttype = "mean"
#   esttype = "median"
#   save_cluster = TRUE
#   # save_cluster = FALSE
#   
#   # Extract ID for simulated dataset (specific to LSF computing cluster)
#   # Note: The LSB_JOBINDEX is specified in the bsub command using the -J
#   # option
#   run_ID = 1
#   
#   # Specify the seed so can repeat simulation later if necessary
#   set.seed(run_ID)
#   
#   # ------------------------------------------------------------------------------
#   # Define of constants (adjust to fit the data generating scenario) -------------
#   # ------------------------------------------------------------------------------
#   
#   # Define number of observations per draw when using Monte Carlo integration to
#   # integrate out confounders
#   D = 1e3
#   
#   # the number of MC integration per iteration for a MCMC chain
#   # gibbs_thin = 10
#   # gibbs_thin = 50
#   gibbs_thin = 1e2
#   # gibbs_num = 1e3
#   # gibbs_thin = floor(gibbs_iter/gibbs_num)
#   
#   # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
#   # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
#   if (N>999){
#     # for big N (e.g. N = 1000), can use less samples (check posterior)
#     gibbs_iter = 2e4
#     gibbs_burnin = 2e4
#   } else if(N>499){
#     # for N = 500
#     gibbs_iter = 2e4
#     gibbs_burnin = 5e4
#   } else{
#     # for small N (e.g. N = 250), need more samples
#     gibbs_iter = 5e4
#     gibbs_burnin = 5e4
#   }
#   
#   # Define number of Monte Carlo draws when using Monte Carlo integration to
#   # integrate out confounders
#   # num_MC = 1e4
#   # num_MC = 1e3
#   num_MC = 1e2
#   
#   # num_MC_prior
#   num_MC_prior = 1e5
#   # num_MC_prior = 1e4
#   
#   # ------------------------------------------------------------------------------
#   # End Definition of Constants --------------------------------------------------
#   # ------------------------------------------------------------------------------
#   
#   # Load data for the specific dataset by run_ID
#   temp_data = generate_data(N,Scenario)
#   
#   # Load outcome
#   Y = temp_data$Y
#   
#   # Load mediation
#   M = temp_data$M
#   
#   # Load covariates
#   V = temp_data$V
#   
#   # Load treatment
#   Z = temp_data$Z
#   
#   # Load covariates
#   C = temp_data$C
#   
#   # Load true causal effect given other variable
#   CE_true = temp_data$CE_true
# }
# N250time = system.time(N250EDPmediation_result <-
#                          EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
#                                       D, num_MC, num_MC_prior, esttype, save_cluster))
# N250time

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Define number of observations for each dataset
# N = 500
# {
#   # esttype = "mean"
#   esttype = "median"
#   save_cluster = TRUE
#   # save_cluster = FALSE
#   
#   # Extract ID for simulated dataset (specific to LSF computing cluster)
#   # Note: The LSB_JOBINDEX is specified in the bsub command using the -J
#   # option
#   run_ID = 1
#   
#   # Specify the seed so can repeat simulation later if necessary
#   set.seed(run_ID)
#   
#   # ------------------------------------------------------------------------------
#   # Define of constants (adjust to fit the data generating scenario) -------------
#   # ------------------------------------------------------------------------------
#   
#   # Define number of observations per draw when using Monte Carlo integration to
#   # integrate out confounders
#   D = 1e3
#   
#   # the number of MC integration per iteration for a MCMC chain
#   # gibbs_thin = 10
#   # gibbs_thin = 50
#   gibbs_thin = 1e2
#   # gibbs_num = 1e3
#   # gibbs_thin = floor(gibbs_iter/gibbs_num)
#   
#   # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
#   # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
#   if (N>999){
#     # for big N (e.g. N = 1000), can use less samples (check posterior)
#     gibbs_iter = 2e4
#     gibbs_burnin = 2e4
#   } else if(N>499){
#     # for N = 500
#     gibbs_iter = 2e4
#     gibbs_burnin = 5e4
#   } else{
#     # for small N (e.g. N = 250), need more samples
#     gibbs_iter = 5e4
#     gibbs_burnin = 5e4
#   }
#   
#   # Define number of Monte Carlo draws when using Monte Carlo integration to
#   # integrate out confounders
#   # num_MC = 1e4
#   # num_MC = 1e3
#   num_MC = 1e2
#   
#   # num_MC_prior
#   num_MC_prior = 1e5
#   # num_MC_prior = 1e4
#   
#   # ------------------------------------------------------------------------------
#   # End Definition of Constants --------------------------------------------------
#   # ------------------------------------------------------------------------------
#   
#   # Load data for the specific dataset by run_ID
#   temp_data = generate_data(N,Scenario)
#   
#   # Load outcome
#   Y = temp_data$Y
#   
#   # Load mediation
#   M = temp_data$M
#   
#   # Load covariates
#   V = temp_data$V
#   
#   # Load treatment
#   Z = temp_data$Z
#   
#   # Load covariates
#   C = temp_data$C
#   
#   # Load true causal effect given other variable
#   CE_true = temp_data$CE_true
# }
# N500time = system.time(N500EDPmediation_result <-
#                          EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
#                                       D, num_MC, num_MC_prior, esttype, save_cluster))
# N500time

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Define number of observations for each dataset
# N = 1000
# {
#   # esttype = "mean"
#   esttype = "median"
#   save_cluster = TRUE
#   # save_cluster = FALSE
#   
#   # Extract ID for simulated dataset (specific to LSF computing cluster)
#   # Note: The LSB_JOBINDEX is specified in the bsub command using the -J
#   # option
#   run_ID = 1
#   
#   # Specify the seed so can repeat simulation later if necessary
#   set.seed(run_ID)
#   
#   # ------------------------------------------------------------------------------
#   # Define of constants (adjust to fit the data generating scenario) -------------
#   # ------------------------------------------------------------------------------
#   
#   # Define number of observations per draw when using Monte Carlo integration to
#   # integrate out confounders
#   D = 1e3
#   
#   # the number of MC integration per iteration for a MCMC chain
#   # gibbs_thin = 10
#   # gibbs_thin = 50
#   gibbs_thin = 1e2
#   # gibbs_num = 1e3
#   # gibbs_thin = floor(gibbs_iter/gibbs_num)
#   
#   # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
#   # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
#   if (N>999){
#     # for big N (e.g. N = 1000), can use less samples (check posterior)
#     gibbs_iter = 2e4
#     gibbs_burnin = 2e4
#   } else if(N>499){
#     # for N = 500
#     gibbs_iter = 2e4
#     gibbs_burnin = 5e4
#   } else{
#     # for small N (e.g. N = 250), need more samples
#     gibbs_iter = 5e4
#     gibbs_burnin = 5e4
#   }
#   
#   # Define number of Monte Carlo draws when using Monte Carlo integration to
#   # integrate out confounders
#   # num_MC = 1e4
#   # num_MC = 1e3
#   num_MC = 1e2
#   
#   # num_MC_prior
#   num_MC_prior = 1e5
#   # num_MC_prior = 1e4
#   
#   # ------------------------------------------------------------------------------
#   # End Definition of Constants --------------------------------------------------
#   # ------------------------------------------------------------------------------
#   
#   # Load data for the specific dataset by run_ID
#   temp_data = generate_data(N,Scenario)
#   
#   # Load outcome
#   Y = temp_data$Y
#   
#   # Load mediation
#   M = temp_data$M
#   
#   # Load covariates
#   V = temp_data$V
#   
#   # Load treatment
#   Z = temp_data$Z
#   
#   # Load covariates
#   C = temp_data$C
#   
#   # Load true causal effect given other variable
#   CE_true = temp_data$CE_true
# }
# N1000time = system.time(N1000EDPmediation_result <-
#                           EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
#                                        D, num_MC, num_MC_prior, esttype, save_cluster))
# N1000time

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Define number of observations for each dataset
# N = 2500
# {
#   # esttype = "mean"
#   esttype = "median"
#   save_cluster = TRUE
#   # save_cluster = FALSE
#   
#   # Extract ID for simulated dataset (specific to LSF computing cluster)
#   # Note: The LSB_JOBINDEX is specified in the bsub command using the -J
#   # option
#   run_ID = 1
#   
#   # Specify the seed so can repeat simulation later if necessary
#   set.seed(run_ID)
#   
#   # ------------------------------------------------------------------------------
#   # Define of constants (adjust to fit the data generating scenario) -------------
#   # ------------------------------------------------------------------------------
#   
#   # Define number of observations per draw when using Monte Carlo integration to
#   # integrate out confounders
#   D = 1e3
#   
#   # the number of MC integration per iteration for a MCMC chain
#   # gibbs_thin = 10
#   # gibbs_thin = 50
#   gibbs_thin = 1e2
#   # gibbs_num = 1e3
#   # gibbs_thin = floor(gibbs_iter/gibbs_num)
#   
#   # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
#   # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
#   if (N>999){
#     # for big N (e.g. N = 1000), can use less samples (check posterior)
#     gibbs_iter = 2e4
#     gibbs_burnin = 2e4
#   } else if(N>499){
#     # for N = 500
#     gibbs_iter = 2e4
#     gibbs_burnin = 5e4
#   } else{
#     # for small N (e.g. N = 250), need more samples
#     gibbs_iter = 5e4
#     gibbs_burnin = 5e4
#   }
#   
#   # Define number of Monte Carlo draws when using Monte Carlo integration to
#   # integrate out confounders
#   # num_MC = 1e4
#   # num_MC = 1e3
#   num_MC = 1e2
#   
#   # num_MC_prior
#   num_MC_prior = 1e5
#   # num_MC_prior = 1e4
#   
#   # ------------------------------------------------------------------------------
#   # End Definition of Constants --------------------------------------------------
#   # ------------------------------------------------------------------------------
#   
#   # Load data for the specific dataset by run_ID
#   temp_data = generate_data(N,Scenario)
#   
#   # Load outcome
#   Y = temp_data$Y
#   
#   # Load mediation
#   M = temp_data$M
#   
#   # Load covariates
#   V = temp_data$V
#   
#   # Load treatment
#   Z = temp_data$Z
#   
#   # Load covariates
#   C = temp_data$C
#   
#   # Load true causal effect given other variable
#   CE_true = temp_data$CE_true
# }
# N2500time = system.time(N2500EDPmediation_result <-
#                           EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin,
#                                        D, num_MC, num_MC_prior, esttype, save_cluster))
# N2500time

N250time
N500time
N1000time
N2500time

# save.image("NtimeEDPM.Rdata")
