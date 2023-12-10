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

Scenario = 1

# Define number of observations for each dataset
Ns = c(250, 500, 1000, 2500)
# Ns = c(500, 1000, 2500, 5000)
for (N in Ns){
  # Extract ID for simulated dataset (specific to LSF computing cluster)
  # Note: The LSB_JOBINDEX is specified in the bsub command using the -J
  # option
  run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  
  # Specify the seed so can repeat simulation later if necessary
  set.seed(run_ID)
  
  # esttype = "mean"
  esttype = "median"
  save_cluster = FALSE
  
  txt.title = paste0("Results/EDP","Scn",Scenario,"N",N,"_result.txt")
  if (run_ID == 1) {
    df = data.frame(matrix(ncol = 13, nrow = 0))
    df_col_names = c("run_ID", 
                     "NIE", "NIE.q025", "NIE.q975", "NIE.CIl95",
                     "NDE", "NDE.q025", "NDE.q975", "NDE.CIl95",
                     "ATE", "ATE.q025", "ATE.q975", "ATE.CIl95")
    colnames(df) = df_col_names
    write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
  
  # ------------------------------------------------------------------------------
  # Define of constants (adjust to fit the data generating scenario) -------------
  # ------------------------------------------------------------------------------
  
  # Define number of observations per draw when using Monte Carlo integration to
  # integrate out confounders
  D = 1e3
  
  # the number of MC integration per iteration for a MCMC chain
  # gibbs_thin = 10
  # gibbs_thin = 50
  gibbs_thin = 1e2
  # gibbs_num = 1e3
  # gibbs_thin = floor(gibbs_iter/gibbs_num)
  
  # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
  # gibbs_iter = 1e3 * gibbs_thin
  # gibbs_burnin = 5e4
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
  
  # Define number of Monte Carlo draws when using Monte Carlo integration to
  # integrate out confounders
  # num_MC = 1e4
  # num_MC = 1e3
  num_MC = 1e2
  
  # num_MC_prior
  num_MC_prior = 1e5
  # num_MC_prior = 1e4
  
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
  
  EDPmediation_result = EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin, 
                                     D, num_MC, num_MC_prior, esttype, save_cluster)
  NIE_table = EDPmediation_result$NIE_result_mc
  NDE_table = EDPmediation_result$NDE_result_mc
  ATE_table = EDPmediation_result$ATE_result_mc
  
  allinfo = data.frame(run_ID,t(c(NIE_table[-2],NDE_table[-2],ATE_table[-2])))
  write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
}
