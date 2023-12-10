# Load R packages
# library(MCMCpack)
# library(miscTools)
# library(mvnfast)
# library(sn)
library(Rcpp)
library(RcppArmadillo)

# library(devtools)
# devtools::install_github("RcppCore/RcppArmadillo")

# setwd("C:/Users/WooJung/Desktop")
# setwd("C:/Users/Jennifer/Desktop/EDPmediationpost/Source")
setwd("C:/Users/WooJung/Documents/Rproject/EDPmediationpost/RealDataAnalysis")

# This code must be in your current directory or you can change the path.
# Load cpp code
sourceCpp("EDPsource_cpp_realdata.cpp")
# Load R code
source("EDPsource_r_realdata.R")

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
#   # txt.title = paste0("Results/D",Scenario,"_result.txt")
#   # if (run_ID == 1) {
#   #   df = data.frame(matrix(ncol = 16, nrow = 0))
#   #   df_col_names = c("run_ID",
#   #                    "simEY00", "simEY01", "simEY11",
#   #                    "EY0", "EY01", "EY1",
#   #                    "NIE_table", "NDE_table", "E_Ym1z1_mc", 
#   #                    "E_Mz0_mc", "E_Mz1_mc",
#   #                    "EM0", "EM1", "simEM00", "simEM11")
#   #   colnames(df) = df_col_names
#   #   write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
#   # }
#   
#   # ------------------------------------------------------------------------------
#   # Define of constants (adjust to fit the data generating scenario) -------------
#   # ------------------------------------------------------------------------------
#   
#   # Define number of observations for each dataset
#   # N = 10000
#   # N = 2000
#   # N = 1000
#   N = 500
#   # N = 240
#   
#   # Define number of observations per draw when using Monte Carlo integration to
#   # integrate out confounders
#   # D = ceiling(2e4/N)
#   D = 2e4
#   # D = 1e3
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
#   # gibbs_iter = 1e3 * gibbs_thin
#   # gibbs_burnin = 5e4
#   # if (N>999){
#   #   # for big N (e.g. N = 1000), can use less samples (check posterior)
#   #   gibbs_iter = 2e4
#   #   gibbs_burnin = 2e4
#   # } else if(N>499){
#   #   # for N = 500
#   #   gibbs_iter = 2e4
#   #   gibbs_burnin = 5e4
#   # } else{
#   #   # for small N (e.g. N = 250), need more samples
#   #   gibbs_iter = 5e4
#   #   gibbs_burnin = 5e4
#   # }
#   gibbs_iter = 5e3
#   gibbs_burnin = 5e3
#   
#   # Define number of Monte Carlo draws when using Monte Carlo integration to
#   # integrate out confounders
#   # num_MC = 1e4
#   num_MC = 1e2
#   
#   # num_MC_prior
#   num_MC_prior = 1e5
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
# 
# EDPmediation_result = EDPmediation(Y, M, V, Z, C, gibbs_iter, gibbs_burnin, gibbs_thin, 
#                                    D, num_MC, num_MC_prior, esttype, save_cluster)
# 
# NIE_table = EDPmediation_result$NIE_result_mc
# NDE_table = EDPmediation_result$NDE_result_mc
# ATE_table = EDPmediation_result$ATE_result_mc
# 
# allinfo = data.frame(run_ID,t(c(NIE_table[-2],NDE_table[-2],ATE_table[-2])))
# 
# allinfomat = matrix(round(allinfo,4)[-1],nrow = 4)
# rownames(allinfomat) = c("Est","q025","q975","CIlength")
# colnames(allinfomat) = c("NIE","NDE","ATE")
# allinfomat
# c(CE_true[1]-CE_true[2],CE_true[2]-CE_true[3],CE_true[1]-CE_true[3])
# c(mean(Y[Z==1]),mean(Y[Z==0]))

# NIE_median_mc - NDE_median_mc
# NDE_median_mc - TE_median_mc
# 
# NIE_mean_mc - NDE_mean_mc
# NDE_mean_mc - TE_mean_mc
# NIE_mean_mc - TE_mean_mc
# 
# mean(Y[Z==1]) - mean(Y[Z==0])
# 
# # write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
# 
# # Make file_name where results will be stored for each dataset.
# # IMPORTANT: Make a folder called 'Results' in the folder where this code is run.
# file_name = paste0("Results/D",Scenario,"_bnp_edp_conOutcome_","n",N,"MCI",D,
#                    "_BurnIn",formatC(gibbs_burnin,format="fg",flag=0,width=5),
#                    "_Gibbs",formatC(gibbs_iter,format="fg",flag=0,width=5),
#                    "run_ID",formatC(run_ID,format="fg",flag=0,width=5),
#                    ".rdata", sep = "")
# 
# # Output results to a file to be read by Gather Results.R
# # while(!file.exists(file_name)) Sys.sleep(3)
# save.image(file = file_name)
# # if(run_ID==500) {save.image(file = file_name)}
