# Load R packages
# library(MCMCpack)
# library(miscTools)
# library(mvnfast)
# library(sn)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("EDPsource_cpp.cpp")
source("EDPsource_r.R")

ff = function(N,Scn){
  temp_data = generate_data(N,Scn)
  CE_true = temp_data$CE_true
}

fff = function(N,sim,Scn){
  ss = sapply(1:sim, function(l) ff(N,Scn))
  apply(ss, 1, mean)
}

run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(run_ID)

txt.title = paste0("Results/EDP_true_result.txt")
if (run_ID == 1) {
  df = data.frame(matrix(ncol = 37, nrow = 0))
  df_col_names = c("run_ID", paste0(rep(c("NIE", "NDE", "ATE"),12),rep(1:12,each=3)))
  colnames(df) = df_col_names
  write.table(df, file = txt.title, sep = "\t", row.names = FALSE, col.names = TRUE)
}

N = 1e4
sim = 2000 # (*500) => 1,000,000
sim_true = sapply(1:12, function(l) fff(N,sim,l))

allinfo = data.frame(run_ID,t(c(sim_true)))
write.table(allinfo, file = txt.title, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)

sim_true[1,]-sim_true[2,]
sim_true[2,]-sim_true[3,]
sim_true[1,]-sim_true[3,]
