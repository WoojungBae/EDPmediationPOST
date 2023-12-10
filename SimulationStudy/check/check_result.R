resulttableprint = function(model){
  # filename = paste0("EDP_true_result.txt")
  # result = as.matrix(read.table(filename, header = TRUE, sep = "", dec = "."))[,-1]
  # t(matrix(apply(result,2,mean),nrow=3))
  
  CE_TRUE = c(6.531447, 5.341431, 4.391350,
              8.844246, 6.954398, 4.391681,
              9.346116, 9.041773, 6.511783,
              6.944848, 5.909748, 5.187499,
              9.356312, 7.712099, 5.187551,
              9.397109, 9.220117, 6.826676,
              5.051179, 3.861218, 2.911183,
              6.539863, 4.649904, 2.910963,
              8.646006, 8.345666, 5.835668,
              5.251504, 4.078904, 3.154426,
              6.788203, 4.925860, 3.154505,
              8.668085, 8.412642, 5.912853)
  
  CE_TRUE = matrix(CE_TRUE, ncol = 3,byrow = T)
  n_senarios = nrow(CE_TRUE)
  
  CE_TRUEtemp = CE_TRUE
  CE_TRUEtemp[,1] = CE_TRUE[,1] - CE_TRUE[,2]
  CE_TRUEtemp[,2] = CE_TRUE[,2] - CE_TRUE[,3]
  CE_TRUEtemp[,3] = CE_TRUE[,1] - CE_TRUE[,3]
  CE_TRUE = CE_TRUEtemp
  
  Ns = c(250, 500, 1000, 2500)
  # Ns = c(500, 1000, 2500, 5000)
  result_list = list(NA)
  result_temp = matrix(nrow=4,ncol=(n_senarios*3))
  result_matrix = matrix(nrow=(n_senarios*3),ncol=(length(Ns)*4))
  for(j in 1:length(Ns)){
    N = Ns[j]
    
    for (Scn in 1:n_senarios) {
      # ------------------------------------------------------------------------------
      # ------------------------------------------------------------------------------
      # ------------------------------------------------------------------------------
      filename = paste0(model,"Scn",Scn,"N",N,"_result.txt")
      result = as.matrix(read.table(filename, header = TRUE, sep = "", dec = "."))[,-1]
      n_sim = nrow(result)
      
      # True values of E_NIE, E_NDE and E_TE
      E_CE_true = CE_TRUE[Scn,]
      
      # NIE_result = result[,1:5]
      # NDE_result = result[,6:10]
      # ATE_result = result[,11:15]
      NIE_result = result[,1:4]
      NDE_result = result[,5:8]
      ATE_result = result[,9:12]
      
      NIE_cover = sum(ifelse((E_CE_true[1]>NIE_result[,2])*(E_CE_true[1]<NIE_result[,3])==1,1,0))
      NDE_cover = sum(ifelse((E_CE_true[2]>NDE_result[,2])*(E_CE_true[2]<NDE_result[,3])==1,1,0))
      ATE_cover = sum(ifelse((E_CE_true[3]>ATE_result[,2])*(E_CE_true[3]<ATE_result[,3])==1,1,0))
      E_CE_cover = c(NIE_cover,NDE_cover,ATE_cover)/n_sim
      
      # CE_cover = cbind(NIE_result[,5],NDE_result[,5],ATE_result[,5]);E_CE_cover = apply(CE_cover,2,mean)
      CE_result = cbind(NIE_result[,c(1,4)],NDE_result[,c(1,4)],ATE_result[,c(1,4)])
      
      # CE_est = matrix(apply(CE_result,2,mean),nrow=2);CE_est
      CE_est = matrix(apply(CE_result,2,median),nrow=2);CE_est
      CE_sd = apply(CE_result,2,sd)[c(1,3,5)];CE_sd
      
      E_CE_est = CE_est[1,]
      E_CE_length = CE_est[2,]
      E_CE_bias = E_CE_est - E_CE_true
      E_CE_mse = E_CE_bias^{2} + CE_sd^{2}
      
      result_table = rbind(E_CE_true,
                           E_CE_est,
                           E_CE_bias,
                           E_CE_mse,
                           E_CE_length,
                           E_CE_cover)
      result_table = round(result_table,4)
      rownames(result_table) = c("TRUE","EDP","Bias","MSE","CIlength","coverage")
      colnames(result_table) = c("E_NIE","E_NDE","E_TE")
      
      result_temp[,((Scn-1)*3+1):((Scn-1)*3+3)] = result_table[-(1:2),]
    }
    
    result_matrix[,((j-1)*4+1):((j-1)*4+4)] = t(result_temp)
  }
  result_matrix = cbind(c(t(CE_TRUE)),result_matrix)
  
  colnames(result_matrix) = c("True",rep(c("Bias","MSE","CIl","Cover"),4))
  rownames(result_matrix) = paste0(rep(c("NIE","NDE","ATE"),n_senarios),rep(1:n_senarios,each=3))
  
  result_matrix
}

wd = paste0("C:/Users/WooJung/Documents/Rproject/EDPpostMediation/SimulationStudy/EDPM/Results")
# wd = paste0("C:/Users/WooJung/Documents/Rproject/EDPpostMediation/SimulationStudy/EDPMwoPOST/Results")
# wd = paste0("C:/Users/WooJung/Documents/Rproject/EDPpostMediation/SimulationStudy/Parametric/Results")
setwd(wd)

model = "EDP"
result_model = resulttableprint("EDP")
result_model = round(result_model,2)
result_model







