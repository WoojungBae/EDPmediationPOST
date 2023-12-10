library(ggplot2)

resultboxplotprint = function(NN){
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
  CE_TRUE = cbind(CE_TRUE[,1] - CE_TRUE[,2],
                  CE_TRUE[,2] - CE_TRUE[,3],
                  CE_TRUE[,1] - CE_TRUE[,3])

  Ns = c(250, 500, 1000, 2500)
  Models = c("EDPM","EDPMwoPOST","Parametric")

  n_Ns = length(Ns)
  n_Models = length(Models)
  n_senarios = nrow(CE_TRUE)
  {
    df = data.frame(matrix(ncol = 8, nrow = 0))
    df_col_names = c("N", "Scn", "Model", "allmin", "allmax", "CE025", "CEmed", "CE975")
    colnames(df) = df_col_names

    wd = paste0("C:/Users/WooJung/Documents/Rproject/EDPpostMediation/SimulationStudy/")
    for(j in 1:n_Ns){
      N = Ns[j]
      for (Scn in 1:n_senarios) {
        for(l in 1:n_Models){
          Model = Models[l]

          # ------------------------------------------------------------------------------
          # ------------------------------------------------------------------------------
          # ------------------------------------------------------------------------------
          filename = paste0(wd,Model,"/Results/EDPScn",Scn,"N",N,"_result.txt")
          result = read.table(filename, header = TRUE, sep = "", dec = ".")[,-c(1,5,9,13,17)]

          # True values of E_NIE, E_NDE and E_ATE
          E_CE_true = CE_TRUE[Scn,]

          # Estimates for NIE, NDE and ATE
          # E_CE_est = apply(result,2,mean)
          E_CE_est = apply(result,2,median)

          CEmin = c(NIEmin = min(result$NIE.q025),
                    NDEmin = min(result$NDE.q025),
                    ATEmin = min(result$ATE.q025))
          CE025 = c(NIE025 = median(result$NIE.q025),
                    NDE025 = median(result$NDE.q025),
                    ATE025 = median(result$ATE.q025))
          CEmed = c(NIEmed = median(result$NIE),
                    NDEmed = median(result$NDE),
                    ATEmed = median(result$ATE))
          CE975 = c(NIE975 = median(result$NIE.q975),
                    NDE975 = median(result$NDE.q975),
                    ATE975 = median(result$ATE.q975))
          CEmax = c(NIEmax = max(result$NIE.q975),
                    NDEmax = max(result$NDE.q975),
                    ATEmax = max(result$ATE.q975))

          # Biases for NIE, NDE and ATE
          CEmin = CEmin - E_CE_true
          CE025 = CE025 - E_CE_true
          CEmed = CEmed - E_CE_true
          CE975 = CE975 - E_CE_true
          CEmax = CEmax - E_CE_true

          df = rbind(df,
                     data.frame(N = N,
                                Scn = paste0("Scn",Scn),
                                Model = paste0(Model),
                                CE = c("NIE","NDE","ATE"),
                                CEmin = CEmin,
                                CEmax = CEmax,
                                CE025 = CE025,
                                CEmed = CEmed,
                                CE975 = CE975))
        }
      }
    }
  }

  df$Scn = factor(df$Scn, levels = paste0("Scn",1:12))
  df$CE = factor(df$CE, levels = c("NIE","NDE","ATE"))
  df = df[which(df$N == NN),]

  wd = paste0("C:/Users/WooJung/Documents/Rproject");setwd(wd)
  ggplot(df, aes(x = CE, y = CEmed, fill = Model)) +
    geom_boxplot(stat = "identity",
                 aes(ymin = CEmin, ymax = CEmax,
                     lower = CE025, middle = CEmed, upper = CE975))  +
    facet_wrap(~Scn, ncol=3, scale="free") + theme(legend.position="top")
}

Ns = c(250, 500, 1000, 2500)

# resultboxplotprint(Ns[1])
# resultboxplotprint(Ns[2])
# resultboxplotprint(Ns[3])
# resultboxplotprint(Ns[4])

N = Ns[1]
cairo_ps(paste0("CEboxplot_","N",N,".eps"),onefile=F,height=12,width=10,fallback_resolution=600)
resultboxplotprint(N)
dev.off()

N = Ns[2]
cairo_ps(paste0("CEboxplot_","N",N,".eps"),onefile=F,height=12,width=10,fallback_resolution=600)
resultboxplotprint(N)
dev.off()

N = Ns[3]
cairo_ps(paste0("CEboxplot_","N",N,".eps"),onefile=F,height=12,width=10,fallback_resolution=600)
resultboxplotprint(N)
dev.off()

N = Ns[4]
cairo_ps(paste0("CEboxplot_","N",N,".eps"),onefile=F,height=12,width=10,fallback_resolution=600)
resultboxplotprint(N)
dev.off()

resultboxplotprint = function(NN){
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
  CE_TRUE = cbind(CE_TRUE[,1] - CE_TRUE[,2],
                  CE_TRUE[,2] - CE_TRUE[,3],
                  CE_TRUE[,1] - CE_TRUE[,3])

  Ns = c(250, 500, 1000, 2500)
  Models = c("EDPM","EDPMwoPOST","Parametric")

  n_Ns = length(Ns)
  n_Models = length(Models)
  n_senarios = nrow(CE_TRUE)
  {
    df = data.frame(matrix(ncol = 6, nrow = 0))
    df_col_names = c("N", "Scn", "Model", "CE")
    colnames(df) = df_col_names

    wd = paste0("C:/Users/WooJung/Documents/Rproject/EDPpostMediation/SimulationStudy/")
    for(j in 1:n_Ns){
      N = Ns[j]
      for (Scn in 1:n_senarios) {
        for(l in 1:n_Models){
          Model = Models[l]

          # ------------------------------------------------------------------------------
          # ------------------------------------------------------------------------------
          # ------------------------------------------------------------------------------
          filename = paste0(wd,Model,"/Results/EDPScn",Scn,"N",N,"_result.txt")
          result = read.table(filename, header = TRUE, sep = "", dec = ".")[,-c(1,5,9,13,17)]

          # True values of E_NIE, E_NDE and E_ATE
          E_CE_true = CE_TRUE[Scn,]

          # Estimates for NIE, NDE and ATE
          E_CE_est = apply(result,2,mean)
          # E_CE_est = apply(result,2,median)

          CE = c(result$NIE - E_CE_true[1],
                 result$NDE - E_CE_true[2],
                 result$ATE - E_CE_true[3])

          df = rbind(df,
                     data.frame(N = N,
                                Scn = paste0("Scn",Scn),
                                Model = paste0(Model),
                                CEnames = rep(c("NIE","NDE", "ATE"),each=nrow(result)),
                                CE = CE))
        }
      }
    }
  }

  df$Scn = factor(df$Scn, levels = paste0("Scn",1:12))
  df = df[which(df$N == NN),]

  wd = paste0("C:/Users/WooJung/Documents/Rproject");setwd(wd)
  ggplot(df, aes(x = CEnames, y = CE, fill = Model)) +
    geom_boxplot(coef = 500) +
    facet_wrap(~Scn, ncol=3, scale="free") + theme(legend.position="top")
}