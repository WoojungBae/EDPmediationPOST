# Load R packages
library(Rcpp)
library(RcppArmadillo)
library(dplyr)
library(mgcv)
library(openxlsx) 
library(lme4)

# Set Directory
setwd("C:/Users/WooJung/Documents/Rproject/EDPpostMediation/RealDataAnalysis")

# -------------------------------------------------------------------------
# ---------------------------------- DATA ---------------------------------
# -------------------------------------------------------------------------
# Here are the details for this analysis:
# - outcome: 24 month weight change (spreadsheet 2)
# - mediator - Phase 2 attendance (measured as a percentage) - spreadsheet 1
# - baseline confounders - mostly Table 1 from the paper and previous spreadsheet I sent
# - potential post treatment confounder - 6 month weight (or 6 month weight change) - spreadsheet 4
# - ok to merge control/low and moderate/high arms (so only two arms)

# Load Data
data1 = read.xlsx('RL_WTS_24mo_MD.xlsx', sheet = 'wts24', startRow = 2)
# data1$ID
# data1$CONDITION
# data1$TRTARM
# data1$WTKG_24

# no edu+income data here
data2 = read.xlsx('final_baseline_data_RL.xlsx')
# data2$ID
# data2$CONDITION
# data2$TRTARM
# data2$sex
# data2$race
# data2$AGE_WK1
# data2$HTSV1
# data2$WTKG_SV2
# data2$BMI_SV2

data3 = read.xlsx('tx data att - 082413 - all summary vars.xlsx')
# data3$ID
# data3$CONDITION
# data3$TRTARM
# data3$p2_attmu_tot_pc

# 
data4 = read.xlsx('RLITE_LABS_WTS.xlsx', sheet = 'WTS_0')
# data4$ID
# data4$CONDITION
# data4$TRTARM
# data4$WTKG_0

# 
data5 = read.xlsx('RLITE_LABS_WTS.xlsx', sheet = 'WTS_6')
# data5$ID
# data5$CONDITION
# data5$TRTARM
# data5$WTKG_6

# 
data6 = read.xlsx('RLITE_LABS_WTS.xlsx', sheet = 'WTS_24')
# data6$ID
# data6$CONDITION
# data6$TRTARM
# data6$WTKG_24

# - mediator - Phase 2 attendance (measured as a percentage) - spreadsheet 1
# - baseline confounders - mostly Table 1 from the paper and previous spreadsheet I sent
sum(abs(data2$ID!=data3$ID))
sum(abs(data2$ID!=data4$ID))
sum(abs(data2$ID!=data5$ID))
sum(abs(data2$ID!=data6$ID))

# we may add education and income 
# Just ignore the warning message! 
LITEdata = data.frame(id = data6$ID, 
                      cond = data6$CONDITION,
                      trt = as.numeric(data6$TRTARM),
                      wgt24 = as.numeric(data6$WTKG_24), 
                      wgt6 = as.numeric(data5$WTKG_6), 
                      wgt0 = as.numeric(data4$WTKG_0),
                      att = as.numeric(data3$p2_attmu_tot_pc), 
                      sex = as.numeric(data2$sex), 
                      race = as.numeric(data2$nih_race_BIQ), # white no not
                      age = as.numeric(data2$AGE_WK1), 
                      bmi = as.numeric(data2$BMI_SV2))
LITEdata$race[-which((LITEdata$race == 4 | LITEdata$race==5)==1)] = NA
LITEdata$race[which(LITEdata$race == 4)] = 0
LITEdata$race[which(LITEdata$race == 5)] = 1

# Drop variables with NAs
N = nrow(LITEdata)
sum(is.na(LITEdata$cond))
sum(is.na(LITEdata$att))
sum(is.na(LITEdata$sex))
sum(is.na(LITEdata$age))
sum(is.na(LITEdata$race))
sum(is.na(LITEdata$bmi))
sum(is.na(LITEdata$wgt0))
sum(is.na(LITEdata$wgt6))
sum(is.na(LITEdata$wgt24))

# wgt6: 56 NAs (among N = 612)
sum(is.na(LITEdata$wgt6))
# wgt24: 121 NAs (among N = 612)
sum(is.na(LITEdata$wgt24))
# both NA: 37 NAs (among N = 612)
sum(is.na(LITEdata$wgt6)*is.na(LITEdata$wgt24))
# Total NA: 140 NAs (among N = 612)
sum((is.na(LITEdata$wgt6) + is.na(LITEdata$wgt24))>=1)

# Outcome: Weight change between wgt 0 and wgt 24
LITEdata$wgtchange0_24 = LITEdata$wgt24 - LITEdata$wgt0
# Post-treatment confounder: Weight change between wgt 0 and wgt 6
LITEdata$wgtchange0_6 = LITEdata$wgt6 - LITEdata$wgt0
# Post-treatment confounder: Weight change between wgt 0 and wgt 6
LITEdata$trtmerge = ifelse(LITEdata$trt<2.5,0,1)

# Delete NAs -> No!
# Delete race which is not Black or White (612-18 = 594)
LITEdata = LITEdata[-which(apply(is.na(LITEdata), 1, sum)>0),]
LITEdata = LITEdata[-which(is.na(LITEdata$race)==1),]

# Y24 = LITEdata$wgt24
# Y06 = LITEdata$wgt6
# Y00 = LITEdata$wgt0
# 
# Y = LITEdata$wgtchange0_24
# M = LITEdata$att
# V = LITEdata$wgtchange0_6
# Z = LITEdata$trtmerge
# C = cbind(LITEdata$sex,
#           LITEdata$race,
#           LITEdata$age,
#           LITEdata$bmi)
# 
# summary(V)
# summary(M)
# summary(Y)
# 
# plot(Y)
# plot(M)
# plot(M,Y)
# 
# Y00Y06Y24M_sum = rbind(c(mean(Y00),mean(Y00[Z==1]),mean(Y00[Z==0])),
#                        c(mean(M),mean(M[Z==1]),mean(M[Z==0])),
#                        c(mean(Y06,na.rm=T),mean(Y06[Z==1],na.rm=T),mean(Y06[Z==0],na.rm=T)),
#                        c(mean(Y24,na.rm=T),mean(Y24[Z==1],na.rm=T),mean(Y24[Z==0],na.rm=T)))
# Y00Y06Y24M_sum = cbind(Y00Y06Y24M_sum, Y00Y06Y24M_sum[,2]-Y00Y06Y24M_sum[,3])
# colnames(Y00Y06Y24M_sum) = c("all","Z==1","Z==0","Z diff")
# rownames(Y00Y06Y24M_sum) = c("Y00","M","Y06","Y24")
# Y00Y06Y24M_sum
# 
# is_naY06 = is.na(Y06)*1
# ind_naY06 = which(is_naY06==1)
# n_naY06 = length(ind_naY06)
# c(mean(Y00[ind_naY06]),mean(Y00[-ind_naY06]))
# c(mean(M[ind_naY06]),mean(M[-ind_naY06]))
# 
# is_naY24 = is.na(Y24)*1
# ind_naY24 = which(is_naY24==1)
# n_naY24 = length(ind_naY24)
# c(mean(Y00[ind_naY24]),mean(Y00[-ind_naY24]))
# c(mean(M[ind_naY24]),mean(M[-ind_naY24]))
# 
# is_naY0624 = is.na(Y24)*is.na(Y06)*1
# ind_naY0624 = which(is_naY0624==1)
# n_naY0624 = length(ind_naY0624)
# c(mean(Y00[ind_naY0624]),mean(Y00[-ind_naY0624]))
# c(mean(M[ind_naY0624]),mean(M[-ind_naY0624]))
# 
# naM0 = cbind(c(sum(M==0),sum(M[ind_naY06]==0),sum(M[ind_naY24]==0),sum(M[ind_naY0624]==0)),
#              c(N, n_naY06, n_naY24, n_naY0624))
# naM0 = cbind(naM0,naM0[,1]/naM0[,2])
# colnames(naM0) = c("NA","Total","ratio")
# rownames(naM0) = c("M","M[naY06]","M[naY24]","M[naY0624]")
# naM0
# 
# hist(M)
# hist(M[ind_naY06])
# hist(M[ind_naY24])
# hist(M[ind_naY0624])
# 
# mean(M[Z==1],na.rm=T)-mean(M[Z==0],na.rm=T)
# mean(Y00[Z==1],na.rm=T)-mean(Y00[Z==0],na.rm=T)
# mean(Y06[Z==1],na.rm=T)-mean(Y06[Z==0],na.rm=T)
# mean(Y24[Z==1],na.rm=T)-mean(Y24[Z==0],na.rm=T)
# 
# # -------------------------------------------------------------------------
# # -------------------------------------------------------------------------
# # -------------------------------------------------------------------------
# round(CE_tableRhoUnif,2)[,c(1,3,4)]
# round(CE_tableRho0,2)[,c(1,3,4)]
# round(CE_tableRho3,2)[,c(1,3,4)]
# round(CE_tableRho5,2)[,c(1,3,4)]
# round(CE_tableRho8,2)[,c(1,3,4)]
# round(Resulttable_noV,2)[,c(1,3,4)]






