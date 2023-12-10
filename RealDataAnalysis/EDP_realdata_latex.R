library(knitr)
library(kableExtra)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
ResulttableRLTmain = rbind(cbind(CE_tableRhoTri01[,c(1,3,4)],
                                 CE_tableRhoUnif01[,c(1,3,4)],
                                 CE_tableRho0[,c(1,3,4)]),
                           cbind(CE_tableRho3[,c(1,3,4)],
                                 CE_tableRho5[,c(1,3,4)],
                                 CE_tableRho8[,c(1,3,4)]))
ResulttableRLTmain = round(ResulttableRLTmain,2);ResulttableRLTmain

ResulttableRLTmain = data.frame(rep(c("NIE", "NDE", "ATE"),2),
                                ResulttableRLTmain)
colnames(ResulttableRLTmain) = NULL
rownames(ResulttableRLTmain) = NULL

kable(ResulttableRLTmain, digits = 3, "latex", booktabs = T, escape = F,
      col.names = c("CE", rep(c("Est.", "95% CI", ""),3))) %>%
  kable_styling() %>%
  collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle")

# -------------------------------------------------------------------------
ResulttableRLTsupp = rbind(cbind(CE_tableRhoTrim10[,c(1,3,4)],
                                 CE_tableRhoUnifm10[,c(1,3,4)],
                                 CE_tableRhoUnifm11[,c(1,3,4)]),
                           cbind(CE_tableRhoM8[,c(1,3,4)],
                                 CE_tableRhoM5[,c(1,3,4)],
                                 CE_tableRhoM3[,c(1,3,4)]))
ResulttableRLTsupp = round(ResulttableRLTsupp,2);ResulttableRLTsupp

ResulttableRLTsupp = data.frame(rep(c("NIE", "NDE", "ATE"),2),
                                ResulttableRLTsupp)
colnames(ResulttableRLTsupp) = NULL
rownames(ResulttableRLTsupp) = NULL

kable(ResulttableRLTsupp, digits = 3, "latex", booktabs = T, escape = F,
      col.names = c("CE", rep(c("Est.", "95% CI", ""),3))) %>%
  kable_styling() %>%
  collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
ResulttableRLTmain_cond = rbind(cbind(CE_tableRhoTri01_black[,c(1,3,4)],
                                      CE_tableRhoUnif01_black[,c(1,3,4)],
                                      CE_tableRho0_black[,c(1,3,4)]),
                                cbind(CE_tableRhoTri01_white[,c(1,3,4)],
                                      CE_tableRhoUnif01_white[,c(1,3,4)],
                                      CE_tableRho0_white[,c(1,3,4)]),
                                cbind(CE_tableRho03_black[,c(1,3,4)],
                                      CE_tableRho05_black[,c(1,3,4)],
                                      CE_tableRho08_black[,c(1,3,4)]),
                                cbind(CE_tableRho03_white[,c(1,3,4)],
                                      CE_tableRho05_white[,c(1,3,4)],
                                      CE_tableRho08_white[,c(1,3,4)]))

ResulttableRLTmain_cond = round(ResulttableRLTmain_cond,2);ResulttableRLTmain_cond
ResulttableRLTmain_cond = data.frame(rep(c(rep("B",3),rep("W",3)),2),
                                     rep(c("NIE", "NDE", "ATE"),4),
                                     ResulttableRLTmain_cond)
colnames(ResulttableRLTmain_cond) = NULL
rownames(ResulttableRLTmain_cond) = NULL

kable(ResulttableRLTmain_cond, digits = 3, "latex", booktabs = T, escape = F,
      col.names = c("race", "CE", rep(c("Est.", "95% CI", ""),3))) %>%
  kable_styling() %>%
  collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle")

# -------------------------------------------------------------------------
ResulttableRLTsupp_cond = rbind(cbind(CE_tableRhoTrim10_black[,c(1,3,4)],
                                      CE_tableRhoUnifm10_black[,c(1,3,4)],
                                      CE_tableRhoUnifm11_black[,c(1,3,4)]),
                                cbind(CE_tableRhoTrim10_white[,c(1,3,4)],
                                      CE_tableRhoUnifm10_white[,c(1,3,4)],
                                      CE_tableRhoUnifm11_white[,c(1,3,4)]),
                                cbind(CE_tableRhoM3_black[,c(1,3,4)],
                                      CE_tableRhoM5_black[,c(1,3,4)],
                                      CE_tableRhoM8_black[,c(1,3,4)]),
                                cbind(CE_tableRhoM3_white[,c(1,3,4)],
                                      CE_tableRhoM5_white[,c(1,3,4)],
                                      CE_tableRhoM8_white[,c(1,3,4)]))
ResulttableRLTsupp_cond = round(ResulttableRLTsupp_cond,2);ResulttableRLTsupp_cond

ResulttableRLTsupp_cond = data.frame(rep(c(rep("B",3),rep("W",3)),2),
                                     rep(c("NIE", "NDE", "ATE"),4),
                                     ResulttableRLTsupp_cond)
colnames(ResulttableRLTsupp_cond) = NULL
rownames(ResulttableRLTsupp_cond) = NULL

kable(ResulttableRLTsupp_cond, digits = 3, "latex", booktabs = T, escape = F,
      col.names = c("race", "CE", rep(c("Est.", "95% CI", ""),3))) %>%
  kable_styling() %>%
  collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle")

















# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Results
# print("EDP with a post-treatment confounder & rho = triangle(-1,0,-1); (lower, upper, mode)")
# CE_tableRhoTrim10
# print("EDP with a post-treatment confounder & rho = unif(-1,0)")
# CE_tableRhoUnifm10
# print("EDP with a post-treatment confounder & rho = unif(-1,1)")
# CE_tableRhoUnifm11
# print("EDP with a post-treatment confounder & rho = unif(0,1)")
# CE_tableRhoUnif01
# print("EDP with a post-treatment confounder & rho = triangle(0,1,1); (lower, upper, mode)")
# CE_tableRhoTri01
# print("EDP with a post-treatment confounder & rho = -0.8")
# CE_tableRhoM8
# print("EDP with a post-treatment confounder & rho = -0.5")
# CE_tableRhoM5
# print("EDP with a post-treatment confounder & rho = -0.3")
# CE_tableRhoM3
# print("EDP with a post-treatment confounder & rho = 0")
# CE_tableRho0
# print("EDP with a post-treatment confounder & rho = 0.3 (true)")
# CE_tableRho3
# print("EDP with a post-treatment confounder & rho = 0.5")
# CE_tableRho5
# print("EDP with a post-treatment confounder & rho = 0.8")
# CE_tableRho8
# print("EDP without a post-treatment confounder")
# Resulttable_noV

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
## Main
# print("EDP with a post-treatment confounder & rho = unif(-1,0)")
# CE_tableRhoUnifm10
# print("EDP with a post-treatment confounder & rho = unif(-1,1)")
# CE_tableRhoUnifm11
# print("EDP with a post-treatment confounder & rho = triangle(0,1,1); (lower, upper, mode)")
# CE_tableRhoTri01
# print("EDP with a post-treatment confounder & rho = -0.8")
# CE_tableRhoM8
# print("EDP with a post-treatment confounder & rho = 0")
# CE_tableRho0
# print("EDP with a post-treatment confounder & rho = 0.5")
# CE_tableRho5

# # -------------------------------------------------------------------------
# ResulttableRLTmain = rbind(cbind(CE_tableRhoUnifm10[,c(1,3,4)],
#                                  CE_tableRhoUnifm11[,c(1,3,4)],
#                                  CE_tableRhoTri01[,c(1,3,4)]),
#                            cbind(CE_tableRhoM8[,c(1,3,4)],
#                                  CE_tableRho0[,c(1,3,4)],
#                                  CE_tableRho5[,c(1,3,4)]))
# ResulttableRLTmain = round(ResulttableRLTmain,2);ResulttableRLTmain
# 
# ResulttableRLTmain = data.frame(rep(c("NIE", "NDE", "ATE"),2),
#                                 ResulttableRLTmain)
# colnames(ResulttableRLTmain) = NULL
# rownames(ResulttableRLTmain) = NULL
# 
# kable(ResulttableRLTmain, digits = 3, "latex", booktabs = T, escape = F,
#       col.names = c("CE", rep(c("Est.", "95% CI", ""),3))) %>%
#   kable_styling() %>%
#   collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle")

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
## Supp
# print("EDP with a post-treatment confounder & rho = triangle(-1,0,-1); (lower, upper, mode)")
# CE_tableRhoTrim10
# print("EDP with a post-treatment confounder & rho = unif(0,1)")
# CE_tableRhoUnif01
# print("EDP with a post-treatment confounder & rho = -0.5")
# CE_tableRhoM5
# print("EDP with a post-treatment confounder & rho = -0.3")
# CE_tableRhoM3
# print("EDP with a post-treatment confounder & rho = 0.3 (true)")
# CE_tableRho3
# print("EDP with a post-treatment confounder & rho = 0.8")
# CE_tableRho8
# print("EDP without a post-treatment confounder")
# Resulttable_noV

# # -------------------------------------------------------------------------
# ResulttableRLTsupp = rbind(cbind(CE_tableRhoTrim10[,c(1,3,4)],
#                                  CE_tableRhoUnif01[,c(1,3,4)],
#                                  CE_tableRhoM5[,c(1,3,4)]),
#                            cbind(CE_tableRhoM3[,c(1,3,4)],
#                                  CE_tableRho3[,c(1,3,4)],
#                                  CE_tableRho8[,c(1,3,4)]))
# ResulttableRLTsupp = round(ResulttableRLTsupp,2);ResulttableRLTsupp
# 
# ResulttableRLTsupp = data.frame(rep(c("NIE", "NDE", "ATE"),2),
#                                 ResulttableRLTsupp)
# colnames(ResulttableRLTsupp) = NULL
# rownames(ResulttableRLTsupp) = NULL
# 
# kable(ResulttableRLTsupp, digits = 3, "latex", booktabs = T, escape = F,
#       col.names = c("CE", rep(c("Est.", "95% CI", ""),3))) %>%
#   kable_styling() %>%
#   collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle")

# # -------------------------------------------------------------------------
# round(Resulttable_noV,2)[,c(1,3,4)]

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# conditional on race: black or white
# print("EDP with a post-treatment confounder & rho = triangle(-1,0,-1); (lower, upper, mode)")
# CE_tableRhoTrim10_black
# CE_tableRhoTrim10_white
# 
# print("EDP with a post-treatment confounder & rho = unif(-1,0)")
# CE_tableRhoUnifm10_black
# CE_tableRhoUnifm10_white
# 
# print("EDP with a post-treatment confounder & rho = unif(-1,1)")
# CE_tableRhoUnifm11_black
# CE_tableRhoUnifm11_white
# 
# print("EDP with a post-treatment confounder & rho = unif(0,1)")
# CE_tableRhoUnif01_black
# CE_tableRhoUnif01_white
# 
# print("EDP with a post-treatment confounder & rho = triangle(0,1,1); (lower, upper, mode)")
# CE_tableRhoTri01_black
# CE_tableRhoTri01_white
# 
# print("EDP with a post-treatment confounder & rho = -0.8")
# CE_tableRhoM8_black
# CE_tableRhoM8_white
# 
# print("EDP with a post-treatment confounder & rho = -0.5")
# CE_tableRhoM5_black
# CE_tableRhoM5_white
# 
# print("EDP with a post-treatment confounder & rho = -0.3")
# CE_tableRhoM3_black
# CE_tableRhoM3_white
# 
# print("EDP with a post-treatment confounder & rho = 0")
# CE_tableRho0_black
# CE_tableRho0_white
# 
# print("EDP with a post-treatment confounder & rho = 0.3 (true)")
# CE_tableRho03_black
# CE_tableRho03_white
# 
# print("EDP with a post-treatment confounder & rho = 0.5")
# CE_tableRho05_black
# CE_tableRho05_white
# 
# print("EDP with a post-treatment confounder & rho = 0.8")
# CE_tableRho08_black
# CE_tableRho08_white

# # -------------------------------------------------------------------------
# ResulttableRLTmain_cond = rbind(cbind(CE_tableRhoUnifm10_black[,c(1,3,4)],
#                                       CE_tableRhoUnifm11_black[,c(1,3,4)],
#                                       CE_tableRhoTri01_black[,c(1,3,4)]),
#                                 cbind(CE_tableRhoUnifm10_white[,c(1,3,4)],
#                                       CE_tableRhoUnifm11_white[,c(1,3,4)],
#                                       CE_tableRhoTri01_white[,c(1,3,4)]),
#                                 cbind(CE_tableRhoM8_black[,c(1,3,4)],
#                                       CE_tableRho0_black[,c(1,3,4)],
#                                       CE_tableRho05_black[,c(1,3,4)]),
#                                 cbind(CE_tableRhoM8_white[,c(1,3,4)],
#                                       CE_tableRho0_white[,c(1,3,4)],
#                                       CE_tableRho05_white[,c(1,3,4)]))
# ResulttableRLTmain_cond = round(ResulttableRLTmain_cond,2);ResulttableRLTmain_cond
# ResulttableRLTmain_cond = data.frame(rep(c(rep("B",3),rep("W",3)),2),
#                                      rep(c("NIE", "NDE", "ATE"),4),
#                                       ResulttableRLTmain_cond)
# colnames(ResulttableRLTmain_cond) = NULL
# rownames(ResulttableRLTmain_cond) = NULL
# 
# kable(ResulttableRLTmain_cond, digits = 3, "latex", booktabs = T, escape = F,
#       col.names = c("race", "CE", rep(c("Est.", "95% CI", ""),3))) %>%
#   kable_styling() %>%
#   collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle")
# 
# # -------------------------------------------------------------------------
# ResulttableRLTsupp_cond = rbind(cbind(CE_tableRhoTrim10_black[,c(1,3,4)],
#                                       CE_tableRhoUnif01_black[,c(1,3,4)],
#                                       CE_tableRhoM5_black[,c(1,3,4)]),
#                                 cbind(CE_tableRhoTrim10_white[,c(1,3,4)],
#                                       CE_tableRhoUnif01_white[,c(1,3,4)],
#                                       CE_tableRhoM5_white[,c(1,3,4)]),
#                                 cbind(CE_tableRhoM3_black[,c(1,3,4)],
#                                       CE_tableRho03_black[,c(1,3,4)],
#                                       CE_tableRho08_black[,c(1,3,4)]),
#                                 cbind(CE_tableRhoM3_white[,c(1,3,4)],
#                                       CE_tableRho03_white[,c(1,3,4)],
#                                       CE_tableRho08_white[,c(1,3,4)]))
# ResulttableRLTsupp_cond = round(ResulttableRLTsupp_cond,2);ResulttableRLTsupp_cond
# 
# ResulttableRLTsupp_cond = data.frame(rep(c(rep("B",3),rep("W",3)),2),
#                                      rep(c("NIE", "NDE", "ATE"),4),
#                                      ResulttableRLTsupp_cond)
# colnames(ResulttableRLTsupp_cond) = NULL
# rownames(ResulttableRLTsupp_cond) = NULL
# 
# kable(ResulttableRLTsupp_cond, digits = 3, "latex", booktabs = T, escape = F,
#       col.names = c("race", "CE", rep(c("Est.", "95% CI", ""),3))) %>%
#   kable_styling() %>%
#   collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle")
