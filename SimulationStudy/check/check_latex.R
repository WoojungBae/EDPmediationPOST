library(knitr)
library(kableExtra)

resultdataframe = data.frame(gamma = rep(paste0("S",1:12), each = 3),
                             rep(c("NIE", "NDE", "ATE"),12),
                             result_EDP)
colnames(resultdataframe) = NULL
rownames(resultdataframe) = NULL

kable(resultdataframe, digits = 3, "latex", booktabs = T, escape = F,
      col.names = c("Scn", "CE", "True", rep(c("Bias","MSE","CIl","Cover"),4)),
      caption = "Simulation results") %>%
  add_header_above(c(" "=3,"250"=4,"500"=4,"1000"=4,"2500"=4)) %>%
  kable_styling() %>%
  collapse_rows(columns = 1:3, latex_hline = "major", valign = "middle") 

