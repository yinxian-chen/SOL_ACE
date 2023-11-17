library(nlme)
library(tidyverse)
library(clubSandwich)

### Convert into the long format for longitudinal analysis
#### Need to replace and dataset and variable names based on your own data
sol.ana.long<-reshape(data = sol.ana.1,
                      varying = list(c("eaa_grim_v1","eaa_grim_v2"),c("eage_grim_v1","eage_grim_v2"),
                                     c("age_v1","age_v2"), c("dunedin_v1","dunedin_v2"),c("NK_1","NK_2"),
                                     c("B_1","B_2"),c("MO_1","MO_2"),c("GR_1","GR_2"),
                                     c("CD4_1","CD4_2"),c("CD8_1","CD8_2")),
                      v.names = c("eaa_grim","eage_grim","age_t","dunedin","NK","B","MO",
                                  "GR","CD4","CD8"),
                      idvar = "ID",
                      timevar = "time",
                      times = c(0,1), 
                      direction = "long")
sol.ana.long<-sol.ana.long[order(sol.ana.long$ID),]
sol.ana.long<-sol.ana.long%>%mutate(age_change = age_t-AGE,
  GENDER = factor(ifelse(GENDER=="M",0,1)))%>%select(-c(age_t,time))


### Outcome definition

outcome_list<-c("eaa_grim","eaa_grim","dunedin")

### Exposure definition

exposure_list<-c("var1","var2","var3"#....
                 )

### GlMM for the primary analysis
glm.ace.ub<-function(data, exposure, outcome) {
  formula.1<-as.formula(paste(outcome, "~", exposure,"+ age_change + AGE + GENDER + fa_edu_3c + ma_edu_3c + age_change:", exposure))
  model.1<-lme(formula.1, 
               data = data,
               random = ~1+age_change|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  #print(summary(model.cell))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  return(robust.1)
}

### Sensitivity analysis for cell proportion
glm.ace.cell<-function(data, exposure, outcome) {
  formula.cell<-as.formula(paste(outcome, "~", exposure,"+ age_change + AGE + GENDER + fa_edu_3c + ma_edu_3c + NK + B + MO + GR + CD4 + CD8 + age_change:", exposure))
  model.cell<-lme(formula.cell, 
                  data = data,
                  random = ~1+age_change|ID,
                  weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
                  method = "REML",
                  na.action = na.omit)
  #print(summary(model.cell))
  robust.cell<-coef_test(model.cell, vcov = "CR0", test = "z")
  return(robust.cell)
}

### Sensitivity analysis for ancestry

glm.ace.ac<-function(data, exposure, outcome) {
  formula.ac<-as.formula(paste(outcome, "~", exposure,"+ age_change + AGE + GENDER + fa_edu_3c + ma_edu_3c + EV1 + EV2 + EV3 + EV4 + EV5 + age_change:", exposure))
  model.ac<-lme(formula.ac, 
                data = data,
                random = ~1+age_change|ID,
                weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
                method = "REML",
                na.action = na.omit)
  #print(summary(model.ac))
  robust.ac<-coef_test(model.ac, vcov = "CR0", test = "z")
  return(robust.ac)
}

### EMM by Nativity

#### Create nativity-specific variable for the exposure of interest
#### The effect size for these variables will represent the nativity-specific estimate
#### given the orginal model of "outcome~exposure + US_BORN + exposure*US_BORN"
data<-data%>%mutate(nusb_ace_c2 = ace_c2*I(US_BORN==0),
                    usb_ace_c2 = ace_c2*I(US_BORN==1),
                    nusb_ACE_TOT = ACE_TOT*I(US_BORN==0),
                    usb_ACE_TOT = ACE_TOT*I(US_BORN==1))

#### Use to compute the nativity-specific estimate
glm.ace.usb<-function(data, exposure.1, exposure.2, outcome) {
  
  formula_1<-as.formula(paste(outcome, "~", exposure.1, "+", exposure.2, "+ age_change + AGE + GENDER + US_BORN + fa_edu_3c + ma_edu_3c + age_change:", exposure.1, "+ age_change:", exposure.2))
  model.1<-lme(formula_1, 
               data = data,
               random = ~1+age_change|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  #print(summary(model.1))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  
  return(robust.1)
}

#### Use to compute the p-value of interaction
glm.ace.int.usb<-function(data, exposure, outcome) {
  
  formula_1<-as.formula(paste(outcome, "~", exposure, "+ age_change + AGE + GENDER + US_BORN + fa_edu_3c + ma_edu_3c + US_BORN:", exposure, "+ age_change:", exposure, "+ age_change:US_BORN:", exposure))
  model.1<-lme(formula_1, 
               data = data,
               random = ~1+age_change|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  #print(summary(model.1))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  
  return(robust.1)
}

### One example of my analysis

cca.list<-lapply(exposure_list,function(x) {
  result.matrix.x<-matrix(NA, nrow = 9, ncol = 5)
  result.matrix.x[,5]<-c(rep("eaa_grim",3),rep("eage_grim",3),rep("dunedin",3))
  for (i in 1:3) {
    res<-glm.ace.ub(sol.ana.long.cca,x,outcome_list[i])
    res<-res[c(2:3,10),]
    result.matrix.x[(3*i-2):(3*i),1]<-res$beta
    result.matrix.x[(3*i-2):(3*i),2]<-res$p_z
    result.matrix.x[(3*i-2):(3*i),3]<-res$beta - 1.96*res$SE
    result.matrix.x[(3*i-2):(3*i),4]<-res$beta + 1.96*res$SE
  }
  colnames(result.matrix.x)<-c("beta","p-value","lower","upper","outcome")
  result.matrix.x<-as.data.frame(result.matrix.x)%>%mutate_at(vars(beta:upper), as.numeric)
  return(result.matrix.x)
})

