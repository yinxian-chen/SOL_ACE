library(geepack)
library(tidyverse)
library(nlme)
#install.packages("clubSandwich")
install.packages('WeMix')
library(clubSandwich)
sol.pheno<-read.csv("pheno_final_983obs_withcells.csv")
sol.v1.eaa<-read.csv("HCHS_epigen_age_V1_updated_100523.csv")
sol.v2.eaa<-read.csv("HCHS_epigen_age_V2_updated_100523.csv")
sol.v1.dun<-read.csv("HCHS_DunedinPACE_Visit1_updated_100523.csv")
sol.v2.dun<-read.csv("HCHS_DunedinPACE_Visit2_updated_100523.csv")
sol.ac<-read.csv("subject_annotation_2017-05-24.csv")

## Data preparation
### Only focus on the GrimAge
sol.v1.grim<-sol.v1.eaa%>%mutate(eaa_grim_v1=AgeAccelGrim,
                                 eage_grim_v1=DNAmGrimAgeBasedOnRealAge,
                                 age_v1=age)%>%
  select(c(ID,eaa_grim_v1,eage_grim_v1,age_v1))


sol.v2.grim<-sol.v2.eaa%>%mutate(eaa_grim_v2=AgeAccelGrim,
                                 eage_grim_v2=DNAmGrimAgeBasedOnRealAge,
                                 age_v2=age)%>%
  select(c(ID,eaa_grim_v2,eage_grim_v2,age_v2))

### Also include the Pace of EAA (DunedinPACE)
sol.v1.dun<-sol.v1.dun%>%mutate(dunedin_v1=DunedinPACE)%>%select(c(ID, dunedin_v1))
sol.v2.dun<-sol.v2.dun%>%mutate(dunedin_v2=DunedinPACE)%>%select(c(ID, dunedin_v2))
sol.ac.pc<-sol.ac%>%select(c(HCHS_ID,EV1,EV2,EV3,EV4,EV5))
sol.list<-list(sol.pheno,sol.v1.grim,sol.v2.grim,sol.v1.dun,sol.v2.dun)

##Merge phenotype + epigenetic age variables
sol.phe.eaa<-sol.list%>%reduce(inner_join,by="ID")

## Create a working dataset
cov_list<-c("ID","AGE","GENDER","CENTER","EDUCATION_C2","US_BORN","YRSUS"
            ,"BKGRD1_C7","NK_1","NK_2","B_1","B_2","MO_1","MO_2",
            "GR_1","GR_2","CD4_1","CD4_2","CD8_1","CD8_2","WEIGHT_NORM_OVERALL_EPIGEN")

ace<-"ACE_TOT"


eaa<-c("eaa_grim_v1","eaa_grim_v2","eage_grim_v1","eage_grim_v2","age_v1","age_v2","dunedin_v1","dunedin_v2")

sol.ana.1<-sol.phe.eaa[,c(cov_list,ace,eaa)]

table((sol.ana.1$age_v1-sol.ana.1$age_v2>0))

sol.ana.1$ID[which((sol.ana.1$age_v1-sol.ana.1$age_v2>0))] #1494646 age_v1 > age_v2

sol.ana.1$age_v1[sol.ana.1$ID==1494646]<-55
sol.ana.1$age_v2[sol.ana.1$ID==1494646]<-60
sol.ana.1$AGE[sol.ana.1$ID==1494646]<-55

sol.ana.long<-reshape(data = sol.ana.1,
                      varying = list(c("eaa_grim_v1","eaa_grim_v2"),c("eage_grim_v1","eage_grim_v2"),
                                     c("age_v1","age_v2"), c("dunedin_v1","dunedin_v2"),c("NK_1","NK_2"),
                                     c("B_1","B_2"),c("MO_1","MO_2"),c("GR_1","GR_2"),
                                     c("CD4_1","CD4_2"),c("CD8_1","CD8_2")),
                      v.names = c("eaa_grim","eage_grim","age_t","dunedin","NK","B","MO",
                                  "GR","CD4","CD8"),
                      idvar = "ID",
                      timevar = "time",
                      times = c(0,6),
                      direction = "long")
sol.ana.long<-sol.ana.long[order(sol.ana.long$ID),]
sol.ana.long$t<-rep(c(1:2),nrow(sol.ana.1))
sol.ana.long$time<-factor(sol.ana.long$time)
sol.ana.long<-sol.ana.long%>%mutate(ace_c2 = ifelse(ACE_TOT<4,0,1),
                                    age_arr_US = factor(ifelse(US_BORN==0|is.na(US_BORN),ifelse(AGE-YRSUS<18,2,3),1)),
                                    BKGRD1_C7 = factor(BKGRD1_C7),
                                    ace_c4 = factor(ifelse(ACE_TOT==0,1,
                                                    ifelse(ACE_TOT<4,2,3))),
                                    age_change = age_t-AGE)

## Function for the GEE aiming for different exposure and outcome
gee.ace<-function(data, exposure, outcome) {
  
  formula_1<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + time:", exposure))
  model.1<-gls(formula_1, 
               data = get(data),
               corr=corCompSymm(form = ~ t | ID),
               weights = varComb(varIdent(form = ~ 1 | t), varFixed(~WEIGHT_NORM_OVERALL_EPIGEN)),
               method = "REML",
               na.action = na.omit)
  
  
  print(summary(model.1))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  
  
  formula_2<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + time:", exposure))
  model.2<-gls(formula_2, 
             data = get(data),
             corr=corCompSymm(form = ~ t | ID),
             weights = varComb(varIdent(form = ~ 1 | t), varFixed(~WEIGHT_NORM_OVERALL_EPIGEN)),
             method = "REML",
             na.action = na.omit)
  print(summary(model.2))
  robust.2<-coef_test(model.2, vcov = "CR0", test = "z")
  return(list(robust.1,robust.2))
}


## Preliminary results
ana.1<-gee.ace("sol.ana.long", "ace_c2", "eaa_grim") #Binary ACE (<4 vs. >=4) and EAA GrimAge
ana.1
ana.2<-gee.ace("sol.ana.long", "ace_c2", "eage_grim") #Binary ACE (<4 vs. >=4) and Epigenetic age change using GrimAge
ana.2
ana.3<-gee.ace("sol.ana.long", "ace_c2", "dunedin") #Binary ACE (<4 vs. >=4) and Epigenetic age Pace
ana.3
ana.4<-gee.ace("sol.ana.long", "ACE_TOT", "eaa_grim") #Continuous ACE and EAA GrimAge
ana.4
ana.5<-gee.ace("sol.ana.long", "ACE_TOT", "eage_grim") #Continuous ACE and Epigenetic age change using GrimAge
ana.5
ana.6<-gee.ace("sol.ana.long", "ACE_TOT", "dunedin") #Continuous ACE and Epigenetic age Pace
ana.6
## Function for the GLM aiming for different exposure and outcome
glm.ace<-function(data, exposure, outcome) {
  
  formula_1<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + time:", exposure))
  model.1<-lme(formula_1, 
               data = get(data),
               random = ~1+time|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  print(summary(model.1))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  
  
  formula_2<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + time:", exposure))
  model.2<-lme(formula_2, 
               data = get(data),
               random = ~1+time|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  print(summary(model.2))
  robust.2<-coef_test(model.2, vcov = "CR0", test = "z")
  return(list(robust.1,robust.2))
}


## Preliminary results
ana.1.ml<-glm.ace("sol.ana.long", "ace_c2", "eaa_grim") #Binary ACE (<4 vs. >=4) and EAA GrimAge
ana.1.ml
ana.2.ml<-glm.ace("sol.ana.long", "ace_c2", "eage_grim") #Binary ACE (<4 vs. >=4) and Epigenetic age change using GrimAge
ana.2.ml
ana.3.ml<-glm.ace("sol.ana.long", "ace_c2", "dunedin") #Binary ACE (<4 vs. >=4) and Epigenetic age Pace
ana.3.ml
ana.4.ml<-glm.ace("sol.ana.long", "ACE_TOT", "eaa_grim") #Continuous ACE and EAA GrimAge
ana.4.ml
ana.5.ml<-glm.ace("sol.ana.long", "ACE_TOT", "eage_grim") #Continuous ACE and Epigenetic age change using GrimAge
ana.5.ml
ana.6.ml<-glm.ace("sol.ana.long", "ACE_TOT", "dunedin") #Continuous ACE and Epigenetic age Pace
ana.6.ml

#test.model<-lme(eaa_grim~ace_c2*time, data = sol.ana.long, random = ~1+time|ID, weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
             #   method = "REML", na.action = na.omit)
#summary(test.model)
#coef_test(test.model, vcov = "CR0", test = "z")

hist(sol.phe.eaa$ACE_TOT)

summary(sol.phe.eaa$ACE_TOT)

## Sensitive analysis
### Adjust for cell type
#### GLM
glm.ace.cell<-function(data, exposure, outcome) {
  formula.cell<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + NK + B + MO + GR + CD4 + CD8 + time:", exposure))
  model.cell<-lme(formula.cell, 
               data = get(data),
               random = ~1+time|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  print(summary(model.cell))
  robust.cell<-coef_test(model.cell, vcov = "CR0", test = "z")
  return(robust.cell)
}

ana.1.ml.cell<-glm.ace.cell("sol.ana.long", "ace_c2", "eaa_grim")
ana.1.ml.cell
ana.2.ml.cell<-glm.ace.cell("sol.ana.long", "ace_c2", "eage_grim")
ana.2.ml.cell
ana.3.ml.cell<-glm.ace.cell("sol.ana.long", "ace_c2", "dunedin")
ana.3.ml.cell
ana.4.ml.cell<-glm.ace.cell("sol.ana.long", "ACE_TOT", "eaa_grim")
ana.4.ml.cell
ana.5.ml.cell<-glm.ace.cell("sol.ana.long", "ACE_TOT", "eage_grim")
ana.5.ml.cell
ana.6.ml.cell<-glm.ace.cell("sol.ana.long", "ACE_TOT", "dunedin")
ana.6.ml.cell

#### GEE
gee.ace.cell<-function(data, exposure, outcome) {
  formula.cell<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + NK + B + MO + GR + CD4 + CD8 + time:", exposure))
  model.cell<-gls(formula.cell, 
                  data = get(data),
                  corr = corCompSymm(form = ~ t | ID),
                  weights = varComb(varIdent(form = ~ 1 | t), varFixed(~WEIGHT_NORM_OVERALL_EPIGEN)) ,
                  method = "REML",
                  na.action = na.omit)
  print(summary(model.cell))
  robust.cell<-coef_test(model.cell, vcov = "CR0", test = "z")
  return(robust.cell)
}

ana.1.cell<-gee.ace.cell("sol.ana.long", "ace_c2", "eaa_grim")
ana.1.cell
ana.2.cell<-gee.ace.cell("sol.ana.long", "ace_c2", "eage_grim")
ana.2.cell
ana.3.cell<-gee.ace.cell("sol.ana.long", "ace_c2", "dunedin")
ana.3.cell
ana.4.cell<-gee.ace.cell("sol.ana.long", "ACE_TOT", "eaa_grim")
ana.4.cell
ana.5.cell<-gee.ace.cell("sol.ana.long", "ACE_TOT", "eage_grim")
ana.5.cell
ana.6.cell<-gee.ace.cell("sol.ana.long", "ACE_TOT", "dunedin")
ana.6.cell

### Use 3-level ACE (0, 1-3, 4+)
#### GLM
ana.1.ml.4c<-glm.ace("sol.ana.long", "ace_c4", "eaa_grim") # EAA GrimAge
ana.1.ml.4c

ana.2.ml.4c<-glm.ace("sol.ana.long", "ace_c4", "eage_grim") #Epigenetic age change using GrimAge
ana.2.ml.4c

ana.3.ml.4c<-glm.ace("sol.ana.long", "ace_c4", "dunedin") # Epigenetic age Pace
ana.3.ml.4c

### GEE
ana.1.4c<-gee.ace("sol.ana.long", "ace_c4", "eaa_grim") # EAA GrimAge
ana.1.4c
ana.1

ana.2.4c<-gee.ace("sol.ana.long", "ace_c4", "eage_grim") #Epigenetic age change using GrimAge
ana.2.4c

ana.3.4c<-gee.ace("sol.ana.long", "ace_c4", "dunedin") # Epigenetic age Pace
ana.3.4c

### Unbalance design
glm.ace.ub<-function(data, exposure, outcome) {
  
  formula_1<-as.formula(paste(outcome, "~", exposure,"+ age_change + AGE + GENDER + CENTER + age_change:", exposure))
  model.1<-lme(formula_1, 
               data = get(data),
               random = ~1+age_change|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  print(summary(model.1))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  
  
  formula_2<-as.formula(paste(outcome, "~", exposure,"+ age_change + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + age_change:", exposure))
  model.2<-lme(formula_2, 
               data = get(data),
               random = ~1+age_change|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  print(summary(model.2))
  robust.2<-coef_test(model.2, vcov = "CR0", test = "z")
  return(list(robust.1,robust.2))
}

ana.1.ub<-glm.ace.ub("sol.ana.long", "ace_c2", "eaa_grim") #Binary ACE (<4 vs. >=4) and EAA GrimAge
ana.1.ub
ana.2.ub<-glm.ace.ub("sol.ana.long", "ace_c2", "eage_grim") #Binary ACE (<4 vs. >=4) and Epigenetic age change using GrimAge
ana.2.ub
ana.2.ml
ana.3.ub<-glm.ace.ub("sol.ana.long", "ace_c2", "dunedin") #Binary ACE (<4 vs. >=4) and Epigenetic age Pace
ana.3.ub
ana.4.ub<-glm.ace.ub("sol.ana.long", "ACE_TOT", "eaa_grim") #Continuous ACE and EAA GrimAge
ana.4.ub
ana.5.ub<-glm.ace.ub("sol.ana.long", "ACE_TOT", "eage_grim") #Continuous ACE and Epigenetic age change using GrimAge
ana.5.ub
ana.6.ub<-glm.ace.ub("sol.ana.long", "ACE_TOT", "dunedin") #Continuous ACE and Epigenetic age Pace
ana.6.ub




