library(geepack)
library(tidyverse)
library(nlme)
library(clubSandwich)
library(kableExtra)

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

summary(sol.ana.1$age_v1-sol.ana.1$AGE)
summary(sol.ana.1$age_v2-sol.ana.1$AGE)
table(sol.ana.1$age_v2-sol.ana.1$AGE)
sol.ana.1$ID[which((sol.ana.1$age_v2-sol.ana.1$AGE<5|sol.ana.1$age_v2-sol.ana.1$AGE>8))]


#table((sol.ana.1$age_v1-sol.ana.1$age_v2>0))

#sol.ana.1$ID[which((sol.ana.1$age_v1-sol.ana.1$age_v2>0))] #1494646 age_v1 > age_v2

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
  
  
  #print(summary(model.1))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  
  
  formula_2<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + time:", exposure))
  model.2<-gls(formula_2, 
             data = get(data),
             corr=corCompSymm(form = ~ t | ID),
             weights = varComb(varIdent(form = ~ 1 | t), varFixed(~WEIGHT_NORM_OVERALL_EPIGEN)),
             method = "REML",
             na.action = na.omit)
  #print(summary(model.2))
  robust.2<-coef_test(model.2, vcov = "CR0", test = "z")
  return(list(robust.1,robust.2))
}

#### Function of summary table
result.summary<-function(data, exposure, covar="normal", model) {
  
  exposure<-exposure
  covar<-covar
  
  if(exposure =="ace_c4"){
    gee.summary.1<-matrix(NA,nrow = 10, ncol = 4)
    for(i in c(1:10)) {
      indices<-c(2,3,4,10,11,2,3,4,13,14)
      if (i<=5) {
        gee.summary.1[i,1]<-data[[1]]$beta[indices[i]]
        gee.summary.1[i,2]<-data[[1]]$beta[indices[i]]-1.96*data[[1]]$SE[indices[i]]
        gee.summary.1[i,3]<-data[[1]]$beta[indices[i]]+1.96*data[[1]]$SE[indices[i]]
        gee.summary.1[i,4]<-data[[1]]$p_z[indices[i]]
      }
      
      else {
        gee.summary.1[i,1]<-data[[2]]$beta[indices[i]]
        gee.summary.1[i,2]<-data[[2]]$beta[indices[i]]-1.96*data[[2]]$SE[indices[i]]
        gee.summary.1[i,3]<-data[[2]]$beta[indices[i]]+1.96*data[[2]]$SE[indices[i]]
        gee.summary.1[i,4]<-data[[2]]$p_z[indices[i]]
      }
    }
    rownames(gee.summary.1)<-rep(c("1-3 ACEs",">=4 ACEs","Time","1-3 ACEs*Time",">=4 ACEs*Time"),2)
  }
  
  else{
    if(covar=="normal") {
      gee.summary.1<-matrix(NA,nrow = 6, ncol = 4)
      for (i in c(1:6)) {
        indices<-c(2,3,9,2,3,12)
        if (i<=3) {
          
          gee.summary.1[i,1]<-data[[1]]$beta[indices[i]]
          gee.summary.1[i,2]<-data[[1]]$beta[indices[i]]-1.96*data[[1]]$SE[indices[i]]
          gee.summary.1[i,3]<-data[[1]]$beta[indices[i]]+1.96*data[[1]]$SE[indices[i]]
          gee.summary.1[i,4]<-data[[1]]$p_z[indices[i]]
        }
        
        else {
          gee.summary.1[i,1]<-data[[2]]$beta[indices[i]]
          gee.summary.1[i,2]<-data[[2]]$beta[indices[i]]-1.96*data[[2]]$SE[indices[i]]
          gee.summary.1[i,3]<-data[[2]]$beta[indices[i]]+1.96*data[[2]]$SE[indices[i]]
          gee.summary.1[i,4]<-data[[2]]$p_z[indices[i]]
        }
      }
      
      if(exposure=="ace_c2"){
        rownames(gee.summary.1)<-rep(c(">=4 ACEs", "Time", ">=4 ACEs*Time"),2)
      }
      else {
        rownames(gee.summary.1)<-rep(c("ACEs total score", "Time", "ACEs total score*Time"),2)
      }
    }
    
    if(covar=="ancenstry") {
      gee.summary.1<-matrix(NA,nrow = 3, ncol = 4)
      for(i in c(1:3)) {
        indices<-c(2,3,17)
        gee.summary.1[i,1]<-data$beta[indices[i]]
        gee.summary.1[i,2]<-data$beta[indices[i]]-1.96*data$SE[indices[i]]
        gee.summary.1[i,3]<-data$beta[indices[i]]+1.96*data$SE[indices[i]]
        gee.summary.1[i,4]<-data$p_z[indices[i]]
      }
      if(exposure=="ace_c2"){
        rownames(gee.summary.1)<-c(">=4 ACEs", "Time", ">=4 ACEs*Time")
      }
      else {
        rownames(gee.summary.1)<-c("ACEs total score", "Time", "ACEs total score*Time")
      }
    }
    
    if(covar=="cell_type") {
      gee.summary.1<-matrix(NA,nrow = 3, ncol = 4)
      for(i in c(1:3)) {
        indices<-c(2,3,18)
        gee.summary.1[i,1]<-data$beta[indices[i]]
        gee.summary.1[i,2]<-data$beta[indices[i]]-1.96*data$SE[indices[i]]
        gee.summary.1[i,3]<-data$beta[indices[i]]+1.96*data$SE[indices[i]]
        gee.summary.1[i,4]<-data$p_z[indices[i]]
      }
      if(exposure=="ace_c2"){
        rownames(gee.summary.1)<-c(">=4 ACEs", "Time", ">=4 ACEs*Time")
      }
      else {
        rownames(gee.summary.1)<-c("ACEs total score", "Time", "ACEs total score*Time")
      }
    }
  }
  if (model=="GEE") {
    colnames(gee.summary.1)<-c("beta.gee","ll.gee","ul.gee","p-value.gee")
  }
  else {
    colnames(gee.summary.1)<-c("beta.lm","ll.lm","ul.lm","p-value.lm")
  }
  return(gee.summary.1)
}

## Preliminary results
ana.1<-gee.ace("sol.ana.long", "ace_c2", "eaa_grim") #Binary ACE (<4 vs. >=4) and EAA GrimAge
#ana.1
ana.2<-gee.ace("sol.ana.long", "ace_c2", "eage_grim") #Binary ACE (<4 vs. >=4) and Epigenetic age change using GrimAge
#ana.2
ana.3<-gee.ace("sol.ana.long", "ace_c2", "dunedin") #Binary ACE (<4 vs. >=4) and Epigenetic age Pace
#ana.3
ana.4<-gee.ace("sol.ana.long", "ACE_TOT", "eaa_grim") #Continuous ACE and EAA GrimAge
#ana.4
ana.5<-gee.ace("sol.ana.long", "ACE_TOT", "eage_grim") #Continuous ACE and Epigenetic age change using GrimAge
#ana.5
ana.6<-gee.ace("sol.ana.long", "ACE_TOT", "dunedin") #Continuous ACE and Epigenetic age Pace
#ana.6


summary.1<-result.summary(ana.1, "ace_c2", model = "GEE")
summary.2<-result.summary(ana.2, "ace_c2", model = "GEE")
summary.3<-result.summary(ana.3, "ace_c2", model = "GEE")
summary.4<-result.summary(ana.4, "ACE_TOT", model = "GEE")
summary.5<-result.summary(ana.5, "ACE_TOT", model = "GEE")
summary.6<-result.summary(ana.6, "ACE_TOT", model = "GEE")

## Function for the GLM aiming for different exposure and outcome
glm.ace<-function(data, exposure, outcome) {
  
  formula_1<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + time:", exposure))
  model.1<-lme(formula_1, 
               data = get(data),
               random = ~1+time|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  #print(summary(model.1))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  
  
  formula_2<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + time:", exposure))
  model.2<-lme(formula_2, 
               data = get(data),
               random = ~1+time|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  #print(summary(model.2))
  robust.2<-coef_test(model.2, vcov = "CR0", test = "z")
  return(list(robust.1,robust.2))
}


## Preliminary results
ana.1.ml<-glm.ace("sol.ana.long", "ace_c2", "eaa_grim") #Binary ACE (<4 vs. >=4) and EAA GrimAge
#ana.1.ml
ana.2.ml<-glm.ace("sol.ana.long", "ace_c2", "eage_grim") #Binary ACE (<4 vs. >=4) and Epigenetic age change using GrimAge
#ana.2.ml
ana.3.ml<-glm.ace("sol.ana.long", "ace_c2", "dunedin") #Binary ACE (<4 vs. >=4) and Epigenetic age Pace
#ana.3.ml
ana.4.ml<-glm.ace("sol.ana.long", "ACE_TOT", "eaa_grim") #Continuous ACE and EAA GrimAge
#ana.4.ml
ana.5.ml<-glm.ace("sol.ana.long", "ACE_TOT", "eage_grim") #Continuous ACE and Epigenetic age change using GrimAge
#ana.5.ml
ana.6.ml<-glm.ace("sol.ana.long", "ACE_TOT", "dunedin") #Continuous ACE and Epigenetic age Pace
#ana.6.ml

summary.1.ml<-result.summary(ana.1.ml, "ace_c2", model="LM")
summary.2.ml<-result.summary(ana.2.ml, "ace_c2",model="LM" )
summary.3.ml<-result.summary(ana.3.ml, "ace_c2", model="LM")
summary.4.ml<-result.summary(ana.4.ml, "ACE_TOT", model="LM")
summary.5.ml<-result.summary(ana.5.ml, "ACE_TOT",model="LM")
summary.6.ml<-result.summary(ana.6.ml, "ACE_TOT",model="LM")

#test.model<-lme(eaa_grim~ace_c2*time, data = sol.ana.long, random = ~1+time|ID, weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
             #   method = "REML", na.action = na.omit)
#summary(test.model)
#coef_test(test.model, vcov = "CR0", test = "z")

#hist(sol.phe.eaa$ACE_TOT)

#summary(sol.phe.eaa$ACE_TOT)

#### Summary table for GEE and LM

test<-rbind(cbind(summary.1,summary.1.ml),
            cbind(summary.4,summary.4.ml),
            cbind(summary.2,summary.2.ml),
            cbind(summary.5,summary.5.ml),
            cbind(summary.3,summary.3.ml),
            cbind(summary.6,summary.6.ml))

kable(test, booktabs = T)


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
  #print(summary(model.cell))
  robust.cell<-coef_test(model.cell, vcov = "CR0", test = "z")
  return(robust.cell)
}

ana.1.ml.cell<-glm.ace.cell("sol.ana.long", "ace_c2", "eaa_grim")
#ana.1.ml.cell
ana.2.ml.cell<-glm.ace.cell("sol.ana.long", "ace_c2", "eage_grim")
#ana.2.ml.cell
ana.3.ml.cell<-glm.ace.cell("sol.ana.long", "ace_c2", "dunedin")
#ana.3.ml.cell
ana.4.ml.cell<-glm.ace.cell("sol.ana.long", "ACE_TOT", "eaa_grim")
#ana.4.ml.cell
ana.5.ml.cell<-glm.ace.cell("sol.ana.long", "ACE_TOT", "eage_grim")
#ana.5.ml.cell
ana.6.ml.cell<-glm.ace.cell("sol.ana.long", "ACE_TOT", "dunedin")
#ana.6.ml.cell

summary.1.ml.cell<-result.summary(ana.1.ml.cell, "ace_c2", "cell_type", "LM")
summary.2.ml.cell<-result.summary(ana.2.ml.cell, "ace_c2", "cell_type", "LM")
summary.3.ml.cell<-result.summary(ana.3.ml.cell, "ace_c2", "cell_type", "LM")
summary.4.ml.cell<-result.summary(ana.4.ml.cell, "ACE_TOT","cell_type", "LM")
summary.5.ml.cell<-result.summary(ana.5.ml.cell, "ACE_TOT","cell_type", "LM")
summary.6.ml.cell<-result.summary(ana.6.ml.cell, "ACE_TOT","cell_type", "LM")

#### GEE
gee.ace.cell<-function(data, exposure, outcome) {
  formula.cell<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + NK + B + MO + GR + CD4 + CD8 + time:", exposure))
  model.cell<-gls(formula.cell, 
                  data = get(data),
                  corr = corCompSymm(form = ~ t | ID),
                  weights = varComb(varIdent(form = ~ 1 | t), varFixed(~WEIGHT_NORM_OVERALL_EPIGEN)) ,
                  method = "REML",
                  na.action = na.omit)
 #print(summary(model.cell))
  robust.cell<-coef_test(model.cell, vcov = "CR0", test = "z")
  return(robust.cell)
}

ana.1.cell<-gee.ace.cell("sol.ana.long", "ace_c2", "eaa_grim")
#ana.1.cell
ana.2.cell<-gee.ace.cell("sol.ana.long", "ace_c2", "eage_grim")
#ana.2.cell
ana.3.cell<-gee.ace.cell("sol.ana.long", "ace_c2", "dunedin")
#ana.3.cell
ana.4.cell<-gee.ace.cell("sol.ana.long", "ACE_TOT", "eaa_grim")
#ana.4.cell
ana.5.cell<-gee.ace.cell("sol.ana.long", "ACE_TOT", "eage_grim")
#ana.5.cell
ana.6.cell<-gee.ace.cell("sol.ana.long", "ACE_TOT", "dunedin")
#ana.6.cell

summary.1.cell<-result.summary(ana.1.cell, "ace_c2", "cell_type", "GEE")
summary.2.cell<-result.summary(ana.2.cell, "ace_c2", "cell_type", "GEE")
summary.3.cell<-result.summary(ana.3.cell, "ace_c2", "cell_type", "GEE")
summary.4.cell<-result.summary(ana.4.cell, "ACE_TOT","cell_type", "GEE")
summary.5.cell<-result.summary(ana.5.cell, "ACE_TOT","cell_type", "GEE")
summary.6.cell<-result.summary(ana.6.cell, "ACE_TOT","cell_type", "GEE")


test.1<-rbind(cbind(summary.1.cell,summary.1.ml.cell),
            cbind(summary.4.cell,summary.4.ml.cell),
            cbind(summary.2.cell,summary.2.ml.cell),
            cbind(summary.5.cell,summary.5.ml.cell),
            cbind(summary.3.cell,summary.3.ml.cell),
            cbind(summary.6.cell,summary.6.ml.cell))

### Use 3-level ACE (0, 1-3, 4+)
#### GLM
ana.1.ml.4c<-glm.ace("sol.ana.long", "ace_c4", "eaa_grim") # EAA GrimAge
ana.1.ml.4c

ana.2.ml.4c<-glm.ace("sol.ana.long", "ace_c4", "eage_grim") #Epigenetic age change using GrimAge
ana.2.ml.4c

ana.3.ml.4c<-glm.ace("sol.ana.long", "ace_c4", "dunedin") # Epigenetic age Pace
ana.3.ml.4c

summary.1.4c.ml<-result.summary(ana.1.ml.4c, "ace_c4", model = "LM")
summary.2.4c.ml<-result.summary(ana.2.ml.4c, "ace_c4", model = "LM")
summary.3.4c.ml<-result.summary(ana.3.ml.4c, "ace_c4", model = "LM")


### GEE
ana.1.4c<-gee.ace("sol.ana.long", "ace_c4", "eaa_grim") # EAA GrimAge
ana.1.4c
ana.1

ana.2.4c<-gee.ace("sol.ana.long", "ace_c4", "eage_grim") #Epigenetic age change using GrimAge
ana.2.4c

ana.3.4c<-gee.ace("sol.ana.long", "ace_c4", "dunedin") # Epigenetic age Pace
ana.3.4c

summary.1.4c<-result.summary(ana.1.4c, "ace_c4", model = "GEE")
summary.2.4c<-result.summary(ana.2.4c, "ace_c4", model = "GEE")
summary.3.4c<-result.summary(ana.3.4c, "ace_c4", model = "GEE")

test.2<-rbind(cbind(summary.1.4c,summary.1.4c.ml),
             cbind(summary.2.4c,summary.2.4c.ml),
             cbind(summary.3.4c,summary.3.4c.ml))

### Unbalance design
glm.ace.ub<-function(data, exposure, outcome) {
  
  formula_1<-as.formula(paste(outcome, "~", exposure,"+ age_change + AGE + GENDER + CENTER + age_change:", exposure))
  model.1<-lme(formula_1, 
               data = get(data),
               random = ~1+age_change|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
  #print(summary(model.1))
  robust.1<-coef_test(model.1, vcov = "CR0", test = "z")
  
  
  formula_2<-as.formula(paste(outcome, "~", exposure,"+ age_change + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + age_change:", exposure))
  model.2<-lme(formula_2, 
               data = get(data),
               random = ~1+age_change|ID,
               weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
               method = "REML",
               na.action = na.omit)
 # print(summary(model.2))
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

summary.1.ub<-result.summary(ana.1.ub, exposure = "ace_c2", model = "LM")
summary.2.ub<-result.summary(ana.2.ub, exposure = "ace_c2", model = "LM")
summary.3.ub<-result.summary(ana.3.ub, exposure = "ace_c2", model = "LM")
summary.4.ub<-result.summary(ana.4.ub, exposure = "ACE_TOT", model = "LM")
summary.5.ub<-result.summary(ana.5.ub, exposure = "ACE_TOT", model = "LM")
summary.6.ub<-result.summary(ana.6.ub, exposure = "ACE_TOT", model = "LM")

test.3<-rbind(summary.1.ub,summary.4.ub,summary.2.ub,summary.5.ub,summary.3.ub,summary.6.ub)

## Adjust for ancestry 
#colnames(sol.ac)[3]<-"ID"
colnames(sol.ac.pc)[1]<-"ID"
sol.ana.ac<-merge(sol.ana.1,sol.ac.pc, by = "ID")

## Long format dataset including ancestry
sol.ana.long.ac<-reshape(data = sol.ana.ac,
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
sol.ana.long.ac<-sol.ana.long.ac[order(sol.ana.long.ac$ID),]
sol.ana.long.ac$t<-rep(c(1:2),nrow(sol.ana.ac))
sol.ana.long.ac$time<-factor(sol.ana.long.ac$time)
sol.ana.long.ac<-sol.ana.long.ac%>%mutate(ace_c2 = ifelse(ACE_TOT<4,0,1),
                                    age_arr_US = factor(ifelse(US_BORN==0|is.na(US_BORN),ifelse(AGE-YRSUS<18,2,3),1)),
                                    BKGRD1_C7 = factor(BKGRD1_C7),
                                    ace_c4 = factor(ifelse(ACE_TOT==0,1,
                                                           ifelse(ACE_TOT<4,2,3))),
                                    age_change = age_t-AGE)


### Function for adjusting for ancestry
####GEE
gee.ace.ac<-function(data, exposure, outcome) {
  
  formula_ac<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + EV1 + EV2 + EV3 + EV4 + EV5 + time:", exposure))
  model.ac<-gls(formula_ac, 
               data = get(data),
               corr=corCompSymm(form = ~ t | ID),
               weights = varComb(varIdent(form = ~ 1 | t), varFixed(~WEIGHT_NORM_OVERALL_EPIGEN)),
               method = "REML",
               na.action = na.omit)
  #print(summary(model.ac))
  robust.ac<-coef_test(model.ac, vcov = "CR0", test = "z")
  return(robust.ac)
}

ana.1.a<-gee.ace.ac("sol.ana.long.ac", "ace_c2", "eaa_grim")
ana.1.a
ana.2.a<-gee.ace.ac("sol.ana.long.ac", "ace_c2", "eage_grim")
ana.2.a
ana.3.a<-gee.ace.ac("sol.ana.long.ac", "ace_c2", "dunedin")
ana.3.a
ana.4.a<-gee.ace.ac("sol.ana.long.ac", "ACE_TOT", "eaa_grim")
ana.4.a
ana.5.a<-gee.ace.ac("sol.ana.long.ac", "ACE_TOT", "eage_grim")
ana.5.a
ana.6.a<-gee.ace.ac("sol.ana.long.ac", "ACE_TOT", "dunedin")
ana.6.a

summary.1.ac<-result.summary(ana.1.a, "ace_c2", "ancenstry", "GEE")
summary.2.ac<-result.summary(ana.2.a, "ace_c2", "ancenstry", "GEE")
summary.3.ac<-result.summary(ana.3.a, "ace_c2", "ancenstry", "GEE")
summary.4.ac<-result.summary(ana.4.a, "ACE_TOT", "ancenstry", "GEE")
summary.5.ac<-result.summary(ana.5.a, "ACE_TOT", "ancenstry", "GEE")
summary.6.ac<-result.summary(ana.6.a, "ACE_TOT", "ancenstry", "GEE")

#### GLM
glm.ace.ac<-function(data, exposure, outcome) {
  formula.ac<-as.formula(paste(outcome, "~", exposure,"+ time + AGE + GENDER + CENTER + EDUCATION_C2 + age_arr_US + EV1 + EV2 + EV3 + EV4 + EV5 + time:", exposure))
  model.ac<-lme(formula.ac, 
                  data = get(data),
                  random = ~1+time|ID,
                  weights = varFixed(~WEIGHT_NORM_OVERALL_EPIGEN),
                  method = "REML",
                  na.action = na.omit)
  #print(summary(model.ac))
  robust.ac<-coef_test(model.ac, vcov = "CR0", test = "z")
  return(robust.ac)
}


ana.1.ac<-glm.ace.ac("sol.ana.long.ac", "ace_c2", "eaa_grim") #Binary ACE (<4 vs. >=4) and EAA GrimAge
#ana.1.ac
ana.2.ac<-glm.ace.ac("sol.ana.long.ac", "ace_c2", "eage_grim") #Binary ACE (<4 vs. >=4) and Epigenetic age change using GrimAge
#ana.2.ac

ana.3.ac<-glm.ace.ac("sol.ana.long.ac", "ace_c2", "dunedin") #Binary ACE (<4 vs. >=4) and Epigenetic age Pace
#ana.3.ac
ana.4.ac<-glm.ace.ac("sol.ana.long.ac", "ACE_TOT", "eaa_grim") #Continuous ACE and EAA GrimAge
#ana.4.ac
ana.5.ac<-glm.ace.ac("sol.ana.long.ac", "ACE_TOT", "eage_grim") #Continuous ACE and Epigenetic age change using GrimAge
#ana.5.ac
ana.6.ac<-glm.ace.ac("sol.ana.long.ac", "ACE_TOT", "dunedin") #Continuous ACE and Epigenetic age Pace
#ana.6.ac

summary.1.ml.ac<-result.summary(ana.1.ac, "ace_c2", "ancenstry", "LM")
summary.2.ml.ac<-result.summary(ana.2.ac, "ace_c2", "ancenstry", "LM")
summary.3.ml.ac<-result.summary(ana.3.ac, "ace_c2", "ancenstry", "LM")
summary.4.ml.ac<-result.summary(ana.4.ac, "ACE_TOT", "ancenstry", "LM")
summary.5.ml.ac<-result.summary(ana.5.ac, "ACE_TOT", "ancenstry", "LM")
summary.6.ml.ac<-result.summary(ana.6.ac, "ACE_TOT", "ancenstry", "LM")

test.4<-rbind(cbind(summary.1.ac,summary.1.ml.ac),
              cbind(summary.4.ac,summary.4.ml.ac),
              cbind(summary.2.ac,summary.2.ml.ac),
              cbind(summary.5.ac,summary.5.ml.ac),
              cbind(summary.3.ac,summary.3.ml.ac),
              cbind(summary.6.ac,summary.6.ml.ac))





