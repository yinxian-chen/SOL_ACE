library(geepack)
library(tidyverse)
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
            "GR_1","GR_2","CD4_1","CD4_2","CD8_1","CD8_2")

ace<-"ACE_TOT"

eaa<-c("eaa_grim_v1","eaa_grim_v2","eage_grim_v1","eage_grim_v2","age_v1","age_v2","dunedin_v1","dunedin_v2")

sol.ana.1<-sol.phe.eaa[,c(cov_list,ace,eaa)]


