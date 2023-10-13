library(geepack)
library(tidyverse)
sol.pheno<-read.csv("pheno_final_983obs_withcells.csv")
sol.v1.eaa<-read.csv("HCHS_epigen_age_V1_updated_100523.csv")
sol.v2.eaa<-read.csv("HCHS_epigen_age_V2_updated_100523.csv")
sol.v1.dun<-read.csv("HCHS_DunedinPACE_Visit1_updated_100523.csv")
sol.v2.dun<-read.csv("HCHS_DunedinPACE_Visit2_updated_100523.csv")
sol.ac<-read.csv("subject_annotation_2017-05-24.csv")

sol.v1.grim<-sol.v1.eaa%>%mutate(eaa_grim_v1=AgeAccelGrim,
                                 eage_grim_v1=DNAmGrimAgeBasedOnRealAge,
                                 age_v1=age)%>%
  select(c(ID,eaa_grim_v1,eage_grim_v1,age_v1))


sol.v2.grim<-sol.v2.eaa%>%mutate(eaa_grim_v2=AgeAccelGrim,
                                 eage_grim_v2=DNAmGrimAgeBasedOnRealAge,
                                 age_v2=age)%>%
  select(c(ID,eaa_grim_v2,eage_grim_v2,age_v2))


sol.v1.dun<-sol.v1.dun%>%mutate(dunedin_v1=DunedinPACE)%>%select(c(ID, dunedin_v1))
sol.v2.dun<-sol.v2.dun%>%mutate(dunedin_v2=DunedinPACE)%>%select(c(ID, dunedin_v2))
sol.ac.pc<-sol.ac%>%select(c(HCHS_ID,EV1,EV2,EV3,EV4,EV5))
sol.list<-list(sol.pheno,sol.v1.grim,sol.v2.grim,sol.v1.dun,sol.v2.dun)

sol.phe.eaa<-sol.list%>%reduce(inner_join,by="ID")




