# visit indexed: BMI (df_vs$visit), smoking (df_smoking$avisit), MMRC (df_adqs$avisit), oxygen (df_adzc$avisit), med (), FEV1 (df_pft$avisit) 

# OUTPUT: 
# data/IV_cohort

# settings ----------------------------------------------------------------
install.packages("tidyverse")
library(tidyverse)
install.packages("here")
library(here)
#install.packages("haven", type="source")
install.packages("haven")
library(haven)
install.packages("janitor")
library(janitor)

#data_1_path <- "N:/data/1/sdtm/data/"
#data_2_path <- "N:/data/2/adam/data/"
data_1_path <- "/Volumes/Untitled/data/1/sdtm/data/"
data_2_path <- "/Volumes/Untitled/data/2/adam/data/"

data_1_files <- list.files(data_1_path)
data_2_files <- list.files(data_2_path)

data_dic <- read_csv("/Volumes/Untitled/data/file_dic.csv")

data_reader <-function(num,filename){
  if(num==1){
    read_sas(paste0(data_1_path,filename,".sas7bdat")) %>% 
      clean_names()
  } else{
    read_sas(paste0(data_2_path,filename,".sas7bdat")) %>% 
      clean_names()
  }
}

df_sub <- data_reader(2,"DEID_ADSL") # subject lvl analysis

df_sub %>% 
  distinct(usubjid) %>% 
  nrow()

df_sub %>% 
  filter(age>=40) %>% 
  nrow()

df_sub %>% 
  filter(age>=40) %>% 
  mutate(FU = rfendt-rfstdt) %>% 
  filter(FU<30) %>% 
  nrow()

df_cohort <- df_sub %>% 
  filter(age>=40) %>% 
  mutate(FU = rfendt-rfstdt) %>% 
  filter(FU>=30)

df_cohort %>% filter(!grp01bl %in% "ASTHMA") -> df_cohort #REMOVE ASTHMA

ids <- df_cohort$usubjid
n_ids <- length(ids)

#####2025.05.27: add spirometry and smoking
df_pft <- data_reader(2,"DEID_ADRE")
df_pft2 <- df_pft %>% filter(usubjid %in% ids & avisitn==1 & paramcd %in% c("FEV1", "FEV1P", "FVC", "FVCP", "FEV1/FVC") & atpt=="PRE PFT")
df_pft3 <- df_pft2 %>% pivot_wider(id_cols=usubjid, names_from=paramcd,  values_from=aval)
df_cohort <- df_cohort %>%  left_join(df_pft3, by = "usubjid")

df_smoking <- data_reader(2,"DEID_ADSU")
df_smoking2 <- df_smoking %>% filter(usubjid %in% ids & avisitn==1 & paramcd %in% c("TOTPACK"))
df_smoking3 <- df_smoking2 %>% pivot_wider(id_cols=usubjid, names_from=paramcd,  values_from=aval)
df_cohort <- df_cohort %>%  left_join(df_smoking3, by = "usubjid")

# SENSITIVITY ANALYSIS: spirometrically confirmed COPD

# df_cohort %>% filter(!is.na(`FEV1/FVC`) & `FEV1/FVC`<=70) -> df_cohort

##### END 2025.05.27

df_analysis <- df_cohort %>% 
  select(usubjid,
         FU,
         rfstdt,
         rfendt,
         covar29:covar32,
         covar02, #Smoking status
         covar14:covar16, # num gp, hosp, exacs at bsl
         covar23, # num exac at yr1
         covar26, # num exac at yr2
         country,geogr1,
         age,sex,
         bmibl,
         race,racegr1,
         ethnic,
         TOTPACK,
         grp01bl, # asth/acopd grouping at bsl
         severity, # disease severity at bsl
         covar07, # treatment setting at bsl
         covar09,hghtbl,
         covar21,
         FEV1, FEV1P, FVC, FVCP, `FEV1/FVC`
  ) %>% # CAAT
  mutate(#FU = pmin(FU,365),
         FU = FU/365) %>% 
  rename(n_exac_y0 = covar16,
         n_exac_y1 = covar23,
         n_exac_y2 = covar26
  )





############ EXACERBATIONS ##########
df_exac <- data_reader(2,"deid_adzc") %>%
  filter(usubjid %in% ids)








### 2025.05.19 (Table 1 stuff)
out <- list()
out$n <-nrow(df_analysis)
out$age <- c(mean=mean(df_analysis$age), sd=sd(df_analysis$age))
out$sex <- c(male=mean(df_analysis$sex=="M"))
out$race <- table(df_analysis$race)
out$race_group <- table(df_analysis$racegr1)
out$ethnicity <- table(df_analysis$ethnic)
out$BMI <- c(mean=mean(df_analysis$bmibl, na.rm=T), sd=sd(df_analysis$bmibl, na.rm=T), n_miss=sum(is.na(df_analysis$bmibl)))
out$smoking_status <- table(df_analysis$covar02)
out$pack_years <- c(mean=mean(df_analysis$TOTPACK, na.rm=T), sd=sd(df_analysis$TOTPACK, na.rm=T), n_miss=sum(is.na(df_analysis$TOTPACK)))
out$FEV1 <- c(mean=mean(df_analysis$FEV1, na.rm=T), sd=sd(df_analysis$FEV1, na.rm=T), n_miss=sum(is.na(df_analysis$FEV1)))
out$FEV1P <- c(mean=mean(df_analysis$FEV1P, na.rm=T), sd=sd(df_analysis$FEV1P, na.rm=T), n_miss=sum(is.na(df_analysis$FEV1P)))
out$FVC <- c(mean=mean(df_analysis$FVC, na.rm=T), sd=sd(df_analysis$FVC, na.rm=T), n_miss=sum(is.na(df_analysis$FVC)))
out$FVCP <- c(mean=mean(df_analysis$FVCP, na.rm=T), sd=sd(df_analysis$FVCP, na.rm=T), n_miss=sum(is.na(df_analysis$FVCP)))
out$`FEV1/FVC` <- c(mean=mean(df_analysis$`FEV1/FVC`, na.rm=T), sd=sd(df_analysis$`FEV1/FVC`, na.rm=T), n_miss=sum(is.na(df_analysis$`FEV1/FVC`)))
out$n_exac <- c(mean=mean(df_analysis$n_exac_y0, na.rm=T), sd=sd(df_analysis$n_exac_y0, na.rm=T), n_miss=sum(is.na(df_analysis$n_exac_y0)))
out$gold <- table(cut(df_analysis$FEV1P, breaks=c(-Inf, 30, 50, 80, 200), labels=c(4,3,2,1)))
out2 <- data.frame(name=character(), param=character(), value=double())
index <- 1
for(i in 1:length(out))
{
  nm <- names(out)[i]
  if(length(out[[i]])==1)
  {
    out2[index,1] <- nm
    if(!is.null(names(out[[i]]))) out2[index,2] <- names(out[[i]])
    out2[index,3] <- out[[i]]
    index <- index+1
  }
  else
  {
    for(j in 1:length(out[[i]]))
    {
      out2[index,1] <- nm
      out2[index,2] <- names(out[[i]])[j]
      out2[index,3] <- out[[i]][j]
      index <- index+1
    }
  }
  cat(nm)
}

write.csv(out2, "table1.csv")