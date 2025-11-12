#### SGRQ ######

# Filter by param and ablfl

library(dplyr)
df_adsgrq <- data_reader(2,"DEID_ADSGRQ")
df_adsgrq %>% filter(param=="SGRQ Total score" & ablfl=="Y" & usubjid %in% ids) %>% 
  select(usubjid, aval, avalc) %>% rename(sgrq = aval, sgrq_cat = avalc) %>% 
  right_join(df_final, by="usubjid") -> df_final

#### Treatment (within 12 months prior to the baseline visit) ######

df_adcm <- data_reader(2,"DEID_ADCM")
df_code_trt <- data_reader(2,"DEID_CODE_TRT")

df_code_trt %>% rename(cmtrt = trt) %>% right_join(df_adcm, by="cmtrt") -> df_adcm
table(df_adcm$class, useNA="ifany")
table(df_adcm$onmlbpbl, useNA="ifany")

df_adcm %>%  
  filter(usubjid %in% ids & onmlbpbl=="Y") %>%  # nrow() #10767
  dplyr::select(usubjid, class, cmtrt) %>% 
  right_join(df_final, by="usubjid") %>% 
  group_by(usubjid) %>% 
  mutate(LAMA = ifelse(any(grepl("LAMA", class, ignore.case=T))==T,1,0),
         LABA = ifelse(any(grepl("LABA", class, ignore.case=T))==T,1,0),
         ICS = ifelse(any(grepl("ICS", class, ignore.case=T))==T,1,0),
         SABA = ifelse(any(grepl("SABA", class, ignore.case=T))==T,1,0), 
         STAT = ifelse(any(grepl("statin", cmtrt, ignore.case=T))==T,1,0), 
         AZIT = ifelse(any(grepl("azitromycin", cmtrt, ignore.case=T)==T | 
                         grepl("azithromycin", cmtrt, ignore.case=T)==T),1,0)) %>%  # nrow() #10961
  select(-class, -cmtrt) %>% slice(1) -> df_final

# CVD

df_admh <- data_reader(2,"DEID_ADMH")

df_admh %>% filter(usubjid %in% df_final$usubjid & mhseq==1) %>% 
  dplyr::select(usubjid, cvdfl) %>% 
  rename(CVD = cvdfl) %>% 
  mutate(CVD=ifelse(CVD=="Y",1,0)) %>%
  right_join(df_final, by="usubjid") -> df_final

# Oxygen therapy

df_adzc <- data_reader(2,"DEID_ADZC")
length(unique(df_adzc$usubjid)) #11170
table(df_adzc$param)
table(df_adzc$avisitn)

df_adzc %>% filter(usubjid %in% ids & param == "Portable Oxygen Therapy?" & zcgrpid==1 & avisitn==1) %>%
  dplyr::select(usubjid, avalc) %>% 
  rename(Oxygen = avalc) %>% 
  mutate(Oxygen=ifelse(Oxygen=="YES",1,ifelse(Oxygen=="NO", 0, NA))) %>%
  right_join(df_final, by="usubjid") -> df_final

########### MRC

df_adqs <- data_reader(2,"DEID_ADQS")
df_adqs %>% filter(usubjid%in%ids & param=="MMRC Dyspnea Scale Grade" & ablfl=="Y") %>%
  dplyr::select(usubjid, aval, avalc) %>% rename(mmrc = aval, mmrc_cat = avalc) %>% 
  right_join(df_final, by="usubjid") -> df_final

###### BEC

df_lb <- data_reader(1,"DEID_LB")
df_lb %>% filter(usubjid %in% ids & lbtest=="Eosinophils" & visit=="BASELINE (VISIT 1)") %>% group_by(usubjid) %>% slice(1) %>%
  select(usubjid, lbstresn, lbstresu) %>% rename(bec = lbstresn, bec_unit = lbstresu) %>% 
  right_join(df_final, by="usubjid") -> df_final

##### Death rate

df_dm <- data_reader(1,"DEID_DM") 
df_dm %>% select(usubjid, dthdtc) %>% right_join(df_final, by="usubjid") -> df_final
class(df_final$dthdtc)
df_final$dthdtc <- as.Date(df_final$dthdtc)
table(df_final$dthdtc, useNA="ifany")
  
###### Exb Date

df_final %>% 
  mutate(event=ifelse(m1s1>0, 1, 0)) %>%
  mutate(event = ifelse(m1s1==0 & !is.na(dthdtc) & dthdtc<=(covar29 %m+% years(1)), 2, 
                        #ifelse(m1==0 & s1==0 & ((!is.na(dthdtc) & dthdtc>(covar29 %m+% years(1))) | is.na(dthdtc)), 0, 1)
                        event),
         event_date = as.Date(ifelse(event==1, d1, 
                                     ifelse(event==2, dthdtc, 
                                            ifelse(FU<1, rfendt, covar29 %m+% years(1))))),
         event_date_severe = as.Date(ifelse(event==1 & s1>0, ds1, 
                                     ifelse(event==2, dthdtc, 
                                            ifelse(FU<1, rfendt, covar29 %m+% years(1)))))
  ) -> df_final

df_final$time2event <- time_length(interval(df_final$covar29, df_final$event_date), "months")
df_final$time2event <- ifelse(df_final$time2event<0, 0, df_final$time2event)
df_final$time2event <- ifelse(df_final$time2event>12, 12, df_final$time2event)
df_final$time2event_severe <- time_length(interval(df_final$covar29, df_final$event_date_severe), "months")
df_final$time2event_severe <- ifelse(df_final$time2event_severe<0, 0, df_final$time2event_severe)
df_final$time2event_severe <- ifelse(df_final$time2event_severe>12, 12, df_final$time2event_severe)

##### Table 1 

df_analysis <- df_final
#df_analysis <- df_imputed

out <- list()
out$n <-nrow(df_analysis)
out$age <- c(mean=mean(df_analysis$age), sd=sd(df_analysis$age), n_miss=sum(is.na(df_analysis$age)))
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
out$n_exac1 <- c(mean=mean(df_analysis$n_exac_y1, na.rm=T), sd=sd(df_analysis$n_exac_y1, na.rm=T), n_miss=sum(is.na(df_analysis$n_exac_y1)))
out$gold <- table(cut(df_analysis$FEV1P, breaks=c(-Inf, 30, 50, 80, 200), labels=c(4,3,2,1)))
out$SGRQ <- c(mean=mean(df_analysis$sgrq, na.rm=T), sd=sd(df_analysis$sgrq, na.rm=T), n_miss=sum(is.na(df_analysis$sgrq)))
out$ICS <- table(df_analysis$ICS, useNA="ifany")
out$LABA <- table(df_analysis$LABA, useNA="ifany")
out$LAMA <- table(df_analysis$LAMA, useNA="ifany")
out$AZIT <- table(df_analysis$AZIT, useNA="ifany")
out$STAT <- table(df_analysis$STAT, useNA="ifany")
out$CVD <- table(df_analysis$CVD, useNA="ifany")
out$Oxygen <- table(df_analysis$Oxygen, useNA="ifany")
out$MMRC <- table(df_analysis$mmrc, useNA="ifany")
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

write.csv(out2, "~/Documents/UBC/accept/table1.csv")
#write.csv(out2, "~/Documents/UBC/accept/table1_imputed.csv")
out2

##### ACCEPT

install.packages("accept")
library(accept)
install.packages("jsonlite", type="source")
library(jsonlite)
install.packages("missRanger")
library(missRanger)

df_final <- df_final %>% mutate(covar02 = ifelse(covar02=="", NA, covar02))
df_final <- df_final %>% mutate(ethnic = ifelse(ethnic=="", NA, ethnic))

# df_imputed <- missRanger(df_final, pmm.k=5, maxit=10, num.trees = 200, seed = 123)

df_final %>% dplyr::select(country, usubjid, time2event, time2event_severe, event, event_date, 
                           age, sex, bmibl, covar02, 
         #covar21, 
         #sgrq,
         mmrc, CVD, 
         ICS, LABA, LAMA, n_exac_y0, m0, s0, FEV1P, FEV1, Oxygen, bec) %>% 
  rename(ID = usubjid, male = sex, BMI = bmibl, smoker = covar02, 
         #SGRQ = sgrq, 
         #CAAT = covar21, 
         mMRC = mmrc, 
         #statin = STAT,
         statin = CVD,
         oxygen = Oxygen, 
         LastYrExacCount = n_exac_y0, 
         LastYrSevExacCount = s0, 
         FEV1 = FEV1P, FEV1_raw = FEV1
         ) %>% 
  mutate(male = as.logical(male=="M"), 
         age = as.double(age),
         BMI = as.double(BMI),
         ICS = as.logical(ICS),
         LAMA = as.logical(LAMA),
         LABA = as.logical(LABA),
         statin = as.logical(statin),
         oxygen = as.logical(oxygen),
         smoker = as.logical(ifelse(smoker=="CURRENT SMOKER",T,F)), #because accept2 coded coefficient for nowsmk
         #SGRQ = as.integer(SGRQ),
         #CAAT = as.integer(CAAT),
         mMRC = as.numeric(mMRC), #not factor because accept2 input it into sgrq using a+b*mmrc
         FEV1 = as.double(ifelse(FEV1>=110, 110, FEV1)),
         LastYrExacCount = as.integer(LastYrExacCount)) -> df_accept
df_imputed <- missRanger(df_accept, formula=~.~-c(ID, event_date, time2event, time2event_severe, event), pmm.k=5, maxit=10, num.trees = 200, seed = 123)
model_accept2 <- accept(newdata=df_imputed, version="flexccept", format = "tibble")
df_analysis <- bind_cols(df_imputed, model_accept2[,c("predicted_exac_probability", 
                                                   "predicted_exac_rate", 
                                                   "predicted_severe_exac_probability", 
                                                   "predicted_severe_exac_rate")])



