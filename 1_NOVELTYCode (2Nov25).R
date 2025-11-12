install.packages("sqldf")
library(sqldf)

#setwd("M:/Projects/2025/Project.COPDGeneNOVELTY/Analysis")

unknown_to_NA <- function(x){
  ifelse(x=="UNKNOWN",NA,x)
}

extract_exac_counts <- function(df_exac, year) #year=0:baseline, 1:year 1, 2:year 2, 3:year 3
{
  exac_vars <- c("ACC_EME","UNSNONEM","DOCTVIS","HLTPROVI","NEBCORT","ORALCORT","INJECORT", "ANTIBIO","HOSPADM", "EXACST") # last - L
  
  i <- c(1,2,4,5)[year+1]
  df_exac_process <- df_exac %>%
    select(usubjid,zcgrpid,paramcd,avalc,avisitn, covar29) %>% #last - L
    filter(paramcd %in% exac_vars) %>%
    filter(avisitn %in% c(i))
  
  #Remove duplicates
  tmp1 <- sqldf("SELECT usubjid, zcgrpid, paramcd, avalc, COUNT(*) AS N FROM df_exac_process GROUP BY usubjid,zcgrpid, paramcd, avalc")
  if(length(unique(tmp1$N))>1)
  {
    df_exac_process$id <- 1:nrow(df_exac_process)
    bads <- tmp1[which(tmp1$N==2),]
    for(j in 1:nrow(bads))
    {
      bad <- bads[j,]
      l1 <- which(df_exac_process$usubjid==bad$usubjid)
      l2 <- which(df_exac_process[l1,]$zcgrpid==bad$zcgrpid)
      l3 <- which(df_exac_process[l1,][l2,]$paramcd==bad$paramcd)
      id <- df_exac_process[l1,][l2,][l3,][1,]$id
      df_exac_process <- df_exac_process[-which(df_exac_process$id==id),]
    }
  }
  
  df_exac_process <- df_exac_process %>% pivot_wider(names_from=paramcd,values_from=avalc) %>% 
    mutate(EXACST = case_when(
        grepl("^\\d{4}$", EXACST) ~ paste0(EXACST, "-12-31"),          # year only → assume Dec 31
        grepl("^\\d{4}-\\d{2}$", EXACST) ~ paste0(EXACST, "-28"),       # year-month → assume 30th of month
        TRUE ~ EXACST                                              # full date, leave as is
        ),
      EXACST = as.Date(EXACST)
    ) %>% 
    mutate(across(c(ACC_EME:ORALCORT),unknown_to_NA)) %>% 
    filter(!(is.na(ACC_EME) & is.na(UNSNONEM) & is.na(DOCTVIS) & is.na(HLTPROVI) & is.na(NEBCORT) &
               is.na(ORALCORT) & is.na(INJECORT) & is.na(ANTIBIO) & is.na(HOSPADM))) %>% 
    mutate(AE_UNSCH = case_when(is.na(ACC_EME) & is.na(UNSNONEM) ~ NA,
                                !is.na(ACC_EME) & ACC_EME=='YES' ~ 'YES',
                                !is.na(UNSNONEM) & UNSNONEM=='YES' ~'YES', 
                                TRUE ~ 'NO'),
           SCS = case_when(is.na(ORALCORT) & is.na(NEBCORT) & is.na(INJECORT) ~ NA,
                           !is.na(ORALCORT) & ORALCORT=='YES' ~ 'YES',
                           !is.na(NEBCORT) & NEBCORT=='YES' ~ 'YES',
                           !is.na(INJECORT) & INJECORT=='YES' ~ 'YES',
                           TRUE ~ 'NO'),
           MODEX = case_when(is.na(SCS) & is.na(ANTIBIO) & 
                               is.na(DOCTVIS) & is.na(HLTPROVI) & 
                               is.na(AE_UNSCH) & is.na(HOSPADM) ~ NA,
                             !is.na(HOSPADM) & HOSPADM=='YES' ~ 'NO',
                             !is.na(ACC_EME) & ACC_EME=='YES' ~ "NO",
                             (!is.na(SCS) & SCS=='YES') |
                               (!is.na(ANTIBIO) & ANTIBIO=='YES' ) |
                               (!is.na(DOCTVIS) & DOCTVIS=='YES') |
                               (!is.na(HLTPROVI) & HLTPROVI=='YES') ~ 'YES',
                             TRUE ~ 'NO'),
           # SEVEX = ifelse(grp01bl=="COPD", 
           #                case_when(is.na(HOSPADM) ~ NA,
           #                  !is.na(HOSPADM) & HOSPADM=="YES" ~ "YES",
           #                TRUE ~ "NO"),
           #                case_when(is.na(SCS) & is.na(AE_UNSCH) & is.na(HOSPADM) ~ NA,
           #                       (!is.na(SCS) & SCS=='YES') |
           #                                (!is.na(HOSPADM) & HOSPADM=='YES') |
           #                                (!is.na(AE_UNSCH) & AE_UNSCH=='YES') ~ 'YES',
           #                       TRUE ~ 'NO'))) %>% 
           SEVEX =    case_when(is.na(HOSPADM) & is.na(ACC_EME) ~ NA,
                                !is.na(HOSPADM) & HOSPADM=="YES" ~ "YES",
                                !is.na(ACC_EME) & ACC_EME=='YES' ~ 'YES',
                                TRUE ~ "NO")) 
  
  df_ex_sum <- df_exac_process %>% 
    # 58 cases where sev ex is NA; treat them as 0
    arrange(usubjid, EXACST) %>% group_by(usubjid) %>% 
    
    # 14-day window for outcome
    filter(MODEX=="YES" | SEVEX=="YES") %>%
    mutate(diff_days = as.numeric(difftime(EXACST, lag(EXACST), units = "days")), 
           within14 = ifelse(!is.na(diff_days) & diff_days <= 14, 1, 0),
           cluster = cumsum(ifelse(is.na(within14) | within14 == 0, 1, 0))) %>%
      ungroup() %>% group_by(usubjid, cluster) %>% 
      summarise(
        EXACST = min(EXACST), # start of cluster
        SEVEX = ifelse(any(SEVEX == "YES"), "YES", "NO"),
        MODEX = ifelse(SEVEX == "YES", "NO", "YES"),
      ) %>%
    summarise(num_ex = n(),
              num_ex_moderate = sum(MODEX=="YES",na.rm=T),
              num_ex_severe = sum(SEVEX == "YES",na.rm=T),
              date_first_msexb = first(EXACST[(MODEX=="YES" | SEVEX == "YES") #& EXACST>=covar29
                                              ]),
              date_first_sexb = first(EXACST[(SEVEX == "YES") #& EXACST>=covar29
              ]),
              #covar29 = first(covar29),
              type_first_msexb = ifelse(is.na(date_first_msexb), 0, ifelse(MODEX[EXACST==date_first_msexb]=="YES", 1, ifelse(SEVEX[EXACST==date_first_msexb]=="YES", 2, NA)))
              )
  
  df_ex_sum
  #df_exac_process
}

df0 <- extract_exac_counts(df_exac, 0) 
df1 <- extract_exac_counts(df_exac, 1)
df2 <- extract_exac_counts(df_exac, 2) 
df3 <- extract_exac_counts(df_exac, 3) 

temp0 <- df0 %>% rename(m0 = num_ex_moderate, s0 = num_ex_severe, d0=date_first_msexb, ds0=date_first_sexb, t0=type_first_msexb) %>% mutate(m0s0 = m0 + s0) %>% right_join(df_analysis, by="usubjid")
temp1 <- df1 %>% rename(m1 = num_ex_moderate, s1 = num_ex_severe, d1=date_first_msexb, ds1=date_first_sexb, t1=type_first_msexb) %>% mutate(m1s1 = m1 + s1) %>% right_join(temp0, by="usubjid")
temp2 <- df2 %>% rename(m2 = num_ex_moderate, s2 = num_ex_severe, d2=date_first_msexb, ds2=date_first_sexb, t2=type_first_msexb) %>% mutate(m2s2 = m2 + s2) %>% right_join(temp1, by="usubjid")
df_final <- df3 %>% rename(m3 = num_ex_moderate, s3 = num_ex_severe, d3=date_first_msexb, ds3=date_first_sexb, t3=type_first_msexb) %>% mutate(m3s3 = m3 + s3) %>% right_join(temp2, by="usubjid")

df_final <- df_final %>% 
  mutate(m0=replace_na(m0,0)) %>%
  mutate(s0=replace_na(s0,0)) %>%
  mutate(m1=replace_na(m1,0)) %>%
  mutate(s1=replace_na(s1,0)) %>%
  mutate(m2=replace_na(m2,0)) %>%
  mutate(s2=replace_na(s2,0)) %>%
  mutate(m3=replace_na(m3,0)) %>%
  mutate(s3=replace_na(s3,0)) %>%
  mutate(m0s0=replace_na(m0s0,0)) %>%
  mutate(m1s1=replace_na(m1s1,0)) %>%
  mutate(m2s2=replace_na(m2s2,0)) %>%
  mutate(m3s3=replace_na(m3s3,0)) 
  
write.csv(df_final[,c('sex','m0','s0','m1','s1','m2','s2','m3','s3')], file = "Exac_patterns.csv")



#Based on FEV1/FVC<0.7

df_SA1 <- df_final %>% filter(`FEV1/FVC`<70) %>%
  mutate(m0=replace_na(m0,0)) %>%
  mutate(s0=replace_na(s0,0)) %>%
  mutate(m1=replace_na(m1,0)) %>%
  mutate(s1=replace_na(s1,0)) %>%
  mutate(m2=replace_na(m2,0)) %>%
  mutate(s2=replace_na(s2,0)) %>%
  mutate(m3=replace_na(m3,0)) %>%
  mutate(s3=replace_na(s3,0))

write.csv(df_SA1[,c('sex','m0','s0','m1','s1','m2','s2','m3','s3')], file = "Exac_patterns_sa1.csv")


### 2025.06.13 (Table 1 for SA)
out <- list()
out$n <-nrow(df_SA1)
out$age <- c(mean=mean(df_SA1$age), sd=sd(df_SA1$age))
out$sex <- c(male=mean(df_SA1$sex=="M"))
out$race <- table(df_SA1$race)
out$race_group <- table(df_SA1$racegr1)
out$ethnicity <- table(df_SA1$ethnic)
out$BMI <- c(mean=mean(df_SA1$bmibl, na.rm=T), sd=sd(df_SA1$bmibl, na.rm=T), n_miss=sum(is.na(df_SA1$bmibl)))
out$smoking_status <- table(df_SA1$covar02)
out$pack_years <- c(mean=mean(df_SA1$TOTPACK, na.rm=T), sd=sd(df_SA1$TOTPACK, na.rm=T), n_miss=sum(is.na(df_SA1$TOTPACK)))
out$FEV1 <- c(mean=mean(df_SA1$FEV1, na.rm=T), sd=sd(df_SA1$FEV1, na.rm=T), n_miss=sum(is.na(df_SA1$FEV1)))
out$FEV1P <- c(mean=mean(df_SA1$FEV1P, na.rm=T), sd=sd(df_SA1$FEV1P, na.rm=T), n_miss=sum(is.na(df_SA1$FEV1P)))
out$FVC <- c(mean=mean(df_SA1$FVC, na.rm=T), sd=sd(df_SA1$FVC, na.rm=T), n_miss=sum(is.na(df_SA1$FVC)))
out$FVCP <- c(mean=mean(df_SA1$FVCP, na.rm=T), sd=sd(df_SA1$FVCP, na.rm=T), n_miss=sum(is.na(df_SA1$FVCP)))
out$`FEV1/FVC` <- c(mean=mean(df_SA1$`FEV1/FVC`, na.rm=T), sd=sd(df_SA1$`FEV1/FVC`, na.rm=T), n_miss=sum(is.na(df_SA1$`FEV1/FVC`)))
out$n_exac <- c(mean=mean(df_SA1$n_exac_y0, na.rm=T), sd=sd(df_SA1$n_exac_y0, na.rm=T), n_miss=sum(is.na(df_SA1$n_exac_y0)))
out$gold <- table(cut(df_SA1$FEV1P, breaks=c(-Inf, 30, 50, 80, 200), labels=c(4,3,2,1)))
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

write.csv(out2, "table1_SA1.csv")