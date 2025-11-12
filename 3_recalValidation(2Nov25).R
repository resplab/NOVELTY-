### Preprocessing

df_analysis$predicted_exac_probability <- ifelse(df_analysis$predicted_exac_probability==0, 0.0001, df_analysis$predicted_exac_probability)
df_analysis$predicted_severe_exac_probability <- ifelse(df_analysis$predicted_severe_exac_probability==0, 0.0001, df_analysis$predicted_severe_exac_probability)
df_analysis <- df_final %>% dplyr::select(t1) %>% cbind(df_analysis) 

df_tdc <- df_analysis %>% rename(event.org = event) %>% mutate(event=ifelse(event.org==2,0,event.org), event_severe = ifelse(t1==2, 1, 0), event_severe = replace_na(event_severe, 0))
df_tdc <- df_tdc %>% rename(CVD=statin, num_history_mod = m0, num_history_sev = LastYrSevExacCount)

### Recalibration of ACCEPT

pkgs <- c("survival", "rms", 
          "timeROC", "riskRegression")
vapply(pkgs, function(pkg) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  require(pkg, character.only = TRUE, quietly = TRUE)
}, FUN.VALUE = logical(length = 1L))

model_cll <- coxph(Surv(time2event, event) ~ log(-log(1-predicted_exac_probability))+frailty(country),data = df_tdc) 
countries <- c("ARG", "AUS", "BRA", "CAN", "COL", "DEU", "DNK", "ESP", "FRA", "GBR", "ITA", "JPN", "KOR", "MEX", "NLD", "NOR", "SWE", "USA")
ints_cll <- data.frame(x=countries, y=model_cll$frail)
for (c in countries) {
  df_tdc[df_tdc$country==c, "int_cll"] <-  ints_cll[ints_cll$x==c,]$y
}
df_new <- df_tdc %>% mutate(time2event=12)
df_new$risk_cll <- 1-exp(-predict(model_cll, newdata=df_new, type="expected")*exp(df_new$int_cll))
df_tdc$risk_cll <- df_new$risk_cll
# Variance of RE
temp <- summary(model_cll)
temp$print2 
# Write country-specific RE and avg predicted risk
df_tdc %>% group_by(country) %>% summarise(obs_risk = mean(event==1)) %>% cbind(ints_cll) %>% dplyr::select(country, obs_risk, y) %>% rename(re = y) -> df_csv
write.csv(df_csv, "predict_re.csv")
library(ggrepel)
fit <- lm(re ~ obs_risk, data = df_csv)
coef_intercept <- round(coef(fit)[1], 3)
coef_slope <- round(coef(fit)[2], 3)
eq <- paste0("RE = ", coef_slope, " x average observed risk - ", abs(coef_intercept))
ggplot(data=df_csv, aes(x=obs_risk, y=re)) + geom_point() +geom_smooth(method="lm", linetype="dashed", col="black", se=F) +  
  geom_text_repel(aes(label = country), size = 3) + 
  labs(x="Avergae observed risk", y="Random effect (RE)", title="Country-specific observed risk and random effects (REs) from recalibrated ACCEPT 2.0") +   
  annotate("text", x = 0.35, y = 0.4, label = eq, hjust = 1.1, vjust = -1.1, size = 4.5, parse = FALSE, fontface="bold") +
  theme_bw()

### Refit (exploratory)

model_re <- coxph(Surv(time2event, event) ~ age+male+BMI+smoker+mMRC+CVD+ICS+LABA+LAMA+num_history_mod+num_history_sev+FEV1+bec+frailty(country),data = df_tdc) 
countries <- c("ARG", "AUS", "BRA", "CAN", "COL", "DEU", "DNK", "ESP", "FRA", "GBR", "ITA", "JPN", "KOR", "MEX", "NLD", "NOR", "SWE", "USA")
ints_re <- data.frame(x=countries, y=model_re$frail)
for (c in countries) {
  df_tdc[df_tdc$country==c, "int_re"] <-  ints_re[ints_re$x==c,]$y
}
df_new$int_re <- df_tdc$int_re
df_new$risk_re <- 1-exp(-predict(model_re, newdata=df_new, type="expected")*exp(df_new$int_re))
df_tdc$risk_re <- df_new$risk_re

model_re_severe <- coxph(Surv(time2event_severe, event_severe) ~ age+male+BMI+smoker+mMRC+CVD+ICS+LABA+LAMA+num_history_mod+num_history_sev+FEV1+bec,data = df_tdc) 
df_new$risk_re_severe <- 1-exp(-predict(model_re_severe, newdata=df_new, type="expected")*exp(df_new$int_re))
df_tdc$risk_re_severe <- df_new$risk_re_severe
summary(df_tdc$risk_re_severe)

### Recalibrtion - severe

model_cll_severe <- coxph(Surv(time2event_severe, event_severe) ~ log(-log(1-predicted_severe_exac_probability)), data = df_tdc) 
df_new <- df_new %>% mutate(time2event_severe=12)
df_new$risk_cll_severe <- 1-exp(-predict(model_cll_severe, newdata=df_new, type="expected")*exp(df_new$int_cll))
temp <- basehaz(model_cll_severe)
temp$hazard[nrow(temp)]
df_tdc$risk_cll_severe <- df_new$risk_cll_severe

######### DCA (Overall)

library(ggplot2)
df_dca <- data.frame(matrix(ncol=3, nrow=100*4))
colnames(df_dca) <- c("Strategy", "Risk Threshold", "NB")
#df_dca$Model <- rep(c("With RE", "Without RE"), each=100)
#df_dca$Model <- rep(c("Recalibrated", "Refitted"), each=100)
#df_dca$Strategy <- rep(c("Treat None", "Refitted ACCEPT 2.0", "Recalibrated ACCEPT 2.0", "Treat All"), each=100)
df_dca$Strategy <- rep(c("Treat None", "Original ACCEPT 2.0", "Recalibrated ACCEPT 2.0", "Treat All"), each=100)
df_dca$`Risk Threshold` <- rep_len(1:100, length.out=400)/100
z <- rep(0,100)
for (i in 1:100) {
  z[i] <- i/100
  #df_dca[df_dca$Strategy=="Without RE",]$NB[i] <- sum(df_tdc$event==1 & df_tdc$risk_wore>z[i])/dim(df_tdc)[1] - (sum(df_tdc$event==0 & df_tdc$risk_wore>z[i])/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  #df_dca[df_dca$Strategy=="With RE",]$NB[i] <- sum(df_tdc$event==1 & df_tdc$risk_re>z[i])/dim(df_tdc)[1] - (sum(df_tdc$event==0 & df_tdc$risk_re>z[i])/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  
  df_dca[df_dca$Strategy=="Original ACCEPT 2.0",]$NB[i] <- sum(df_tdc$event==1 & df_tdc$predicted_exac_probability>z[i])/dim(df_tdc)[1] - (sum(df_tdc$event==0 & df_tdc$predicted_exac_probability>z[i])/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  df_dca[df_dca$Strategy=="Recalibrated ACCEPT 2.0",]$NB[i] <- sum(df_tdc$event==1 & df_tdc$risk_cll>z[i])/dim(df_tdc)[1] - (sum(df_tdc$event==0 & df_tdc$risk_cll>z[i])/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  #df_dca[df_dca$Strategy=="Original ACCEPT 2.0",]$NB[i] <- sum(df_tdc$event_severe==1 & df_tdc$predicted_severe_exac_probability>z[i])/dim(df_tdc)[1] - (sum(df_tdc$event_severe==0 & df_tdc$predicted_severe_exac_probability>z[i])/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  #df_dca[df_dca$Strategy=="Recalibrated ACCEPT 2.0",]$NB[i] <- sum(df_tdc$event_severe==1 & df_tdc$risk_cll_severe>z[i])/dim(df_tdc)[1] - (sum(df_tdc$event_severe==0 & df_tdc$risk_cll_severe>z[i])/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  #df_dca[df_dca$Strategy=="Refitted ACCEPT 2.0",]$NB[i] <- sum(df_tdc$event==1 & df_tdc$risk_re>z[i])/dim(df_tdc)[1] - (sum(df_tdc$event==0 & df_tdc$risk_re>z[i])/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  #df_dca[df_dca$Strategy=="Refitted ACCEPT 2.0",]$NB[i] <- sum(df_tdc$event_severe==1 & df_tdc$risk_re_severe>z[i])/dim(df_tdc)[1] - (sum(df_tdc$event_severe==0 & df_tdc$risk_re_severe>z[i])/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  
  df_dca[df_dca$Strategy=="Treat None",]$NB[i] <- 0
  df_dca[df_dca$Strategy=="Treat All",]$NB[i] <- sum(df_tdc$event==1)/dim(df_tdc)[1] - (sum(df_tdc$event==0)/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  #df_dca[df_dca$Strategy=="Treat All",]$NB[i] <- sum(df_tdc$event_severe==1)/dim(df_tdc)[1] - (sum(df_tdc$event_severe==0)/dim(df_tdc)[1])*(z[i] / (1-z[i]))
  
}
ggplot(df_dca, aes(x=`Risk Threshold`, y=round(NB,2), group=Strategy, col=Strategy))+
  geom_line()+
  scale_color_manual(values=c("blue", "orange", "black", "grey")) +
  #scale_y_continuous(limits=c(-0.4,0.3), breaks=seq(-0.4, 0.3, 0.1)) +
  scale_y_continuous(limits=c(-0.4,0.1), breaks=seq(-0.4, 0.1, 0.1)) +
  scale_x_continuous(limits=c(0,1), breaks=seq(0, 1, 0.2)) +
  geom_vline(xintercept=c(0.05, 0.1, 0.2), linetype="dashed", col="red") +
  annotate("text", x = c(0.01, 0.13, 0.23), y = -0.4, label = as.character(c(0.05, 0.1, 0.2)),
           vjust = -0.5, color = "red",size = 3) +
  #ggtitle("Decision Curve Analysis") +
  labs(y="Net Benefit") +
  theme_bw()

df_dca[df_dca$`Risk Threshold`%in% c(0.32, 0.33, 0.34),]

######### DCA (Country)

library(ggplot2)
library(tidyr)
df_dca <- data.frame(matrix(ncol=4, nrow=100*2*19))
colnames(df_dca) <- c("Country", "Model", "Risk Threshold", "NB")
df_dca$Country <- rep(c("Total", unique(df_tdc$country)), each=2*100)
df_dca$Model <- rep_len(c(rep("With RE", times=100), rep("Without RE", times=100)), length.out = 100*2*19)
#df_dca$Model <- rep_len(c(rep("Recalibrated", times=100), rep("Refitted", times=100)), length.out = 100*2*19)
df_dca$`Risk Threshold` <- rep_len(1:100, length.out=100*2*19)/100
z <- rep(0,100)
for (c in c("Total", unique(df_dca$Country))) {
  for (i in 1:100) {
    if (c=="Total") {
      df_temp <- df_tdc
    } else {
      df_temp <- df_tdc[df_tdc$country==c,]
    }
    z[i] <- i/100
    df_dca[df_dca$Country==c & df_dca$Model=="Without RE",]$NB[i] <- round(sum(df_temp$event==1 & df_temp$risk_wore>z[i])/dim(df_temp)[1] - (sum(df_temp$event==0 & df_temp$risk_wore>z[i])/dim(df_temp)[1])*(z[i] / (1-z[i])),3)
    df_dca[df_dca$Country==c & df_dca$Model=="With RE",]$NB[i] <- round(sum(df_temp$event==1 & df_temp$risk_re>z[i])/dim(df_temp)[1] - (sum(df_temp$event==0 & df_temp$risk_re>z[i])/dim(df_temp)[1])*(z[i] / (1-z[i])),3)
    #df_dca[df_dca$Country==c & df_dca$Model=="Recalibrated",]$NB[i] <- round(sum(df_temp$event==1 & df_temp$risk_cll>z[i])/dim(df_temp)[1] - (sum(df_temp$event==0 & df_temp$risk_cll>z[i])/dim(df_temp)[1])*(z[i] / (1-z[i])),3)
    #df_dca[df_dca$Country==c & df_dca$Model=="Refitted",]$NB[i] <- round(sum(df_temp$event==1 & df_temp$risk_re>z[i])/dim(df_temp)[1] - (sum(df_temp$event==0 & df_temp$risk_re>z[i])/dim(df_temp)[1])*(z[i] / (1-z[i])),3)
  }
}
df_dca <- df_dca %>% pivot_wider(names_from = "Model", values_from = "NB")
df_dca[df_dca$`Risk Threshold`==0.05,]
df_dca[df_dca$`Risk Threshold`==0.1,]
df_dca[df_dca$`Risk Threshold`==0.2,]

############# Validation #############

# Discrimination ---------------------------------------

alpha <- .05

# Uno's time dependent AUC

Uno_gbsg5 <-
  timeROC::timeROC(
    #T = df_tdc$time2event, 
    T = df_tdc$time2event_severe, 
    
    #delta = df_tdc$event,
    delta = df_tdc$event_severe,
    
    #marker = df_tdc$predicted_exac_probability,
    #marker = df_tdc$predicted_severe_exac_probability,
    #marker = df_tdc$risk_cll,
    marker = df_tdc$risk_cll_severe,
    #marker = df_tdc$risk_re,
    #marker = df_tdc$risk_re_severe,
    
    cause = 1, 
    weighting = "marginal", 
    times = 11.99,
    iid = TRUE
  )

Uno_AUC_res <- c(
  "Uno AUC" = unname(Uno_gbsg5$AUC[2]),
  "2.5 %" = unname(Uno_gbsg5$AUC["t=11.99"] -
                     qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=11.99"]),
  "97. 5 %" = unname(Uno_gbsg5$AUC["t=11.99"] +
                       qnorm(1 - alpha / 2) * Uno_gbsg5$inference$vect_sd_1["t=11.99"])
)

Uno_AUC_res

# Calibration -----------------------------------------


# Observed / Expected ratio
t_horizon <- 12

# Observed
for (c in unique(df_tdc$country)) {
df_byc <- df_tdc %>% filter(country==c)
  obj <- summary(survfit(
  Surv(time2event, event) ~ 1, 
  #Surv(time2event_severe, event_severe) ~ 1, 
  #data = df_tdc),
  data = df_byc),
  times = t_horizon)
  print(c)
print(obs_t <- 1 - obj$surv)
}
# Expected
#exp_t <- mean(df_tdc$predicted_exac_probability)
#exp_t <- mean(df_tdc$predicted_severe_exac_probability)
#exp_t <- mean(df_tdc$risk_cll)
exp_t <- mean(df_tdc$risk_cll_severe)

OE_t <- obs_t / exp_t

alpha <- .05
OE_summary <- c(
  "OE" = OE_t,
  "2.5 %" = OE_t * exp(-qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event)),
  "97.5 %" = OE_t * exp(+qnorm(1 - alpha / 2) * sqrt(1 / obj$n.event))
)

OE_summary


# Calibration plot ----------------------------------
#df_tdc <- df_new
#df_tdc$pred.cll <- log(-log(1 - df_tdc$predicted_exac_probability))
#df_tdc$pred.cll <- log(-log(1 - df_tdc$predicted_severe_exac_probability))
df_tdc$pred.cll <- log(-log(1 - df_tdc$risk_cll))
#df_tdc$pred.cll <- log(-log(1 - df_tdc$risk_cll_severe))
#df_tdc$pred.cll <- log(-log(1 - df_tdc$risk_re))
#df_tdc$pred.cll <- log(-log(1 - df_tdc$risk_re_severe))

#df_temp <- df_tdc
#df_tdc <- df_temp[df_temp$country=="ARG",]
#df_tdc <- df_temp

# Estimate actual risk
dd <- rms::datadist(df_tdc)
options(datadist = "dd")
vcal <- rms::cph(Surv(time2event, event) ~ rcs(pred.cll, 3),
#vcal <- rms::cph(Surv(time2event_severe, event_severe) ~ rcs(pred.cll, 3),
                 x = T,
                 y = T,
                 surv = T,
                 data = df_tdc) 

dat_cal <- cbind.data.frame(
  "obs" = 1 - rms::survest(vcal,
                           times = 12, 
                           newdata = df_tdc)$surv,
  
  "lower" = 1 - rms::survest(vcal, 
                             times = 12, 
                             newdata = df_tdc)$upper,
  
  "upper" = 1 - rms::survest(vcal, 
                             times = 12, 
                             newdata = df_tdc)$lower,
  #"pred" = df_tdc$predicted_exac_probability,
  #"pred" = df_tdc$predicted_severe_exac_probability,
  "pred" = df_tdc$risk_cll,
  #"pred" = df_tdc$risk_cll_severe,
  #"pred" = df_tdc$risk_re,
  #"pred" = df_tdc$risk_re_severe,
  
  "cll" = df_tdc$pred.cll,
  
  "id" = 1:nrow(df_tdc)
)

dat_cal <- dat_cal[order(dat_cal$pred), ]

dev.new() 
par(xaxs = "i", yaxs = "i", las = 1)
# Row 1 = calibration plot, Row 2 = histogram (with x-axis), Row 3 = empty
layout(matrix(c(1, 2), ncol = 1), heights = c(3, 1))

## 1️⃣ Calibration plot (no x-axis)
par(mar = c(0, 4, 2, 2)) # no bottom margin
plot(
  dat_cal$pred, 
  dat_cal$obs,
  type = "l", 
  lty = 1, 
  xlim = c(0, 1),
  ylim = c(0, 1), 
  #xlim = c(0, 0.5),
  #ylim = c(0, 0.5),
  lwd = 2,
  xaxt = "n", # hide x-axis
  xlab = "",
  ylab = "Predicted risk from proxy model", 
  bty = "n"
)
lines(dat_cal$pred, dat_cal$lower, lty = 2, lwd = 2)
lines(dat_cal$pred, dat_cal$upper, lty = 2, lwd = 2)
abline(0, 1, lwd = 2, lty = 2, col = 2)

legend("bottomright",
       c("Ideal calibration",
         "Calibration curve based on secondary Cox model",
         "95% confidence interval"),
       col = c(2, 1, 1),
       lty = c(2, 1, 2),
       lwd = c(2, 2, 2),
       bty = "n",
       cex = 0.85)

## 2️⃣ Histogram with x-axis
par(mar = c(4, 4, 0, 2)) # normal bottom margin
hist(dat_cal$pred,
     breaks = 20,
     col = "grey80",
     border = "white",
     xlim = c(0, 1),
     #xlim = c(0, 0.5),
     main = "",
     #xlab = "Predicted risk from ACCEPT 2.0",
     xlab = "Predicted risk from ACCEPT 3.0",
     #xlab = "Predicted risk from refitted model",
     ylab = "Count")

# calibratoin slope (fixed time point)-------------------------------------

df_tdc$pred.cll <- log(-log(1 - df_tdc$predicted_exac_probability))
(gval <- coxph(Surv(time2event, event) ~ pred.cll, data=df_tdc)) 
#(gval <- glm(pseudo_obs~offset(pred.cll)+pred.cll, data=df_tdc_cen, family=gaussian))

calslope_summary <- c(
  "calibration slope" = gval$coef,
  "2.5 %"  = gval$coef - qnorm(1 - alpha / 2) * sqrt(gval$var),
  "97.5 %" = gval$coef + qnorm(1 - alpha / 2) * sqrt(gval$var)
)
calslope_summary

# KM Plot

library(survival)
install.packages("survminer")
library(survminer)
library(dplyr)

km_fit <- survfit(Surv(time2event, df_tdc$event) ~ 1, data = df_tdc)
#km_fit <- survfit(Surv(time2event_severe, df_tdc$event_severe) ~ 1, data = df_tdc)
ggsurvplot(km_fit, data = df_tdc,
           conf.int = FALSE,       # show confidence interval
           risk.table = T,     # add risk table
           palette = "black",
           xlab = "Time (months)",
           ylab = "Survival Probability",
           title = "Kaplan-Meier Survival Curve")
