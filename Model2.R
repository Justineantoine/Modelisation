rm(list=ls())
library(deSolve)
library(FME)
library(plotrix)
library(Hmisc)
library(corrplot)

dat <- read.table(file = "CLINICAL.csv", sep= ',', header = TRUE)
######################################
# BODYWEIGHT, FAT MASS AND LEAN MASS #
######################################
BW0 <- dat$BW0
BW2 <- dat$BW2
BW5 <- dat$BW5

FM0 <- dat$TBFD0
FM2 <- dat$TBFD2
FM5 <- dat$TBFD5


adjust = subset(FM0/dat$TBFB0, is.na(FM0/dat$TBFB0)==F)
mean_adjust = mean(adjust)

FM0[18] <- mean_adjust*dat$TBFB0[18]

#On a BW = LM + FM

LM0 <- BW0-FM0
LM2 <- BW2-FM2
LM5 <- BW5-FM5

##################
# AGE AND HEIGHT #
##################
age0 <- dat$Age0
age2 <- dat$Age2
age5 <- dat$Age5
H <- dat$Height

BMI0 <- dat$BMI0
BMI1 <- dat$BMI1
BMI2 <- dat$BMI2
BMI5 <- dat$BMI5

adjust2 = subset(BMI1-BMI2, is.na(BMI1)==F)
mean_adjust2 = mean(adjust2)
BMI1[9] <- BMI2[9] - mean_adjust2
BMI1[20] <- BMI2[20] - mean_adjust2
BMI1[21] <- BMI2[21] - mean_adjust2

########
# TIME #
########
T0 <- rep(0,41)
T2 <- dat$T2
T1 <- T2/2
T5 <- dat$T5

##########################
# LIPID AGE AND TURNOVER #
##########################
LA0 <- dat$LA0*365
LA2 <- dat$LA2*365
LA5 <- dat$LA5*365

kout0 <- dat$TO0
kout2 <- dat$TO2
kout5 <- dat$TO5

kin0 <- kout0 * FM0
kin2 <- kout2 * FM2
kin5 <- kout5 * FM5

##############
# PARAMETERS #
##############
pf <- 39500 #kJ/kg
pl <- 7600 #kJ/kg
gf <- 13 #kJ/kg/day
gl <- 92 #kJ/kg/day
nf <- 750 #kJ/kg 
nl <- 960 #kJ/kg
Btef <- 0.1
Bat <- 0.14
C <- 10.4*pl/pf #kg
tau <- 14 #days

#####################
# PHYSICAL ACTIVITY #
#####################
d0 <- ((1-Btef)*1.5-1)*(10*BW0+625*H-5*age0-161)*4.184
d2 <- ((1-Btef)*1.5-1)*(10*BW2+625*H-5*age2-161)*4.184
d5 <- ((1-Btef)*1.5-1)*(10*BW5+625*H-5*age5-161)*4.184

#######################
# ENERGY PARTITIONING #
#######################
p2 <- C/(C+FM2) 
p5 <- C/(C+FM5)

#########################
# INITIAL ENERGY INTAKE #
#########################
EI0 <- (10*BW0 + 625*H -5*age0-161)*1.5*4.184
EIsurg <- c()
EIfinal <- c()

EI <- function(t, Ts, Tf, EIi, EIs, EIf) #i pour l'individu auquel on s'interesse
{
  s = (EIf-EIs)/(Tf-Ts) #garde le regime pendant 150 jours (~5 mois) puis reprise progressive de nourriture jusqu'au 900eme jour (~3 ans)
  if (t<=0){
    EIi
    
  }
  else if (t<=Ts){
    EIs
    
  }
  else{
    intake = (t-500)*s+EIs
    if (intake<EIf){
      intake
    }
    else{
      EIf
      
    }
  }
}

##############################
# INITIAL ENERGY EXPENDITURE #
##############################
EE0 <- EI0

EE <- function(t, Ts, Tf, EIi, EIs, EIf, k, h, a, FM, LM){
  B <- Btef+Bat*(1-exp(-t/tau))
  PA <- ((1-Btef)*1.5-1)*4.184*(10*(FM+LM)+625*h-5*(a+t/365)-161)
  partF <- FM/(C+FM)
  partL <- 1-partF
  result <- (k + gf*FM + gl*LM + PA - EIi*B + EI(t, Ts, Tf, EIi, EIs, EIf)*(B + partF*nf/pf + partL*nl/pl))/(1 + partF*nf/pf + partL*nl/pl)
  result
  }

##############################
# K CONSTANT AND K2 CONSTANT #
##############################
K <- EI0 - gf*FM0 - gl*LM0 - d0
K2 <- 5e-6

##########################
# EDO SYSTEM : BW AND LA #
##########################

  # # # #
  # FIT #
  # # # #

EqBW <- function(t, y, parameters){

  EIi <- parameters[1]
  k <- parameters[2]
  h <- parameters[3]
  a <- parameters[4]
  Ts <- parameters[5]
  Tf <- parameters[6]
  EIs <- parameters[7]
  EIf <- parameters[8]
  partF <- y[1]/(C+y[1])
  partL <- 1-partF

  dF <- partF/pf * (EI(t, Ts, Tf, EIi, EIs, EIf) - EE(t, Ts, Tf, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dL <- partL/pl * (EI(t, Ts, Tf, EIi, EIs, EIf) - EE(t, Ts, Tf, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  list(c(dF, dL)) 
}

soltime <- seq(0, 2200)
for (i in 1:41){
  parameters <- c(EI0[i], K[i], H[i], age0[i], 500, 900)
  init <- c(fatmass = FM0[i],leanmass = LM0[i])
  Data <- data.frame(time = c(T2[i],T5[i]), fatmass = c(FM2[i],FM5[i]) , leanmass= c(LM2[i],LM5[i]))
  
  modelcost <- function(P) {
    sol <- lsoda(y=init, times=soltime, func = EqBW, parms = c(parameters,P[1],P[2]))
    return(modCost(sol,Data))
  }
  
  Fit <- modFit(f = modelcost, p = c(7500,9800))

  EIsurg[i] <- Fit$par[1]
  EIfinal[i] <- Fit$par[2]
}
  
  # # # # # # #
  # LIPID AGE #
  # # # # # # #
EqBWLA <- function(t, y, parameters){
  
  EIi <- parameters[1]
  k <- parameters[2]
  h <- parameters[3]
  a <- parameters[4]
  Ts <- parameters[5]
  Tf <- parameters[6]
  EIs <- parameters[7]
  EIf <- parameters[8]
  partF <- y[1]/(C+y[1])
  partL <- 1-partF
  kin <- K2 * EI(t, Ts, Tf, EIi, EIs, EIf)
  
  dF <- partF/pf * (EI(t, Ts, Tf, EIi, EIs, EIf) - EE(t, Ts, Tf, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dL <- partL/pl * (EI(t, Ts, Tf, EIi, EIs, EIf) - EE(t, Ts, Tf, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dA <- 1 - kin*y[3]/y[1]

  list(c(dF, dL, dA)) 
}

#################################################
# PLOT FM AND BODY WEIGHT, EE AND EI, LIPID AGE #
#################################################
for (k in 1:41){
  graphtime <- seq(-100,2200)
  parameters <- c(EI0[k], K[k], H[k], age0[k], 500, 900, EIsurg[k], EIfinal[k])
  init <- c(fatmass = FM0[k],leanmass = LM0[k], age =LA0[k])
  bestfit <- lsoda(y=init, times=soltime, func = EqBWLA, parms = c(parameters))
  
  # # # # # # #
  # EI AND EE #
  # # # # # # #
  
  graphEI= rep(EI0[k], 100)
  graphEE = rep(EE0[k], 100)
  for (i in 101:2301){
    graphEI[i] <- EI(graphtime[i], 500, 900, EI0[k], EIsurg[k], EIfinal[k])
    graphEE[i] <- EE(graphtime[i], 500, 900, EI0[k], EIsurg[k], EIfinal[k], K[k], H[k], age0[k], bestfit[i-100,2], bestfit[i-100,3])
  }
  
  plot(graphtime, graphEI, type="l", xlab="Days", ylab="Energy rate ", col=1, ylim=c(1000, 15000))
  lines(graphtime, graphEE, type ="l", lty = 1, col=2)
  title(main=c("Energy rates for patient : ", k))
  legend("bottomright",lty=c(1,1), cex=0.7, col=c(1,2), legend=c("Energy Intake rate", "Energy Expenditure rate"))
  
  # # # # # # # # #
  # BW, LM AND FM #
  # # # # # # # # #
  
  graphFM <- c(rep(FM0[k], 100), bestfit[,2])
  graphLM <- c(rep(LM0[k], 100), bestfit[,3])
  graphBW <- c(rep(BW0[k], 100), bestfit[,2]+bestfit[,3])
  
  plot(graphtime, graphFM, type="l", ylim=c(0,160), xlab="Days", ylab="Weight in kg")
  lines(graphtime, graphBW, type = "l", lty =1, col=4)
  lines(graphtime, graphLM, type = "l", lty =1, col="dodgerblue4")
  points(T0[k], FM0[k], pch=19)
  points(T2[k], FM2[k], pch=19)
  points(T5[k], FM5[k], pch=19)
  points(T0[k], BW0[k], pch=19, col=4)
  points(T2[k], BW2[k], pch=19, col=4)
  points(T5[k], BW5[k], pch=19, col=4)
  points(T0[k], LM0[k], pch=3, col="dodgerblue4")
  points(T2[k], LM2[k], pch=3, col="dodgerblue4")
  points(T5[k], LM5[k], pch=3, col="dodgerblue4")
  title(main=c("BodyWeight time course for patient : ", k))
  legend("bottomright", lty=c(1,1,1), legend=c("BW", "FM","LM"), col=c(4,1,"dodgerblue4"), cex=0.7)
  
  # # # # # # #
  # LIPID AGE #
  # # # # # # #
  
  graphA <- c(rep(LA0[k], 100), bestfit[,4])
  plot(graphtime, graphA, type="l", xlab="Days", ylab="Lipid Age", ylim=c(0, 2700))
  title(main=c("Lipid age for patient : ", k))
}

####################
# AVERAGE PATIENTS #
####################

dBMI <- BMI0 - BMI5
tertile<- unname(quantile(dBMI, probs = c(0.33, 0.67)))

gr <- c()
gs <- c()

for (i in 1:41){
  if (dBMI[i] <= tertile[1]){
    gr <- c(gr, i)
  }
  else{
    gs <- c(gs, i)
  }
}

st_error<- function(x){
  result <- sd(x)/sqrt(length(x))
  result
}

Tg <- c(0,1,2,5)
BMIgr <- c(mean(BMI0[gr]), mean(BMI1[gr]), mean(BMI2[gr]), mean(BMI5[gr]))
errgr <- c(st_error(BMI0[gr]), st_error(BMI1[gr]), st_error(BMI2[gr]), st_error(BMI5[gr]))
BMIgs <- c(mean(BMI0[gs]), mean(BMI1[gs]), mean(BMI2[gs]), mean(BMI5[gs]))
errgs <- c(st_error(BMI0[gs]), st_error(BMI1[gs]), st_error(BMI2[gs]), st_error(BMI5[gs]))


plot(Tg, BMIgr, type = "l", ylim=c(25, 47), xlab="Year of investigation", ylab="BMI (kg/m²)")
plotCI(Tg, BMIgr, uiw = errgr, lwd =2, col = 1, add =T)
lines(Tg, BMIgs, col=2)
plotCI(Tg, BMIgs, uiw = errgs, lwd =2, col = 2, add =T)
legend("topright", cex=0.7, lty=c(1,1), col=c(1,2), legend=c("Tertile 1 : rebounders", "Tertile 2-3 : weight stable"))

                  # # # # # # #
                  # REBOUNDER #
                  # # # # # # #

#PARAMETERS AND DATA
R_n <- length(gr)
R_age0 <- mean(age0[gr])
R_BW0 <- mean(BW0[gr])
R_FM0 <- mean(FM0[gr])
R_LM0 <- mean(LM0[gr])
R_LA0 <- mean(LA0[gr])
R_H <- mean(H[gr])
R_d0 <- ((1-Btef)*1.5-1)*(10*R_BW0+625*R_H-5*R_age0-161)*4.184
R_EI0 <- (10*R_BW0 + 625*R_H -5*R_age0-161)*1.5*4.184
R_EE0 <- R_EI0
R_K <- R_EI0 - gf*R_FM0 - gl*R_LM0 - R_d0
R_parameters <- c(R_EI0, R_K, R_H, R_age0, 380, 1500)
R_init <- c(fatmass = R_FM0,leanmass = R_LM0)
R_Data <- data.frame(time = c(T2[gr],T5[gr]), fatmass = c(FM2[gr],FM5[gr]) , leanmass= c(LM2[gr],LM5[gr]))
soltime = seq(0,2200)


#FIT AND STATISTICS
R_modelcost <- function(P) {
  sol <- lsoda(y=R_init, times=soltime, func = EqBW, parms = c(R_parameters,P[1],P[2]))
  return(modCost(sol,R_Data))
}
R_Fit <- modFit(f = R_modelcost, p = c(mean(EIsurg[gr]),mean(EIfinal[gr])))
R_sssurg <- 0
R_ssfinal <- 0
for (i in 1:length(gr)){
  R_sssurg <- R_sssurg + (EIsurg[i]-mean(EIsurg[gr]))^2
  R_ssfinal <- R_ssfinal + (EIfinal[i]-mean(EIfinal[gr]))^2
}
R_mse <- mean(R_Fit$residuals^2)
R_sesurg <- sqrt(R_mse * (1/R_n + 1 + (7000-mean(EIsurg[gr])^2/R_sssurg)))
R_sefinal <- sqrt(R_mse * (1/R_n + 1 + (7000-mean(EIfinal[gr])^2/R_ssfinal)))


#BESTFIT AND PREDICTION INTERVAL
R_init <- c(fatmass = R_FM0,leanmass = R_LM0, age = R_LA0)
R_EIsurg <- R_Fit$par[1]
R_EIfinal <- R_Fit$par[2]
R_EIsurgCI <- unname(c(R_EIsurg + R_sesurg*qt(0.05, R_n-2), R_EIsurg + R_sesurg*qt(0.95, R_n-2)))
R_EIfinalCI <- unname(c(R_EIfinal + R_sefinal*qt(0.05, R_n-2), R_EIfinal + R_sefinal*qt(0.95, R_n-2)))
R_bestfit <- lsoda(y=R_init, times=soltime, func = EqBWLA, parms = c(R_parameters,R_EIsurg, R_EIfinal))
R_fitinf <- lsoda(y=R_init, times=soltime, func = EqBWLA, parms = c(R_parameters,R_EIsurgCI[1], R_EIfinalCI[1]))
R_fitsup <- lsoda(y=R_init, times=soltime, func = EqBWLA, parms = c(R_parameters,R_EIsurgCI[2], R_EIfinalCI[2]))

    #bmi
R_BMI <- (R_bestfit[,2]+R_bestfit[,3])/R_H^2
R_BMIinf <- (R_fitinf[,2]+R_fitinf[,3])/R_H^2
R_BMIsup <- (R_fitsup[,2]+R_fitsup[,3])/R_H^2

    #EI and EE
R_EI <- c() ; R_EIinf <- c() ; R_EIsup <- c()
R_EE <- c() ; R_EEinf <- c() ; R_EEsup <- c()
for (i in soltime){
  R_EI <- c(R_EI, EI(i, 500, 900, R_EI0, R_EIsurg, R_EIfinal))
  R_EIinf <- c(R_EIinf, EI(i, 500, 900, R_EI0, R_EIsurgCI[1], R_EIfinalCI[1]))
  R_EIsup <- c(R_EIsup, EI(i, 500, 900, R_EI0, R_EIsurgCI[2], R_EIfinalCI[2]))
  
  R_EE <- c(R_EE, EE(i, 500, 900, R_EI0, R_EIsurg, R_EIfinal, R_K, R_H, R_age0, R_bestfit[i+1,2], R_bestfit[i+1,3]))
  R_EEinf <- c(R_EEinf, unname(EE(i, 500, 900, R_EI0, R_EIsurgCI[1], R_EIfinalCI[1], R_K, R_H, R_age0, R_fitinf[i+1,2], R_fitinf[i+1,3])))
  R_EEsup <- c(R_EEsup, unname(EE(i, 500, 900, R_EI0, R_EIsurgCI[2], R_EIfinalCI[2], R_K, R_H, R_age0, R_fitsup[i+1,2], R_fitsup[i+1,3])))
}


#BMI GRAPHS
plot(soltime, R_BMI, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50))
lines(soltime, R_BMIinf, type='l', lty=2)
lines(soltime, R_BMIsup, type='l', lty=2)
points(T2[gr], BMI2[gr])
points(T5[gr], BMI5[gr])
title(main="BMI")
legend("topleft", cex=0.7, lty=c(1,2), col=c(1,1), legend=c("Average BMI", "Expected inter-individual BMI variability"))


#ENERGY RATES GRAPHS
plot(soltime, R_EI, type="l", ylim=c(7000, 14000), xlab="Days", ylab="Energy rates (kJ/d)")
lines(soltime, R_EE, type="l", col=2)
title(main="Energy rates")
legend("bottomright", cex = 0.7, lty=c(1,1), col=c(1,2), legend=c("Average energy intake rate", "Average energy expenditure rate"))

plot(soltime, R_EI, type="l", ylim=c(7000, 14000), xlab="Days", ylab="Energy Intake rate (kJ/d)")
lines(soltime, R_EIinf, lty=2)
lines(soltime, R_EIsup, lty=2)
title(main="Energy rates")
legend("bottomright", cex = 0.7, lty=c(1,2), legend=c("Average energy intake rate", "Expected inter-individual EI variability"))

plot(soltime, R_EE, type="l", col=2, ylim=c(8500, 13000), xlab="Days", ylab="Energy Expenditure rate (kJ/d)")
lines(soltime, R_EEinf, col=2, lty=2)
lines(soltime, R_EEsup, col=2, lty=2)
title(main="Energy rates")
legend("bottomright", cex = 0.7, lty=c(1,2), col=c(2,2), legend=c("Average energy expenditure rate", "Expected inter-individual EE variability"))

#LIPID AGE GRAPH
plot(soltime, R_bestfit[,4], type ="l", xlab = "Days", ylab="Lipid Age (d)", ylim=c(500, 1500))
lines(soltime, R_fitinf[,4], lty=2)
lines(soltime, R_fitsup[,4], lty=2)
title(main="Lipid Age")


                  # # # # # # # # #
                  # WEIGHT STABLE #
                  # # # # # # # # #

#PARAMETERS AND DATA
S_n <- length(gs)
S_age0 <- mean(age0[gs])
S_BW0 <- mean(BW0[gs])
S_FM0 <- mean(FM0[gs])
S_LM0 <- mean(LM0[gs])
S_LA0 <- mean(LA0[gs])
S_H <- mean(H[gs])
S_d0 <- ((1-Btef)*1.5-1)*(10*S_BW0+625*S_H-5*S_age0-161)*4.184
S_EI0 <- (10*S_BW0 + 625*S_H -5*S_age0-161)*1.5*4.184
S_EE0 <- S_EI0
S_K <- S_EI0 - gf*S_FM0 - gl*S_LM0 - S_d0
S_parameters <- c(S_EI0, S_K, S_H, S_age0, 500, 900)
S_init <- c(fatmass = S_FM0,leanmass = S_LM0)
S_Data <- data.frame(time = c(T2[gs],T5[gs]), fatmass = c(FM2[gs],FM5[gs]) , leanmass= c(LM2[gs],LM5[gs]))
soltime = seq(0,2200)

#FIT AND STATISTICS
S_modelcost <- function(P) {
  sol <- lsoda(y=S_init, times=soltime, func = EqBW, parms = c(S_parameters,P[1],P[2]))
  return(modCost(sol,S_Data))
}
S_Fit <- modFit(f = S_modelcost, p = c(mean(EIsurg[gs]),mean(EIfinal[gs])))
S_EIsurg <- S_Fit$par[1]
S_EIfinal <- S_Fit$par[2]
S_sssurg <- 0
S_ssfinal <- 0
for (i in 1:length(gs)){
  S_sssurg <- S_sssurg + (EIsurg[i]-mean(EIsurg[gs]))^2
  S_ssfinal <- S_ssfinal + (EIfinal[i]-mean(EIfinal[gs]))^2
}
S_mse <- mean(S_Fit$residuals^2)
S_sesurg <- sqrt(S_mse * (1/S_n + 1 + (7000-mean(EIsurg[gs])^2/S_sssurg)))
S_sefinal <- sqrt(S_mse * (1/S_n + 1 + (7000-mean(EIfinal[gs])^2/S_ssfinal)))

#BESTFIT AND PREDICTION INTERVAL

S_init <- c(fatmass = S_FM0,leanmass = S_LM0, age = S_LA0)
S_EIsurgCI <- unname(c(S_EIsurg + S_sesurg*qt(0.05, S_n-2), S_EIsurg + S_sesurg*qt(0.95, S_n-2)))
S_EIfinalCI <- unname(c(S_EIfinal + S_sefinal*qt(0.05, S_n-2), S_EIfinal + S_sefinal*qt(0.95, S_n-2)))
S_bestfit <- lsoda(y=S_init, times=soltime, func = EqBWLA, parms = c(S_parameters,S_EIsurg, S_EIfinal))
S_fitinf <- lsoda(y=S_init, times=soltime, func = EqBWLA, parms = c(S_parameters,S_EIsurgCI[1], S_EIfinalCI[1]))
S_fitsup <- lsoda(y=S_init, times=soltime, func = EqBWLA, parms = c(S_parameters,S_EIsurgCI[2], S_EIfinalCI[2]))

    #bmi
S_BMI <- (S_bestfit[,2]+S_bestfit[,3])/S_H^2
S_BMIinf <- (S_fitinf[,2]+S_fitinf[,3])/S_H^2
S_BMIsup <- (S_fitsup[,2]+S_fitsup[,3])/S_H^2

    #EI and EE
S_EI <- c() ; S_EIinf <- c() ; S_EIsup <- c()
S_EE <- c() ; S_EEinf <- c() ; S_EEsup <- c()
for (i in soltime){
  S_EI <- c(S_EI, EI(i, 500, 900, S_EI0, S_EIsurg, S_EIfinal))
  S_EIinf <- c(S_EIinf, EI(i, 500, 900, S_EI0, S_EIsurgCI[1], S_EIfinalCI[1]))
  S_EIsup <- c(S_EIsup, EI(i, 500, 900, S_EI0, S_EIsurgCI[2], S_EIfinalCI[2]))
  
  S_EE <- c(S_EE, EE(i, 500, 900, S_EI0, S_EIsurg, S_EIfinal, S_K, S_H, S_age0, S_bestfit[i+1,2], S_bestfit[i+1,3]))
  S_EEinf <- c(S_EEinf, unname(EE(i, 500, 900, S_EI0, S_EIsurgCI[1], S_EIfinalCI[1], S_K, S_H, S_age0, S_fitinf[i+1,2], S_fitinf[i+1,3])))
  S_EEsup <- c(S_EEsup, unname(EE(i, 500, 900, S_EI0, S_EIsurgCI[2], S_EIfinalCI[2], S_K, S_H, S_age0, S_fitsup[i+1,2], S_fitsup[i+1,3])))
}

#BMI GRAPHS
plot(soltime, S_BMI, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50), col=4)
lines(soltime, S_BMIinf, type='l', lty=2, col=4)
lines(soltime, S_BMIsup, type='l', lty=2, col=4)
points(T2[gs], BMI2[gs], col=4)
points(T5[gs], BMI5[gs], col=4)
title(main="BMI")
legend("topleft", cex=0.7, lty=c(1,2), col=c(4,4), legend=c("Average BMI", "Expected inter-individual BMI variability"))

#ENERGY RATES GRAPHS
plot(soltime, S_EI, type="l", ylim=c(7000, 14000), xlab="Days", ylab="Energy rates (kJ/d)")
lines(soltime, S_EE, type="l", col=2)
title(main="Energy rates")
legend("bottomright", cex = 0.7, lty=c(1,1), col=c(1,2), legend=c("Average energy intake rate", "Average energy expenditure rate"))

plot(soltime, S_EI, type="l", ylim=c(5800, 12000), xlab="Days", ylab="Energy intake rate (kJ/d)")
lines(soltime, S_EIinf, lty=2)
lines(soltime, S_EIsup, lty=2)
title(main="Energy rates")
legend("bottomright", cex = 0.7, lty=c(1,2), legend=c("Average energy intake rate", "Expected inter-individual EI variability"))

plot(soltime, S_EE, type="l", col=2, ylim=c(7500, 12000), xlab="Days", ylab="Energy expenditure rate (kJ/d)")
lines(soltime, S_EEinf, col=2, lty=2)
lines(soltime, S_EEsup, col=2, lty=2)
title(main="Energy rates")
legend("bottomright", cex = 0.7, lty=c(1,2), col=c(2,2), legend=c("Average energy expenditure rate", "Expected inter-individual EE variability"))


#LIPID AGE GRAPH
plot(soltime, S_bestfit[,4] , type ="l", xlab = "Days", ylab="Lipid Age (d)", ylim=c(500, 1500))
lines(soltime, S_fitinf[,4], lty=2)
lines(soltime, S_fitsup[,4], lty=2)
title(main="Lipid Age")


############################
# MEAN COMPAR TEST R AND S #
############################
comp <- t.test(BMI0[g1], BMI0[g4])


#######################
# COMPARISON BOXPLOTS #
#######################
boxplot()

####################
# CORRELATION TEST #
####################
EE <- cor
for (i in 1:41){
  EEcor[i] <- EI0[i] - EI(T2[i], 500, 900, EI0[i], EIsurg[i], EIfinal[i])
  EI5[i] <- EI0[i] - EI(T5[i], 500, 900, EI0[i], EIsurg[i], EIfinal[i])
}

cordat <- matrix(c(EI0, EI2, EI5, kin0, kin0-kin2, kin0-kin5), nrow = 41, ncol= 6, byrow=F)
colnames(cordat) <- c("EI0", "EI2", "EI5","kin0", "kin2", "kin5")
mcor <- rcorr(cordat)$r[1:3,4:6]
corrplot(mcor, type="upper", order="hclust", tl.col="black", tl.srt=45)
pairs(cordat)
