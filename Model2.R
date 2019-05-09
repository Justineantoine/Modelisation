rm(list=ls())
library(deSolve)
library(FME)
library(plotrix)

dat <- read.table(file = "CLINICAL(copietravail).csv", sep= ',', header = TRUE)
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
adjust2 = subset(BW1-BW2, is.na(BW1)==F)
mean_adjust2 = mean(adjust2)

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

EI <- function(t, EIi, EIs, EIf) #i pour l'individu auquel on s'interesse
{
  s = (EIf-EIs)/(900-500) #garde le regime pendant 150 jours (~5 mois) puis reprise progressive de nourriture jusqu'au 900eme jour (~3 ans)
  if (t<=0){
    EIi
    
  }
  else if (t<=500){
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

EE <- function(t,EIi, EIs, EIf, k, h, a, FM, LM){
  B <- Btef+Bat*(1-exp(-t/tau))
  PA <- ((1-Btef)*1.5-1)*4.184*(10*(FM+LM)+625*h-5*(a+t/365)-161)
  partF <- FM/(C+FM)
  partL <- 1-partF
  result <- (k + gf*FM + gl*LM + PA - EIi*B + EI(t, EIi, EIs, EIf)*(B + partF*nf/pf + partL*nl/pl))/(1 + partF*nf/pf + partL*nl/pl)
  result
  }

##############
# K CONSTANT #
##############
K <- EI0 - gf*FM0 - gl*LM0 - d0

##############
# EDO SYSTEM #
##############

EqBW <- function(t, y, parameters){

  EIi <- parameters[1]
  k <- parameters[2]
  h <- parameters[3]
  a <- parameters[4]
  EIs <- parameters[5]
  EIf <- parameters[6]
  partF <- y[1]/(C+y[1])
  partL <- 1-partF
  
  dF <- partF/pf * (EI(t, EIi, EIs, EIf) - EE(t,EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dL <- partL/pl * (EI(t, EIi, EIs, EIf) - EE(t,EIi, EIs, EIf, k, h, a, y[1], y[2]))
  list(c(dF, dL))
}

soltime <- seq(0, 2200)

for (i in 1:41){
  parameters <- c(EI0[i], K[i], H[i], age0[i])
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

######################################
# PLOT FM AND BODY WEIGHT, EE AND EI #
######################################
for (k in 1:41){
  ind <- k
  graphtime <- seq(-100,2200)
  parameters <- c(EI0[ind], K[ind], H[ind], age0[ind], EIsurg[ind], EIfinal[ind])
  init <- c(fatmass = FM0[ind],leanmass = LM0[ind])
  bestfit <- lsoda(y=init, times=soltime, func = EqBW, parms = c(parameters))
  
  graphEI= rep(EI0[ind], 100)
  graphEE = rep(EE0[ind], 100)
  for (i in 101:2301){
    graphEI[i] <- EI(graphtime[i], EI0[ind], EIsurg[ind], EIfinal[ind])
    graphEE[i] <- EE(graphtime[i], EI0[ind], EIsurg[ind], EIfinal[ind], K[ind], H[ind], age0[ind], bestfit[i-100,2], bestfit[i-100,3])
  }
  
  plot(graphtime, graphEI, type="l", xlab="Days", ylab="Energy rate ", col=1, ylim=c(1000, 15000))
  lines(graphtime, graphEE, type ="l", lty = 1, col=2)
  legend("bottomright",lty=c(1,1), cex=0.7, col=c(1,2), legend=c("Energy Intake rate", "Energy Expenditure rate"))
  
  
  graphFM <- c(rep(FM0[ind], 100), bestfit[,2])
  graphLM <- c(rep(LM0[ind], 100), bestfit[,3])
  graphBW <- c(rep(BW0[ind], 100), bestfit[,2]+bestfit[,3])
  
  plot(graphtime, graphFM, type="l", ylim=c(0,160), xlab="Days", ylab="Weight in kg")
  lines(graphtime, graphBW, type = "l", lty =1, col=4)
  lines(graphtime, graphLM, type = "l", lty =1, col="dodgerblue4")
  points(T0[ind], FM0[ind], pch=19)
  points(T2[ind], FM2[ind], pch=19)
  points(T5[ind], FM5[ind], pch=19)
  points(T0[ind], BW0[ind], pch=19, col=4)
  points(T2[ind], BW2[ind], pch=19, col=4)
  points(T5[ind], BW5[ind], pch=19, col=4)
  points(T0[ind], LM0[ind], pch=3, col="dodgerblue4")
  points(T2[ind], LM2[ind], pch=3, col="dodgerblue4")
  points(T5[ind], LM5[ind], pch=3, col="dodgerblue4")
  title(main= k)
  legend("bottomright", lty=c(1,1,1), legend=c("BW", "FM","LM"), col=c(4,1,"dodgerblue4"), cex=0.7)
  
}

####################
# AVERAGE PATIENTS #
####################

dBMI <- BMI0 - BMI5
tertile<- unname(quantile(dBMI, probs = c(0.33, 0.67)))

g1 <- c()
g2 <- c()
g3 <- c()

for (i in 1:41){
  if (dBMI[i] <= tertile[1]){
    g1 <- c(g1, i)
  }
  else if (dBMI[i] <= tertile[2]){
    g2 <- c(g2, i)
  }
  else {
    g3 <- c(g3, i)
  }
}

st_error<- function(x){
  result <- sd(x)/sqrt(length(x))
  result
}

Tg <- c(0,1,2,5)
BMIg1 <- c(mean(BMI0[g1]), mean(BMI1[g1]), mean(BMI2[g1]), mean(BMI5[g1]))
errg1 <- c(st_error(BMI0[g1]), st_error(BMI1[g1]), st_error(BMI2[g1]), st_error(BMI5[g1]))
BMIg2 <- c(mean(BMI0[g2]), mean(BMI1[g2]), mean(BMI2[g2]), mean(BMI5[g2]))
errg2 <- c(st_error(BMI0[g2]), st_error(BMI1[g2]), st_error(BMI2[g2]), st_error(BMI5[g2]))
BMIg3 <- c(mean(BMI0[g3]), mean(BMI1[g3]), mean(BMI2[g3]), mean(BMI5[g3]))
errg3 <- c(st_error(BMI0[g3]), st_error(BMI1[g3]), st_error(BMI2[g3]), st_error(BMI5[g3]))

plot(Tg, BMIg1, type = "l", ylim=c(25, 50), xlab="Year of investigation", ylab="BMI (kg/m²)")
plotCI(Tg, BMIg1, uiw = errg1, lwd =2, col = 1, add =T)
lines(Tg, BMIg2, col=2)
plotCI(Tg, BMIg2, uiw = errg2, lwd =2, col = 2, add =T)
lines(Tg, BMIg3, col=3)
plotCI(Tg, BMIg3, uiw = errg3, lwd =2, col = 3, add =T)

# # # # # #
# GROUP 1 #
# # # # # #
A1_age0 <- mean(age0[g1])
A1_BW0 <- mean(BW0[g1])
A1_FM0 <- mean(FM0[g1])
A1_LM0 <- mean(LM0[g1])
A1_H <- mean(H[g1])
A1_d0 <- ((1-Btef)*1.5-1)*(10*A1_BW0+625*A1_H-5*A1_age0-161)*4.184
A1_EI0 <- (655 + 9.56*A1_BW0 + 186*A1_H -4.68*A1_age0)*1.33*4.184
A1_EE0 <- A1_EI0
A1_K <- A1_EI0 - gf*A1_FM0 - gl*A1_LM0 - A1_d0

A1_parameters <- c(A1_EI0, A1_K, A1_H, A1_age0)
A1_init <- c(fatmass = A1_FM0,leanmass = A1_LM0)
A1_Data <- data.frame(time = c(T2[g1],T5[g1]), fatmass = c(FM2[g1],FM5[g1]) , leanmass= c(LM2[g1],LM5[g1]))
soltime = seq(0,2200)

A1_modelcost <- function(P) {
  sol <- lsoda(y=A1_init, times=soltime, func = EqBW, parms = c(A1_parameters,P[1],P[2]))
  return(modCost(sol,A1_Data))
}

A1_Fit <- modFit(f = A1_modelcost, p = c(mean(EIsurg[g1]),mean(EIfinal[g1])))

A1_EIsurg <- A1_Fit$par[1]
A1_EIfinal <- A1_Fit$par[2]
A1_EIsurgCI <- unname(c(A1_EIsurg + summary(A1_Fit)$par[1,2]*qt(0.05, 2), A1_EIsurg + summary(A1_Fit)$par[1,2]*qt(0.95, 2)))
A1_EIfinalCI <- unname(c(A1_EIfinal + summary(A1_Fit)$par[2,2]*qt(0.05, 2), A1_EIfinal + summary(A1_Fit)$par[2,2]*qt(0.95, 2)))
A1_bestfit <- lsoda(y=A1_init, times=soltime, func = EqBW, parms = c(A1_parameters,A1_EIsurg, A1_EIfinal))
A1_fitinf <- lsoda(y=A1_init, times=soltime, func = EqBW, parms = c(A1_parameters,A1_EIsurgCI[1], A1_EIfinalCI[1]))
A1_fitsup <- lsoda(y=A1_init, times=soltime, func = EqBW, parms = c(A1_parameters,A1_EIsurgCI[2], A1_EIfinalCI[2]))

plot(soltime, (A1_bestfit[,2]+A1_bestfit[,3])/A1_H^2, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50))
lines(soltime, (A1_fitinf[,2]+A1_fitinf[,3])/A1_H^2, type='l', lty=2)
lines(soltime, (A1_fitsup[,2]+A1_fitsup[,3])/A1_H^2, type='l', lty=2)
points(T2[g1], BMI2[g1])
points(T5[g1], BMI5[g1])
title(main="BMI change overtime for the rebounder group")
legend("topleft", cex=0.7, lty=c(1,2), col=c(1,1), legend=c("Mean BMI time course", "Expected inter-individual BMI variability"))

# # # # # #
# GROUP 2 #
# # # # # #

A2_age0 <- mean(age0[g2])
A2_BW0 <- mean(BW0[g2])
A2_FM0 <- mean(FM0[g2])
A2_LM0 <- mean(LM0[g2])
A2_H <- mean(H[g2])
A2_d0 <- ((1-Btef)*1.5-1)*(10*A2_BW0+625*A2_H-5*A2_age0-161)*4.184
A2_EI0 <- (655 + 9.56*A2_BW0 + 186*A2_H -4.68*A2_age0)*1.33*4.184
A2_EE0 <- A2_EI0
A2_K <- A2_EI0 - gf*A2_FM0 - gl*A2_LM0 - A2_d0

A2_parameters <- c(A2_EI0, A2_K, A2_H, A2_age0)
A2_init <- c(fatmass = A2_FM0,leanmass = A2_LM0)
A2_Data <- data.frame(time = c(T2[g2],T5[g2]), fatmass = c(FM2[g2],FM5[g2]) , leanmass= c(LM2[g2],LM5[g2]))
soltime = seq(0,2200)

A2_modelcost <- function(P) {
  sol <- lsoda(y=A2_init, times=soltime, func = EqBW, parms = c(A2_parameters,P[1],P[2]))
  return(modCost(sol,A2_Data))
}

A2_Fit <- modFit(f = A2_modelcost, p = c(mean(EIsurg[g2]),mean(EIfinal[g2])))

A2_EIsurg <- A2_Fit$par[1]
A2_EIfinal <- A2_Fit$par[2]
A2_EIsurgCI <- unname(c(A2_EIsurg + summary(A2_Fit)$par[1,2]*qt(0.05, 2), A2_EIsurg + summary(A2_Fit)$par[1,2]*qt(0.95, 2)))
A2_EIfinalCI <- unname(c(A2_EIfinal + summary(A2_Fit)$par[2,2]*qt(0.05, 2), A2_EIfinal + summary(A2_Fit)$par[2,2]*qt(0.95, 2)))
A2_bestfit <- lsoda(y=A2_init, times=soltime, func = EqBW, parms = c(A2_parameters,A2_EIsurg, A2_EIfinal))
A2_fitinf <- lsoda(y=A2_init, times=soltime, func = EqBW, parms = c(A2_parameters,A2_EIsurgCI[1], A2_EIfinalCI[1]))
A2_fitsup <- lsoda(y=A2_init, times=soltime, func = EqBW, parms = c(A2_parameters,A2_EIsurgCI[2], A2_EIfinalCI[2]))

plot(soltime, (A2_bestfit[,2]+A2_bestfit[,3])/A2_H^2, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50), col=2)
lines(soltime, (A2_fitinf[,2]+A2_fitinf[,3])/A2_H^2, type='l', lty=2, col=2)
lines(soltime, (A2_fitsup[,2]+A2_fitsup[,3])/A2_H^2, type='l', lty=2, col=2)
points(T2[g2], BMI2[g2], col=2)
points(T5[g2], BMI5[g2], col=2)
title(main="BMI change overtime for the first weight stable group")
legend("topleft", cex=0.7, lty=c(1,2), col=c(2,2), legend=c("Mean BMI time course", "Expected inter-individual BMI variability"))

# # # # # #
# GROUP 3 #
# # # # # #

A3_age0 <- mean(age0[g3])
A3_BW0 <- mean(BW0[g3])
A3_FM0 <- mean(FM0[g3])
A3_LM0 <- mean(LM0[g3])
A3_H <- mean(H[g3])
A3_d0 <- ((1-Btef)*1.5-1)*(10*A3_BW0+625*A3_H-5*A3_age0-161)*4.184
A3_EI0 <- (655 + 9.56*A3_BW0 + 186*A3_H -4.68*A3_age0)*1.33*4.184
A3_EE0 <- A3_EI0
A3_K <- A3_EI0 - gf*A3_FM0 - gl*A3_LM0 - A3_d0

A3_parameters <- c(A3_EI0, A3_K, A3_H, A3_age0)
A3_init <- c(fatmass = A3_FM0,leanmass = A3_LM0)
A3_Data <- data.frame(time = c(T2[g3],T5[g3]), fatmass = c(FM2[g3],FM5[g3]) , leanmass= c(LM2[g3],LM5[g3]))
soltime = seq(0,2200)

A3_modelcost <- function(P) {
  sol <- lsoda(y=A3_init, times=soltime, func = EqBW, parms = c(A3_parameters,P[1],P[2]))
  return(modCost(sol,A3_Data))
}

A3_Fit <- modFit(f = A3_modelcost, p = c(mean(EIsurg[g3]),mean(EIfinal[g3])))

A3_EIsurg <- A3_Fit$par[1]
A3_EIfinal <- A3_Fit$par[2]
A3_EIsurgCI <- unname(c(A3_EIsurg + summary(A3_Fit)$par[1,2]*qt(0.05, 2), A3_EIsurg + summary(A3_Fit)$par[1,2]*qt(0.95, 2)))
A3_EIfinalCI <- unname(c(A3_EIfinal + summary(A3_Fit)$par[2,2]*qt(0.05, 2), A3_EIfinal + summary(A3_Fit)$par[2,2]*qt(0.95, 2)))
A3_bestfit <- lsoda(y=A3_init, times=soltime, func = EqBW, parms = c(A3_parameters,A3_EIsurg, A3_EIfinal))
A3_fitinf <- lsoda(y=A3_init, times=soltime, func = EqBW, parms = c(A3_parameters,A3_EIsurgCI[1], A3_EIfinalCI[1]))
A3_fitsup <- lsoda(y=A3_init, times=soltime, func = EqBW, parms = c(A3_parameters,A3_EIsurgCI[2], A3_EIfinalCI[2]))

plot(soltime, (A3_bestfit[,2]+A3_bestfit[,3])/A3_H^2, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50), col=3)
lines(soltime, (A3_fitinf[,2]+A3_fitinf[,3])/A3_H^2, type='l', lty=2, col=3)
lines(soltime, (A3_fitsup[,2]+A3_fitsup[,3])/A3_H^2, type='l', lty=2, col=3)
points(T2[g3], BMI2[g3], col=3)
points(T5[g3], BMI5[g3], col=3)
title(main="BMI change overtime for the second weight stable group")
legend("topleft", cex=0.7, lty=c(1,2), col=c(3,3), legend=c("Mean BMI time course", "Expected inter-individual BMI variability"))