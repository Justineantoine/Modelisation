rm(list=ls())
library(deSolve)
library(FME)
library(plotrix)

dat <- read.table(file = "CLINICAL(copietravail).csv", sep= ',', header = TRUE)
######################################
# BODYWEIGHT, FAT MASS AND LEAN MASS #
######################################
BW0 <- dat$BW0
BW1 <- dat$BW1
BW2 <- dat$BW2
BW5 <- dat$BW5

FM0 <- dat$TBFD0
FM2 <- dat$TBFD2
FM5 <- dat$TBFD5

adjust = subset(FM0/dat$TBFB0, is.na(FM0/dat$TBFB0)==F)
mean_adjust = mean(adjust)

FM0[18]=mean_adjust*dat$TBFB0[18]

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
tau <- 14

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
EI0 <- (10*BW0 + 625*H -5*age0)*1.5*4.184
EIsurg <- c()
EIfinal <- c()
Tsurg <- c()

EI <- function(t, Ts, EIi, EIs, EIf)
{
  # if (t<=0){
  #   EIi
  # }
  # else if (t<Ts){
  #   EIs
  # }
  # else{
  #   EIf
  # }
  s = (EIf-EIs)/(900-Ts) #garde le regime pendant 150 jours (~5 mois) puis reprise progressive de nourriture jusqu'au 900eme jour (~3 ans)
  if (t<=0){
    EIi

  }
  else if (t<Ts){
    EIs

  }
  else{
    intake = (t-Ts)*s+EIs
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
EE <- function(t, Ts, EIi, EIs, EIf, k, h, a, FM, LM){
  B <- Btef+Bat*(1-exp(-t/tau))
  PA <- ((1-Btef)*1.5-1)*4.184*(10*(FM+LM)+625*h-5*(a+t/365)-161)
  partF <- FM/(C+FM)
  partL <- 1-partF
  result <- (k + gf*FM + gl*LM + PA - EIi*B + EI(t, Ts, EIi, EIs, EIf)*(B + partF*nf/pf + partL*nl/pl))/(1 + partF*nf/pf + partL*nl/pl)
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
  Ts <- parameters[7]
  partF <- y[1]/(C+y[1])
  partL <- 1-partF
  
  dF <- partF/pf * (EI(t, Ts, EIi, EIs, EIf) - EE(t, Ts, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dL <- partL/pl * (EI(t, Ts, EIi, EIs, EIf) - EE(t, Ts, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  list(c(dF, dL))
}

soltime <- seq(0, 2000)

# for (i in 1:41){
#   parameters <- c(EI0[i], K[i], H[i], age0[i])
#   init <- c(fatmass = FM0[i],leanmass = LM0[i])
#   Data <- data.frame(time = c(T2[i],T5[i]), fatmass = c(FM2[i],FM5[i]) , leanmass= c(LM2[i],LM5[i]))
# 
#   modelcost <- function(P) {
#     sol <- lsoda(y=init, times=soltime, func = EqBW, parms = c(parameters,P[1],P[2],P[3]))
#     return(modCost(sol,Data))
#   }
# 
#   Fit <- modMCMC(f = modelcost, p = c(7300,9300,500))
#   print(summary(Fit))
#   EIsurg[i] <- summary(Fit)$p1[1]
#   EIfinal[i] <- summary(Fit)$p2[1]
#   Tsurg[i] <- summary(Fit)$p3[1]
# 
# }

EIsurg <- c(5885,6897,7103,7533,6188,7129,6451,3100,865,4789,9226,4285,5752,1012,6531,9446,634,1836,8634,6124,6641,9300,5636,1884,6000,9445,1449,6341,6106,5268,9101,5909,5492,79,6747,6127,3676,5206,9967,8349,13209)
EIfinal <- c(12110,10049,11914,10421,10834,10184,8365,10863,12604,12080,11863,8507,10178,11202,9549,13434,11299,13191,10014,8782,9172,11087,12933,11773,11248,787,9550,9965,13093,11357,12141,10862,9122,12820,11513,9531,9794,14593,11286,15986,11593)
Ts <- c(521,556,278,731,383,563,780,287,561,412,739,436,273,247,425,381,437,308,1043,895,1053,876,400,384,282,840,312,394,184,443,780,402,571,247,364,660,349,359,364,695,991)

####################
# AVERAGE PATIENTS #
####################

dBMI <- BMI2 - BMI5
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

gr <- c(g1,g2)
gws <- g3

st_error<- function(x){
  result <- sd(x)/sqrt(length(x))
  result
}

Tg <- c(0,1,2,5)
BMIgr <- c(mean(BMI0[gr]), mean(BMI1[gr]), mean(BMI2[gr]), mean(BMI5[gr]))
errgr <- c(st_error(BMI0[gr]), st_error(BMI1[gr]), st_error(BMI2[gr]), st_error(BMI5[gr]))
BMIgws <- c(mean(BMI0[gws]), mean(BMI1[gws]), mean(BMI2[gws]), mean(BMI5[gws]))
errgws <- c(st_error(BMI0[gws]), st_error(BMI1[gws]), st_error(BMI2[gws]), st_error(BMI5[gws]))

plot(Tg, BMIgr, type = "l", ylim=c(25, 50), xlab="Year of investigation", ylab="BMI (kg/m²)")
plotCI(Tg, BMIgr, uiw = errgr, lwd =2, col = 1, add =T)
lines(Tg, BMIg2, col=2)
plotCI(Tg, BMIg2, uiw = errg2, lwd =2, col = 2, add =T)
lines(Tg, BMIgws, col=3)
plotCI(Tg, BMIgws, uiw = errgws, lwd =2, col = 3, add =T)


      # # # # # #
      # GROUP 1 #
      # # # # # #

A1_age0 <- mean(age0[gr])
A1_BW0 <- mean(BW0[gr])
A1_FM0 <- mean(FM0[gr])
A1_LM0 <- mean(LM0[gr])
A1_H <- mean(H[gr])
A1_d0 <- ((1-Btef)*1.5-1)*(10*A1_BW0+625*A1_H-5*A1_age0-161)*4.184
A1_EI0 <- (10*A1_BW0 + 625*A1_H -5*A1_age0-161)*1.5*4.184
A1_EE0 <- A1_EI0
A1_K <- A1_EI0 - gf*A1_FM0 - gl*A1_LM0 - A1_d0

A1_parameters <- c(A1_EI0, A1_K, A1_H, A1_age0)
A1_init <- c(fatmass = A1_FM0,leanmass = A1_LM0)
A1_Data <- data.frame(time = c(T2[gr],T5[gr]), fatmass = c(FM2[gr],FM5[gr]) , leanmass= c(LM2[gr],LM5[gr]))


A1_modelcost <- function(P) {
  sol <- lsoda(y=A1_init, times=soltime, func = EqBW, parms = c(A1_parameters,P[1],P[2],P[3]))
  return(modCost(sol,A1_Data))
}

P1 <- modFit(f = A1_modelcost, p =c(mean(EIsurg[gr]),mean(EIfinal[gr]), mean(Ts[gr])), control=list(epsfcn=0.1, factor=0.001))
Covar <- summary(P1)$cov.scaled * 0.5^2/2

A1_Fit <- modMCMC(f = A1_modelcost, p = P1$par, jump=Covar, lower = c(0,0,0))

A1_EIsurg <- summary(A1_Fit)$p1[1]
A1_EIsurgsd <- summary(A1_Fit)$p1[2]
A1_EIfinal <- summary(A1_Fit)$p2[1]
A1_EIfinalsd <- summary(A1_Fit)$p2[2]
A1_Ts <- summary(A1_Fit)$p3[1]
A1_Tssd <- summary(A1_Fit)$p3[2]
A1_EIsurgCI <- unname(c(A1_EIsurg + A1_EIsurgsd*qt(0.025, 2), A1_EIsurg +A1_EIsurgsd*qt(0.975, 2)))
A1_EIfinalCI <- unname(c(A1_EIfinal + A1_EIfinalsd*qt(0.025, 2), A1_EIfinal + A1_EIfinalsd*qt(0.975, 2)))
A1_TsCI <- unname(c(A1_Ts + A1_Tssd*qt(0.025, 2), A1_Ts + A1_Tssd*qt(0.975, 2)))
A1_bestfit <- lsoda(y=A1_init, times=soltime, func = EqBW, parms = c(A1_parameters,A1_EIsurg, A1_EIfinal,A1_Ts))
A1_fitinf <- lsoda(y=A1_init, times=soltime, func = EqBW, parms = c(A1_parameters,A1_EIsurgCI[1], A1_EIfinalCI[1], A1_TsCI[2]))
A1_fitsup <- lsoda(y=A1_init, times=soltime, func = EqBW, parms = c(A1_parameters,A1_EIsurgCI[2], A1_EIfinalCI[2], A1_TsCI[1]))

plot(soltime, (A1_bestfit[,2]+A1_bestfit[,3])/A1_H^2, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(0,50))
lines(soltime, (A1_fitinf[,2]+A1_fitinf[,3])/A1_H^2, type='l', lty=2)
lines(soltime, (A1_fitsup[,2]+A1_fitsup[,3])/A1_H^2, type='l', lty=2)
points(T2[gr], BMI2[gr])
points(T5[gr], BMI5[gr])
title(main="BMI change overtime for the rebounder group")
legend("topleft", cex=0.7, lty=c(1,2), col=c(1,1), legend=c("Mean BMI time course", "Expected inter-individual BMI variability"))


      # # # # # #
      # GROUP 3 #
      # # # # # #

A3_age0 <- mean(age0[gws])
A3_BW0 <- mean(BW0[gws])
A3_FM0 <- mean(FM0[gws])
A3_LM0 <- mean(LM0[gws])
A3_H <- mean(H[gws])
A3_d0 <- ((1-Btef)*1.5-1)*(10*A3_BW0+625*A3_H-5*A3_age0-161)*4.184
A3_EI0 <- (655 + 9.56*A3_BW0 + 186*A3_H -4.68*A3_age0)*1.33*4.184
A3_EE0 <- A3_EI0
A3_K <- A3_EI0 - gf*A3_FM0 - gl*A3_LM0 - A3_d0

A3_parameters <- c(A3_EI0, A3_K, A3_H, A3_age0)
A3_init <- c(fatmass = A3_FM0,leanmass = A3_LM0)
A3_Data <- data.frame(time = c(T2[gws],T5[gws]), fatmass = c(FM2[gws],FM5[gws]) , leanmass= c(LM2[gws],LM5[gws]))
soltime = seq(0,2200)

A3_modelcost <- function(P) {
  sol <- lsoda(y=A3_init, times=soltime, func = EqBW, parms = c(A3_parameters,P[1],P[2],P[3]))
  return(modCost(sol,A3_Data))
}

A3_Fit <- modMCMC(f = A3_modelcost, p = c(mean(EIsurg[gws]),mean(EIfinal[gws]), mean(Ts[gws])))

A3_EIsurg <- summary(A3_Fit)$p1[1]
A3_EIsurgsd <- summary(A3_Fit)$p1[2]
A3_EIfinal <- summary(A3_Fit)$p2[1]
A3_EIfinalsd <- summary(A3_Fit)$p2[2]
A3_Ts <- summary(A3_Fit)$p3[1]
A3_Tssd <- summary(A3_Fit)$p3[2]
A3_EIsurgCI <- unname(c(A3_EIsurg + A3_EIsurgsd*qt(0.025, 2), A3_EIsurg +A3_EIsurgsd*qt(0.975, 2)))
A3_EIfinalCI <- unname(c(A3_EIfinal + A3_EIfinalsd*qt(0.025, 2), A3_EIfinal + A3_EIfinalsd*qt(0.975, 2)))
A3_TsCI <- unname(c(A3_Ts + A3_Tssd*qt(0.025, 2), A3_Ts + A3_Tssd*qt(0.975, 2)))
A3_bestfit <- lsoda(y=A3_init, times=soltime, func = EqBW, parms = c(A3_parameters,A3_EIsurg, A3_EIfinal,A3_Ts))
A3_fitinf <- lsoda(y=A3_init, times=soltime, func = EqBW, parms = c(A3_parameters,A3_EIsurgCI[1], A3_EIfinalCI[1], A3_TsCI[2]))
A3_fitsup <- lsoda(y=A3_init, times=soltime, func = EqBW, parms = c(A3_parameters,A3_EIsurgCI[2], A3_EIfinalCI[2], A3_TsCI[1]))

plot(soltime, (A3_bestfit[,2]+A3_bestfit[,3])/A3_H^2, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(0,50), col=3)
lines(soltime, (A3_fitinf[,2]+A3_fitinf[,3])/A3_H^2, type='l', lty=2, col=3)
lines(soltime, (A3_fitsup[,2]+A3_fitsup[,3])/A3_H^2, type='l', lty=2, col=3)
points(T2[gws], BMI2[gws], col=3)
points(T5[gws], BMI5[gws], col=3)
title(main="BMI change overtime for the second weight stable group")
legend("topleft", cex=0.7, lty=c(1,2), col=c(3,3), legend=c("Mean BMI time course", "Expected inter-individual BMI variability"))
