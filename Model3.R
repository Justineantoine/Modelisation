rm(list=ls())
library(deSolve)
library(FME)
library(fields)

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
Ts <- c()

EI <- function(t, Ts, EI0, EIsurg, EIfinal) #i pour l'individu auquel on s'interesse
{
  s = (EIfinal-EIsurg)/(900-Ts) #garde le regime pendant 150 jours (~5 mois) puis reprise progressive de nourriture jusqu'au 900eme jour (~3 ans)
  if (t<=0){
    EI0
  }
  else if (t<=Ts){
    EIsurg
  }
  else{
    EIfinal
  }
}

##############################
# INITIAL ENERGY EXPENDITURE #
##############################
EE0 <- EI0

##############
# K CONSTANT #
##############
K <- EI0 - gf*FM0 - gl*LM0 - d0

##############
# EDO SYSTEM #
##############

EqBW <- function(t, y, parameters){
  B <- Btef+Bat
  PA <- ((1-Btef)*1.5-1)*4.184*(10*(y[1]+y[2])+625*parameters[3]-5*(parameters[4]+t/365)-161)
  partF <- y[1]/(C+y[1])
  partL <- 1-partF
  P1 <- parameters[5]
  P2 <- parameters[6]
  P3 <- parameters[7]
  EE <- (parameters[2] + gf*y[1] + gl*y[2] + PA + parameters[1]*B + EI(t, P3, parameters[1], P1, P2)*(-B + partF*nf/pf + partL*nl/pl))/(1 + partF*nf/pf + partL*nl/pl)
  
  dF <- partF/pf * (EI(t, P3, parameters[1], P1, P2) - EE)
  dL <- partL/pl * (EI(t, P3, parameters[1], P1, P2) - EE)
  list(c(dF, dL))
}

soltime <- seq(0, 2200)

for (i in 1:41){
  parameters <- c(EI0[i], K[i], H[i], age0[i])
  init <- c(fatmass = FM0[i],leanmass = LM0[i])
  Data <- data.frame(time = c(T2[i],T5[i]), fatmass = c(FM2[i],FM5[i]) , leanmass= c(LM2[i],LM5[i]))
  
  modelcost <- function(P) {
    sol <- lsoda(y=init, times=soltime, func = EqBW, parms = c(parameters,P[1],P[2],P[3]))
    return(modCost(sol,Data))
  }
  
  Fit <- modFit(f = modelcost, p = c(7300,9300,500))
  
  EIsurg[i] <- Fit$par[1]
  EIfinal[i] <- Fit$par[2]
  Ts[i] <- Fit$par[3]
  
}

######################################
# PLOT FM AND BODY WEIGHT, EE AND EI #
######################################
for (k in 1:41){
  ind <- k
  graphtime <- seq(-100,2200)
  parameters <- c(EI0[ind], K[ind], H[ind], age0[ind], EIsurg[ind], EIfinal[ind], Ts[ind])
  init <- c(fatmass = FM0[ind],leanmass = LM0[ind])
  bestfit <- lsoda(y=init, times=soltime, func = EqBW, parms = c(parameters))
  
  graphEI= rep(EI0[ind], 100)
  graphEE = rep(EE0[ind], 100)
  S = c()
  for (i in 101:2301){
    graphEI[i] <- EI(graphtime[i], Ts[ind], EI0[ind], EIsurg[ind], EIfinal[ind])
    graphEE[i] <- (K[ind] + gf*bestfit[i-100,2] + gl*bestfit[i-100,3] + ((1-Btef)*1.5-1)*4.184*(10*(bestfit[i-100,2]+bestfit[i-100,3])+625*H[ind]-5*(age0[ind]+(i-100)/365)-161) + EI0[ind]*(Bat+Btef) + EI(graphtime[i], Ts[ind], EI0[ind], EIsurg[ind], EIfinal[ind])*(-Btef -Bat + bestfit[i-100,2]/(C+bestfit[i-100,2])*nf/pf + C/(C+bestfit[i-100,2])*nl/pl))/(1 + bestfit[i-100,2]/(C+bestfit[i-100,2])*nf/pf + C/(C+bestfit[i-100,2])*nl/pl)
    if (graphEI[i] < 1.001*graphEE[i] & graphEI[i] >0.999*graphEE[i] & i >101) {S = c(S, i-100)}
  }
  
  plot(graphtime, graphEI, type="l", xlab="Days", ylab="Energy rate ", col=1, ylim=c(4000, 15000))
  lines(graphtime, graphEE, type ="l", lty = 1, col=2)
  abline(v=S[1], lty=4, col="brown4")
  legend("bottomright",lty=c(1,1,3), cex=0.7, col=c(1,2,"brown4"), legend=c("Energy Intake rate", "Energy Expenditure rate", "Stable state"))
  
  
  graphFM <- c(rep(FM0[ind], 100), bestfit[,2])
  graphLM <- c(rep(LM0[ind], 100), bestfit[,3])
  graphBW <- c(rep(BW0[ind], 100), bestfit[,2]+bestfit[,3])
  
  plot(graphtime, graphFM, type="l", ylim=c(0,130), xlab="Days", ylab="Weight in kg")
  #lines(graphtime, graphBW, type = "l", lty =1, col=4)
  #lines(graphtime, graphLM, type = "l", lty =1, col="dodgerblue4")
  points(T0[ind], FM0[ind], pch=19)
  points(T2[ind], FM2[ind], pch=19)
  points(T5[ind], FM5[ind], pch=19)
  # points(T0[ind], BW0[ind], pch=19, col=4)
  # points(365, BW1[ind], pch=19, col=4)
  # points(T2[ind], BW2[ind], pch=19, col=4)
  # points(T5[ind], BW5[ind], pch=19, col=4)
  # points(T0[ind], LM0[ind], pch=3, col="dodgerblue4")
  # points(T2[ind], LM2[ind], pch=3, col="dodgerblue4")
  # points(T5[ind], LM5[ind], pch=3, col="dodgerblue4")
  abline(v=S[1], lty=4, col="brown4")
  title(main= k)
  legend("bottomright", lty=c(1,1,1,4), legend=c("BW", "FM","LM", "Stable state"), col=c(4,1,"dodgerblue4","brown4"), cex=0.7)
  
}
