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

FM0[17] <- mean_adjust*dat$TBFB0[17]

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
BMI1[20] <- BMI2[20] - mean_adjust2
BMI1[19] <- BMI2[19] - mean_adjust2

########
# TIME #
########
T0 <- rep(0,39)
T2 <- dat$T2
T1 <- T2/2
T5 <- dat$T5

##########################
# LIPID AGE AND TURNOVER #
##########################
LA0 <- dat$LA0*365
LA2 <- dat$LA2*365
LA5 <- dat$LA5*365

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

EI <- function(t, EIi, EIs, EIf)
{
  #' Return the Energy Intake course of a patient. Adherence to the diet during 500 days, then a linear intake recovery until stabilisation at 900 days
  s = (EIf-EIs)/400
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

EI2 <- function(t, EIi, EIs, EIf)
{
  #' Return the Energy Intake course of a patient. Adherence to the diet during 500 days and then a stable intake recovery.
  if (t<=0){
    EIi
  }
  else if (t<=500){
    EIs
  }
  else{
    EIf
  }
}

##############################
# INITIAL ENERGY EXPENDITURE #
##############################
EE0 <- EI0

EE <- function(t, EIi, EIs, EIf, k, h, a, FM, LM){
  B <- Btef+Bat*(1-exp(-t/tau))
  PA <- ((1-Btef)*1.5-1)*4.184*(10*(FM+LM)+625*h-5*(a+t/365)-161)
  partF <- FM/(C+FM)
  partL <- 1-partF
  result <- (k + gf*FM + gl*LM + PA - EIi*B + EI(t, EIi, EIs, EIf)*(B + partF*nf/pf + partL*nl/pl))/(1 + partF*nf/pf + partL*nl/pl)
  result
}

EE2 <- function(t, EIi, EIs, EIf, k, h, a, FM, LM){
  B <- Btef+Bat*(1-exp(-t/tau))
  PA <- ((1-Btef)*1.5-1)*4.184*(10*(FM+LM)+625*h-5*(a+t/365)-161)
  partF <- FM/(C+FM)
  partL <- 1-partF
  result <- (k + gf*FM + gl*LM + PA - EIi*B + EI2(t, EIi, EIs, EIf)*(B + partF*nf/pf + partL*nl/pl))/(1 + partF*nf/pf + partL*nl/pl)
  result
}

##############################
# K CONSTANT AND K2 CONSTANT #
##############################
K <- EI0 - gf*FM0 - gl*LM0 - d0
K2 <- 5e-6

################################
# EDO SYSTEM : BW, LA AND Kout #
################################

  # # # #
  # FIT #
  # # # #
soltime <- seq(0,2200)

EqBW <- function(t, y, parameters){

  EIi <- parameters[1]
  k <- parameters[2]
  h <- parameters[3]
  a <- parameters[4]
  EIs <- parameters[5]
  EIf <- parameters[6]
  partF <- y[1]/(C+y[1])
  partL <- 1-partF

  dF <- partF/pf * (EI(t, EIi, EIs, EIf) - EE(t, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dL <- partL/pl * (EI(t, EIi, EIs, EIf) - EE(t, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  list(c(dF, dL)) 
}

EqBW2 <- function(t, y, parameters){
  
  EIi <- parameters[1]
  k <- parameters[2]
  h <- parameters[3]
  a <- parameters[4]
  EIs <- parameters[5]
  EIf <- parameters[6]
  partF <- y[1]/(C+y[1])
  partL <- 1-partF
  
  dF <- partF/pf * (EI2(t, EIi, EIs, EIf) - EE2(t, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dL <- partL/pl * (EI2(t, EIi, EIs, EIf) - EE2(t, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  list(c(dF, dL)) 
}

for (i in 1:39){
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
  
  # # # # # # #
  # LIPID AGE #
  # # # # # # #
EqBWLA <- function(t, y, parameters){
  
  EIi <- parameters[1]
  k <- parameters[2]
  h <- parameters[3]
  a <- parameters[4]
  EIs <- parameters[5]
  EIf <- parameters[6]
  partF <- y[1]/(C+y[1])
  partL <- 1-partF
  
  kin <- K2 * EI(t, EIi, EIs, EIf)
  
  dF <- partF/pf * (EI(t, EIi, EIs, EIf) - EE(t, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dL <- partL/pl * (EI(t, EIi, EIs, EIf) - EE(t, EIi, EIs, EIf, k, h, a, y[1], y[2]))
  dA <- 1 - kin*y[3]/y[1]

  list(c(dF, dL, dA)) 
}

  # # # # #
  # K OUT #
  # # # # #

EqKout <- function(parameters,fit){
  Kout <- c()
  EIi <- parameters[1]
  k <- parameters[2]
  h <- parameters[3]
  a <- parameters[4]
  EIs <- parameters[5]
  EIf <- parameters[6]
  A0 <- parameters[7]
  
  for (i in soltime){
    partF <- fit[i+1,2]/(C+fit[i+1,2])
    partL <- 1-partF
    
    kin <- K2 * EI(i, EIi, EIs, EIf)
    
    # if (EI(i, EIi, EIs, EIf) - EE(i, EIi, EIs, EIf, k, h, a, fit[i+1,2], fit[i+1,3])< 0){kin <- 0}
    # else{kin <- EI(i, EIi, EIs, EIf) - EE(i, EIi, EIs, EIf, k, h, a, fit[i+1,2], fit[i+1,3])}
    
    # kin <- fit[i+1,2]/A0
    
    dF <- partF/pf * (EI(i, EIi, EIs, EIf) - EE(i, EIi, EIs, EIf, k, h, a, fit[i+1,2], fit[i+1,3]))
    Kout <- c(Kout, unname((kin - dF)/fit[i+1,2]))
  }
  
  Kout
}

##################
# PATIENT GROUPS #
##################
dBMI <- BMI0 - BMI5
tertile<- unname(quantile(dBMI, probs = c(0.33, 0.67)))

gr <- c()
gs <- c()

for (i in 1:39){
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

groupsBMI <- function(){
  #' Plot the BMI time course of each group of patients
  plot(Tg, BMIgr, type = "l", ylim=c(25, 47), xlab="Year of investigation", ylab="BMI (kg/m²)")
  plotCI(Tg, BMIgr, uiw = errgr, lwd =2, col = 1, add =T)
  lines(Tg, BMIgs, col=2)
  plotCI(Tg, BMIgs, uiw = errgs, lwd =2, col = 2, add =T)
  legend("topright", cex=0.7, lty=c(1,1), col=c(1,2), legend=c("Tertile 1 : rebounders", "Tertile 2-3 : weight stable"))
}


##########################################################
# PLOT FM AND BODY WEIGHT, EE AND EI, LIPID AGE AND KOUT #
##########################################################
graphtime <- seq(-100,2200)
individuals <- c(gr,gs)

  # # # # # # # # #
  # BW, LM AND FM #
  # # # # # # # # #
indBW <- function(k, bestfit)
{
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
}

  # # # # # # #
  # EI AND EE #
  # # # # # # #
indE <- function(k, bestfit)
{
  graphEI <- rep(EI0[k], 100)
  graphEE <- rep(EE0[k], 100)
  for (i in 101:2301){
    graphEI[i] <- EI(graphtime[i], EI0[k], EIsurg[k], EIfinal[k])
    graphEE[i] <- EE(graphtime[i], EI0[k], EIsurg[k], EIfinal[k], K[k], H[k], age0[k], bestfit[i-100,2], bestfit[i-100,3])
  }
  
  plot(graphtime, graphEI, type="l", xlab="Days", ylab="Energy rate ", col=1, ylim=c(1000, 15000))
  lines(graphtime, graphEE, type ="l", lty = 1, col=2)
  title(main=c("Energy rates for patient : ", k))
  legend("bottomright",lty=c(1,1), cex=0.7, col=c(1,2), legend=c("Energy Intake rate", "Energy Expenditure rate"))
}
  
  # # # # # # #
  # LIPID AGE #
  # # # # # # #
indLA <- function(k, bestfit) 
{
  graphA <- c(rep(LA0[k], 100), bestfit[,4])
  
  plot(graphtime, graphA, type="l", xlab="Days", ylab="Lipid Age", ylim=c(0, 2700))
  title(main=c("Lipid age for patient : ", k))
}

  # # # # #
  # K OUT #
  # # # # #
indKout <- function(k, parms, bestfit) 
{
  table <- EqKout(parms, bestfit)
  graphKout <- c(rep(table[1], 100), table)
  
  plot(graphtime, graphKout, type ="l", xlab="Days", ylab="Kout (/d)", ylim=c(0, 0.01))
  title(main=c("Turnover for patient : ", k))

} 

  # # # # # # # # #  #
  # INDIVIDUAL PLOTS #
  # # # # # # # # #  #
indPlots <- function(BW=T, E=F, LA=F, Kout=F)
{
  for (k in c(gr,gs)){
    if (k <= length(gr)){parameters <- c(EI0[k], K[k], H[k], age0[k], EIsurg[k], EIfinal[k])}
    else{parameters <- c(EI0[k], K[k], H[k], age0[k], EIsurg[k], EIfinal[k])}
    init <- c(fatmass = FM0[k],leanmass = LM0[k], age =LA0[k])
    bestfit <- lsoda(y=init, times=soltime, func = EqBWLA, parms = c(parameters)) 
    
    if(BW){
      indBW(k, bestfit)
    }
    
    if(E){
      indE(k, bestfit)
    }
    
    if(LA){
      indLA(k, bestfit)
    }
    
    if(Kout){
      indKout(k, parameters, bestfit)
    }
  }
}

  # # # # # #
  # DATASET #
  # # # # # #

indTables <- function()
{
  #' Creates
  LAdata <- c()
  Koutdata <- c()
  allKout <- c()
  allLA <- c()
  for (i in 1:39){
    k <- individuals[i]
    parameters <- c(EI0[k], K[k], H[k], age0[k], EIsurg[k], EIfinal[k])
    init <- c(fatmass = FM0[k],leanmass = LM0[k], age = LA0[k])
    bestfit <- lsoda(y=init, times=soltime, func = EqBWLA, parms = c(parameters))
    
    allLA <- cbind(allLA, unname(bestfit[,4]))
    LAdata <- rbind(LAdata, c(patient = k, T0 = allLA[1,i], T2 = allLA[T2[i],i], T5 = allLA[T5[i],i]))
    allKout <- cbind(allKout, EqKout(parameters, bestfit))
    Koutdata <- rbind(Koutdata, c(patient = k, T0 = allKout[1,i], T2 = allKout[T2[i],i], T5 = allKout[T5[i],i]))
  }
  
  colnames(allLA) <- individuals
  colnames(allKout) <- individuals
  list(LAdata, allLA, Koutdata, allKout)
}

tables <- indTables()
LAdata <- tables[[1]]
allLA <- tables[[2]]
Koutdata <- tables[[3]]
allKout <- tables[[4]]

  # # # # # # # # # # # # #
  # INTERINDIVIDUAL PLOTS #
  # # # # # # # # # # # # #

interPlots <- function(){
  #' Plots the Kout and Lipid Age of each patient.
  plot(soltime, allKout[,1], type ="l", xlab="Days", ylab="Kout (/d)", ylim=c(0.0002, 0.004))
  for (i in 1:39){
    if (i > length(gr)){
        lines(soltime, allKout[,i], col=2)
    }
    else{
      lines(soltime, allKout[,i])
    }
    
  }
  title(main="Kout for all patients")
  legend("topright", legend=c("Rebounders", "Weight stable"), lty=c(1,1), col=c(1,2))
  
  plot(soltime, allLA[,1], type ="l", xlab="Days", ylab="Lipid Age (d)", ylim=c(0, 2700))
  for (i in 1:39){
    if (i > length(gr)){
      lines(soltime, allLA[,i], col=2)
    }
    else{
      lines(soltime, allLA[,i])
    }
    
  }
  title(main="Lipid Age for all patients")
  legend("topright", legend=c("Rebounders", "Weight stable"), lty=c(1,1), col=c(1,2))
}

    
##################################
# INDIVIDUAL INFORMATION STORAGE #
##################################
write.csv(cbind(EI0, EIsurg, EIfinal, EI0-EIfinal), file = "Intakes.csv")
write.csv(cbind(LAdata, Delta = LAdata[,1]-LAdata[,3]), file = "LipidAge.csv")
write.csv(cbind(Koutdata, Delta = Koutdata[,1]-Koutdata[,3]), file = "Kout.csv")
write.csv(cbind(K2*EI0, K2*EIsurg, K2*EIfinal, K2*(EI0-EIfinal)), file = "Kin.csv")
####################
# AVERAGE PATIENTS #
####################

                  # # # # # # #
                  # REBOUNDER #
                  # # # # # # #

#PARAMETERS AND DATA
R_n <- length(gr)
R_age0 <- mean(age0[gr])
R_BW0 <- mean(BW0[gr])
R_FM0 <- mean(FM0[gr])
R_FM0CI <- c(R_FM0 + qnorm(0.05)*sd(FM0[gr]/sqrt(length(gr))), R_FM0 + qnorm(0.95)*sd(FM0[gr]/sqrt(length(gr))))
R_LM0 <- mean(LM0[gr])
R_LM0CI <- c(R_LM0 + qnorm(0.05)*sd(LM0[gr]/sqrt(length(gr))), R_LM0 + qnorm(0.95)*sd(LM0[gr]/sqrt(length(gr))))
R_LA0 <- mean(LA0[gr])
R_LA0CI <- c(R_LA0 + qnorm(0.05)*sd(LA0[gr]/sqrt(length(gr))), R_LA0 + qnorm(0.95)*sd(LA0[gr]/sqrt(length(gr))))
R_H <- mean(H[gr])
R_d0 <- ((1-Btef)*1.5-1)*(10*R_BW0+625*R_H-5*R_age0-161)*4.184
R_EI0 <- (10*R_BW0 + 625*R_H -5*R_age0-161)*1.5*4.184
R_EE0 <- R_EI0
R_K <- R_EI0 - gf*R_FM0 - gl*R_LM0 - R_d0
R_parameters <- c(R_EI0, R_K, R_H, R_age0)
R_init <- c(fatmass = R_FM0,leanmass = R_LM0)
R_Data <- data.frame(time = c(T2[gr],T5[gr]), fatmass = c(FM2[gr],FM5[gr]) , leanmass= c(LM2[gr],LM5[gr]))


#FIT AND STATISTICS
R_modelcost <- function(P) {
  sol <- lsoda(y=R_init, times=soltime, func = EqBW, parms = c(R_parameters,P[1],P[2]))
  return(modCost(sol,R_Data))
}
R_Fit <- modFit(f = R_modelcost, p = c(mean(EIsurg[gr]),mean(EIfinal[gr])))
R_EIsurg <- R_Fit$par[1]
R_EIfinal <- R_Fit$par[2]

    # error sum of squares #
R_sssurg <- sum((EIsurg-mean(EIsurg[gr]))^2)
R_ssfinal <- sum((EIfinal-mean(EIfinal[gr]))^2)

    # mse #
R_mse <- mean(R_Fit$residuals^2)

    # standard error #
R_sesurg <- sqrt(R_mse * (1/R_n + 1 + (7000-mean(EIsurg[gr])^2/R_sssurg)))
R_sefinal <- sqrt(R_mse * (1/R_n + 1 + (7000-mean(EIfinal[gr])^2/R_ssfinal)))


#BESTFIT AND PREDICTION INTERVAL

    # initialisation #
R_init <- c(fatmass = R_FM0,leanmass = R_LM0, age = R_LA0)
R_initinf <- c(fatmass = R_FM0CI[1],leanmass = R_LM0CI[1], age = R_LA0CI[1])
R_initsup <- c(fatmass = R_FM0CI[2],leanmass = R_LM0CI[2], age = R_LA0CI[2]) #R_initinf and sup takes into consderation the intial variability between individuals

    # prediction intervals #
R_EIsurgCI <- unname(c(R_EIsurg + R_sesurg*qt(0.05, R_n-2), R_EIsurg + R_sesurg*qt(0.95, R_n-2)))
R_EIfinalCI <- unname(c(R_EIfinal + R_sefinal*qt(0.05, R_n-2), R_EIfinal + R_sefinal*qt(0.95, R_n-2)))
R_bestfit <- lsoda(y=R_init, times=soltime, func = EqBWLA, parms = c(R_parameters,R_EIsurg, R_EIfinal))
R_fitinf <- lsoda(y=R_init, times=soltime, func = EqBWLA, parms = c(R_parameters,R_EIsurgCI[1], R_EIfinalCI[1]))
R_fitsup <- lsoda(y=R_init, times=soltime, func = EqBWLA, parms = c(R_parameters,R_EIsurgCI[2], R_EIfinalCI[2]))
R_fitinf2 <- lsoda(y=R_initinf, times=soltime, func = EqBWLA, parms = c(R_parameters,R_EIsurgCI[1], R_EIfinalCI[1]))
R_fitsup2 <- lsoda(y=R_initsup, times=soltime, func = EqBWLA, parms = c(R_parameters,R_EIsurgCI[2], R_EIfinalCI[2]))

    # bmi #
R_BMI <- (R_bestfit[,2]+R_bestfit[,3])/R_H^2
R_BMIinf <- (R_fitinf[,2]+R_fitinf[,3])/R_H^2
R_BMIsup <- (R_fitsup[,2]+R_fitsup[,3])/R_H^2
R_BMIinf2 <- (R_fitinf2[,2]+R_fitinf2[,3])/R_H^2
R_BMIsup2 <- (R_fitsup2[,2]+R_fitsup2[,3])/R_H^2

    # EI and EE #
R_EI <- c() ; R_EIinf <- c() ; R_EIsup <- c()
R_EE <- c() ; R_EEinf <- c() ; R_EEsup <- c()
for (i in soltime){
  R_EI <- c(R_EI, EI(i, R_EI0, R_EIsurg, R_EIfinal))
  R_EIinf <- c(R_EIinf, EI(i, R_EI0, R_EIsurgCI[1], R_EIfinalCI[1]))
  R_EIsup <- c(R_EIsup, EI(i, R_EI0, R_EIsurgCI[2], R_EIfinalCI[2]))
  
  R_EE <- c(R_EE, EE(i, R_EI0, R_EIsurg, R_EIfinal, R_K, R_H, R_age0, R_bestfit[i+1,2], R_bestfit[i+1,3]))
  R_EEinf <- c(R_EEinf, unname(EE(i, R_EI0, R_EIsurgCI[1], R_EIfinalCI[1], R_K, R_H, R_age0, R_fitinf[i+1,2], R_fitinf[i+1,3])))
  R_EEsup <- c(R_EEsup, unname(EE(i, R_EI0, R_EIsurgCI[2], R_EIfinalCI[2], R_K, R_H, R_age0, R_fitsup[i+1,2], R_fitsup[i+1,3])))
}

    # kout #
R_Kout <- EqKout(c(R_parameters, R_EIsurg, R_EIfinal, R_LA0), R_bestfit)
R_Koutinf <- EqKout(c(R_parameters, R_EIsurgCI[1], R_EIfinalCI[1], R_LA0CI[1]), R_fitinf2)
R_Koutsup <- EqKout(c(R_parameters, R_EIsurgCI[2], R_EIfinalCI[2], R_LA0CI[2]), R_fitsup2)


                  # # # # # # # # #
                  # WEIGHT STABLE #
                  # # # # # # # # #

#PARAMETERS AND DATA
S_n <- length(gs)
S_age0 <- mean(age0[gs])
S_BW0 <- mean(BW0[gs])
S_FM0 <- mean(FM0[gs])
S_FM0CI <- c(S_FM0 + qnorm(0.05)*sd(FM0[gs]/sqrt(length(gs))), S_FM0 + qnorm(0.95)*sd(FM0[gs]/sqrt(length(gs))))
S_LM0 <- mean(LM0[gs])
S_LM0CI <- c(S_LM0 + qnorm(0.05)*sd(LM0[gs]/sqrt(length(gs))), S_LM0 + qnorm(0.95)*sd(LM0[gs]/sqrt(length(gs))))
S_LA0 <- mean(LA0[gs])
S_LA0CI <- c(S_LA0 + qnorm(0.05)*sd(LA0[gs]/sqrt(length(gs))), S_LA0 + qnorm(0.95)*sd(LA0[gs]/sqrt(length(gs))))
S_H <- mean(H[gs])
S_d0 <- ((1-Btef)*1.5-1)*(10*S_BW0+625*S_H-5*S_age0-161)*4.184
S_EI0 <- (10*S_BW0 + 625*S_H -5*S_age0-161)*1.5*4.184
S_EE0 <- S_EI0
S_K <- S_EI0 - gf*S_FM0 - gl*S_LM0 - S_d0
S_parameters <- c(S_EI0, S_K, S_H, S_age0)
S_init <- c(fatmass = S_FM0,leanmass = S_LM0)
S_Data <- data.frame(time = c(T2[gs],T5[gs]), fatmass = c(FM2[gs],FM5[gs]) , leanmass= c(LM2[gs],LM5[gs]))

#FIT AND STATISTICS
S_modelcost <- function(P) {
  sol <- lsoda(y=S_init, times=soltime, func = EqBW, parms = c(S_parameters,P[1],P[2]))
  return(modCost(sol,S_Data))
}
S_Fit <- modFit(f = S_modelcost, p = c(mean(EIsurg[gs]),mean(EIfinal[gs])))
S_EIsurg <- S_Fit$par[1]
S_EIfinal <- S_Fit$par[2]

    # error sum of squares #
S_sssurg <- sum((EIsurg-mean(EIsurg[gs]))^2)
S_ssfinal <- sum((EIfinal-mean(EIfinal[gs]))^2)

    # mse #
S_mse <- mean(S_Fit$residuals^2)

    # standard error #
S_sesurg <- sqrt(S_mse * (1/S_n + 1 + (7000-mean(EIsurg[gs])^2/S_sssurg)))
S_sefinal <- sqrt(S_mse * (1/S_n + 1 + (7000-mean(EIfinal[gs])^2/S_ssfinal)))

#BESTFIT AND PREDICTION INTERVAL

    # initialisation #
S_init <- c(fatmass = S_FM0,leanmass = S_LM0, age = S_LA0)
S_initinf <- c(fatmass = S_FM0CI[1],leanmass = S_LM0CI[1], age = S_LA0CI[1])
S_initsup <- c(fatmass = S_FM0CI[2],leanmass = S_LM0CI[2], age = S_LA0CI[2])
    
    # prediction interval #
S_EIsurgCI <- unname(c(S_EIsurg + S_sesurg*qt(0.05, S_n-2), S_EIsurg + S_sesurg*qt(0.95, S_n-2)))
S_EIfinalCI <- unname(c(S_EIfinal + S_sefinal*qt(0.05, S_n-2), S_EIfinal + S_sefinal*qt(0.95, S_n-2)))
S_bestfit <- lsoda(y=S_init, times=soltime, func = EqBWLA, parms = c(S_parameters,S_EIsurg, S_EIfinal))
S_fitinf <- lsoda(y=S_init, times=soltime, func = EqBWLA, parms = c(S_parameters,S_EIsurgCI[1], S_EIfinalCI[1]))
S_fitsup <- lsoda(y=S_init, times=soltime, func = EqBWLA, parms = c(S_parameters,S_EIsurgCI[2], S_EIfinalCI[2]))
S_fitinf2 <- lsoda(y=S_initinf, times=soltime, func = EqBWLA, parms = c(S_parameters,S_EIsurgCI[1], S_EIfinalCI[1]))
S_fitsup2 <- lsoda(y=S_initsup, times=soltime, func = EqBWLA, parms = c(S_parameters,S_EIsurgCI[2], S_EIfinalCI[2]))

    # bmi #
S_BMI <- (S_bestfit[,2]+S_bestfit[,3])/S_H^2
S_BMIinf <- (S_fitinf[,2]+S_fitinf[,3])/S_H^2
S_BMIsup <- (S_fitsup[,2]+S_fitsup[,3])/S_H^2
S_BMIinf2 <- (S_fitinf2[,2]+S_fitinf2[,3])/S_H^2
S_BMIsup2 <- (S_fitsup2[,2]+S_fitsup2[,3])/S_H^2

    # EI and EE #
S_EI <- c() ; S_EIinf <- c() ; S_EIsup <- c()
S_EE <- c() ; S_EEinf <- c() ; S_EEsup <- c()
for (i in soltime){
  S_EI <- c(S_EI, EI(i, S_EI0, S_EIsurg, S_EIfinal))
  S_EIinf <- c(S_EIinf, EI(i, S_EI0, S_EIsurgCI[1], S_EIfinalCI[1]))
  S_EIsup <- c(S_EIsup, EI(i, S_EI0, S_EIsurgCI[2], S_EIfinalCI[2]))
  
  S_EE <- c(S_EE, EE(i, S_EI0, S_EIsurg, S_EIfinal, S_K, S_H, S_age0, S_bestfit[i+1,2], S_bestfit[i+1,3]))
  S_EEinf <- c(S_EEinf, unname(EE(i, S_EI0, S_EIsurgCI[1], S_EIfinalCI[1], S_K, S_H, S_age0, S_fitinf[i+1,2], S_fitinf[i+1,3])))
  S_EEsup <- c(S_EEsup, unname(EE(i, S_EI0, S_EIsurgCI[2], S_EIfinalCI[2], S_K, S_H, S_age0, S_fitsup[i+1,2], S_fitsup[i+1,3])))
}

    # kout #
S_Kout <- EqKout(c(S_parameters, S_EIsurg, S_EIfinal, S_LA0), S_bestfit)
S_Koutinf <- EqKout(c(S_parameters, S_EIsurgCI[1], S_EIfinalCI[1], S_LA0CI[1]), S_fitinf2)
S_Koutsup <- EqKout(c(S_parameters, S_EIsurgCI[2], S_EIfinalCI[2], S_LA0CI[2]), S_fitsup2)





#BMI GRAPHS

avBMI <- function(initialvar = T, superposition = T)
{
  if (initialvar){
    plot(soltime, R_BMI, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50))
    lines(soltime, R_BMIinf2, lty=2)
    lines(soltime, R_BMIsup2, lty=2)
    points(T2[gr], BMI2[gr])
    points(T5[gr], BMI5[gr])
    title(main="BMI for Weight Rebounders (WR)")
    legend("topleft", cex=0.7, lty=c(1,2), col=c(1,1), legend=c("Average BMI", "Expected inter-individual BMI variability"))
    
    plot(soltime, S_BMI, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50), col=2)
    lines(soltime, S_BMIinf2, lty=2, col=2)
    lines(soltime, S_BMIsup2, lty=2, col=2)
    points(T2[gs], BMI2[gs], col=2)
    points(T5[gs], BMI5[gs], col=2)
    title(main="BMI for Weight Stable (WS)")
    legend("topleft", cex=0.7, lty=c(1,2), col=c(2,2), legend=c("Average BMI", "Expected inter-individual BMI variability"))
  }
  
  else{
    plot(soltime, R_BMI, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50))
    lines(soltime, R_BMIinf, lty=2)
    lines(soltime, R_BMIsup, lty=2)
    points(T2[gr], BMI2[gr])
    points(T5[gr], BMI5[gr])
    title(main="BMI for WR")
    legend("topleft", cex=0.7, lty=c(1,2), col=c(1,1), legend=c("Average BMI", "Expected inter-individual BMI variability"))
    
    plot(soltime, S_BMI, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50), col=2)
    lines(soltime, S_BMIinf, lty=2, col=2)
    lines(soltime, S_BMIsup, lty=2, col=2)
    points(T2[gs], BMI2[gs], col=2)
    points(T5[gs], BMI5[gs], col=2)
    title(main="BMI for WS")
    legend("topleft", cex=0.7, lty=c(1,2), col=c(2,2), legend=c("Average BMI", "Expected inter-individual BMI variability"))
  }
  
  if (superposition){
    plot(soltime, R_BMI, type='l', xlab = "Days", ylab="BMI (kg/m²)", ylim=c(20,50))
    lines(soltime, S_BMI, col=2)
    title(main="BMI Comparison")
    legend("topright", lty=c(1,1), col=c(1,2), legend=c("WR", "WS"))
  }
}

#ENERGY RATES GRAPHS

avE <- function(superposition = T, details = F)
{
  plot(soltime, R_EI, type="l", ylim=c(7000, 14000), xlab="Days", ylab="Energy rates (kJ/d)")
  lines(soltime, R_EE, type="l", lty=4)
  title(main="Energy rates for WR")
  legend("bottomright", cex=0.7, lty=c(1,4), legend=c("Average energy intake (EI) rate", "Average energy expenditure (EE) rate"))
  
  plot(soltime, S_EI, type="l", ylim=c(7000, 14000), xlab="Days", ylab="Energy rates (kJ/d)", col=2)
  lines(soltime, S_EE, type="l", lty=4, col=2)
  title(main="Energy rates for WS")
  legend("bottomright", cex=0.7, lty=c(1,4), col=c(2,2), legend=c("Average EI rate", "Average EE rate"))
  
  if(superposition){
    plot(soltime, R_EI, type="l", ylim=c(7000, 14000), xlab="Days", ylab="Energy rates (kJ/d)")
    lines(soltime, R_EE, lty=4)
    lines(soltime, S_EI, col=2)
    lines(soltime, S_EE, lty=4, col=2)
    title(main="Energy rates Comparison")
    legend("bottomright", cex=0.7, lty=c(1,4,1,4), col=c(1,1,2,2), legend=c("Average EI rate for WR", "Average EE rate for WR", "Average EI rate for WS", "Average EE rate for WS"))
  }
  
  if (details){
    plot(soltime, R_EI, type="l", ylim=c(7000, 14000), xlab="Days", ylab="Energy Intake rate (kJ/d)", col=6)
    lines(soltime, R_EIinf, lty=2, col=6)
    lines(soltime, R_EIsup, lty=2, col=6)
    title(main="EI for WR")
    legend("bottomright", cex = 0.7, lty=c(1,2), col=c(6,6), legend=c("Average EI rate", "Prediction Interval"))
    
    plot(soltime, R_EE, type="l", ylim=c(8500, 13000), xlab="Days", ylab="Energy Expenditure rate (kJ/d)", col=4)
    lines(soltime, R_EEinf, col=4, lty=2)
    lines(soltime, R_EEsup, col=4, lty=2)
    title(main="EE for WR")
    legend("bottomright", cex = 0.7, lty=c(1,2), col=c(4,4), legend=c("Average EE rate", "Prediction Interval"))
    
    plot(soltime, S_EI, type="l", ylim=c(5800, 12000), xlab="Days", ylab="Energy intake rate (kJ/d)", col=6)
    lines(soltime, S_EIinf, lty=2, col=6)
    lines(soltime, S_EIsup, lty=2, col=6)
    title(main="EI for WS")
    legend("bottomright", cex = 0.7, lty=c(1,2), col=c(6,6), legend=c("Average EI rate", "Prediction Interval"))
    
    plot(soltime, S_EE, type="l", ylim=c(7500, 12000), xlab="Days", ylab="Energy expenditure rate (kJ/d)", col=4)
    lines(soltime, S_EEinf, col=4, lty=2)
    lines(soltime, S_EEsup, col=4, lty=2)
    title(main="EE for WS")
    legend("bottomright", cex = 0.7, lty=c(1,2), col=c(4,4), legend=c("Average EE rate", "Prediction Interval"))
  }
}

#LIPID AGE GRAPH

avLA <- function(superposition = T)
{
  plot(soltime, R_bestfit[,4], type ="l", xlab = "Days", ylab="Lipid Age (d)", ylim=c(min(R_fitinf2[,4]), max(R_fitsup2[,4])))
  lines(soltime, R_fitinf2[,4], lty=2)
  lines(soltime, R_fitsup2[,4], lty=2)
  points(T2[gr], LA2[gr])
  points(T5[gr], LA5[gr])
  title(main="Lipid Age for WR")
  legend("topright", cex = 0.7, lty=c(1,2), legend=c("Average LA", "Prediction Interval"))
  
  plot(soltime, S_bestfit[,4] , type ="l", xlab = "Days", ylab="Lipid Age (d)", ylim=c(min(S_fitinf2[,4]), max(S_fitsup2[,4])), col=2)
  lines(soltime, S_fitinf2[,4], lty=2, col=2)
  lines(soltime, S_fitsup2[,4], lty=2, col=2)
  points(T2[gs], LA2[gs], col=2)
  points(T5[gs], LA5[gs], col=2)
  title(main="Lipid Age for WS")
  legend("topright", cex = 0.7, lty=c(1,2), col=c(2,2), legend=c("Average LA", "Prediction Interval"))
  
  if (superposition){
    plot(soltime, S_bestfit[,4] , type ="l", xlab = "Days", ylab="Lipid Age (d)", ylim=c(min(S_fitinf2[,4]), max(S_fitsup2[,4])), col=2)
    lines(soltime, R_bestfit[,4]) 
    title(main="Lipid Age Comparison")
    legend("topright", lty=c(1,1), col=c(1,2), legend=c("WR", "WS"))
  }
  
}

#KOUT GRAPH

avKout <- function(superposition = T)
{
  plot(soltime, R_Kout, type="l", xlab="Days", ylab="Kout (/d)", ylim=c(min(R_Koutsup), max(R_Koutinf)))
  lines(soltime, R_Koutinf, lty=2)
  lines(soltime, R_Koutsup, lty=2)
  title(main="Kout for WR")
  legend("topright", cex = 0.7, lty=c(1,2), legend=c("Average Kout", "Prediction Interval"))
  
  plot(soltime, S_Kout, type="l", xlab="Days", ylab="Kout (/d)", ylim=c(min(S_Koutsup), max(S_Koutinf)), col=2)
  lines(soltime, S_Koutinf, lty=2, col=2)
  lines(soltime, S_Koutsup, lty=2, col=2)
  title(main="Kout for WS")
  legend("topright", cex = 0.7, lty=c(1,2), col=c(2,2), legend=c("Average Kout", "Prediction Interval"))
  
  if (superposition){
    plot(soltime, S_Kout, type="l", xlab="Days", ylab="Kout (/d)", ylim=c(min(S_Koutsup), max(S_Koutinf)), col=2)
    lines(soltime, R_Kout)
    legend("topright", legend=c("WR", "WS"), lty=c(1,1), col=c(1,2))
    title(main="Kout Comparison")
  }
}

############################
# COMPARISON TO BW PLANNER #
############################
ind <- 38
dat2 <- read.table(file = "planner.csv", sep= ',', header = TRUE)
ptime <- dat2$Day
pBW <- dat2$Weight
pFM <- dat2$FM
pLM <- pBW - pFM
pEI <- dat2$Intake
pEE <- dat2$Expenditure
pEIsurg <- pEI[1]
pEIfinal <- pEI[length(pEI)]

parametersc <- c(EI0[ind], K[ind], H[ind], age0[ind], pEIsurg, pEIfinal)
initc <- c(fatmass = FM0[ind],leanmass = LM0[ind])
bestfitc <- lsoda(y=initc, times=soltime, func = EqBW2, parms = c(parametersc))

graphEIc <- c()
graphEEc <- c()
for (i in 1:2201){
  graphEIc[i] <- EI2(soltime[i], EI0[ind], pEIsurg, pEIfinal)
  graphEEc[i] <- EE2(soltime[i], EI0[ind], pEIsurg, pEIfinal, K[ind], H[ind], age0[ind], bestfitc[i,2], bestfitc[i,3])
}

plot(soltime, graphEIc, type="l", xlab="Days", ylab="Energy rate ", col=1, xlim=c(0, 500), ylim=c(1000, 15000))
lines(ptime, pEI, lty =2)
lines(soltime, graphEEc, type ="l", lty = 1, col=2)
lines(ptime, pEE, type ="l", lty = 2, col=2)
legend("bottomright",lty=c(1,1,2,2), cex=0.7, col=c(1,2,1,2), legend=c("Energy Intake rate", "Energy Expenditure rate", "Planner EI", "Planner EE"))

graphFMc <- bestfitc[,2]
graphLMc <- bestfitc[,3]
graphBWc <- bestfitc[,2]+bestfitc[,3]

plot(soltime, graphFMc, type="l", xlim=c(0, 500), ylim=c(0,140), xlab="Days", ylab="Weight in kg")
lines(ptime, pFM, type = "l", lty =2)
lines(soltime, graphBWc, type = "l", lty =1, col=4)
lines(ptime, pBW, type = "l", lty =2, col=4)
lines(soltime, graphLMc, type = "l", lty =1, col="dodgerblue4")
lines(ptime, pLM, type = "l", lty =2, col="dodgerblue4")
legend("bottomright", lty=c(1,1,1,2,2,2), legend=c("BW", "FM","LM", "Planner BW", "Planner FM", "Planner LM"), col=c(4,1,"dodgerblue4",4,1,"dodgerblue4"), cex=0.7)

############################
# MEAN COMPAR TEST R AND S #
############################
comp <- t.test(BMI0[gr], BMI0[gs])
comp2 <- t.test(LA0[gr], LA0[gs])
comp3 <- t.test(allLA[2201,1:13], allLA[2201,13:39])
comp4 <- t.test(LA0[gr], allLA[2201,1:13])
comp5 <- t.test(LA0[gs], allLA[2201,14:39])

comp6 <- t.test(allLA[1,1:13], allLA[2201,1:13])
comp7 <- t.test(allLA[1,14:39], allLA[2201,14:39])


# sdr <- sd(R_bestfit[T0[gr]+1,4]-R_bestfit[T5[gr],4])
sdr <- sd(LA0[gr]-LA5[gr])
sds <- sd(LA0[gs]-LA5[gs])
sdr2 <- sd(LA0[gr]-allLA[1:13,4])
sds2 <- sd(LA0[gs]-allLA[14:39,4])
deltaR <- mean(LA0[gr]-LA5[gr])
deltaR2 <- mean(R_bestfit[1,4]-R_bestfit[2201,4])
deltaS <- mean(LA0[gs]-LA5[gs])
deltaS2 <- mean(S_bestfit[1,4]-S_bestfit[2201,4])

deltaRS <- abs(deltaR-deltaS)
deltaRS2 <- abs(deltaR2-deltaS2)

power.t.test(n=13, delta=deltaR, sd=sdr)
power.t.test(delta=deltaR, sd=sdr, power=0.8)
power.t.test(n=13, sd=sdr, power=0.8)
power.t.test(delta=deltaR2, sd=sdr2, power=0.8)
power.t.test(n=13,delta=deltaR2, sd=sdr)

power.t.test(n=26, delta=deltaS, sd=sds)
power.t.test(delta=deltaS, sd=sds, power=0.8)
power.t.test(n=26, sd=sds, power=0.8)
power.t.test(delta=deltaS2, sd=sdr, power=0.8)
power.t.test(n=26,delta=deltaS2, sd=sdr)


power.t.test(n=13, delta=deltaRS, sd=sds, type="paired")
power.t.test(n=13, delta=deltaRS, sd=sds)
power.t.test(delta=deltaRS, sd=sds, p=0.8)
power.t.test(n=13, sd=sds, p=0.8)

power.t.test(n=13, delta=deltaRS2, sd=sds)
power.t.test(delta=deltaRS2, sd=sds, p=0.8)
power.t.test(n=13, sd=sds, p=0.8)









