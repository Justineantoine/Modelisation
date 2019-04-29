rm(list=ls())
library(deSolve)

dat <- read.table(file = "CLINICAL(copietravail).csv", sep= ',', header = TRUE)
 ######################################
 # BODYWEIGHT, FAT MASS AND LEAN MASS #
 ######################################
BW0 <- dat$BW0
BW1 <- dat$BW1
BW2 <- dat$BW2
BW5 <- dat$BW5

FM0 <- dat$TBFI0
FM2 <- dat$TBFI2
FM5 <- dat$TBFI5

adjust = subset(FM2/dat$TBFB2, is.na(FM2/dat$TBFB2)==F)
mean_adjust = mean(adjust)

for (i in 1:41){
  if (is.na(FM0[i])){
    FM0[i] = mean_adjust*dat$TBFB0[i]
  }
}
FM5[32]=mean_adjust*dat$TBFB5[32]

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
T5 <- dat$T5

 ##############
 # PARAMETERS #
 ##############
pf <- 39.5 #kJ/g
pl <- 7.6 #kJ/g
gf <- 13 #kJ/kg/day
gl <- 92 #kJ/kg/day
nf <- 0.75 #kJ/g 
nl <- 0.96 #kJ/g
Btef <- 0.1
Bat <- 0.14
PAL <- 1.5
C <- 10.4

 #####################
 # PHYSICAL ACTIVITY #
 #####################
d0 <- ((1-Btef)*1.5-1)*(10*BW0+6.25*H-5*age0-161)*4.184
d2 <- ((1-Btef)*1.5-1)*(10*BW2+6.25*H-5*age2-161)*4.184
d5 <- ((1-Btef)*1.5-1)*(10*BW5+6.25*H-5*age5-161)*4.184

 #######################
 # ENERGY PARTITIONING #
 #######################
p0 <- pl/pf*C/(pl/pf*C+FM0) 
p2 <- pl/pf*C/(pl/pf*C+FM2) 
p5 <- pl/pf*C/(pl/pf*C+FM5) 

 #########################
 # INITIAL ENERGY INTAKE #
 #########################
EI0 <- (655 + 9.56*BW0 + 186*H -4.68*age0)*1.33*4.184

 ##############################
 # INITIAL ENERGY EXPENDITURE #
 ##############################
EE0 <- EI0

 ##############
 # K CONSTANT #
 ##############
K <- EI0 - gf*FM0 - gl*LM0 - d0

 #################
 # ENERGY INTAKE #
 #################
EIsurg <- rep(7000, 41)

EI <- function(i, t) #i pour l'individu auquel on s'interesse
{
  s = (EI5[i]-EIsurg[i])/(900-150) #garde le regime pendant 150 jours (~5 mois) puis reprise progressive de nourriture jusqu'au 900eme jour (~3 ans)
  if (t<=0){
    EI0[i]
  }
  else if (t<=150){
    EIsurg[i]
  }
  else{
    intake = (t-150)*s+EIsurg[i]
    if (intake<EI5[i]){
      intake
      }
    else{
      EI5[i]
    }
  }
}

EI2 <- c()
EI5 <- (K + gf*FM5 + gl*LM5 + d5 +  EI0*(Btef+Bat))/(1+Btef+Bat)

for (i in 1:41){
  EI2[i] <- EI(i, T2[i])
}

######################
# ENERGY EXPENDITURE #
######################
EE0 <- EI0
EE2 <- (K + gf*FM2 + gl*LM2 + (Bat+Btef)*(EI2-EI0) + EI2*(p2*nl/pl + (1-p2)*nf/pf))/(1 + p2*nl/pl + (1-p2)*nf/pf)
EE5 <- EI5


test = c()
test2 = c()

temps <- seq(-100, 2200)
for (c in 1:2301){
  test[c] <- EI(2,temps[c])
  test2[c] <- EI(1,temps[c])
}

plot(temps, test, type="l", xlab="Temps en jours", ylab="Energy rate ", col=2)
points(0, EE0[2], col=2)
points(T2[2], EE2[2], col=2)
points(T5[2], EE5[2], col=2)
lines(temps, test2, type="l", col=1)
points(0, EE0[1], col=1)
points(T2[1], EE2[1], col=1)
points(T5[1], EE5[1], col=1)
