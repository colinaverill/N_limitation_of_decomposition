#Building a linked C-N model examining N fertilization effects on POM and MAOM formation. 
###FIXED CUE AND NUE###

###ASSUMPTIONS###
# C:N stoichiometry of MAOM, POM, and microbial biomass are fixed
# microbes take up organic N (from decomposition of POM) before taking up inorganic N
# microbes aim to have the C:N ratio of uptake match their biomass stoichiometry
     # when microbes are C-limited, they mineralize excess N
     # when microbes are N-limited, they immobilize sufficient inorganic N to match stoichiometric demand
     # if there is not enough inorganic N to match demand, then they respire excess C such that the C:N ratio of uptake matches    biomass stoichiometry

#clear R environment.
rm(list=ls())

#parameters: C cycling
I    <- 1.0      #C input rate                        (mg time-1)
CUE  <- 0.3      #carbon use efficiency               (unitless)
NUE <- 0.8 		 #N use efficiency					  (unitless)
v1   <- .4       #biomass-specific decay multiplier   (mg time-1)
v2   <- 1        #Vmax of sorption                    (mg time-1)
k1   <- 100      #half saturation of decomp           (mg)
k2   <- 10       #half saturation of sorption         (mg)
h1   <- 0.02    #biomass turnover rate               (1/time)
h2   <- 0.002    #C-specific desorption rate          (1/time)
h3   <- 0.001    #fraction of POM that potentially sorbs (1/time)
h4 <- 0.01       #inorganic N loss rate

#parameters: N cycling
IN <- 50 #C:N ratio of inputs
CN <- 50 #C:N ratio of particulate organic matter
MN <- 25 #C:N ratio of mineral-associated organic matter
BN <- 7 #C:N ratio of the microbial biomass

#initial pool sizes
C <- 100       #C pool in POM               (mg / g)
M <-  20       #C pool in MAOM              (mg / g)
B <- C*0.1    #microbial biomass C 			(mg / g)
R <- 0        # respired C					(mg / g)
N1 <- C/CN    #N pool in POM				(mg / g)
N2 <- M/MN	  #N pool in MAOM 				(mg / g)
N3 <- B/BN    #microbial biomass N			(mg / g)
N4 <- 0.001	  #inorganic N pool size 		(mg / g)

#number of days to step the dynamic simulation through time.
t <- 1500

#choose outputs you want to follow. 
outputs <- c('B','C','M','R','N1','N2','N3','N4')
out <- matrix(rep(0,t*length(outputs)),nrow=t,dimnames=list(NULL,outputs))


for(i in 1:t){
  #defining fluxes
  DEATH.C      <- h1*B^1.5
  DECOMP.C     <- v1*B*C / (k1 + C)
  POM2MOM.C    <- h3*C
  SORPTION.C   <- v2*(POM2MOM.C+DEATH.C) / (k2 + POM2MOM.C+DEATH.C)
  SORPTION.B.C <- (DEATH.C  /(POM2MOM.C+DEATH.C))*SORPTION.C
  SORPTION.P.C <- (POM2MOM.C/(POM2MOM.C+DEATH.C))*SORPTION.C
  DESORPTION.C <- h2*M
  
  DEATH.N      <- h1*N3^1.5
  DECOMP.N     <- (v1*B*C / (k1 + C))/(C/N1)
  POM2MOM.N    <- (h3*C)/(C/N1)
  SORPTION.N   <- v2*(POM2MOM.N+DEATH.N) / (k2 + POM2MOM.N+DEATH.N)
  SORPTION.B.N <- (DEATH.N  /(POM2MOM.N+DEATH.N))*SORPTION.N
  SORPTION.P.N <- (POM2MOM.N/(POM2MOM.N+DEATH.N))*SORPTION.N
  DESORPTION.N <- (h2*M)/(M/N2)
  

  #ODEs
  dBdt <- (CUE * DECOMP.C) - DEATH.C
  dCdt <- I + (DEATH.C - SORPTION.B.C) + DESORPTION.C - DECOMP.C - SORPTION.P.C
  dMdt <- SORPTION.C - DESORPTION.C
  dRdt <- (1-CUE)*DECOMP.C
  
  dN1dt <- (I/IN) + (DEATH.N - SORPTION.B.N) + DESORPTION.N - DECOMP.N - SORPTION.P.N
  dN2dt <- SORPTION.N - DESORPTION.N

  #microbial N uptake is dependent on nutrient limitation status
  N.uptake.initial <- (NUE * DECOMP.N) - DEATH.N
  if ((dBdt/N.uptake.initial) < BN) {
  	x <- -((dBdt/BN)-N.uptake.initial)
  	dN4dt <- x
  	dN3dt <- N.uptake.initial - dN4dt
  } 
  else {
  	x <- -((dBdt/BN)-N.uptake.initial)
  	if (N4 >= x) {immob = x} else {immob = N4}
  	dN4dt <- -immob
  	dN3dt <- N.uptake.initial + immob
  	if (dBdt/dN3dt > BN) {
  		y <- (dN3dt*BN) - dBdt
  		dBdt <- dBdt - y
  		dRdt <- dRdt + y
  	} else {dBdt = dBdt}
  }

  #update pools
  B <- B + dBdt
  C <- C + dCdt
  M <- M + dMdt
  R <- R + dRdt
  N1 <- N1 + dN1dt
  N2 <- N2 + dN2dt
  N3 <- N3 + dN3dt
  N4 <- N4 + dN4dt - (h4*N4)
  if (N4 > 0) {N4=N4} else {N4 = 0} #inorganic N pool cannot be negative
  
  
  #save outputs
  out[i,] <- c(B,C,M,R,N1,N2,N3,N4)
}

par(mfrow=c(3,3))
plot(out[,1], ylab="microbial biomass C")
plot(out[,2], ylab="POM C")
plot(out[,3], ylab= "MAOM C")
plot(out[,4], ylab="respired C")
plot(out[,5], ylab="POM N")
plot(out[,6], ylab= "MAOM N")
plot(out[,7], ylab="microbial biomass N")
plot(out[,8], ylab="inorganic N")

#use stode to get numerical solution for state variables.
require(rootSolve)
model<- function(t,y,pars){
  with(as.list(c(y,pars)),{
  dBdt <- (CUE * (v1*B*C / (k1 + C))) - (h1*B^1.5)
  dCdt <- I + ((h1*B^1.5) - (((h1*B^1.5)  /((h3*C)+(h1*B^1.5)))*(v2*((h3*C)+(h1*B^1.5)) / (k2 + (h3*C)+(h1*B^1.5))))) + h2*M - (v1*B*C / (k1 + C)) - (((h3*C)/((h3*C)+(h1*B^1.5)))*(v2*((h3*C)+(h1*B^1.5)) / (k2 + (h3*C)+(h1*B^1.5))))
  dMdt <- (v2*((h3*C)+(h1*B^1.5)) / (k2 + (h3*C)+(h1*B^1.5))) - h2*M
  dRdt <- (1-CUE)*(v1*B*C / (k1 + C))
  dN1dt <- (I/IN) + ((h1*N3^1.5) - (((h1*N3^1.5)/(((h3*C)/(C/N1))+(h1*N3^1.5)))*(v2*(((h3*C)/(C/N1))+(h1*N3^1.5)) / (k2 + ((h3*C)/(C/N1)) + (h1*N3^1.5))))) + ((h2*M)/(M/N2)) - ((v1*B*C / (k1 + C))/(C/N1)) - ((((h3*C)/(C/N1))/(((h3*C)/(C/N1))+(h1*N3^1.5)))*(v2*(((h3*C)/(C/N1))+(h1*N3^1.5)) / (k2 + (h3*C)/(C/N1)+(h1*N3^1.5))))
 dN2dt <- (v2*(((h3*C)/(C/N1))+(h1*N3^1.5)) / (k2 + ((h3*C)/(C/N1))+(h1*N3^1.5))) - ((h2*M)/(M/N2))
   list(c(dBdt, dCdt, dMdt, dN1dt, dN2dt))
  })
}

#define parameters and initial pool sizes.
pars <- c(CUE=CUE, NUE=NUE, BN=BN, I=I, IN=IN, v1=v1, v2=v2, k1=k1, k2=k2, h1=h1, h2=h2, h3=h3)
y <- c(B=B, C=C, M=M, N1=N1, N2=N2)

#numerically solve the model.
ST <- stode(y=y,func=model,parms=pars,pos=T)
ST$y; sum(ST$y)
