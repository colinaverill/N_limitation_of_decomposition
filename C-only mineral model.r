#making the simplest mineral POM and MAOM model I can.
#No input sensitivity of POM, as suspected. Increasing inputs increasing MBC. 
#This also increases MAOM, as we increase the input rate at steady state, even with non-linear desorption kinetics.
#sorption follows langmuir isotherm adsorption model, only microbial biomass turnover enters the sorbed pool. 
#There is also potential for some fraction of POM to enter the MAOM pool directly, rather than passing through MB.

#clear R environment.
rm(list=ls())

#parameters
I    <- 1.0      #C input rate                        (mg time-1)
CUE  <- 0.3      #carbon use efficiency               (unitless)
v1   <- .5       #biomass-specific decay multiplier   (mg time-1)
v2   <- 1        #Vmax of sorption                    (mg time-1)
k1   <- 100      #half saturation of decomp           (mg)
k2   <- 10       #half saturation of sorption         (mg)
h1   <- 0.02     #biomass turnover rate               (1/time)
h2   <- 0.002    #C-specific desorption rate          (1/time)
h3   <- 0.001    #fraction of POM that potentially sorbs (1/time)

#initial pool sizes
C <- 71      #POM                (mg / g)
M <- 24      #MAOM               (mg / g)
B <-  8.6   #microbial biomass  (mg / g)

#number of days to step the dynamic simulation through time.
t <- 1500

#choose outputs you want to follow. In this case, the 3 state variables. 
outputs <- c('B','C','M','test')
out <- matrix(rep(0,t*length(outputs)),nrow=t,dimnames=list(NULL,outputs))


for(i in 1:t){
  #defining fluxes
  DEATH      <- h1*B^1.5
  DECOMP     <- v1*B*C / (k1 + C)
  POM2MOM    <- h3*C
  SORPTION   <- v2*(POM2MOM+DEATH) / (k2 + POM2MOM+DEATH)
  SORPTION.B <- (DEATH  /(POM2MOM+DEATH))*SORPTION
  SORPTION.P <- (POM2MOM/(POM2MOM+DEATH))*SORPTION
  DESORPTION <- h2*M
  
  #ODEs
  dBdt <- CUE * DECOMP - DEATH
  dCdt <- I + (DEATH - SORPTION.B) + DESORPTION - DECOMP - SORPTION.P
  dMdt <- SORPTION - DESORPTION
  
  #update pools
  B <- B + dBdt
  C <- C + dCdt
  M <- M + dMdt
  test <- dCdt
  
  #save outputs
  out[i,] <- c(B,C,M,test)
}


par(mfrow=c(2,2))
plot(out[,1], ylab="Biomass C")
plot(out[,2], ylab="POM C")
plot(out[,3], ylab="MAOM C")
plot(out[,4], ylab="dCdt")


#use stode to get numerical solution for state variables.
require(rootSolve)
model<- function(t,y,pars){
  with(as.list(c(y,pars)),{
    dBdt <- CUE * v1 * B * C / (k1 + C) - h1*B^1.5
    dCdt <- I + (h1*B^1.5 - (h1*B^1.5  /(h3*C+h1*B^1.5))*v2*(h3*C+h1*B^1.5) / (k2 + h3*C+h1*B^1.5)) + h2*M - v1 * B * C / (k1 + C) - (h3*C/(h3*C+h1*B^1.5))*v2*(h3*C+h1*B^1.5) / (k2 + h3*C+h1*B^1.5)
    dMdt <- v2*(h3*C+h1*B^1.5) / (k2 + h3*C+h1*B^1.5) - h2*M
    list(c(dBdt,dCdt,dMdt))
  })
}

#define parameters and initial pool sizes.
pars <- c(CUE=CUE, v1=v1, v2=v2, k1=k1, k2=k2, I=I, h1=h1, h2=h2, h3=h3)
y <- c(B=B,C=C,M=M)

#numerically solve the model.
ST <- stode(y=y,func=model,parms=pars,pos=T)
ST$y; sum(ST$y)
