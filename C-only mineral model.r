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
v1   <- .4       #biomass-specific decay multiplier   (mg time-1)
v2   <- 1        #Vmax of sorption                    (mg time-1)
k1   <- 100      #half saturation of decomp           (mg)
k2   <- 10       #half saturation of sorption         (mg)
h1   <- 0.05     #biomass turnover rate               (1/time)
h2   <- 0.002    #C-specific desorption rate          (1/time)
h3   <- 0.001    #fraction of POM that potentially sorbs (1/time)

#initial pool sizes
C <- 100       #POM               (mg / g)
M <-  20       #MAOM              (mg / g)
B <- C*0.1    #microbial biomass  (mg / g)

#number of days to step the dynamic simulation through time.
t <- 1500

#choose outputs you want to follow. In this case, the 3 state variables. 
outputs <- c('B','C','M')
out <- matrix(rep(0,t*length(outputs)),nrow=t,dimnames=list(NULL,outputs))


for(i in 1:t){
  #defining fluxes
  DEATH      <- h1*B
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
  
  #save outputs
  out[i,] <- c(B,C,M)
}

par(mfrow=c(1,3))
plot(out[,1])
plot(out[,2])
plot(out[,3])


#use stode to get numerical solution for state variables.
require(rootSolve)
model<- function(t,y,pars){
  with(as.list(c(y,pars)),{
    dBdt <- CUE * v1 * B * C / (k1 + C) - h1*B
    dCdt <- I + (h1*B - (h1*B  /(h3*C+h1*B))*v2*(h3*C+h1*B) / (k2 + h3*C+h1*B)) + h2*M - v1 * B * C / (k1 + C) - (h3*C/(h3*C+h1*B))*v2*(h3*C+h1*B) / (k2 + h3*C+h1*B)
    dMdt <- v2*(h3*C+h1*B) / (k2 + h3*C+h1*B) - h2*M
    list(c(dBdt,dCdt,dMdt))
  })
}

#define parameters and initial pool sizes.
pars <- c(CUE=CUE, v1=v1, v2=v2, k1=k1, k2=k2, I=I, h1=h1, h2=h2)
y <- c(B=B,C=C,M=M)

#numerically solve the model.
ST <- stode(y=y,func=model,parms=pars,pos=T)
ST$y; sum(ST$y)
