#colin had a thought about representing organic N vs. organic C pools explicitly
#requires invoking 'shielding'. C-only must be degraded to access organic N that is shielded.
#C:N of return on org N acquisition is a function of C:N of whole soil.
#When investing in N decomp, C-only comes with.
#Solving optimal allocation must account for this. 
#Effective C:N of org N must = MB_CN / CUE to induce N limitation.
#This needs to happen at a whole soil C:N of ~ 30 (expert knowledge...)

#Status:
#1. we have correct allocation solution, but stoichiometry is violated.
#all my microbes are dying... which is probably linked to problem 1.
#anddd somehow I just broke the whole damn thing. 

#clear R environment.
rm(list=ls())

#parameters: C cycling
I    <- 0.05     #C input rate                           (mg time-1)
CUE  <- 0.3      #carbon use efficiency                  (unitless)
v1   <- .3       #biomass-specific decay multiplier      (mg time-1)
k1   <- 100      #half saturation of decomp              (mg)
k3   <- 0.001    #half saturation of inorganic N uptake  (mg)
h1   <- 0.02     #biomass turnover rate                  (1/time)
h3   <- 0.9      #inorganic N loss rate

#parameters: N cycling
IN      <- 50 #C:N ratio of inputs (constant)
CN.i    <- 30 #C:N ratio of initial particulate organic matter (will change)
m.CN    <-  7 #C:N ratio of the microbial biomass (constant)
CN.org  <-  4 #organic N C:N
soil.CN <- 40 #initial soil CN

#initial pool sizes
C1 <- 100                     #C pool in C1          (mg / g)
C2 <- C1/(soil.CN/CN.org - 1) #C pool in C2          (mg / g)
B <- C1*0.1                   #microbial biomass C 	 (mg / g)
N1 <- C2/CN.org               #initial N pool in C2	 (mg / g)
N2 <- B/m.CN                  #microbial biomass N   (mg / g)
N3 <- 0.001	                  #inorganic N pool size (mg / g)
eco.C <- C1 + C2 + B          #total soil C          (mg / g)
eco.N <- N1 + N2 + N3         #total soil N          (mg / g)

#number of days to step the dynamic simulation through time.
t <- 3000

#choose outputs you want to follow. 
outputs <- c('B','C1','C2','R','N1','N2','N3','eco.C','eco.N')
out <- matrix(rep(0,t*length(outputs)),nrow=t,dimnames=list(NULL,outputs))


for(i in 1:t){
  #tradeoff solution to have uptake C:N match your TER
  x = (C2*(k1 + C1)*(m.CN - CUE*CN.org))/(m.CN*k1*C2 + m.CN*C1*C2 + CUE*k1*CN.org*C1 - CUE*k1*CN.org*C2)
  #set x=.5 to turn off tradeoff for testing.
  x = .5
  #C fluxes
  I.C1          <- I - (I/IN) * CN.org
  I.C2          <-     (I/IN) * CN.org
  DECOMP.C1     <-    x *v1*B*C1 / (k1 + C1)
  DECOMP.C2     <- (1-x)*v1*B*C2 / (k1 + C2)
  uptake.C      <- DECOMP.C1 + DECOMP.C2
  DEATH.C       <- h1*B^1.5
  DEATH.C1      <- DEATH.C - (DEATH.C/m.CN) * CN.org
  DEATH.C2      <-           (DEATH.C/m.CN) * CN.org
  respiration   <- (1-CUE)*uptake.C + overflow.R
  
  #N fluxes
  DEATH.N      <- DEATH.C/m.CN
  DECOMP.N     <- DECOMP.C2/ (C2 / N1)
  uptake.N     <- DECOMP.N
  
  #MINERALIZATION-IMMOBILIZATION.
  #N uptake minus demand. If in excess, positive. If limited, negative.
  mineralization.immobilization <- (uptake.N) - (CUE * uptake.C) / m.CN
  #mineralization statement.
  if(mineralization.immobilization  > 0) {
    #if in excess (+ value), mineralize. 
    mineralization <- mineralization.immobilization
    immobilization = 0
  }
  #immobilization statement.
  if(mineralization.immobilization <= 0) {
    #if in debt (- value), immobilize, w/ non-linear uptake kinetics.
    mineralization = 0
    immobilization.potential <- (N3 / (k3 + N3)) * N3
    immobilization <- immobilization.potential
    #however, don't take up more N than you need.
    if(immobilization > abs(mineralization.immobilization)){
      immobilization = abs(mineralization.immobilization)
    }
  }
  #leaching/plant uptake happens on post-mineralization/immobilization inorganic N pool.
  INORG.N.LOSS <- h3*(N3 + mineralization - immobilization)
  
  #OVERFLOW RESPIRATION
  #excess C uptake must go somewhere because mass balance.
  if(mineralization.immobilization > 0){
    overflow.R = 0
  }
  if(mineralization.immobilization <= 0){  
    overflow.R <- uptake.C*CUE - (uptake.N + immobilization)*m.CN
  }
  
  #ODEs
  dBdt  <- (CUE * uptake.C) - DEATH.C - overflow.R
  dC1dt <- I.C1 + DEATH.C1 - DECOMP.C1
  dC2dt <- I.C2 + DEATH.C2 - DECOMP.C2
  
  dN1dt <- I.C2/IN + DEATH.N - DECOMP.C2/CN.org
  dN2dt <- uptake.N + immobilization - mineralization - DEATH.N
  dN3dt <- mineralization - immobilization - INORG.N.LOSS

  #update pools
  B  <- B  + dBdt
  C1 <- C1 + dC1dt
  C2 <- C1 + dC2dt
  N1 <- N1 + dN1dt
  N2 <- N2 + dN2dt
  N3 <- N3 + dN3dt
  
  #ecosystem sacle C and N
  eco.C <-  B + C1 + C2
  eco.N <- N1 + N2 + N3
  
  #respiration
  R  <- respiration
  
  #save outputs
  out[i,] <- c(B,C1,C2,R,N1,N2,N3,eco.C,eco.N)
}
par(mfrow=c(4,3))
plot(out[,1], ylab="microbial biomass C")
plot(out[,2], ylab="C1")
plot(out[,3], ylab= "C2")
plot(out[,6], ylab="microbial biomass N")
plot(out[,5], ylab="N1 (from C2 pool)")
plot(out[,7], ylab= "inorganic N")
plot(out[,4], ylab="CO2 respiration")
plot(out[,8], ylab="ecosystem C")
plot(out[,9], ylab='ecosystem N')
#check stoichiometry of pools
plot(out[,1]/out[, 6], cex = 0.1, col='blue' , ylim = c(0,40))
par(new=T)
plot(out[,3]/out[, 5], cex = 0.1, col='green', ylim = c(0,40))
par(new=T)
plot(out[,8]/out[,9], cex = 0.1, col='purple'  , ylim = c(0,40))


cat(paste('Ecosystem C:N = ',signif(out[nrow(out),9]/out[nrow(out),10], digits = 4)))
cat(paste('microbial C:N',out[nrow(out),1] / out[nrow(out),6]))
cat(paste(       'C2 C:N',out[nrow(out),3] / out[nrow(out),5]))
