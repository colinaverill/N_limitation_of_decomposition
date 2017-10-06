#This script is used to check model benchmarks at end of N-fertilization experiment.
rm(list=ls())

#specify model, parameter, and output paths
parameters.path <- 'av_parameters.r'

#number of days to step the dynamic simulation through time.
source(parameters.path)
t <- 400000
fert.level <- 0.05 / 365
fert.day   <- 300000      #which day to begin fertilization.
day <- 0
fert = 0

#Input N has be at 30, 60, 80
IN <- 80
#v2- sorption - is at 0.5 and 2.5 for low and high, respectively.
v2 <- 2.5
#whether or not we are in pH mode. If pH mode set to F, pH constant at 5.5.
pH_mode <- F

#choose outputs you want to follow. 
outputs <- c('C','M','B','R','N1','N2','N3','N4','pH','day')
out <- matrix(rep(0,t*length(outputs)),nrow=t,dimnames=list(NULL,outputs))

#model loop
for(i in 1:t){
  #C fluxes
  DEATH.C      <- h1*B^1.5
  DECOMP.C     <- v1*B*C / (k1 + C) * (1 - ((7-pH)*ph.mod)) #adding pH mod to decomp flux.
  exo.loss.C   <- h5*C
  POM2MOM.C    <- h3*C
  SORPTION.C   <- v2*(POM2MOM.C+DEATH.C) / (k2 + POM2MOM.C+DEATH.C)
  SORPTION.B.C <- (DEATH.C  /(POM2MOM.C+DEATH.C))*SORPTION.C
  SORPTION.P.C <- (POM2MOM.C/(POM2MOM.C+DEATH.C))*SORPTION.C
  DESORPTION.C <- h2*M

  #N fluxes
  DEATH.N      <- DEATH.C / BN
  DECOMP.N     <- DECOMP.C/ (C / N1)
  exo.loss.N   <- exo.loss.C / (C/N1)
  POM2MOM.N    <- POM2MOM.C / (C/N1)
  SORPTION.N   <- v2*(POM2MOM.N+DEATH.N) / (k2 + POM2MOM.N+DEATH.N)
  SORPTION.B.N <- (DEATH.N  /(POM2MOM.N+DEATH.N))*SORPTION.N
  SORPTION.P.N <- (POM2MOM.N/(POM2MOM.N+DEATH.N))*SORPTION.N
  DESORPTION.N <- DESORPTION.C/(M/N2)
  
  
  #MINERALIZATION-IMMOBILIZATION.
  #N uptake minus demand. If in excess, positive. If limited, negative.
  mineralization.immobilization <- (NUE * DECOMP.N) - (CUE * DECOMP.C) / BN
  if(mineralization.immobilization  > 0) {
    #if in excess (+ value), mineralize. 
    mineralization <- mineralization.immobilization
    immobilization = 0
  }
  if(mineralization.immobilization <= 0) {
    #if in debt (- value), immobilize, w/ non-linear uptake kinetics.
    mineralization = 0
    immobilization.potential <- ((N4 / (k3 + N4)) * N4)
    immobilization <- immobilization.potential
    #however, don't take up more N than you need.
    if(immobilization > abs(mineralization.immobilization)){
       immobilization = abs(mineralization.immobilization)
    }
  }
  
  #Drop in N fertilization when day = fert.day
  if(day > fert.day){
    fert = fert.level
  }
  
  #leaching/plant uptake happens on post-mineralization/immobilization inorganic N pool.
  #Make sure to include inputs to N4 from imperfect NUE values
  INORG.N.LOSS <- h4*(N4 + DECOMP.N*(1-NUE) + fert + mineralization - immobilization)

  #OVERFLOW RESPIRATION
  #excess C uptake must go somewhere because mass balance.
  if(mineralization.immobilization > 0){
    overflow.R = 0
  }
  if(mineralization.immobilization <= 0){  
    overflow.R <- DECOMP.C*CUE - (DECOMP.N*NUE + immobilization)*BN
  }
  
  #pH sub-routine.
  #acidification is linked to nitrification.
  #nitrification is linear with respect to the sum of all inputs to the inorganic N pool.
  #nitrification decreases as pH decreases.
       nitrif <- (pH / pH_opt)*(DECOMP.N*(1-NUE) + mineralization + fert)
  proton_gain <- nitrif * nitr.acid
  proton_loss <- protons* acid.loss
  #Soil Buffering - this never kicks in with our parameter set.
  if (protons < 1e-07) {proton_loss = 0} else {proton_loss = proton_loss} #pH cannot rise above 7
  if (protons > 1e-04) {proton_gain = 0} else {proton_gain = proton_gain} #pH cannot drop below 4
  
  #ODEs
  dBdt  <- (CUE * DECOMP.C) - DEATH.C - overflow.R
  dCdt  <- I + (DEATH.C - SORPTION.B.C) + DESORPTION.C - DECOMP.C - SORPTION.P.C - exo.loss.C
  dMdt  <- SORPTION.C - DESORPTION.C
  dN1dt <- (I/IN) + (DEATH.N - SORPTION.B.N) + DESORPTION.N - DECOMP.N - SORPTION.P.N - exo.loss.N
  dN2dt <- SORPTION.N - DESORPTION.N
  dN3dt <- (NUE * DECOMP.N) + immobilization - mineralization - DEATH.N
  dN4dt <- DECOMP.N*(1-NUE) + mineralization + fert - immobilization - INORG.N.LOSS
  dPROTdt <- proton_gain - proton_loss
  
  #respiration
  R  <- (1-CUE)*DECOMP.C + overflow.R
  
  #update pools
  C  <- C  + dCdt
  M  <- M  + dMdt
  B  <- B  + dBdt
  N1 <- N1 + dN1dt
  N2 <- N2 + dN2dt
  N3 <- N3 + dN3dt
  N4 <- N4 + dN4dt
  day <- day + 1
  protons <- protons + dPROTdt
  
  #pH mode toggles whether or not pH is constant (5.5), or emergent from proton pool.
  if(pH_mode == T){pH <- -log10(protons)}
  if(pH_mode == F){pH = 5.5}
  
  #save outputs
  out[i,] <- c(C,M,B,R,N1,N2,N3,N4,pH,day)
} #end model loop.

out <- subset(out, out[,10] > 290800)
#out <- subset(out, out[,10] < 350000)   #most simulation responses stable by this point.
out <- out[seq(1,nrow(out),by = 100),] #only need every 100th day for plotting.
out <- data.frame(out)

#plot output.
par(mfrow=c(2,3))
plot(out$C, ylim = c(0,200)); mtext('POM', side = 3)
plot(out$M, ylim = c(0,200)); mtext('MAOM', side = 3)
plot(out$B, ylim = c(0,10)); mtext('biomass', side = 3)
plot((out$C + out$M + out$B), ylim = c(0,200)); mtext('total C', side = 3)
plot(out$pH, ylim = c(4,7)); mtext('pH', side = 3)
plot(out$R); mtext('respiration', side = 3)

cat('BENCHMARKS',
    '\npH =',pH,
    '\nmicrobial biomass C:N =',B/N3,'(better be 7)',
    '\ntotal C is',(C+B+M),'mg per g soil, (target 10-200mg)',
    '\nmicrobial biomass C is',100*B/(C+M+B),'% of total C (target 0.5-10%)',
    '\ninorganic N is ',1000*N4,'ug per g soil (target 0-10ug)',
    '\nrespiration is',(R*1000) / ((C+M+B)/1000),'ug C per g soil C per day (target 70-1500ug)')
