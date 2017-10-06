#POM_MAON_CN_model as a stode function.
CN_NUE_model<- function(t,y,pars){
  with(as.list(c(y,pars)),{
    
    #C fluxes
    DEATH.C      <- h1*B^1.5
    DECOMP.C     <- v1*B*C / (k1 + C) * (1 - ((7-pH)*ph.mod))
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
    #if(day > fert.day){
    #  fert = fert.level
    #}
    
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
    
    #pH mode toggles whether or not pH is constant (5.5), or emergent from proton pool.
    #if(pH_mode == T){pH <- -log10(protons)}
    #if(pH_mode == F){pH = 5.5}
    
    #list of differential equations to solve for.
    list(c(dCdt,dMdt,dBdt,dN1dt,dN2dt,dN3dt,dN4dt,dPROTdt))
  })
}