#POM_MAON_CN_model as a stode function.
CN_NUE_model<- function(t,y,pars){
  with(as.list(c(y,pars)),{
    
    #C fluxes
    DEATH.C      <- h1*B^1.5
    DECOMP.C     <- v1*B*C / (k1 + C)
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
    
    #leaching/plant uptake happens on post-mineralization/immobilization inorganic N pool.
    #Make sure to include inputs to N4 from imperfect NUE values
    INORG.N.LOSS <- h4*(N4 + DECOMP.N*(1-NUE) + mineralization - immobilization)
    
    #OVERFLOW RESPIRATION
    #excess C uptake must go somewhere because mass balance.
    if(mineralization.immobilization > 0){
      overflow.R = 0
    }
    if(mineralization.immobilization <= 0){  
      overflow.R <- DECOMP.C*CUE - (DECOMP.N*NUE + immobilization)*BN
    }
    
    #ODEs
    dCdt  <- I + (DEATH.C - SORPTION.B.C) + DESORPTION.C - DECOMP.C - SORPTION.P.C - exo.loss.C
    dMdt  <- SORPTION.C - DESORPTION.C
    dBdt  <- (CUE * DECOMP.C) - DEATH.C - overflow.R
    dN1dt <- (I/IN) + (DEATH.N - SORPTION.B.N) + DESORPTION.N - DECOMP.N - SORPTION.P.N - exo.loss.N
    dN2dt <- SORPTION.N - DESORPTION.N
    dN3dt <- DECOMP.N*NUE + immobilization - mineralization - DEATH.N
    dN4dt <- DECOMP.N*(1-NUE) + mineralization - immobilization - INORG.N.LOSS
    
    #list of shit to solve for.
    list(c(dCdt,dMdt,dBdt,dN1dt,dN2dt,dN3dt,dN4dt))
  })
}