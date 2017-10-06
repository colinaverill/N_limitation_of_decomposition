CN_NUE_pH_model <- function(t,y,pars){
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
      
      #pH
      pH_mod <- pH/pH_opt 
	    N_turnover <- dN4dt
	    Nitr <- pH_mod*N_turnover #nitrification is proportional to pH (further from optimum, the slower it goes) and rate of N turnover sensu Parton et al 1996
	    proton_change <- Nitr*1e-3
	    protons <- protons + proton_change - proton_loss
	    if (protons < 1e-08) {protons = 1e-08} else {protons = protons} #pH cannot rise above a certain threshold due to soil buffering capacity
	    pH <- -log10(protons)
	    if (pH < 4) {pH = 4} else {pH=pH} #pH cannot drop beneath a certain threshold due to soil buffering capacity
      
      #pH effects on biomass - growth is reduced by 18% for each unit change in pH sensu Rousk et al. 2010
      pH_offset <- 7 - pH
      growth_decrease <- pH_offset * 0.18
      
      #ODEs
      dBdt  <- (1-growth_decrease)*(CUE * DECOMP.C) - DEATH.C - overflow.R
      dCdt  <- I + (DEATH.C - SORPTION.B.C) + DESORPTION.C + growth_decrease*(CUE*DECOMP.C) - DECOMP.C - SORPTION.P.C - exo.loss.C
      dMdt  <- SORPTION.C - DESORPTION.C
      dN1dt <- (I/IN) + (DEATH.N - SORPTION.B.N) + DESORPTION.N + growth_decrease*(NUE*DECOMP.N) - DECOMP.N - SORPTION.P.N - exo.loss.N
      dN2dt <- SORPTION.N - DESORPTION.N
      dN3dt <- (1-growth_decrease)*(NUE * DECOMP.N) + immobilization - mineralization - DEATH.N
      dN4dt <- DECOMP.N*(1-NUE) + mineralization - immobilization - INORG.N.LOSS
      
      #respiration
      resp  <- (1-CUE)*DECOMP.C + overflow.R
      
      #update pools
      C  <- C  + dCdt 
      M  <- M  + dMdt
      B  <- B  + dBdt
      N1 <- N1 + dN1dt 
      N2 <- N2 + dN2dt
      N3 <- N3 + dN3dt
      N4 <- N4 + dN4dt
      pH <- -log10(protons)
      
      #solver
      list(c(dCdt,dMdt,dBdt,dN1dt,dN2dt,dN3dt,dN4dt,N_turnover,proton_change))
      })
}
