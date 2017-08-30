#fertilize with a 0.05mg inorganic N dose, equivalent to a 1 time application of 50 kg N / ha / yr

#clear R environment, load packages

rm(list=ls())




#specify model, parameter, and output paths

model.path <- 'model_pH.r'

parameters.path <- 'parameters_pH.r'

output.path <- 'experiment_output/N_pH_feedbacks.experiment.rds'

#Specify vector of input C:N values, and pl (proton loss, or soil buffering) values
cn.range <- c(30,60,80)
pl.range <- c(5e-5,5e-30)

#number of days to step the dynamic simulation through time.
t <- 400000

#addition rate of 0.05 mg / g soil / yr in our model, ~50kg N / ha / yr.
#This assumes an effect soil depth of 10cm, and a bulk density of 1 g / cm3.
fert.level <- 0.05 / 365
fert.day   <- 300000      #which day to begin fertilization.


#create empty meta.list for storing multiple lists of matries
meta.list <- list()

#create outermost loop for 3 levels of proton loss
for(k in 1:length(pl.range)){
  #grab parameters and starting values
  source('parameters_pH.r')
  #set proton loss rate (simulating soil buffering capacity)
  proton_loss <- pl.range[k]
  #create list for saving outputs across 3 levels of CN.
  out.list <- list()
  
    for(j in 1:length(cn.range)){
    #choose outputs you want to follow. 
    outputs <- c('C','M','B','R','N1','N2','N3','N4','pH','day')
    out <- matrix(rep(0,t*length(outputs)),nrow=t,dimnames=list(NULL,outputs))
    
    #set current loop input CN value
    IN <- cn.range[j]
    day <- 0
    
    for(i in 1:t){
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
      
      #Drop in N fertilization at t=80000 (index 4000 w/ thinning at 20)
      if(day > fert.day){
        dN4dt = dN4dt + fert.level
      }
      
      #respiration
      resp  <- (1-CUE)*DECOMP.C + overflow.R
      
      #pH
      pH_mod <- pH/pH_opt 
	  N_turnover <- dN4dt
	  Nitr <- pH_mod*N_turnover #nitrification is proportional to pH (further from optimum, the slower it goes) and rate of N turnover sensu Parton et al 1996
	 proton_change <- Nitr*1e-3
	 protons <- protons + proton_change - proton_loss
	 if (protons < 1e-07) {protons = 1e-07} else {protons = protons} #pH cannot rise above 7
	 pH <- -log10(protons)
	 if (pH < 4) {pH = 4} else {pH=pH} #pH cannot drop beneath a certain threshold due to soil buffering capacity
      
      #pH effects on biomass - growth is reduced by 18% for each unit change in pH sensu Rousk et al. 2010
      pH_offset <- 7 - pH
      growth_decrease <- pH_offset * 0.18
      dBdt  <- (1-growth_decrease)*(CUE * DECOMP.C) - DEATH.C - overflow.R
      

      #update pools
      C  <- C  + dCdt
      M  <- M  + dMdt
      B  <- B  + dBdt
      N1 <- N1 + dN1dt
      N2 <- N2 + dN2dt
      N3 <- N3 + dN3dt
      N4 <- N4 + dN4dt
      R  <- resp
      day <- day + 1
      
      #save outputs
      out[i,] <- c(C,M,B,R,N1,N2,N3,N4,pH,day)
    } #end model loop.
    
    #thin output matrix for plotting. Grabbing days just before fertilization for plotting.
    out <- subset(out, out[,10] > 290800)
    out <- subset(out, out[,10] < 350000)   #most simulation responses stable by this point.
    out <- out[seq(1,nrow(out),by = 100),] #only need every 100th day for plotting.
    out <- data.frame(out)
    out$tot <- out$C + out$M + out$B         #calculate total C pool
    out$year <- (out$day - min(out$day))/365 #calcualte time in years, setting first day to zero. 
    #save each level of CN within a soil buffering level as a matrix within a list. 
    out.list[[j]] <- out
    
  } #end CN.range loop.
  
  #Name each list of C:N ratios within a level of soil buffering
  meta.list[[k]] <- out.list
  
} #end v.range loop.

#3 levels of input CN by 2 levels of soil buffering capacity are now stored in a nested list (meta.list).
#save output. 
saveRDS(meta.list, output.path)
