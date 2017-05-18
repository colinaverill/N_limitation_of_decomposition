#CN by clay experiment
#model parameters for zero clay.

I    <- 0.015    #C input rate                            (mg time-1)     
CUE  <- 0.3      #carbon use efficiency                   (unitless)
NUE  <- .8 		   #N use efficiency					              (unitless)
v1   <- .017     #biomass-specific decay multiplier       (mgC mgB-1 time-1)
v2   <- 1        #Vmax of sorption                        (mgC time-1)
k1   <- 33       #half saturation of decomp               (mg)              
k2   <- 10       #half saturation of sorption             (mg)
k3   <- 0.01     #half saturation of inorganic N uptake   (mg)
h1   <- 0.003    #biomass turnover rate                   (1/time)    
h2   <- 0.001    #C-specific desorption rate              (1/time)
h3   <- 0.0002   #fraction of POM that potentially sorbs  (1/time)
h4   <- 0.2      #inorganic N loss rate                   (1/time)
h5   <- 0.00005  #exogenous losses of POM                 (1/time)

#parameters: N cycling
IN <- 20         #C:N ratio of inputs (constant).
BN <-  7         #C:N ratio of the microbial biomass (constant)
CN <- 50         #C:N ratio of initial particulate organic matter (will change)
MN <- 25         #C:N ratio of initial mineral-associated organic matter (will change)

#initial values of state varaibles
C <- 100      #C pool in POM            (mg / g)
M <-  20      #C pool in MAOM           (mg / g)
B <- C*0.1    #microbial biomass C 			(mg / g)
N1 <- C/CN    #initial N pool in POM		(mg / g)
N2 <- M/MN	  #N pool in MAOM 			  	(mg / g)
N3 <- B/BN    #microbial biomass N			(mg / g)
N4 <- 0.001	  #inorganic N pool size 		(mg / g)

y       <- c(C=C,M=M,B=B,N1=N1,N2=N2,N3=N3,N4=N4)
