#CN by clay experiment
#model parameters for zero clay.

I    <- 0.03     #C input rate                            (mgC time-1)
CUE  <- 0.4      #C use efficiency                        (unitless)
NUE  <- .8 		   #N use efficiency					              (unitless)
v1   <- 0.0299   #biomass-specific decay multiplier       (mgC mgB-1 time-1) 
v2   <- 2        #Vmax of sorption                        (mgC time-1)
k1   <- 30       #half saturation of decomp               (mg)
k2   <- 1        #half saturation of sorption             (mg)
k3   <- 0.01     #half saturation of inorganic N uptake   (mg)
h1   <- 0.002    #biomass turnover rate multiplier        (1/time)
h2   <- 0.00075  #C-specific desorption rate              (1/time)
h3   <- 0.00001  #fraction of POM that potentially sorbs  (1/time)
h4   <- 0.2      #inorganic N loss rate                   (1/time)
h5   <- 0.00012  #exogenous losses of POM                 (1/time)

#pH parameters
   pH_opt <- 7
nitr.acid <- 1e-3      
acid.loss <- 1.5e-2
   ph.mod <- 0.2       #this reduces biomass growth by 20% per unit pH decline relative to 7
  pH_mode <- F         #pH mode is default off. can be turned on once parameters are sourced.

#parameters: N cycling
IN <- 20         #C:N ratio of inputs                                     (constant)
BN <-  7         #C:N ratio of the microbial biomass                      (constant)
CN <- 50         #C:N ratio of initial particulate organic matter         (will change)
MN <- 25         #C:N ratio of initial mineral-associated organic matter  (will change)
fert <- 0        #default inorganic N fertilization rate is 0.

#initial values of state varaibles
pH <- 7
protons <- 10^(-pH)
C <- 100      #C pool in POM            (mg / g)
M <-  20      #C pool in MAOM           (mg / g)
B <- C*0.1    #microbial biomass C 			(mg / g)
N1 <- C/CN    #initial N pool in POM		(mg / g)
N2 <- M/MN	  #N pool in MAOM 			  	(mg / g)
N3 <- B/BN    #microbial biomass N			(mg / g)
N4 <- 0.001	  #inorganic N pool size 		(mg / g)

#vector of state variables to track or solve model for.
y       <- c(C=C,M=M,B=B,N1=N1,N2=N2,N3=N3,N4=N4,protons=protons)