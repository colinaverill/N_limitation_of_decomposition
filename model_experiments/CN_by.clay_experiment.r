#CN model experiment
#Model experiment 4.05.17.
#Experiment 1: vary CN of inputs from from 10-60
#track effects on all 3 C pools and 4 N pools at zero, lo and hi clay.

#clear R environment, load packages
rm(list=ls())
require(rootSolve)
require(wesanderson)

#grab the model as a function stored in another R script. 
source('model_functions/CN_model.r')

#Specify vector of input C:N values.
cn.range <- seq(20,100, by = 1)
#Specify vector of clay values.
v2.range <- c(0,1,5)

#create a matrix to store outputs of loop.
out <- matrix(rep(0,length(cn.range)*7),nrow=length(cn.range),dimnames=list(NULL,c('C','M','B','N1','N2','N3','N4')))

#Experiment 1a: no clay
I    <- 1.0      #C input rate                           (mg time-1)
CUE  <- 0.3      #carbon use efficiency                  (unitless)
NUE  <- 1 		   #N use efficiency					             (unitless)
v1   <- .4       #biomass-specific decay multiplier      (mg time-1)
v2   <- 1        #Vmax of sorption                       (mg time-1)
k1   <- 200      #half saturation of decomp              (mg)
k2   <- 10       #half saturation of sorption            (mg)
k3   <- 0.001    #half saturation of inorganic N uptake  (mg)
h1   <- 0.01     #biomass turnover rate                  (1/time)
h2   <- 0.002    #C-specific desorption rate             (1/time)
h3   <- 0.001    #fraction of POM that potentially sorbs (1/time)
h4   <- 0.9      #inorganic N loss rate
h5   <- 0.01     #exogenous losses of POM                (1/time)

#parameters: N cycling
#IN- input C:N, is in the for loop. 
CN <- 50 #C:N ratio of initial particulate organic matter (will change)
MN <- 25 #C:N ratio of initial mineral-associated organic matter (will change)
BN <-  7 #C:N ratio of the microbial biomass (constant)

#initial pool sizes
C <- 100      #C pool in POM            (mg / g)
M <-  20      #C pool in MAOM           (mg / g)
B <- C*0.1    #microbial biomass C 			(mg / g)
R <- 0        # respired C				    	(mg / g)
N1 <- C/CN    #initial N pool in POM		(mg / g)
N2 <- M/MN	  #N pool in MAOM 			  	(mg / g)
N3 <- B/BN    #microbial biomass N			(mg / g)
N4 <- 0.001	  #inorganic N pool size 		(mg / g)


#no clay
for(i in 1:length(cn.range)){
  #set parameters,across respective ranges.
  pars    <- c(IN=cn.range[i], BN=BN,                #parameters
               I=I,
               v1=v1, v2=v2.range[1], 
               k1=k1, k2=k2, k3=k3, 
               h1=h1, h2=h2, h3=h3, h4=h4, h5=h5
               )
  y       <- c(C=C,M=M,B=B,N1=N1,N2=N2,N3=N3,N4=N4)  #state variable starting values for numeric solver.
  ST      <- runsteady(y=y,func=CN_model,parms=pars) #solve model
  out[i,] <- c(ST$y)                                 #save output of each iteration of i.
}

out    <- data.frame(out)
out$CN <- cn.range
out$tot <- out$C + out$M + out$B
saveRDS(out,'experiment_output/no_clay.CN_model.rds')

#lo clay
#create a matrix to store outputs of loop.
out <- matrix(rep(0,length(cn.range)*7),nrow=length(cn.range),dimnames=list(NULL,c('C','M','B','N1','N2','N3','N4')))

for(i in 1:length(cn.range)){
  #set parameters,across respective ranges.
  pars    <- c(IN=cn.range[i], BN=BN,                #parameters
               I=I,
               v1=v1, v2=v2.range[2], 
               k1=k1, k2=k2, k3=k3, 
               h1=h1, h2=h2, h3=h3, h4=h4, h5=h5
  )
  y       <- c(C=C,M=M,B=B,N1=N1,N2=N2,N3=N3,N4=N4)  #state variable starting values for numeric solver.
  ST      <- runsteady(y=y,func=CN_model,parms=pars) #solve model
  out[i,] <- c(ST$y)                                 #save output of each iteration of i.
}

out    <- data.frame(out)
out$CN <- cn.range
out$tot <- out$C + out$M + out$B
saveRDS(out,'experiment_output/lo_clay.CN_model.rds')

#hi clay
#create a matrix to store outputs of loop.
out <- matrix(rep(0,length(cn.range)*7),nrow=length(cn.range),dimnames=list(NULL,c('C','M','B','N1','N2','N3','N4')))

for(i in 1:length(cn.range)){
  #set parameters,across respective ranges.
  pars    <- c(IN=cn.range[i], BN=BN,                #parameters
               I=I,
               v1=v1, v2=v2.range[3], 
               k1=k1, k2=k2, k3=k3, 
               h1=h1, h2=h2, h3=h3, h4=h4, h5=h5
  )
  y       <- c(C=C,M=M,B=B,N1=N1,N2=N2,N3=N3,N4=N4)  #state variable starting values for numeric solver.
  ST      <- runsteady(y=y,func=CN_model,parms=pars) #solve model
  out[i,] <- c(ST$y)                                 #save output of each iteration of i.
}

out    <- data.frame(out)
out$CN <- cn.range
out$tot <- out$C + out$M + out$B
saveRDS(out,'experiment_output/hi_clay.CN_model.rds')