###Sensitivity experiments###
####Examine sensitivity of model to range of parameter values###
rm(list=ls())
require(rootSolve)

#grab the model as a function stored in another R script. 
source('model.r')

#grab parameters
source('parameters.r')

#Specify vector of input C:N values, and 'parameter of interest' values (for assessing model sensitivity to that parameter)
cn.range <- c(30,60)
param.range <- c(0.00001,0.001)

out.list <- list()

#create outermost loop for parameter-of-interest values
for(j in 1:length(param.range)){
  h5 <- param.range[j] ###alter which parameter is modified
  out <- matrix(rep(0,length(cn.range)*7),nrow=length(cn.range),dimnames=list(NULL,c('C','M','B','N1','N2','N3','N4')))

  for(i in 1:length(cn.range)){
    #examine model performance across range of C:N and intermediate soil sorptive capacity. If assessing sensitivity to v2 parameter, set 'v2=v2' 
    pars    <- c(IN=cn.range[i], BN=BN,                    #parameters
                 I=I,
                 v1=v1, v2=1.25, 
                 k1=k1, k2=k2, k3=k3, 
                 h1=h1, h2=h2, h3=h3, h4=h4, h5=h5)
    y       <- y                                           #state variable starting values for numeric solver from parameters file. 
    ST      <- runsteady(y=y,func=CN_NUE_model,parms=pars) #solve model
    out[i,] <- c(ST$y)                                     #save output of each iteration of i.
  }#end model loop.
    out <- data.frame(out)
    out$CN  <- cn.range
    out$tot <- out$C + out$M + out$B         #calculate total C pool
  #save each level of CN within a parameter-of-interest level as a matrix within a list. 
    out.list[[j]] <- out  
}   

##calculate model sensitivity to parameter of interest, under C limitation and N limitation
low <- out.list[[1]]
high <- out.list[[2]]


param.range <- as.data.frame(param.range)
C.limitation.sensitivity <- ((log10(high[1,9])) - (log10(low[1,9])))/((log10(param.range[2,])) - (log10(param.range[1,])))
N.limitation.sensitivity <- ((log10(high[2,9])) - (log10(low[2,9])))/((log10(param.range[2,])) - (log10(param.range[1,])))

C.limitation.sensitivity
N.limitation.sensitivity