#C-only model experiment
#Model experiment 3.22.17.
#Experiment 1: vary CUE from 0.1 - 0.6, track effects on all 3 C pools and total C, at both high and low clay.
#Experiment 2: vary I   from 0.1 - 1.5, track effects on all 3 C pools and total C, at both high and low clay. (this has not been done yet)
#clear R environment, load packages
rm(list=ls())
require(rootSolve)
require(wesanderson)

#grab the model as a function stored in another R script. 
source('model_functions/C_model.r')

#Specify range of CUE values. Lower than CUE=0.13 we don't get a steady state. 
cue.range <- seq(0.13,0.6, by = 0.01)
#create a matrix to store outputs of loop.
out <- matrix(rep(0,length(cue.range)*4),nrow=length(cue.range),dimnames=list(NULL,c('B','C','M','tot')))

#Experiment 1a: no clay
#loop it (probably a way to do this with lapply...)
for(i in 1:length(cue.range)){
  #set a parameter range
    pars    <- c(CUE=cue.range[i], I=1, v1=8, v2=0, k1=800, k2=10, h1=0.11, h2=0.003, h3=0.001)
  y       <- c(B=10, C=100, M=20) #starting values for numeric solver.
  ST      <- stode(y=y,func=C_model,parms=pars,pos=T) #solve model
  out[i,] <- c(ST$y, sum(ST$y)) #save output of each iteration of i.
}
out     <- data.frame(out)
out$CUE <- cue.range
saveRDS(out,'experiment_output/no_clay.C_only.rds')


#test plot.
cols <- wes_palette("Zissou", 5)
plot(tot ~ CUE, data = out, ylim = c(0,100), cex = 0, ylab = 'mg C per g soil')
lines(smooth.spline(out$tot~out$CUE), lwd=2)
lines(smooth.spline(out$C  ~out$CUE), lwd=2, col=cols[1])
lines(smooth.spline(out$M  ~out$CUE), lwd=2, col=cols[3])
lines(smooth.spline(out$B  ~out$CUE), lwd=2, col=cols[5])
mtext('no clay', line = 0.2)

#legend
legend(.4, 95,c('total C','POM','MAOM','microbial'), lwd=2, col=c('black',cols[1],cols[3],cols[5]), bty='n', y.intersp = 1, x.intersp = 0.75, cex=1.2, seg.len=1)


#Experiment 1b: lo clay
#create a matrix to store outputs of loop.
out <- matrix(rep(0,length(cue.range)*4),nrow=length(cue.range),dimnames=list(NULL,c('B','C','M','tot')))

#loop it (probably a way to do this with lapply...)
for(i in 1:length(cue.range)){
  #set a parameter range
  pars    <- c(CUE=cue.range[i], I=1, v1=8, v2=0.1, k1=800, k2=10, h1=0.11, h2=0.003, h3=0.001)
  y       <- c(B=10, C=100, M=20) #starting values for numeric solver.
  ST      <- stode(y=y,func=C_model,parms=pars,pos=T) #solve model
  out[i,] <- c(ST$y, sum(ST$y)) #save output of each iteration of i.
}
out     <- data.frame(out)
out$CUE <- cue.range
saveRDS(out,'experiment_output/lo_clay.C_only.rds')


#Experiment 1c: hi clay
#create a matrix to store outputs of loop.
out <- matrix(rep(0,length(cue.range)*4),nrow=length(cue.range),dimnames=list(NULL,c('B','C','M','tot')))
#loop it (probably a way to do this with lapply...)
for(i in 1:length(cue.range)){
  #set a parameter range
  pars    <- c(CUE=cue.range[i], I=1, v1=8, v2=1, k1=800, k2=10, h1=0.11, h2=0.003, h3=0.001)
  y       <- c(B=10, C=100, M=20) #starting values for numeric solver.
  ST      <- stode(y=y,func=C_model,parms=pars,pos=T) #solve model
  out[i,] <- c(ST$y, sum(ST$y)) #save output of each iteration of i.
}
out     <- data.frame(out)
out$CUE <- cue.range
saveRDS(out,'experiment_output/hi_clay.C_only.rds')