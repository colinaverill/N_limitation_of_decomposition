#C-only model.
C_model<- function(t,y,pars){
  with(as.list(c(y,pars)),{
    dBdt <- CUE * v1 * B * C / (k1 + C) - h1*B^1.5
    dCdt <- I + (h1*B^1.5 - (h1*B^1.5  /(h3*C+h1*B^1.5))*v2*(h3*C+h1*B^1.5) / (k2 + h3*C+h1*B^1.5)) + h2*M - v1 * B * C / (k1 + C) - (h3*C/(h3*C+h1*B^1.5))*v2*(h3*C+h1*B^1.5) / (k2 + h3*C+h1*B^1.5)
    dMdt <- v2*(h3*C+h1*B^1.5) / (k2 + h3*C+h1*B^1.5) - h2*M
    list(c(dBdt,dCdt,dMdt))
  })
}
