#C-only model.
C_model <- function(t,y,pars){
  with(as.list(c(y,pars)),{
    dBdt <- CUE * v1 * B * C / (k1 + C) - h1*B
    dCdt <- I + (h1*B - (h1*B  /(h3*C+h1*B))*v2*(h3*C+h1*B) / (k2 + h3*C+h1*B)) + h2*M - v1 * B * C / (k1 + C) - (h3*C/(h3*C+h1*B))*v2*(h3*C+h1*B) / (k2 + h3*C+h1*B)
    dMdt <- v2*(h3*C+h1*B) / (k2 + h3*C+h1*B) - h2*M
    list(c(dBdt,dCdt,dMdt))
  })
}