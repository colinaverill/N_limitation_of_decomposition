#derive analytical solution of alpha parameter.
m.CN / CUE = x*DECOMP.C1 + (1-x)*DECOMP.C2 / ((1-x)*DECOMP.C2 / CN.org)
m.CN / CUE = (x*(v1*B*C1 / (k1 + C1)) + (1-x)*(v1*B*C2 / (k1 + C2))) / ((1-x)*(v1*B*C2 / (k1 + C2)) / CN.org)
m.CN   = a
CUE    = b
v1     = c
k1     = d
CN.org = e
B      = q
C1     = r
C2     = s
#re-express
a / b = (x*(c*q*r / (d + r)) + (1-x)*(c*q*s / (d + s))) / ((1-x)*(c*q*s / (d + s)) / e)

#Solve for x using wolframAlpha
x = (s*(d + r)*(a - b*e))/(a*d*s + a*r*s + b*d*e*r - b*d*e*s)
x = (C2*(k1 + C1)*(m.CN - CUE*CN.org))/(m.CN*k1*C2 + m.CN*C1*C2 + CUE*k1*CN.org*C1 - CUE*k1*CN.org*C2)