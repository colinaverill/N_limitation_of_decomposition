#C-fertilization figure.
#clear R environment, load packages
rm(list=ls())
library(wesanderson)

#load experiment data.
d <- readRDS('experiment_output/AV_C_fert.experiment.rds')

#set output path
output.path <- 'figures/AV_Experiment 3. C-fertilization figure.png'

#we have 3 levels of clay: 0, 1, 2
#we have 3 levels of CN  :30,60,80

v0.30 <- d[[1]][[1]]
v0.60 <- d[[1]][[2]]
v0.80 <- d[[1]][[3]]
v1.30 <- d[[2]][[1]]
v1.60 <- d[[2]][[2]]
v1.80 <- d[[2]][[3]]
v2.30 <- d[[3]][[1]]
v2.60 <- d[[3]][[2]]
v2.80 <- d[[3]][[3]]

#choose some colors
cols <- wes_palette("Zissou", 5)

#set ylimits
limy <- c(0,260)

#outer label cex value
o.lab = 1.5

#save dimensions, destination
png(filename=output.path,width=8.1,height=7.2,units='in',res=300)


#setup panels
par(mfrow=c(3,3), oma=c(5,8,5,5), mar=c(0,0,0,0),las=1)

#make the first panel that all other panels will be modeled off of
#v0.30 panel
dd <- v0.30
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
Axis(side=2)
mtext('C:N = 30'    , side = 3, line = .5, cex = o.lab)
mtext('C limitation', side = 3, line = 3,  cex = 1.2)
mtext('a.', side =3, adj = 0.975, line = -1.25)
legend(30,250,c('total C','POM','MAOM','microbial','begin +C'), 
       lwd=c(2,2,2,2,1), col=c('black',cols[1],cols[3],cols[5],'black'), lty=c(1,1,1,1,3), 
       bty='n', y.intersp = 1, x.intersp = 0.75, cex=1.5, seg.len=1.5)


#v0.60 panel
dd <- v0.60
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
mtext('C:N = 60', side = 3, line = .5,         cex = o.lab)
mtext('mild N limitation', side = 3, line = 3, cex = 1.2)

mtext('b.', side =3, adj = 0.975, line = -1.25)

#v0.80 panel
dd <- v0.80
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
mtext('C:N = 80', side = 3, line = .5,          cex = o.lab)
mtext('severe N limitation', side = 3, line = 3, cex = 1.2)
mtext('no sorption',side = 4, line =  1, las = 0, cex = o.lab)
mtext('c.', side =3, adj = 0.975, line = -1.25)

#v1.30 panel
dd <- v1.30
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3], lty=2) #dashed line because POM and MAOM lines on top of eachother.
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
Axis(side=2)
mtext('d.', side =3, adj = 0.975, line = -1.25)

#v1.60 panel
dd <- v1.60
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
mtext('e.', side =3, adj = 0.975, line = -1.25)

#v1.80 panel
dd <- v1.80
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
mtext('low sorption',side = 4, line = 1, las = 0, cex = o.lab)
mtext('f.', side =3, adj = 0.975, line = -1.25)

#v2.30 panel
dd <- v2.30
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
Axis(side=1, labels=T)
Axis(side=2)
mtext('g.', side =3, adj = 0.975, line = -1.25)

#v2.60 panel
dd <- v2.60
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
Axis(side=1, labels=T)
mtext('h.', side =3, adj = 0.975, line = -1.25)

#v2.80 panel
dd <- v2.80
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
Axis(side=1, labels=T)
mtext('i.', side =3, adj = 0.975, line = -1.25)
mtext('high sorption',side = 4, line = 1, las = 0, cex = o.lab)

#outer labels
mtext('time (years)'                         , side=1, out=T, cex=o.lab, line = 3)
mtext(expression(paste('mg C (g soil)'^'-1')), side=2, out=T, cex=o.lab, line = 3, las = 0)

#end plot
dev.off()