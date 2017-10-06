#Plotting N fertilization by C:N of inputs by soil sorptive capcity with pH effects on.
#Clear environment, load packages.
rm(list=ls())
library(wesanderson)

#load experiment data.
d <- readRDS('experiment_output/N_pH_feedbacks.experiment.rds')

#set output path for saving figure
output.path <- "figures/Experiment 4. N by pH figure.png"

#we have 2 levels of soil clay content
#we have 3 levels of CN: 30, 60, 80

pl1.30 <- d[[1]][[1]]
pl1.60 <- d[[1]][[2]]
pl1.80 <- d[[1]][[3]]
pl2.30 <- d[[2]][[1]]
pl2.60 <- d[[2]][[2]]
pl2.80 <- d[[2]][[3]]


#choose some colors
cols <- wes_palette("Zissou", 5)

#set ylimits
limy <- c(0,160)

#outer label cex value
o.lab = 1.5

#save dimensions, destination
png(filename=output.path, width=8.1, height=7.2, units='in',res=300)

#setup panels
par(mfrow=c(2,3), oma=c(5,8,5,5), mar=c(0,0,0,0),las=1)

#make the first panel that all other panels will be modeled off of

dd <- pl1.30
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
legend(30,150,c('total C','POM','MAOM','microbial','begin +N'), 
       lwd=c(2,2,2,2,1), col=c('black',cols[1],cols[3],cols[5],'black'), lty=c(1,1,1,1,3), 
       bty='n', y.intersp = 1, x.intersp = 0.75, cex=1.5, seg.len=1.5)



dd <- pl1.60
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
mtext('C:N = 60', side = 3, line = .5,         cex = o.lab)
mtext('mild N limitation', side = 3, line = 3, cex = 1.2)

mtext('b.', side =3, adj = 0.975, line = -1.25)


dd <- pl1.80
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
mtext('C:N = 80', side = 3, line = .5,          cex = o.lab)
mtext('severe N limitation', side = 3, line = 3, cex = 1.2)
mtext('low sorption',side = 4, line =  1, las = 0, cex = 1.3)
mtext('c.', side =3, adj = 0.975, line = -1.25)


dd <- pl2.30
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
Axis(side=1, labels=T)
Axis(side=2)
mtext('d.', side =3, adj = 0.975, line = -1.25)


dd <- pl2.60
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
Axis(side=1, labels=T)
mtext('e.', side =3, adj = 0.975, line = -1.25)


dd <- pl2.80
plot(tot~year, data = dd, cex = 0, xlab=NA, ylab=NA, xaxt='n', yaxt='n', ylim=limy)
lines(smooth.spline(dd$tot ~ dd$year), lwd=2)
lines(smooth.spline(dd$C   ~ dd$year), lwd=2, col=cols[1])
lines(smooth.spline(dd$M   ~ dd$year), lwd=2, col=cols[3])
lines(smooth.spline(dd$B   ~ dd$year), lwd=2, col=cols[5])
abline(v=25.21, lty = 3)
Axis(side=1, labels=T)
mtext('high sorption',side = 4, line = 1, las = 0, cex = 1.3)
mtext('f.', side =3, adj = 0.975, line = -1.25)


#outer labels
mtext('time (years)'                         , side=1, out=T, cex=o.lab, line = 3)
mtext(expression(paste('mg C (g soil)'^'-1')), side=2, out=T, cex=o.lab, line = 3, las = 0)

#end plot
dev.off()