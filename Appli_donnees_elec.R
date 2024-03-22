source('functions.R')
# Data
energy_data <- read.csv("energydata_complete.csv")
View(energy_data)

energy_data$date <- strptime(as.character(energy_data$date),format="%Y-%m-%d %H:%M:%S")
energy_data$date <- as.POSIXct(energy_data$date,tz = "UTC")

dim(energy_data)
names(energy_data)

n = 137
t = seq(from=0,by=24/144,length=144)

# From the temporal series energy_data[,2] to a matrix of functional data
energy_data_day = matrix(energy_data[1:19728,2],137,144,byrow=TRUE)

plot(range(t),range(energy_data_day),type='n',xlab='time(hour)',ylab='Electric consumption (kWh)')
for (i in 1:n){
  points(t,energy_data_day[i,],type='b',col=i,pch=16,cex=0.5)
}

# Kernel smoothing



h = c(0.01,1,10)
tt = seq(0,max(t),length=500)
tildeXwn1 = mNW(tt,t,energy_data_day[1,],h[1])
tildeXwn2 = mNW(tt,t,energy_data_day[1,],h[2])
tildeXwn3 = mNW(tt,t,energy_data_day[1,],h[3])

plot(t,energy_data_day[1,],type='b',xlab='time (hour)',ylab='Electric consumption (kWh)',pch=16,cex=0.5)
points(tt,tildeXwn1,type='l',col='red')
points(tt,tildeXwn2,type='l',col='chartreuse3')
points(tt,tildeXwn3,type='l',col='blue')
legend('top',c('observed data','h=0.01','h=1','h=10'),lty=rep(1,4),col=c('black','red','chartreuse3','blue'),pch=c(16,NA,NA,NA))


# Approximation by histograms

D = 1
tildeXhisto1 = approx.histo(t/max(t),energy_data_day[1,],D)
plot(t,energy_data_day[1,],type='b',xlab='time (hour)',ylab='Electric consumption (kWh)',pch=16,cex=0.5)
points(tt,tildeXhisto1,type='l',col='blue')
D=round(length(t)/10)
tildeXhisto2 = approx.histo(t/max(t),energy_data_day[1,],D)
points(tt,tildeXhisto2,type='l',col='chartreuse3')
D = length(t)
tildeXhisto3 = approx.histo(t/max(t),energy_data_day[1,],D)
points(tt,tildeXhisto3,type='l',col='red')
legend('top',c('observed data','D=p=144','D=14','D=1'),lty=rep(1,4),col=c('black','red','chartreuse3','blue'),pch=c(16,NA,NA,NA))

# Approximation by Fourier basis

D = 1
tildeXfourier1 = approx.fourier(t/max(t),energy_data_day[1,],D)
plot(t,energy_data_day[1,],type='b',xlab='time (hour)',ylab='Electric consumption (kWh)',pch=16,cex=0.5)
points(tt,tildeXfourier1,type='l',col='blue')
D=13
tildeXfourier2 = approx.fourier(t/max(t),energy_data_day[1,],D)
points(tt,tildeXfourier2,type='l',col='chartreuse3')
D = 51
tildeXfourier3 = approx.fourier(t/max(t),energy_data_day[1,],D)
points(tt,tildeXfourier3,type='l',col='red')
legend('top',c('observed data','D=51','D=13','D=1'),lty=rep(1,4),col=c('black','red','chartreuse3','blue'),pch=c(16,NA,NA,NA))

# Splines reconstruction
library(fda)

nbasis = 144
basisobj = create.bspline.basis(c(0,1),nbasis)
# Yfunc = fdata(Y,t[itd])
# Ysmooth.Fourier = fdata2fd(Yfunc,type.basis="fourier")
# plot(Ysmooth.Fourier,col='darkgreen',lwd=2)
# points(t[itd],Y,col='darkred',pch=4,lwd=2)
# points(t,X,type='l',col='blue')
# legend('bottomleft',c('true X','reconstruction','observations'),lty=c(1,1,0),col=c('blue','black','darkred'),pch=c(NA,NA,4),lwd=c(1,1,2))

Ysmooth.splines = smooth.basis(argvals=t/max(t),y=energy_data_day[1,],fdParobj=basisobj)
plot(Ysmooth.splines,col='brown',lwd=3,ylim=range(energy_data_day[1,]),xlab='time (hour)',ylab='Electric consumption (kWh)')
points(t/max(t),energy_data_day[1,],type='b',pch=16,cex=0.5)
legend('top',c('splines reconstruction','observations'),lty=c(1,0),col=c('black','black'),pch=c(NA,16),lwd=c(1,2),cex=0.5)

