#####################
####  Load data  ####
#####################

setwd('~/Documents/Research/ALPHA/prevalence/')
source('prevalence-functions.R')

prepare.output <- function(sites, sexes, min.time=NULL, max.time=NULL, nk.time=NULL, fixcoef.time=NULL,
                           min.age=15.0, max.age=55.0, nk.age=11L, fixcoef.age=4L, dt=0.1){

  dat <- prepare.data(subset(hiv, site %in% sites & sex %in% sexes), min.time, max.time, nk.time, fixcoef.time, min.age, max.age, nk.age, fixcoef.age)

  pred.time <- seq(dt*round(dat$k_time[4]/ dt), dt*round(tail(dat$k_time, 4)[1] / dt), dt)
  pred.age <- seq(dt*round(min.age / dt), dt*round(max.age / dt), dt)

  dat$dt <- dt
  
  dat$pred_time <- pred.time
  dat$X_time <- splineDesign(dat$k_time, pred.time, outer.ok=TRUE)
  dat$Xmid_time <- splineDesign(dat$k_time, pred.time[-1]-dt/2)

  dat$pred_age <- pred.age
  dat$X_age <- splineDesign(dat$k_age, pred.age, outer.ok=TRUE)
  dat$Xmid_age <- splineDesign(dat$k_age, pred.age[-1]-dt/2)

  return(dat)
}

karongam.dat <- prepare.output('Karonga', 'Male')
karongaf.dat <- prepare.output('Karonga', 'Female')
kisesam.dat <- prepare.output('Kisesa', 'Male')
kisesaf.dat <- prepare.output('Kisesa', 'Female')
kisumum.dat <- prepare.output('Kisumu', 'Male')
kisumuf.dat <- prepare.output('Kisumu', 'Female', max.time=2013.6)
manicm.dat <- prepare.output('Manicaland', 'Male')
manicf.dat <- prepare.output('Manicaland', 'Female')
masakam.dat <- prepare.output('Masaka', 'Male')
masakaf.dat <- prepare.output('Masaka', 'Female')
rakaim.dat <- prepare.output('Rakai', 'Male')
rakaif.dat <- prepare.output('Rakai', 'Female')
umkhanm.dat <- prepare.output('uMkhanyakude', 'Male', min.time=2002.0)
umkhanf.dat <- prepare.output('uMkhanyakude', 'Female', min.time=2002.0)



##############################
####  Analysis functions  ####
##############################

load("alpha-prevalence-fits_2015-02-28.RData")

karongam.samp <- extract(karongam.fit)
karongaf.samp <- extract(karongaf.fit)
kisesam.samp <- extract(kisesam.fit)
kisesaf.samp <- extract(kisesaf.fit)
kisumum.samp <- extract(kisumum.fit)
kisumuf.samp <- extract(kisumuf.fit)
manicm.samp <- extract(manicm.fit)
manicf.samp <- extract(manicf.fit)
masakam.samp <- extract(masakam.fit)
masakaf.samp <- extract(masakaf.fit)
rakaim.samp <- extract(rakaim.fit)
rakaif.samp <- extract(rakaif.fit)
umkhanm.samp <- extract(umkhanm.fit)
umkhanf.samp <- extract(umkhanf.fit)


invlogit <- function(x) return(exp(x)/(1+exp(x)))

calc.prev.age <- function(samp, dat, ages=c(20, 30, 40, 50)){
  prev.age <- array(NA, c(length(ages), length(dat$pred_time), length(samp$lp__)))
  for(idx in seq_along(ages))
    prev.age[idx,,] <- invlogit(dat$X_time %*% apply(samp$coef_time_age, 1, "%*%", dat$X_age[dat$pred_age==ages[idx],]))
  dimnames(prev.age) <- list(ages, dat$pred_time, NULL)
  return(prev.age)
}

karongam.prev.age <- calc.prev.age(karongam.samp, karongam.dat)
karongaf.prev.age <- calc.prev.age(karongaf.samp, karongaf.dat)
kisesam.prev.age <- calc.prev.age(kisesam.samp, kisesam.dat)
kisesaf.prev.age <- calc.prev.age(kisesaf.samp, kisesaf.dat)
kisumum.prev.age <- calc.prev.age(kisumum.samp, kisumum.dat)
kisumuf.prev.age <- calc.prev.age(kisumuf.samp, kisumuf.dat)
manicm.prev.age <- calc.prev.age(manicm.samp, manicm.dat)
manicf.prev.age <- calc.prev.age(manicf.samp, manicf.dat)
masakam.prev.age <- calc.prev.age(masakam.samp, masakam.dat)
masakaf.prev.age <- calc.prev.age(masakaf.samp, masakaf.dat)
rakaim.prev.age <- calc.prev.age(rakaim.samp, rakaim.dat)
rakaif.prev.age <- calc.prev.age(rakaif.samp, rakaif.dat)
umkhanm.prev.age <- calc.prev.age(umkhanm.samp, umkhanm.dat)
umkhanf.prev.age <- calc.prev.age(umkhanf.samp, umkhanf.dat)

calc.prev.time <- function(samp, dat, times=c(1999, 2003, 2007, 2011, 2014)){
  prev.time <- array(NA, c(length(times), length(dat$pred_age), length(samp$lp__)))
  for(idx in seq_along(times))
    prev.time[idx,,] <- invlogit(dat$X_age %*% apply(samp$coef_time_age, 1, function(Bmat) dat$X_time[dat$pred_time==times[idx],] %*% Bmat))
  dimnames(prev.time) <- list(times, dat$pred_age, NULL)
  return(prev.time)
}

karongam.prev.time <- calc.prev.time(karongam.samp, karongam.dat, c(2007, 2011, 2014))
karongaf.prev.time <- calc.prev.time(karongaf.samp, karongaf.dat, c(2007, 2011, 2014))
kisesam.prev.time <- calc.prev.time(kisesam.samp, kisesam.dat, c(1995, 1999, 2003, 2007, 2011, 2014))
kisesaf.prev.time <- calc.prev.time(kisesaf.samp, kisesaf.dat, c(1995, 1999, 2003, 2007, 2011, 2014))
kisumum.prev.time <- calc.prev.time(kisumum.samp, kisumum.dat, c(2011, 2014))
kisumuf.prev.time <- calc.prev.time(kisumuf.samp, kisumuf.dat, c(2011, 2014))
manicm.prev.time <- calc.prev.time(manicm.samp, manicm.dat, c(1999, 2003, 2007, 2011))
manicf.prev.time <- calc.prev.time(manicf.samp, manicf.dat, c( 1999, 2003, 2007, 2011))
masakam.prev.time <- calc.prev.time(masakam.samp, masakam.dat, c(1991, 1995, 1999, 2003, 2007, 2011, 2014))
masakaf.prev.time <- calc.prev.time(masakaf.samp, masakaf.dat, c(1991, 1995, 1999, 2003, 2007, 2011, 2014))
rakaim.prev.time <- calc.prev.time(rakaim.samp, rakaim.dat, c(1999, 2003, 2007, 2011, 2014))
rakaif.prev.time <- calc.prev.time(rakaif.samp, rakaif.dat, c(1999, 2003, 2007, 2011, 2014))
umkhanm.prev.time <- calc.prev.time(umkhanm.samp, umkhanm.dat, c(2003, 2007, 2011, 2014))
umkhanf.prev.time <- calc.prev.time(umkhanf.samp, umkhanf.dat, c(2003, 2007, 2011, 2014))

save(list=c(grep(".prev.age", ls(), value=TRUE), grep(".prev.time", ls(), value=TRUE)),
     file="ALPHA-prev-posterior_2015-03-02.RData")


########################
####  Plot results  ####
########################

load("ALPHA-prev-posterior_2015-03-02.RData")

library(adegenet)
library(RColorBrewer)

cred.region <- function(x, y, ...)
  polygon(c(x, rev(x)), c(y[1,], rev(y[2,])), ...)

## prevalence

plot.fit <- function(prev.age, prev.time, name, prev.range=c(0, 0.50), xlim.time=c(1992, 2015)){
  par(tcl=-0.25, mgp=c(2, 0.5, 0), mar=c(3.1, 3.1, 2.6, 0.6))
  ##
  plot(NA, NA, xlim=xlim.time, ylim=prev.range, xlab="Year", ylab="HIV prevalence", main=name)
  for(idx in 1:4){
    col <- brewer.pal(4, "Dark2")[idx]
    cred.region(as.numeric(dimnames(prev.age)[[2]]), apply(prev.age[idx,,], 1, quantile, c(0.1, 0.9)), col=transp(col, 0.2), border=NA)
    ## cred.region(as.numeric(dimnames(prev.age)[[2]]), apply(prev.age[idx,,], 1, quantile, c(0.025, 0.975)), col=transp(col, 0.2), border=NA)
    lines(as.numeric(dimnames(prev.age)[[2]]), apply(prev.age[idx,,], 1, mean), lty=1, lwd=2, col=col)
    lines(as.numeric(dimnames(prev.age)[[2]]), apply(prev.age[idx,,], 1, median), lty=2, lwd=2, col=col)
  }
  legend("topleft", paste("age", rev(dimnames(prev.age)[[1]])), lwd=2, col=rev(brewer.pal(4, "Dark2")))
  ##
  plot(NA, NA, xlim=c(15, 55), ylim=prev.range, xlab="Age", ylab="HIV prevalence", main=name)
  for(idx in 1:dim(prev.time)[1]){
    col <- brewer.pal(dim(prev.time)[1], "Dark2")[idx]
    cred.region(as.numeric(dimnames(prev.time)[[2]]), apply(prev.time[idx,,], 1, quantile, c(0.1, 0.9)), col=transp(col, 0.2), border=NA)
    ## cred.region(as.numeric(dimnames(prev.time)[[2]]), apply(prev.time[idx,,], 1, quantile, c(0.025, 0.975)), col=transp(col, 0.2), border=NA)
    lines(as.numeric(dimnames(prev.time)[[2]]), apply(prev.time[idx,,], 1, mean), lty=1, lwd=2, col=col)
    lines(as.numeric(dimnames(prev.time)[[2]]), apply(prev.time[idx,,], 1, median), lty=2, lwd=2, col=col)
  }
  legend("topleft", paste(dimnames(prev.time)[[1]]), lwd=2, col=brewer.pal(dim(prev.time)[1], "Dark2"))
  ##
  return(NULL)
}


##
X11(h=6, w=6, pointsize=8)


pdf("ALPHA-prevalence-fits_2015-03-02.pdf", h=6, w=6, pointsize=8)
##
par(mfrow=c(2,2), cex=1)
##
plot.fit(karongam.prev.age, karongam.prev.time, "Karonga men")
plot.fit(karongaf.prev.age, karongaf.prev.time, "Karonga women")
plot.fit(kisesam.prev.age, kisesam.prev.time, "Kisesa men")
plot.fit(kisesaf.prev.age, kisesaf.prev.time, "Kisesa women")
plot.fit(kisumum.prev.age, kisumum.prev.time, "Kisumu men")
plot.fit(kisumuf.prev.age, kisumuf.prev.time, "Kisumu women")
plot.fit(manicm.prev.age, manicm.prev.time, "Manicaland men")
plot.fit(manicf.prev.age, manicf.prev.time, "Manicaland women")
plot.fit(masakam.prev.age, masakam.prev.time, "Masaka men")
plot.fit(masakaf.prev.age, masakaf.prev.time, "Masaka women")
plot.fit(rakaim.prev.age, rakaim.prev.time, "Rakai men")
plot.fit(rakaif.prev.age, rakaif.prev.time, "Rakai women")
plot.fit(umkhanm.prev.age, umkhanm.prev.time, "uMkhanyakude men", prev.range=c(0, 0.7))
plot.fit(umkhanf.prev.age, umkhanf.prev.time, "uMkhanyakude women", prev.range=c(0, 0.7))
##
dev.off()
