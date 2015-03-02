setwd(paste(Sys.getenv("WORK"), "/ALPHA-prevalence/", sep=""))
## setwd("~/Documents/Research/ALPHA/mortality/cluster/")

library(rstan)

lapply(paste("karongam-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("karongaf-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("kisesam-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("kisesaf-fit", c(1,3), ".RData", sep=""), load, globalenv())
lapply(paste("kisumum-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("kisumuf-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("manicm-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("manicf-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("masakam-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("masakaf-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("rakaim-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("rakaif-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("umkhanm-fit", 1:6, ".RData", sep=""), load, globalenv())
lapply(paste("umkhanf-fit", 1:6, ".RData", sep=""), load, globalenv())

rm(karongam.fit, karongaf.fit, kisesam.fit, kisesaf.fit, kisesam.fit, kisesaf.fit, manicm.fit, manicf.fit,
   masakam.fit, masakaf.fit, rakaim.fit, rakaif.fit, umkhanm.fit, umkhanf.fit)

karongam.fit <- sflist2stanfit(lapply(grep("karongam.fit", ls(), value=TRUE), get))
karongaf.fit <- sflist2stanfit(lapply(grep("karongaf.fit", ls(), value=TRUE), get))
kisesam.fit <- sflist2stanfit(lapply(grep("kisesam.fit", ls(), value=TRUE), get))
kisesaf.fit <- sflist2stanfit(lapply(grep("kisesaf.fit", ls(), value=TRUE), get))
kisumum.fit <- sflist2stanfit(lapply(grep("kisumum.fit", ls(), value=TRUE), get))
kisumuf.fit <- sflist2stanfit(lapply(grep("kisumuf.fit", ls(), value=TRUE), get))
manicm.fit <- sflist2stanfit(lapply(grep("manicm.fit", ls(), value=TRUE), get))
manicf.fit <- sflist2stanfit(lapply(grep("manicf.fit", ls(), value=TRUE), get))
masakam.fit <- sflist2stanfit(lapply(grep("masakam.fit", ls(), value=TRUE), get))
masakaf.fit <- sflist2stanfit(lapply(grep("masakaf.fit", ls(), value=TRUE), get))
rakaim.fit <- sflist2stanfit(lapply(grep("rakaim.fit", ls(), value=TRUE), get))
rakaif.fit <- sflist2stanfit(lapply(grep("rakaif.fit", ls(), value=TRUE), get))
umkhanm.fit <- sflist2stanfit(lapply(grep("umkhanm.fit", ls(), value=TRUE), get))
umkhanf.fit <- sflist2stanfit(lapply(grep("umkhanf.fit", ls(), value=TRUE), get))

save(karongam.fit, karongaf.fit, kisesam.fit, kisesaf.fit, kisumum.fit, kisumuf.fit,
     manicm.fit, manicf.fit, masakam.fit, masakaf.fit,
     rakaim.fit, rakaif.fit, umkhanm.fit, umkhanf.fit,
     file="alpha-prevalence-fits_2015-02-28.RData")
