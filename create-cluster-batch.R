create.cluster.scripts <- function(label, site, sex, chains, iter=1000, min.time="NULL", max.time="NULL"){

  for(i in chains){
    fileConn <- file(paste("cluster/", label, "-fit", i, ".R", sep=""))
    writeLines(c("setwd(paste(Sys.getenv('HOME'), '/ALPHA-prevalence/', sep=''))",
                 "source('prevalence-functions.R')",
                 paste(label, ".dat <- prepare.data(subset(hiv, site == '", site, "' & sex == '", sex, "'), min.time=", min.time, ", max.time=", max.time, ")", sep=""),
                 paste(label, ".fit", i, " <- stan(model_code = stan.code, data = ", label, ".dat, iter = ", iter, ", chains=1, chain_id=", i, ", refresh=10)", sep=""),
                 paste("save(", label, ".fit", i, ", file=paste(Sys.getenv('WORK'), '/ALPHA-prevalence/", label, "-fit", i, ".RData', sep=''))", sep="")), fileConn)
    close(fileConn)

    fileConn <- file(paste("cluster/", label, "-fit", i, ".pbs", sep=""))
    writeLines(c("#!/bin/sh", 
                 "#PBS -l walltime=48:00:00",
                 "#PBS -l select=01:ncpus=1:mem=4gb",
                 "#PBS -j oe",
                 "",
                 "module load R/3.1.0",
                 "module load intel-suite",
                 "",
                 paste("R CMD BATCH --no-restore --no-save $HOME/ALPHA-prevalence/", label, "-fit", i, ".R $WORK/ALPHA-prevalence/", label, "-fit", i, ".Rout", sep=""),
                 "",
                 "qstat -f $PBS_JOBID"), fileConn)
    close(fileConn)
  }

  return(NULL)
}

create.cluster.scripts("karongam", "Karonga", "Male", 1:6)
create.cluster.scripts("karongaf", "Karonga", "Female", 1:6)

create.cluster.scripts("kisesam", "Kisesa", "Male", 1:6)
create.cluster.scripts("kisesaf", "Kisesa", "Female", 1:6)

create.cluster.scripts("kisumum", "Kisumu", "Male", 1:6)
create.cluster.scripts("kisumuf", "Kisumu", "Female", 1:6, max.time=2013.6)

create.cluster.scripts("manicm", "Manicaland", "Male", 1:6)
create.cluster.scripts("manicf", "Manicaland", "Male", 1:6)

create.cluster.scripts("masakam", "Masaka", "Male", 1:6)
create.cluster.scripts("masakaf", "Masaka", "Female", 1:6)

create.cluster.scripts("rakaim", "Rakai", "Male", 1:6)
create.cluster.scripts("rakaif", "Rakai", "Female", 1:6)

create.cluster.scripts("umkhanm", "uMkhanyakude", "Male", 1:6, min.time=2002)
create.cluster.scripts("umkhanf", "uMkhanyakude", "Female", 1:6, min.time=2002)
