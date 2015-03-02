library(reshape2)
library(splines)
library(mgcv)
library(rstan)
set_cppo("fast")  # for best running speed

## Load dataset
load("../prepare-data/ALPHA-hiv-data_v2015-02-23.RData")
hiv <- subset(hiv, type %in% c("survey") & !is.na(result) & !is.na(testdate))

prepare.data <- function(dat, min.time=NULL, max.time=NULL, nk.time=NULL, fixcoef.time=NULL,
                         min.age=15.0, max.age=55.0, nk.age=11L, fixcoef.age=4L){
  
  ## ################## ##
  ##  Prepare the data  ##
  ## ################## ##

  dat$time <- dat$testdate
  dat$age <- dat$time - dat$dob

  if(is.null(min.time)) min.time <- min(dat$time)
  if(is.null(max.time)) max.time <- max(dat$time)

  dat <- subset(dat, time >= min.time & time <= max.time & age >= min.age & age <= max.age)
  

  ## ###################### ##
  ##  Prepare spline model  ##
  ## ###################### ##

  if(is.null(nk.time)){ # place knots every year for span of data
    k.time <- seq(floor(min.time) - 3.0, ceiling(max.time) + 3.0, 1.0)
    nk.time <- length(k.time) - 4L
  } else {
    time.dur <- max.time - min.time
    k.time <- seq(min.time - 3*time.dur/(nk.time-3), max.time + 3*time.dur/(nk.time-3), time.dur/(nk.time-3))
  }

  age.dur <- max.age - min.age
  k.age <- seq(min.age - 3*age.dur/(nk.age-3), max.age + 3*age.dur/(nk.age-3), age.dur/(nk.age-3))

  P.time <- diff(diag(nk.time), diff=1)
  P.age <- diff(diag(nk.age), diff=1)

  if(is.null(fixcoef.time))
    fixcoef.time <- as.integer(nk.time / 2)


  ## ###################### ##
  ##  Create design matrix  ##
  ## ###################### ##

  Xmat <- tensor.prod.model.matrix(list(splineDesign(k.age, dat$age),
                                        splineDesign(k.time, dat$time)))

  ## ################### ##
  ##  List of Stan data  ##
  ## ################### ##

  stan.data <- list(nk_time=nk.time,
                    k_time=k.time,
                    P_time=P.time,
                    fixcoef_time_idx=fixcoef.time,
                    nk_age=nk.age,
                    k_age=k.age,
                    P_age=P.age,
                    fixcoef_age_idx=fixcoef.age,
                    N=nrow(dat),
                    x_time = dat$time,
                    x_age  = dat$age,
                    X_mat = Xmat,
                    hivpos = as.integer(dat$result))

  return(stan.data)
}

stan.code <- "
functions{
  vector inv_logit_vec(vector x){

    vector[rows(x)] val;
    
    for(i in 1:rows(x))
      val[i] <- inv_logit(x[i]);

    return val;
  }
}
data {
        // spline model setup
        int<lower=5> nk_time;
        matrix[nk_time-1, nk_time] P_time;
        int fixcoef_time_idx;

        int<lower=5> nk_age;
        matrix[nk_age-1, nk_age] P_age;
        int fixcoef_age_idx;

        // data
        int<lower=1> N;
        matrix[N, nk_time*nk_age] X_mat;
        int<lower=0> hivpos[N];
}
transformed data{
}
parameters {
        matrix[nk_time, nk_age] coef_time_age;
        real<lower=0> sigma2_time;
        real<lower=0> sigma2_age;
        real<lower=0> sigma2_time_age;
}
transformed parameters{
}
model {
        vector[nk_time] coef_time;
        row_vector[nk_age] coef_age;

        1/sigma2_time ~ gamma(1.0, 0.0005);
        1/sigma2_age ~ gamma(1.0, 0.0005);
        1/sigma2_time_age ~ gamma(1.0, 0.0005);

        coef_time <- col(coef_time_age, fixcoef_age_idx);
        coef_age <- row(coef_time_age, fixcoef_time_idx) - coef_time_age[fixcoef_time_idx, fixcoef_age_idx] ;

        P_time * coef_time ~ normal(0, sqrt(sigma2_time));
        P_age * coef_age' ~ normal(0, sqrt(sigma2_age));
        to_vector(P_time * (coef_time_age - rep_matrix(coef_time, nk_age) - rep_matrix(coef_age, nk_time)) * P_age') ~ normal(0, sqrt(sigma2_time_age));

        hivpos ~ bernoulli_logit(X_mat * to_vector(coef_time_age));
}
"

## dat <- subset(hiv, site=="Manicaland" & sex == "Male")
## stan.data <- prepare.data(dat)
## fit <- stan(model_code = stan.code, data=stan.data, iter=100)
