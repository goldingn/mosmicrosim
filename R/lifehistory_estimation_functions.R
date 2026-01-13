# estimate lifehistory parameters as a function of microclimate conditions, for
# use in population dynamic and transmission modelling

# functions copied over from
# https://github.com/idem-lab/anopheles_stephensi_expansion
# to centralise them in this package


# First: refit the models in Villena et al. to get and interpolate the posterior
# predicted relationship against temperature of key life history parameters.
# This requires adapting the R and JAGS code provided in Villena et al. (the
# basis of the below code), and defining here a number of functions that code
# appears to use from an R package that was formerly hosted on github:
# lorecatta/DENVclimate, which appears to no longer exist anymore, but the
# source code for which is still available on rdrr.io

# These functions are defined below, sourced, reformatted, and slightly modified
# from: https://rdrr.io/github/lorecatta/DENVclimate/src/R/mcmc_utils_all.R
make.briere.samps <- function(coda.samps,
                              nchains = 2,
                              samp.lims = c(151, 5000),
                              sig = TRUE) {
  T0 <- Tm <- cc <- sigma <- NULL
  l1 <- samp.lims[1]
  l2 <- samp.lims[2]
  for (i in 1:nchains) {
    T0 <- c(T0, coda.samps[[i]][l1:l2, 1])
    Tm <- c(Tm, coda.samps[[i]][l1:l2, 2])
    cc <- c(cc, coda.samps[[i]][l1:l2, 3])
    if (sig) sigma <- c(sigma, coda.samps[[i]][l1:l2, 4])
  }
  if (sig) {
    samps <- data.frame(matrix(c(T0, Tm, cc, sigma),
                               ncol = 4,
                               byrow = FALSE))
    names(samps) <- c("T0", "Tm", "c", "sigma")
  }
  else{
    samps <- data.frame(matrix(c(T0, Tm, cc), ncol = 3, byrow=FALSE))
    names(samps) <- c("T0", "Tm", "c")
  }

  return(samps)
}

make.quad.samps<-function(coda.samps, nchains=2, samp.lims=c(151, 5000), sig=TRUE){
  T0<-Tm<-qd<-sigma<-NULL
  l1<-samp.lims[1]
  l2<-samp.lims[2]
  for(i in 1:nchains){
    T0<-c(T0, coda.samps[[i]][l1:l2,1])
    Tm<-c(Tm, coda.samps[[i]][l1:l2,2])
    qd<-c(qd, coda.samps[[i]][l1:l2,3])
    if(sig) sigma<-c(sigma, coda.samps[[i]][l1:l2,4])
  }
  if(sig){
    samps<-data.frame(matrix(c(T0, Tm, qd, sigma), ncol=4, byrow=FALSE))
    names(samps)<-c("T0", "Tm", "qd", "sigma")
  }
  else{
    samps<-data.frame(matrix(c(T0, Tm, qd), ncol=3, byrow=FALSE))
    names(samps)<-c("T0", "Tm", "qd")

  }
  return(samps)
}

briere<-function(t, c, Tm, T0){
  b=c()
  for (i in 1:length(t))
  {
    if(t[i]>T0 && t[i]<Tm){  b[i]<-(c*t[i]*(t[i]-T0)*sqrt(Tm-t[i]))  }
    else {b[i]<-0}
  }
  b
}
quad.2<-function(t, T0, Tm, qd){
  b=c()
  for (i in 1:length(t)){
    if(t[i]>T0 && t[i]<Tm) {b[i]<--qd*(t[i]-T0)*(t[i]-Tm)}
    else {b[i]<-0}
  }
  b
}

make.sims.temp.resp<-function(sim, samps, Temps, thinned, p.name="PDR", trunc.num=0){

  out<-data.sim<-NULL
  out<-list()
  data.sim<-matrix(NA, nrow=length(Temps), ncol=length(thinned))
  for(i in 1:length(thinned)){

    if(sim == "briere"){
      c<-as.numeric(samps$c[thinned[i]])
      Tm<-as.numeric(samps$Tm[thinned[i]])
      T0<-as.numeric(samps$T0[thinned[i]])
      w0<-which(Temps<=T0)
      wM<-which(Temps>=Tm)
      data.sim[,i]<-briere(Temps, c, Tm, T0)
      data.sim[c(w0,wM),i]<-0
    }

    if(sim == "briere.trunc"){
      c<-as.numeric(samps$c[thinned[i]])
      Tm<-as.numeric(samps$Tm[thinned[i]])
      T0<-as.numeric(samps$T0[thinned[i]])
      w0<-which(Temps<=T0)
      wM<-which(Temps>=Tm)
      data.sim[,i]<-briere.trunc(Temps, c, Tm, T0)
      data.sim[c(w0,wM),i]<-0
    }

    if(sim == "quad"){
      T0<-as.numeric(samps$T0[i])
      Tm<-as.numeric(samps$Tm[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad.2(Temps, T0=T0, Tm=Tm, qd=qd)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }

    if(sim == "quad.pos.trunc"){    # added by EAM on 8/31/15
      inter<-as.numeric(samps$inter[i])
      n.slope<- as.numeric(samps$n.slope[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad(Temps, inter=inter, n.slope=n.slope, qd=qd)
      w<-which(data.sim[,i]<trunc.num)
      data.sim[w,i]<-trunc.num
    }

    if(sim == "quad.trunc"){
      T0<-as.numeric(samps$T0[i])
      Tm<-as.numeric(samps$Tm[i])
      qd<-as.numeric(samps$qd[i])
      data.sim[,i]<-quad.trunc(Temps, T0=T0, Tm=Tm, qd=qd)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }

    if(sim == "linear"){
      n.inter<-as.numeric(samps$n.inter[i])
      slope<-as.numeric(samps$slope[i])
      data.sim[,i]<-linear(Temps, inter=-n.inter, slope=slope)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }

    if(sim == "nlinear"){
      inter<-as.numeric(samps$inter[i])
      n.slope<-as.numeric(samps$n.slope[i])
      data.sim[,i]<-linear(Temps, inter=inter, slope=-n.slope)
      w<-which(data.sim[,i]<0)
      data.sim[w,i]<-0
    }

  }

  list(param = p.name,
       T   = Temps,
       fits = data.sim)

}

temp.sim.quants<-function(sim.data, l.temps, byCol=FALSE,
                          probs=c(0.025, 0.975)){

  q<-matrix(NA, nrow=length(probs), ncol=l.temps)
  if(byCol) for(i in 1:l.temps) q[,i]<-quantile(sim.data[,i], probs, na.rm=TRUE)
  else for(i in 1:l.temps) q[,i]<-quantile(sim.data[i,], probs, na.rm=TRUE)

  return(q)
}

# the following code is given by Villena et al. for fitting JAGS models of the
# different types:

# Briere model (MDR, PDR, a)
jags_briere.bug <- "model {

for (i in 1:N) {
Y[i] ~ dnorm(mu[i], tau)T(0,)
mu.temp[i] <- c*T[i]*(T[i]-T0)*sqrt((Tm-T[i])*(Tm>T[i]))
mu[i] <- 0*(mu.temp[i]<0) + mu.temp[i]*(mu.temp[i]>0)

}

c ~ dgamma(1,10)
Tm ~ dunif(25,45)
T0  ~ dunif(0, 24)
sigma<-1/tau
tau ~ dgamma(0.0001, 0.0001)

}"


# Concave down quadratic model (PEA, EFD, bc)
jags_quad.bug <- "model {

for (i in 1:N) {
Y[i] ~ dnorm(mu[i], tau)
# Y[i] ~ dnorm(mu[i], tau)T(0,)
mu[i] <- -qd*(T[i]-T0)*(T[i]-Tm)*((T[i]>T0))*((T[i]<Tm))
}

Tm  ~ dunif(25,45)
T0 ~ dunif(0,24)
qd  ~ dgamma(1,1)
sigma<-1/tau
tau ~ dgamma(0.0001, 0.0001)

}"

define_jags_model <- function(data,
                              model_text,
                              mcmc_params,
                              inits) {
  jags.model(textConnection(model_text),
             data = list(
               Y = data$trait,
               T = data$T,
               N = length(data$T)
             ),
             n.chains = mcmc_params$n_chains,
             inits = inits,
             n.adapt = mcmc_params$n_adapt)

}

# this is much faster than pmax(0, x)
ensure_positive <- function(x) {
  x * as.numeric(x > 0)
}

# spline through the prediction and temperatures, and then constrain to be
# non-negative (as spline can induce negatives)
positive_spline <- function(pred, temps_out) {

  # duplicate this here due to weird lexical scoping issues
  ensure_positive <- function(x) {
    x * as.numeric(x > 0)
  }

  function_raw <- splinefun(temps_out, pred)
  function_positive <- function(temperature) {
    ensure_positive(function_raw(temperature))
  }

  function_positive

}

# make a positive-constrained function from MCMC samples, for the required model
# type
make_function_jags <- function(coda_samples,
                               mcmc_params,
                               model_type = c("briere", "quad"),
                               temps_out = seq(-20, 80, by = 0.1)) {

  model_type <- match.arg(model_type)

  sample_function <- switch(model_type,
                            briere = make.briere.samps,
                            quad = make.quad.samps)

  # This command combines the samples from the chains into a format that we can
  # use for further analyses. Use appropriate model for specific traits
  samps <- sample_function(coda_samples,
                           nchains = mcmc_params$n_chains,
                           samp.lims = c(1, mcmc_params$n_samps))

  # use the parameter samples to get posterior samples of the temperature
  # response curve, and compute the posterior mean curve
  out <- make.sims.temp.resp(sim = model_type,
                             samps,
                             temps_out,
                             thinned = seq(1,
                                           mcmc_params$n_samps,
                                           length = 1000))

  post_mean <- rowMeans(out$fits)

  positive_spline(post_mean, temps_out)

}

# Use the Villena et al JAGS code to estimate the relationship between
# temperature and aquatic development rate for the required species
fit_mdr_temp <- function(data,
                         species = c("An. stephensi", "An. gambiae"),
                         # specify the parameters that control the MCMC
                         mcmc_params = list(
                           n_chains = 5,
                           n_adapt = 10000,
                           n_samps = 20000,
                           n_burn = 10000
                         ),
                         plot_fit = TRUE
) {
  species <- match.arg(species)

  data_sub <- data %>%
    filter(
      trait.name == "mdr",
      specie == species
    )

  # Use the corresponding JAGS model for each trait
  model <- define_jags_model(data_sub,
                             jags_briere.bug,
                             mcmc_params,
                             inits = list(
                               Tm = 31,
                               T0 = 5,
                               c = 0.00007
                             ))

  update(model, mcmc_params$n_samps)
  model_samps_coda <- coda.samples(model,
                                   c("c", "Tm", "T0", "sigma"),
                                   mcmc_params$n_samps)
  if (plot_fit) {
    # visual check for convergence
    plot(model_samps_coda, ask = TRUE)
  }

  make_function_jags(model_samps_coda,
                     mcmc_params,
                     model_type = "briere")

}

# Use the Villena et al JAGS code to estimate the relationship between
# temperature and eggs per female per day for the required species

# note that we model the EFD curve as a down quadratic, as described in the
# paper and SI, but the one plotted in the MS is clearly a Gaussian. I also had
# to remove the observation truncation in the model definition, since a bunch of
# 0s were observed and these were being thrown out, and also to flip the sign on
# the 'quad.2' function above to match the downward quadratic
fit_efd_temp <- function(data,
                         species = c("An. stephensi", "An. gambiae"),
                         # specify the parameters that control the MCMC
                         mcmc_params = list(
                           n_chains = 5,
                           n_adapt = 10000,
                           n_samps = 20000,
                           n_burn = 10000
                         ),
                         plot_fit = TRUE
) {
  species <- match.arg(species)

  data_sub <- data %>%
    filter(
      trait.name == "efd",
      specie == species
    )

  # do EFD as convex down
  model <- define_jags_model(data_sub,
                             jags_quad.bug,
                             mcmc_params,
                             inits = list(
                               Tm = 31,
                               T0 = 5,
                               qd = 0.00007
                             ))

  update(model, mcmc_params$n_samps)
  model_samps_coda <- coda.samples(model,
                                   c("qd", "Tm", "T0", "sigma"),
                                   mcmc_params$n_samps)
  if (plot_fit) {
    # visual check for convergence
    plot(model_samps_coda, ask = TRUE)
  }

  make_function_jags(model_samps_coda,
                     mcmc_params,
                     model_type = "quad")

}

# given an mgcv mortality model, return the function for probability of survival
make_function_mgcv <- function(mortality_model,
                               temps_out = seq(-20, 80, by = 0.1)) {
  pred_df <- data.frame(temperature = temps_out)
  daily_preds <- 1 - predict(mortality_model, pred_df, type = "response")
  positive_spline(daily_preds, temps_out)
}

# DAS: Daily probability of survival during aquatic stages, computed from PEA:
# probability of surviving from egg to adult as a function of temperature,
# refitting the PEA data as a daily rate, using a cox proportional hazards
# model, with the exposure time given by the expected time to emergence from the
# MDR model. Ie. the daily survival probability may be lower at temperatures
# where the MDR is high, since they need to survive for less long. This enables
# us to model larval survival on a much shorter timestep
fit_das_temp <- function(data,
                         mdr_temp_fun,
                         species = c("An. stephensi", "An. gambiae"),
                         plot_fit = TRUE) {

  data_sub <- data %>%
    filter(
      trait.name == "e2a",
      specie == species,
      !is.na(initial)
    )

  data_survival <- data_sub %>%
    mutate(
      survived = round(initial * trait),
      died = initial - survived,
      exposure_period = 1 / mdr_temp_fun(temperature = T),
      exposure_offset = log(exposure_period)
    ) %>%
    rename(
      temperature = T
    )

  # model daily survival probabilities at these temperatures via a proportional
  # hazards model, with the exposure period given by the MDR model
  mortality_model <- mgcv::gam(
    cbind(died, survived) ~ s(temperature),
    family = stats::binomial("cloglog"),
    offset = exposure_offset,
    method = "REML",
    data = data_survival,
    # enforce extra smoothing to make this consistent with the others
    gamma = ifelse(species == "An. gambiae", 2, 20)
  )

  if (plot_fit) {
    plot(mortality_model)
  }

  # function to return function of temperature
  make_function_mgcv(mortality_model)

}

# given the daily survival rate in aquatic stages (DAS) and the aquatic stage
# development rate (MDR), return the probability of surviving from an egg to an
# adult at all (PEA)
make_pea_temp <- function(das_temp, mdr_temp) {
  function(temperature) {
    das_temp(temperature) ^ (1 / mdr_temp(temperature))
  }
}

# density dependence effects on daily aquatic survival (from egg to
# adult) for An stephensi from Evans figure 1 panels D&F
# https://doi.org/10.1002/eap.2334
load_stephensi_survival_data <- function(){

  read_csv("data-raw/life_history_params/evans/data/clean/CSVs/survival.csv",
           show_col_types = FALSE) %>%
    filter(
      Species == "Stephensi",
      AeDens == 0
    ) %>%
    mutate(
      NumInitial = StDens / 2
    ) %>%
    select(
      sex = Sex,
      temperature = Temp,
      density = StDens,
      replicate = Replicate,
      initial = NumInitial,
      survived = NumSurvived
    ) %>%
    mutate(
      # add a random ID for the experimental trial (M and F in together, but
      # multiple replicates per trial)
      trial = paste(temperature, density, replicate, sep = "_"),
      trial_id = as.numeric(factor(trial)),
      # remove the intercept term, and just add a parameter for the males
      is_male = as.numeric(sex == "Male")
    )

}

# load data on temperatire-dependence of life history traits, provided in the
# supplemental information to Villena et al., and clean andaugment it
load_villena_data <- function() {

  # https://github.com/oswaldov/Malaria_Temperature/blob/main/data/traits.csv
  read.csv("data-raw/life_history_params/oswaldov-Malaria_Temperature-16c9d29/data/traits.csv",
           header = TRUE,
           row.names = 1) %>%
    # I can't find a study with this name and year that does aquatic stage
    # survival, only adult survival, so I am assuming this is a mistake and removing
    # it (6 observations)
    filter(
      !(trait.name == "e2a" & ref == "Murdock et al. 2016")
    ) %>%
    # Ditto this study, there is a 2004 study by these authors which does adult
    # survival (3 observations) https://doi.org/10.1079/ber2004316
    filter(
      !(trait.name == "e2a" & ref == "Kirby and Lindsay 2009")
    ) %>%
    # add on the initial number of eggs in the aquatic survival experiements (from
    # going back to the literature)
    mutate(
      initial = case_when(
        # https://doi.org/10.1111/gcb.12240
        ref == "Paaijmans et al. 2013" & trait.name == "e2a" ~ 50,
        # https://doi.org/10.1079/BER2003259
        ref == "Bayoh and Lindsay 2003" & trait.name == "e2a" ~ 30,
        # J VBDs, no doi, initial number not given:
        # https://www.mrcindia.org/journal/issues/464295.pdf
        ref == "Olayemi and Ande 2009" & trait.name == "e2a" ~ NA,
        .default = NA
      )
    )

}

# given a fitted density dependence parameter (relative to some absolute number
# of individuals in a dish), the surface area of that dish, and the type of
# relationship fitted, return a function to *modify* a survival probability
# (e.g. computed just based on temperature) to account for density
make_surv_densmod_function <- function(dd_effect,
                                        surface_area_cm2,
                                        type = c("cox_ph", "logit")) {

  type <- match.arg(type)

  # use faster logit/inverse logit for minor speedups (vs qlogis/plogis)
  logit <- function(p) {
    log(p / (1 - p))
  }

  ilogit <- function(x) {
    1 / (1 + exp(-x))
  }

  if(type == "cox_ph") {
    # if the parameter for the density-dependent effect on survival is estimated
    # as in a Cox proportional hazards model, use a complementary log-log function
    # mapping:
    fun <- function(surv_prob, density) {
      # scale density to the experimental dish size used in estimating the
      # equation
      density_dish <- density * surface_area_cm2
      # get daily *mortality probability* from the survival probability
      daily_mortality_zero_density <- 1 - surv_prob
      # convert to the log hazard for a single day
      loghaz_mortality_zero_density <- log(-log(1 - daily_mortality_zero_density))
      # add on the density effect linear in log density (per dish used to
      # estimate the parameter)
      loghaz_mortality <- loghaz_mortality_zero_density +
        dd_effect * density_dish
      # convert back to a daily *mortality probability*, including the density
      # effect
      daily_mortality <- 1 - exp(-exp(loghaz_mortality))
      # and return as daily survival
      1 - daily_mortality
    }
  } else if (type == "logit") {
    # if the parameter for the density-dependent effect on survival is estimated
    # as the slope in a logit model of log densities use that to transform with
    # complementary log-log function mapping
    fun <- function(surv_prob, density) {
      # scale density to the experimental dish size used in estimating the
      # equation
      density_dish <- density * surface_area_cm2
      # get logit of daily survival probability at zero/low density
      logit_daily_survival_zero_density <- logit(surv_prob)
      # add on the density effect linear in log density (per dish used to
      # estimate the parameter)
      logit_daily_survival <- logit_daily_survival_zero_density +
        dd_effect * density_dish
      # convert back to a daily survival probability, including the density
      # effect, and return
      ilogit(logit_daily_survival)
    }
  }

  # add an attribute to the function denoting the type of prediction, and return
  # the appropriate function
  attr(fun, "type") <- type
  fun

}

# given a survival temperature function, a fitted density dependence parameter
# (relative to some absolute number of individuals in a dish), the surface area
# of that dish, and the type of relationship fitted, return a survival function
# of temperature and density (retained for consistency)
make_surv_temp_dens_function <- function(surv_temp_function,
                                         surv_densmod_function) {

  function(temperature, density) {
    surv_prob_raw <- surv_temp_function(temperature)
    surv_prob_mod <- surv_densmod_function(surv_prob_raw, density)
    surv_prob_mod
  }

}


#
# # given a survival temperature function, a fitted density dependence parameter
# # (relative to some absolute number of individuals in a dish), the surface area
# # of that dish, and the type of relationship fitted,
# # return a survival function of temperature and density
# make_surv_temp_dens_function <- function(surv_temp_function,
#                                          dd_effect,
#                                          surface_area_cm2,
#                                          type = c("cox_ph", "logit")) {
#
#   type <- match.arg(type)
#
#   if(type == "cox_ph") {
#     # if the parameter for the density-dependent effect on survival is estimated
#     # as in a Cox proportional hazards model, use a complementary log-log function
#     # mapping:
#     fun <- function(temperature, density) {
#       # scale density to the experimental dish size used in estimating the
#       # equation
#       density_dish <- density * surface_area_cm2
#       # get daily *mortality probability* at zero/low density at this
#       # temperature
#       daily_mortality_zero_density <- 1 - surv_temp_function(temperature)
#       # convert to the log hazard for a single day
#       loghaz_mortality_zero_density <- log(-log(1 - daily_mortality_zero_density))
#       # add on the density effect linear in log density (per dish used to
#       # estimate the parameter)
#       loghaz_mortality <- loghaz_mortality_zero_density +
#         dd_effect * density_dish
#       # convert back to a daily *mortality probability*, including the density
#       # effect
#       daily_mortality <- 1 - exp(-exp(loghaz_mortality))
#       # and return as daily survival
#       1 - daily_mortality
#     }
#   } else if (type == "logit") {
#     # if the parameter for the density-dependent effect on survival is estimated
#     # as the slope in a logit model of log densities use that to transform with
#     # complementary log-log function mapping
#     fun <- function(temperature, density) {
#       # scale density to the experimental dish size used in estimating the
#       # equation
#       density_dish <- density * surface_area_cm2
#       # get logit of daily survival probability at zero/low density at this
#       # temperature
#       logit_daily_survival_zero_density <- qlogis(surv_temp_function(temperature))
#       # add on the density effect linear in log density (per dish used to
#       # estimate the parameter)
#       logit_daily_survival <- logit_daily_survival_zero_density +
#         dd_effect * density_dish
#       # convert back to a daily survival probability, including the density
#       # effect, and return
#       plogis(logit_daily_survival)
#     }
#   }
#
#   # add an attribute to the function denoting the type of prediction, and return
#   # the appropriate function
#   attr(fun, "type") <- type
#   fun
#
# }

dehydrate_lifehistory_function <- function(fun, path_to_object) {

  arguments <- formals(fun)
  e <- environment(fun)

  # determine how to store the required components
  if(identical(names(arguments), "temperature")) {
    object <- list(
      arguments = arguments,
      temperature_function_raw = e$function_raw,
      rectifier = ensure_positive,
      dummy_function = function() {
        object$rectifier(
          object$temperature_function_raw(
            temperature
          )
        )
      }
    )
  } else if (identical(names(arguments), c("temperature", "density"))) {
    # temperature and density-dependent effects: we need to handle the fact that
    # there are two different functions being used
    type <- attr(fun, "type")
    if(type == "cox_ph") {

      dummy_function <- function() {
        density_dish <- density * object$surface_area_cm2
        daily_mortality_zero_density <- 1 - object$temperature_function_raw(temperature)
        loghaz_mortality_zero_density <- log(-log(1 - daily_mortality_zero_density))
        loghaz_mortality <- loghaz_mortality_zero_density +
          object$dd_effect * density_dish
        daily_mortality <- 1 - exp(-exp(loghaz_mortality))
        1 - daily_mortality
      }

    } else if(type == "logit") {

      dummy_function <- function(temperature, density) {
        density_dish <- density * object$surface_area_cm2
        logit_daily_survival_zero_density <- qlogis(object$temperature_function_raw(temperature))
        logit_daily_survival <- logit_daily_survival_zero_density +
          object$dd_effect * density_dish
        plogis(logit_daily_survival)
      }

    } else {
      stop("unknown type of model")
    }

    object <- list(
      arguments = arguments,
      surface_area_cm2 = e$surface_area_cm2,
      temperature_function_raw = e$surv_temp_function,
      dd_effect = e$dd_effect,
      dummy_function = dummy_function
    )

  } else if (identical(names(arguments), c("temperature", "humidity", "species"))) {
    object <- list(
      arguments = arguments,
      temp_humid_function_raw = e$ds_temp_humid_raw,
      model = e$adult_mortality_model,
      log_hazard_correction = e$log_hazard_correction,
      dummy_function = function() {
        object$temp_humid_function_raw(
          temperature = temperature,
          humidity = humidity,
          species = species,
          model = object$model,
          log_hazard_correction = object$log_hazard_correction)
      }
    )
  } else {
    stop("cannot dehydrate this function")
  }

  saveRDS(object, path_to_object)

}

# load An. gambiae adult survival data under temperature and humidity treatments
# from Bayoh's thesis
load_bayoh_data <- function() {
  read_csv(
    "data-raw/life_history_params/adult_survival/bayoh/bayoh_an_gambiae_adult_survival.csv",
    col_types = cols(
      temperature = col_double(),
      humidity = col_double(),
      sex = col_character(),
      time = col_double(),
      died_cumulative = col_double(),
      alive = col_double()
    )) %>%
    mutate(
      sex = case_when(
        is.na(sex) ~ "mixed",
        .default = sex),
      replicate = 1,
      species = "An. gambiae",
      study = "bayoh"
    )
}

# load data from Krajacich et al. 2020
# https://doi.org/10.1186/s13071-020-04276-y on An gambiae adult survival under
# temperature and humidity combinations, when attempting to induce aestivation,
# and when not.
load_krajacich_data <- function() {
  krajacich <- read_csv("data-raw/life_history_params/adult_survival/krajacich/aestivation.manu.files.scripts/22-Jan-18-R.formatted.masterlist.csv",
                        col_types = cols(
                          Experiment = col_character(),
                          Primed.as = col_character(),
                          Primed = col_character(),
                          Temp = col_character(),
                          `Date of death` = col_double(),
                          Censor = col_double()
                        )) %>%
    # temperatures and humidities were variable for some of these, so remove and
    # keep only the constant and clearly recorded temperatures
    filter(
      Temp != "SE",
      Temp != "18.male"
    ) %>%
    mutate(
      Temp = as.numeric(Temp)
    ) %>%
    # combine priming and experiments to get different replicates
    mutate(
      replicate = paste(Experiment, Primed.as, Primed, sep = "_"),
      replicate = match(replicate, unique(replicate))
    ) %>%
    select(
      -Experiment,
      -Primed.as,
      -Primed
    ) %>%
    rename(
      temperature = Temp,
      time = `Date of death`,
      status = Censor
    )


  # each row is an individual mosquito, Date of death is actually the last day
  # in the timeseries for that mosquito, and status is whether they were alive
  # (ie. when the experiment stopped) or dead (ie. day is the day they died)
  # then. We need to get the number alive at the start of each day (per
  # experiment and temp and humidity), and the number dying on that day

  # get all possible days, for all experiments
  all_combos <- expand_grid(
    replicate = unique(krajacich$replicate),
    temperature = unique(krajacich$temperature),
    time = seq(min(krajacich$time), max(krajacich$time))
  )

  # get the number of mosquitos at the start of each of these experiments
  starting <- krajacich %>%
    group_by(
      replicate,
      temperature
    ) %>%
    summarise(
      starting = n(),
      .groups = "drop"
    )

  # and the numbers dying on each day that one or more died on
  died <- krajacich %>%
    group_by(
      replicate,
      temperature,
      time
    ) %>%
    summarise(
      died = sum(status == 1),
      .groups = "drop"
    )

  # pull these all together, and add on other info
  all_combos %>%
    left_join(
      starting,
      by = join_by(replicate, temperature)
    ) %>%
    left_join(
      died,
      by = join_by(replicate, temperature, time)
    ) %>%
    mutate(
      died = replace_na(died, 0)
    ) %>%
    arrange(
      replicate,
      temperature,
      time
    ) %>%
    group_by(
      replicate,
      temperature
    ) %>%
    mutate(
      died_cumulative = cumsum(died)
    ) %>%
    ungroup() %>%
    mutate(
      alive = starting - died_cumulative
    ) %>%
    select(
      -starting,
      -died
    ) %>%
    mutate(
      sex = "F",
      species = "An. gambiae",
      humidity = 85,
      study = "krajacich"
    )
}

# load An stephensi adult survival data under temperature treatments from
# Miazgowicz et al. 2020 https://doi.org/10.1098/rspb.2020.1093, from Data
# dryad: https://doi.org/10.5061/dryad.8cz8w9gmd and prepare for modelling
load_miazgowicz_data <- function() {
  miazgowicz <- read_csv("data-raw/life_history_params/adult_survival/miazgowicz/constant_master.csv",
                         col_types = cols(
                           Date = col_character(),
                           Time = col_time(format = ""),
                           Treatment = col_double(),
                           Block = col_double(),
                           Donor = col_double(),
                           Female = col_double(),
                           Feed = col_double(),
                           Size = col_character(),
                           Laid = col_double(),
                           Count = col_double(),
                           Dead = col_double(),
                           Day = col_double()
                         )) %>%
    group_by(
      Treatment,
      Block,
      Day
    ) %>%
    # Dead column is coded as 1 for dead, 0 for alive, and 2 for censored
    # (alive, but study discontinued), so recode these censored observations as alive, and count the number of observations
    # mutate(
    #   dead = as.numeric(Dead != 1)
    # ) %>%
    summarise(
      alive = n_distinct(Female[Dead %in% c(0, 2)]),
      died = n_distinct(Female[Dead %in% c(1)]),
      .groups = "drop"
    ) %>%
    group_by(
      Treatment,
      Block
    ) %>%
    mutate(
      died_cumulative = cumsum(died),
      humidity = 80,
      sex = "F",
      species = "An. stephensi",
      study = "miazgowicz"
    ) %>%
    ungroup() %>%
    select(-died) %>%
    rename(
      replicate = Block,
      temperature = Treatment,
      time = Day
    )

}

# load temperature-dependent data from Shapiro et al. 2017
# https://doi.org/10.1371/journal.pbio.2003489 from the data uploaded to data
# dryad https://doi.org/10.5061/dryad.74839
load_shapiro_data <- function() {

  shapiro <- read_csv("data-raw/life_history_params/adult_survival/shapiro/temp.surv.csv",
                      col_types = cols(
                        id = col_double(),
                        expt = col_double(),
                        temp = col_double(),
                        cup = col_double(),
                        day = col_double(),
                        status = col_double()
                      )) %>%
    # combine experimental blocks and cups into a single replicate effect
    mutate(
      replicate = paste(expt, cup, sep = "_"),
      replicate = match(replicate, unique(replicate))
    ) %>%
    rename(
      temperature = temp,
      time = day
    )

  # each row is an individual mosquito, day is the last day in the timeseries
  # for that mosquito, and status is whether they were alive (ie. when the
  # experiment stopped) or dead (ie. day is the day they died) then. We need to
  # get the number alive at the start of each day (per experiment/cup and temp),
  # and the number dying on that day

  # get all possible days, for all experiments
  all_combos <- expand_grid(
    replicate = unique(shapiro$replicate),
    temperature = unique(shapiro$temperature),
    time = seq(min(shapiro$time), max(shapiro$time))
  )

  # get the number of mosquitos at the start of each of these experiments
  starting <- shapiro %>%
    group_by(
      replicate,
      temperature
    ) %>%
    summarise(
      starting = n(),
      .groups = "drop"
    )

  # and the numbers dying on each day that one or more died on
  died <- shapiro %>%
    group_by(
      replicate,
      temperature,
      time
    ) %>%
    summarise(
      died = sum(status == 1),
      .groups = "drop"
    )

  # pull these all together, and add on other info
  all_combos %>%
    left_join(
      starting,
      by = join_by(replicate, temperature)
    ) %>%
    left_join(
      died,
      by = join_by(replicate, temperature, time)
    ) %>%
    mutate(
      died = replace_na(died, 0)
    ) %>%
    arrange(
      replicate,
      temperature,
      time
    ) %>%
    group_by(
      replicate,
      temperature
    ) %>%
    mutate(
      died_cumulative = cumsum(died)
    ) %>%
    ungroup() %>%
    mutate(
      alive = starting - died_cumulative
    ) %>%
    select(
      -starting,
      -died
    ) %>%
    mutate(
      sex = "F",
      species = "An. stephensi",
      humidity = 80,
      study = "shapiro"
    )

}

# Read off the modelled probabilities of An. gambiae surviving to pupation at
# different larval densities, from Figure 1A of Muriu et al. (2013)
# https://doi.org/10.1111/1365-2656.12002, and convert to intercept and slope
# parameters for analysis
load_muriu_dd_survival_parameters <- function() {
  tibble::tribble(
    ~group, ~survival_low_density, ~survival_high_density,
    "yellow", 0.93, 0.37,
    "black", 0.9, 0.6,
    "cyan", 0.87, 0.35,
    "magenta", 0.84, 0.48,
    "red", 0.68, 0.45,
    "blue", 0.59, 0.35,
    "green", 0.44, 0.1,
    "main", 0.82, 0.37
  ) %>%
    mutate(
      density_low = 32,
      density_high = 512,
      species = "An. gambiae",
      temperature = mean(c(22, 30))
    ) %>%
    mutate(
      # transform to linear model scale
      across(starts_with("survival"), qlogis),
      # note Muriu et al. modelled as a log-linear relationship - we skip this
      # and remodel based on the estimates at two densities
      # across(starts_with("density"), log),
      # compute slope for each group
      survival_diff = survival_high_density - survival_low_density,
      density_diff = density_high - density_low,
      slope = survival_diff / density_diff,
      # and intercept
      intercept = survival_low_density - slope * density_low
    ) %>%
    # average these
    group_by(species, temperature) %>%
    summarise(
      intercept = mean(intercept),
      slope = mean(slope),
      .groups = "drop"
    )
}

mortality_prob_to_log_hazard <- function(mortality_prob) {
  log(-log(1 - mortality_prob))
}
# log_hazard_to_mortality_prob <- function(log_hazard) 1 - exp(-exp(log_hazard))

# given the adult mortality model, make a function for daily survival (of
# adults) for the given species, as a function of temperature and humidity,
# accounting for field hazards (log_hazard_correction). This is defined here to
# assist in lexical scoping to make all objects and functions visible later when
# the resulting function is saved in the package
make_surv_temp_humid_function <- function(adult_mortality_model,
                                          log_hazard_correction = 0,
                                          species = c("An. gambiae",
                                                      "An. stephensi")) {
  species <- match.arg(species)
  epsilon <- sqrt(.Machine$double.eps)

  # return a function to make these predictions
  function(temperature, humidity) {

    df <- data.frame(
      temperature = pmax(epsilon, temperature),
      humidity = pmax(epsilon, humidity),
      time = epsilon,
      sex = "F",
      id = 1,
      species = species,
      non_preferred = 0,
      study = ifelse(species == "An. gambiae",
                     "krajacich",
                     "miazgowicz"),
      replicate = 1
    )

    # get the (log) daily hazard of the Cox survival model, based on temperature
    # and humidity
    lab_daily_log_hazard <- predict(adult_mortality_model,
                                    newdata = df,
                                    type = "link")

    # add on the constant hazard due to field conditions
    field_daily_log_hazard <- lab_daily_log_hazard + log_hazard_correction

    # convert to a probability of dying in a single day, under this hazard
    field_mortality_prob <- 1 - exp(-exp(field_daily_log_hazard))

    # convert to probability of surviving
    1 - field_mortality_prob

  }

}

# given a function of temperature and humidity, return a new function to emulate
# this by bilinear interpolation, at a much lower computational cost.
make_temp_humid_interpolator <- function(original_function,
                                         temperature_limits = c(0, 45),
                                         humidity_limits = c(0, 100),
                                         grid_resolution = 250) {

  # set up the grid
  temperature_steps <- seq(temperature_limits[1],
                           temperature_limits[2],
                           length.out = grid_resolution)
  humidity_steps <- seq(humidity_limits[1],
                        humidity_limits[2],
                        length.out = grid_resolution)
  interpolation_grid <- expand.grid(
    temperature = temperature_steps,
    humidity = humidity_steps
  )

  # populate the values
  interpolation_grid$value <- original_function(
    temperature = interpolation_grid$temperature,
    humidity = interpolation_grid$humidity
  )

  # convert the values to a matrix
  value_matrix <- matrix(interpolation_grid$value,
                         nrow = grid_resolution,
                         ncol = grid_resolution)

  # create the list of. interpolation information needed by
  # fields::interp.surface()
  interpolation_info <- list(
    x = temperature_steps,
    y = humidity_steps,
    z = value_matrix
  )

  # make and return the function
  new_function <- function(temperature, humidity) {
    fields::interp.surface(
      interpolation_info,
      cbind(temperature, humidity)
    )
  }

  new_function

}

