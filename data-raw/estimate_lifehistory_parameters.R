# Estimate life history functions (life history parameters as a function of
# environmental conditions) of Anopheles gambiae s.l. and Anopheles stephensi
# from a variety of sources, in particular: Villena et al.
# https://doi.org/10.1002/ecy.3685

# This must be run after data-raw/ye_gads.R to build data for the package!

# We use slight modifications of the models defined in that publication, with
# some modification to ensure greater agreement between e.g. daily survival
# probabilities and development times, and to incorporate the effects of
# humidity

# code copied over from
# https://github.com/idem-lab/anopheles_stephensi_expansion
# to centralise it in this package

# load the required R packages

# for modelling
library(rjags)
library(mgcv)
library(lme4)
library(MASS)
library(DHARMa)

# for data handling and plotting
library(tidyverse)

# for contour labelling
library(metR)

# functions in this package
pkgload::load_all(".")

# Note: I had to get jags running on my M2 mac, so I did the following before
# loading rjags:
# Install jags at terminal with: `brew install jags`
# then install rjags pointing to this jags (modify path to wherever jags is,
# found with `which jags` at terminal)
# devtools::install_url("http://sourceforge.net/projects/mcmc-jags/files/rjags/4/rjags_4-4.tar.gz",
#                       args="--configure-args='--with-jags-include=/opt/homebrew/bin/jags/include/JAGS
#                                               --with-jags-lib=/opt/homebrew/bin/jags/lib'")
# Alternative method may be needed
# devtools::install_url("https://sourceforge.net/projects/mcmc-jags/files/rjags/4/rjags_4-4.tar.gz",
#                       configure.args = c("--with-jags-include=/opt/homebrew/bin/jags/include/JAGS
#                                               --with-jags-lib=/opt/homebrew/bin/jags/lib"))


set.seed(2025-12-18)

# load data, provided in the supplemental information to Villena et al., and clean/augment it
data_villena <- load_villena_data()

# use Villena et al. code to fit functions for MDR and EFD, against temperature
# for An. stephensi, and refit PEA in a way that enables daily survival
# probabilities to be computed

# MDR: mosquito development rate (time to move through aquatic stages from egg
# to adult) as a function of temperature
mdr_temp_As <- fit_mdr_temp(data_villena, species = "An. stephensi")
mdr_temp_Ag <- fit_mdr_temp(data_villena, species = "An. gambiae")

# EFD: eggs per female per day as a function of temperature
efd_temp_As <- fit_efd_temp(data_villena, species = "An. stephensi")
efd_temp_Ag <- fit_efd_temp(data_villena, species = "An. gambiae")

# DAS: daily probability of survival in the aquatic stages (eggs, larvae, pupae)
das_temp_As <- fit_das_temp(data_villena,
                            mdr_temp_fun = mdr_temp_As,
                            species = "An. stephensi")
das_temp_Ag <- fit_das_temp(data_villena,
                            mdr_temp_fun = mdr_temp_As,
                            species = "An. gambiae")

# PEA: overall probability of surviving from an egg to an adult (used only for
# plotting)
pea_temp_As <- make_pea_temp(das_temp_As, mdr_temp_As)
pea_temp_Ag <- make_pea_temp(das_temp_Ag, mdr_temp_Ag)

# plot fitted curves and data
curves <- expand_grid(
  temperature = seq(0, 60, by = 1),
  trait_name = c("MDR", "PEA", "EFD"),
  species = c("An. stephensi", "An. gambiae")
) %>%
  mutate(
    trait = case_when(
      trait_name == "MDR" & species == "An. stephensi" ~ mdr_temp_As(temperature),
      trait_name == "MDR" & species == "An. gambiae" ~ mdr_temp_Ag(temperature),
      trait_name == "PEA" & species == "An. stephensi" ~ pea_temp_As(temperature),
      trait_name == "PEA" & species == "An. gambiae" ~ pea_temp_Ag(temperature),
      trait_name == "EFD" & species == "An. stephensi" ~ efd_temp_As(temperature),
      trait_name == "EFD" & species == "An. gambiae" ~ efd_temp_Ag(temperature)
    ),
    trait_name = case_when(
      trait_name == "MDR" ~ "Aquatic development",
      trait_name == "PEA" ~ "Aquatic survival",
      trait_name == "EFD" ~ "Adult eggs/day"
    )
  )

data_villena_plot <- data_villena %>%
  mutate(
    trait_name = toupper(trait.name),
    trait_name = case_when(
      trait_name == "E2A" ~ "PEA",
      .default = trait_name
    ),
    species = specie,
    temperature = T,
    trait_name = case_when(
      trait_name == "MDR" ~ "Aquatic development",
      trait_name == "PEA" ~ "Aquatic survival",
      trait_name == "EFD" ~ "Adult eggs/day"
    )
  ) %>%
  filter(
    trait_name %in% curves$trait_name,
    species %in% curves$species
  )

ggplot(
  curves,
  aes(y = trait,
      x = temperature)
) +
  geom_line() +
  facet_grid(trait_name ~ species,
             scales = "free_y",
             switch = "y") +
  geom_point(
    data = data_villena_plot,
    alpha = 0.3
  ) +
  theme_minimal() +
  theme(strip.placement = "outside") +
  xlab("Temperature (C)") +
  ylab("")

ggsave("data-raw/figures/lifehistory_temperature.png",
       bg = "white",
       width = 5,
       height = 5)

ggplot(
  filter(curves,
         species == "An. stephensi"),
  aes(y = trait,
      x = temperature)
) +
  geom_line() +
  facet_grid(trait_name ~ species,
             scales = "free_y",
             switch = "y") +
  geom_point(
    data = filter(data_villena_plot,
                  species == "An. stephensi"),
    alpha = 0.3
  ) +
  theme_minimal() +
  theme(strip.placement = "outside") +
  xlab("Temperature (C)") +
  ylab("")

ggsave("data-raw/figures/lifehistory_temperature_stephensi.png",
       bg = "white",
       width = 3,
       height = 5)

ggplot(
  filter(curves,
         species == "An. gambiae"),
  aes(y = trait,
      x = temperature)
) +
  geom_line() +
  facet_grid(trait_name ~ species,
             scales = "free_y",
             switch = "y") +
  geom_point(
    data = filter(data_villena_plot,
                  species == "An. gambiae"),
    alpha = 0.3
  ) +
  theme_minimal() +
  theme(strip.placement = "outside") +
  xlab("Temperature (C)") +
  ylab("")

ggsave("data-raw/figures/lifehistory_temperature_gambiae.png",
       bg = "white",
       width = 3,
       height = 5)

# now add density dependence effect to daily aquatic survival (from egg to
# adult) too for An stephensi

# load density and temperature treatment survival data from Evans figure 1
# panels D&F and use it to compute a density coefficent on the log hazard scale
# for the probability of survival. Note that Evans don't give the data in
# survival curve form and don't specify the durations of the experiments so we
# again use the modelled MDR as the exposure period and a Cox PH model to model
# that. Note that this assumes density does not affect development rate, but any
# impact of that on fraction reaching adulthood will be accounted for in the
# survival effect. Evan et al.: https://doi.org/10.1002/eap.2334
stephensi_survival <- load_stephensi_survival_data() %>%
  mutate(
    died = initial - survived,
    exposure_period = 1 / mdr_temp_As(temperature = temperature),
    exposure_offset = log(exposure_period),
    mortality_zero_density = 1 - das_temp_As(temperature),
    cloglog_mortality_zero_density = log(-log(1 - mortality_zero_density))
  )

mortality_model <- lme4::glmer(
  cbind(died, survived) ~ 1 +
    # dummy for sex differences
    is_male +
    # linear effect of density on log hazard
    density +
    # random effect for the trial id (unique per replicate/condition
    # interaction)
    (1|trial_id),
  family = stats::binomial("cloglog"),
  offset = exposure_offset + cloglog_mortality_zero_density,
  data = stephensi_survival
)

# residuals look fine (not surprising, given all the random effects)
plot(mortality_model)

# this is the density effect on *mortality* in the COX PH model
# female is the default, and we don't care about the scaling on the survival
# probability (the daily egg to adult survival prob from Villena is more
# useful), so just scale that by the density term on the logit scale
dd_effect_As <- fixef(mortality_model)[["density"]]

# Note that we could roll all of this into the GAM temperature survival model,
# but that would be hard to model with multiple datasets and noise terms, and
# would require a 2D spline which would be more computationally costly to
# predict from when simulating the model

# The experiment used (Evans et al.) has 250ml water in a 'quart size mason
# jar', which is rather quaint, but not particularly specific. I'm assuming it's
# a Ball brand 'regular mouth canning jar'. According to masonjarlifestyle.com,
# that has a 2 3/8" internal diameter. In real money, that's 6.0325cm, for an
# area of 28.5814687428cm^2, or 0.00285814687m2.
evans_area_cm2 <- pi * (6.0325 / 2) ^ 2

# create the final function, using the Cox proportional hazards mapping
das_temp_dens_As <- make_surv_temp_dens_function(
  surv_temp_function = das_temp_As,
  dd_effect = dd_effect_As,
  surface_area_cm2 = evans_area_cm2,
  type = "cox_ph"
)

# model larval density dependence for An. gambiae by extracting the slope of the
# density dependence relationsship from Muriu et al. (2013) and combining it
# with the parameter for daily aquatic survival as a funciton of temperature
dd_effect_Ag <- load_muriu_dd_survival_parameters() %>%
  pull(slope)

# the dish in muriu is 35cm diameter, so ~~962cm2 surface area
muriu_area_cm2 <- pi * (35 / 2) ^ 2

# this was modelled a logistic function of log density, so we need to model it
# differently for for stephensi
das_temp_dens_Ag <- make_surv_temp_dens_function(
  surv_temp_function = das_temp_Ag,
  dd_effect = dd_effect_Ag,
  surface_area_cm2 =  muriu_area_cm2,
  type = "logit"
)

# plot the survival model
density_plotting <- expand_grid(
  temperature = seq(5, 40, length.out = 100),
  density = seq(0, 7, length.out = 100)
) %>%
  rowwise() %>%
  mutate(
    prob_As = das_temp_dens_As(temperature, density),
    prob_Ag = das_temp_dens_Ag(temperature, density),
  ) %>%
  pivot_longer(
    cols = starts_with("prob"),
    names_to = "species",
    names_prefix = "prob_",
    values_to = "prob"
  ) %>%
  mutate(
    species = case_when(
      species == "As" ~ "An. stephensi",
      species == "Ag" ~ "An. gambiae"
    )
  )

max_grey <- 0.6
cols <- grey(1 - max_grey * seq(0, 1, by = 0.1))
cols[1] <- "transparent"

das_As_plot <- density_plotting %>%
  filter(species == "An. stephensi") %>%
  ggplot(
    aes(y = density,
        x = temperature,
        z = prob)
  ) +
  geom_contour_filled(
    binwidth = 0.1
  ) +
  geom_contour(
    binwidth = 0.1,
    colour = grey(0.2),
    linewidth = 0.5
  ) +
  geom_text_contour(
    binwidth = 0.1,
    nudge_y = -0.2,
    # rotate = FALSE,
    skip = 0
    # label.placer = label_placer_n(3)
  ) +
  coord_cartesian(
    ylim = c(0, 7),
  ) +
  scale_fill_discrete(type = cols) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle(
    "Daily survival probability of aquatic stages",
    expression(italic("Anopheles stephensi"))
  ) +
  ylab(
    expression(paste("Individuals per cm"^"2"))
  ) +
  xlab("Temperature (C)")

ggsave("data-raw/figures/aquatic_survival_stephensi.png",
       plot = das_As_plot,
       bg = "white",
       width = 5,
       height = 5)

# and for An gambiae
das_Ag_plot <- density_plotting %>%
  filter(species == "An. gambiae") %>%
  ggplot(
    aes(y = density,
        x = temperature,
        z = prob)
  ) +
  geom_contour_filled(
    binwidth = 0.1
  ) +
  geom_contour(
    binwidth = 0.1,
    colour = grey(0.2),
    linewidth = 0.5
  ) +
  geom_text_contour(
    binwidth = 0.1,
    nudge_y = -0.05,
    # rotate = FALSE,
    skip = 0
    # label.placer = label_placer_n(3)
  ) +
  coord_cartesian(
    ylim = c(0, 1.8)
  ) +
  scale_fill_discrete(type = cols) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle(
    "Daily survival probability of aquatic stages",
    expression(italic("Anopheles gambiae"))
  ) +
  ylab(
    expression(paste("Individuals per cm"^"2"))
  ) +
  xlab("Temperature (C)")


# need to fix the low-temperature survival for An. gambiae
ggsave("data-raw/figures/aquatic_survival_gambiae.png",
       plot = das_Ag_plot,
       bg = "white",
       width = 5,
       height = 5)

# model adult survival as a function of air temperature and humidity. Reanalyse
# Bayoh data on an gambiae, include other data on longer-lived An gambiae, and
# jointly model An stephensi (with temperature effects, but fixed humidity.)


# fit temperature- and humidity-dependent survival curve

# load Mohammed Bayoh's thesis data for An gambiae with temperature and humidity
# treatments, 16CSS strain originally colonised from Lagos, Nigeria, in 1974
bayoh_Ag <- load_bayoh_data()

# load Krajacich et al 2020's data on inducing aestivation in An gambiae (a
# younger colony than Bayoh; coluzzi from Mali in 2012 and ss from Cameroon in
# 2008) strangely, the species is not distinguished in the results or data
# https://doi.org/10.1186/s13071-020-04276-y
krajacich_Ag <- load_krajacich_data()

# Villena has two datasets for temperature dependent An stephensi adult
# mortality rates, Shapiro et al. 2017
# https://doi.org/10.1371/journal.pbio.2003489 and Miazgowicz et al. 2020
# https://doi.org/10.1098/rspb.2020.1093 (listed as Kerri 2019, presumably
# referring to the preprint). The only adult survival data is from Bayoh's
# thesis, which we already have.
# from the 'long-standing colony' at Walt Reed.
shapiro_As <- load_shapiro_data()

# load data from Miazgowicz et al. 2020 https://doi.org/10.1098/rspb.2020.1093,
# from Data dryad: https://doi.org/10.5061/dryad.8cz8w9gmd from a 'long-standing
# colony (~40 years) of An. stephensi mosquitoes from Pennsylvania State
# University which were originally obtained from the Walter Reed Army Institute
# of Research'
miazgowicz_As <- load_miazgowicz_data()

# fit a proportional hazards model on probability of *mortality* via cloglog
# trick: offset handles the cumulative effect, first term is baseline hazard,
# temperature and humidity have a multiplicative effect, fixed intercepts for
# each sex group

adult_survival_data <- bind_rows(
  bayoh_Ag,
  krajacich_Ag,
  miazgowicz_As,
  shapiro_As
) %>%
  mutate(
    species = factor(species),
    study = factor(study),
    # Define a preferred study for each species, to account for study
    # differences, without confounding the species differences. For An gambiae,
    # the preferred study is the one with the youngest colony. For An. stephensi
    # they are from the same colony (v. old), so just picking one at random
    non_preferred = case_when(
      study %in% c("krajacich", "miazgowicz") ~ 0,
      .default = 1
    ),
    id = row_number()
  ) %>%
  # treat each observation period as an independent observation, so account for
  # the number starting each time period (alive_start) and compute the number
  # that died (died_end) and number that lived (alive_end) in that period
  group_by(temperature, humidity, sex, replicate, species, study) %>%
  mutate(
    alive_start = c(0, alive[-n()]),
    died_end = c(0, diff(died_cumulative)),
    alive_end = alive_start - died_end,
    # how long was this survival interval?
    duration = c(0, diff(time)),
    .after = alive
  ) %>%
  ungroup() %>%
  # use the survival interval times the number of trials as the offset (NB approximation to the betabinomial)
  mutate(
    offset = log(duration * alive_start)
  ) %>%
  # remove rows for time zero starting values (anywhere there weren't mossies
  # alive at the start)
  filter(
    alive_start != 0,
    sex == "F"
  )

# we could model this as binomial with cloglog (survival model with proportional
# hazards), but there is additional dispersion even after fitting smooth terms,
# so we need extra observation-level variance. Betabinomial is not possible in.
# mgcv, so we use a negative binomial approximation to the BB

adult_mortality_model <- mgcv::gam(died_end ~ 1 +
                                     # intercept for species, and for the preferred study type (to
                                     # account for differences in colony age)
                                     species +
                                     # nonlinear interaction between them
                                     s(temperature, log(humidity), k = 15) +
                                     # linear effect of mosquito age (duration of exposure)
                                     time * study,
                                   # penalise all effects to zero (lasso-like regulariser)
                                   select = TRUE,
                                   data = adult_survival_data,
                                   # apply extra smoothing, to incorporate a
                                   # priori information that this should be a
                                   # smooth function
                                   gamma = 10,
                                   # fit with REML
                                   method = "REML",
                                   # account for the differing durations of the observation periods
                                   # using a cloglog and offset (proportional hazards model with
                                   # double censoring)
                                   offset = adult_survival_data$offset,
                                   family = mgcv::nb(link = "log")
)
summary(adult_mortality_model)
gratia::draw(adult_mortality_model)
mgcv::gam.check(adult_mortality_model)

# use randomised quantile residuals to check model fit


sims <- simulateResiduals(adult_mortality_model)
testOutliers(sims, type = "bootstrap")
plot(sims)

# mostly it's fine, but there are still some outliers

# get target-normally-distributed residuals, to explore lack of fit
resid_checks <- adult_survival_data %>%
  select(
    time,
    temperature,
    humidity,
    species,
    study
  ) %>%
  mutate(
    resid = qnorm(sims$scaledResiduals)
  )

# this should look uniform
hist(pnorm(resid_checks$resid))

# what proportion of datapoints are outliers? # percentage?
round(100 * mean(!is.finite(resid_checks$resid)), 1)

# plot where these fall in time vs temp
ggplot(resid_checks,
       aes(y = jitter(time),
           x = jitter(temperature),
           color = is.finite(resid),
           size = !is.finite(resid),
           group = study)) +
  geom_point(
    alpha = 0.3,
    pch = 16
  ) +
  facet_grid(~ study)


# This is in a lab; field mortality is much greater. Charlwood (1997) estimates
# daily adult survival for An gambiae at around 0.83 in Ifakara with temperature
# = 25.6, humidity = 80%. Compute a correction on the hazard scale for the
# lab-based survival estimates.
df_ifakara <- data.frame(temperature = 25.6,
                         humidity = 80,
                         sex = "F",
                         time = sqrt(.Machine$double.eps),
                         species = "An. gambiae",
                         replicate = 1,
                         id = 1,
                         non_preferred = 0,
                         # freshest An. gambiae population we have data for
                         study = "krajacich",
                         off = 0)

lab_daily_log_hazard <- predict(adult_mortality_model, df_ifakara, type = "link")[1]
field_daily_log_hazard <- mortality_prob_to_log_hazard(1 - 0.83)
log_hazard_correction <- field_daily_log_hazard - lab_daily_log_hazard

# construct function of temperature and humidity
ds_temp_humid_raw <- function(temperature, humidity, species, model, log_hazard_correction) {
  epsilon <- sqrt(.Machine$double.eps)
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
  lab_daily_log_hazard <- predict(model, df, type = "link")
  field_daily_log_hazard <- lab_daily_log_hazard + log_hazard_correction
  field_mortality_prob <- 1 - exp(-exp(field_daily_log_hazard))
  1 - field_mortality_prob
}

ds_temp_humid <- function(temperature, humidity, species) {
  ds_temp_humid_raw(temperature, humidity, species,
                    model = adult_mortality_model,
                    log_hazard_correction = log_hazard_correction)
}

ds_temp_humid_Ag <- function(temperature, humidity) {
  ds_temp_humid(temperature, humidity, species = "An. gambiae")
}

ds_temp_humid_As <- function(temperature, humidity) {
  ds_temp_humid(temperature, humidity, species = "An. stephensi")
}

# plot the survival model

# contour bin width and colours
bw <- 0.05
max_grey <- 0.6
cols <- grey(1 - max_grey * seq(0, 1, by = bw))
cols[1] <- "transparent"

density_plotting <- expand_grid(
  temperature = seq(5, 40, length.out = 100),
  humidity = seq(0, 100, length.out = 100),
  species = c("An. stephensi", "An. gambiae")
) %>%
  mutate(
    prob = case_when(
      species == "An. stephensi" ~ ds_temp_humid_As(temperature, humidity),
      species == "An. gambiae" ~ ds_temp_humid_Ag(temperature, humidity),
      .default = NA
    ),
    prob_7d = prob ^ 7
  )


# summarise survivial data to overplot
survival_summary <- adult_survival_data %>%
  group_by(
    species,
    temperature,
    humidity,
    replicate,
    sex
  ) %>%
  ungroup() %>%
  group_by(
    species,
    temperature,
    humidity
  ) %>%
  summarise(
    alive = sum(alive_end),
    total = sum(alive_start),
    .groups = "drop"
  ) %>%
  mutate(
    survival = alive / total,
    survival_7d = survival ^ 7,
    survival_7d_cut = cut(survival_7d, seq(0, 1, by = bw))
  )


ggplot(density_plotting,
       aes(y = humidity,
           x = temperature)) +
  facet_wrap(~species) +
  geom_contour_filled(
    aes(z = prob_7d),
    binwidth = bw
  ) +
  geom_contour(
    aes(z = prob_7d),
    binwidth = bw,
    colour = grey(0.2),
    linewidth = 0.5
  ) +
  geom_text_contour(
    aes(z = prob_7d),
    binwidth = bw,
    nudge_y = -5,
    skip = 0
  ) +
  geom_point(
    data = survival_summary,
    mapping = aes(
      colour = survival_7d_cut
    ),
    shape = 16,
    size = 4
  ) +
  geom_point(
    data = survival_summary,
    shape = 21,
    colour = grey(0.4),
    size = 4
  ) +
  scale_fill_discrete(type = cols) +
  scale_color_discrete(type = cols) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Adult survival vs air temperature and humidity",
          "Probability of surviving 1 week") +
  ylab("Humidity(%)") +
  xlab("Temperature (C)")

ggsave("data-raw/figures/lifehistory_adult_survival.png",
       bg = "white",
       width = 7,
       height = 5)


ggplot(
  density_plotting %>%
    filter(species == "An. stephensi"),
  aes(
    y = humidity,
    x = temperature
  )
) +
  facet_wrap(~species) +
  geom_contour_filled(
    aes(z = prob_7d),
    binwidth = bw
  ) +
  geom_contour(
    aes(z = prob_7d),
    binwidth = bw,
    colour = grey(0.2),
    linewidth = 0.5
  ) +
  geom_text_contour(
    aes(z = prob_7d),
    binwidth = bw,
    nudge_y = -5,
    skip = 0
  ) +
  geom_point(
    data = survival_summary %>%
      filter(species == "An. stephensi"),
    mapping = aes(
      colour = survival_7d_cut
    ),
    shape = 16,
    size = 4
  ) +
  geom_point(
    data = survival_summary %>%
      filter(species == "An. stephensi"),
    shape = 21,
    colour = grey(0.4),
    size = 4
  ) +
  scale_fill_discrete(type = cols) +
  scale_color_discrete(type = cols) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Adult survival vs air temperature and humidity",
          "Probability of surviving 1 week") +
  ylab("Humidity(%)") +
  xlab("Temperature (C)")

ggsave("data-raw/figures/lifehistory_adult_survival.png",
       bg = "white",
       width = 7,
       height = 5)


ggplot(density_plotting,
       aes(y = humidity,
           x = temperature)) +
  facet_wrap(~species) +
  geom_contour_filled(
    aes(z = prob_7d),
    binwidth = bw
  ) +
  geom_contour(
    aes(z = prob_7d),
    binwidth = bw,
    colour = grey(0.2),
    linewidth = 0.5
  ) +
  geom_text_contour(
    aes(z = prob_7d),
    binwidth = bw,
    nudge_y = -5,
    skip = 0
  ) +
  geom_point(
    data = survival_summary,
    mapping = aes(
      colour = survival_7d_cut
    ),
    shape = 16,
    size = 4
  ) +
  geom_point(
    data = survival_summary,
    shape = 21,
    colour = grey(0.4),
    size = 4
  ) +
  scale_fill_discrete(type = cols) +
  scale_color_discrete(type = cols) +
  theme_minimal() +
  theme(legend.position = "none") +
  ggtitle("Adult survival vs air temperature and humidity",
          "Probability of surviving 1 week") +
  ylab("Humidity(%)") +
  xlab("Temperature (C)")

ggsave("data-raw/figures/lifehistory_adult_survival.png",
       bg = "white",
       width = 7,
       height = 5)


# save these in the package to be used later, by loading, extending, and
# resaving sysdata.rda to influde the lifehistory functions

# load the current sysdata.R into a contained environment
e <- new.env()
load("R/sysdata.rda",envir = e)

# add the lifehistory functions to the environment
e$lifehistory_functions <- list(
  mdr_temp_As = mdr_temp_As,
  mdr_temp_Ag = mdr_temp_Ag,
  efd_temp_As = efd_temp_As,
  efd_temp_Ag = efd_temp_Ag,
  das_temp_As = das_temp_As,
  das_temp_Ag = das_temp_Ag,
  pea_temp_As = pea_temp_As,
  pea_temp_Ag = pea_temp_Ag,
  das_temp_dens_As = das_temp_dens_As,
  das_temp_dens_Ag = das_temp_dens_Ag,
  ds_temp_humid = ds_temp_humid
)

# save this again
save(list = names(e),
     file = "R/sysdata.rda",
     envir = e,
     compress = "xz")
