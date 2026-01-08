# simulate mosquito population dynamics over time, from microclimate and water
# availability data

# functions largely copied over from
# https://github.com/idem-lab/anopheles_stephensi_expansion
# to centralise them in this package

# Build a simple two stage model (aquatic stages and adults), with effects of
# density dependence (daily aquatic survival; DAS), water temperature (DAS and
# aquatic development rate; MDR), air temperature (adult survival; DS, and egg
# laying; EFD), and humidity (DS).

# iterate the state of the model
iterate_state <- function(state,
                          t,
                          das_densmod_function,
                          larval_habitat_area,
                          water_temperature,
                          mdr,
                          efd,
                          das_zerodensity,
                          ds,
                          timestep = 1 / 24) {

  # extract the state
  aquatic <- state$aquatic
  adult <- state$adult

  # extract all of the environmental conditions, for this time

  # given the current aquatic density and current survival at zero density (ie.
  # as a function of water temperature), update daily aquatic survival to be
  # density- (and temperature-) dependent
  aquatic_density <- aquatic / larval_habitat_area[t]
  das_t <- das_densmod_function(surv_prob = das_zerodensity[t],
                                density = aquatic_density)

  # convert all of these to the required timestep

  # survival probabilities cumulative (bernoulli process of surviving)
  aquatic_survival_prob <- das_t ^ timestep
  adult_survival_prob <- ds[t] ^ timestep

  # set fraction emerging to the emergence probability from on emergence hazard
  # model. MDR is the mosquito development rate (rate of aquatic stages
  # developing into adults, per day) so 1/MDR is the expected duration of an
  # aquatic lifestage in days, and 1 - exp(-MDR) is the expected fraction
  # emerging in a single day. We convert MDR to emergence rate in units
  # of 'timestep', and compute the expected fraction emerging in this interval
  emergence_rate <- timestep * mdr[t]
  emergence_fraction <- 1 - exp(-emergence_rate)

  # egg laying rate linear in time
  egg_laying_rate <- efd[t] * timestep

  # iterate the states

  # surviving adults and surviving larvae
  surviving_adult <- adult * adult_survival_prob
  surviving_aquatic <- aquatic * aquatic_survival_prob

  # new aquatic stages (surviving adults times egg laying rate)
  new_aquatic <- surviving_adult * egg_laying_rate

  # remaining larvae (surviving and non-developing larvae)
  remaining_aquatic <- surviving_aquatic * (1 - emergence_fraction)

  # new adults (surviving and developing larvae)
  new_adult <- surviving_aquatic * emergence_fraction

  # update and return the state
  list(
    aquatic = remaining_aquatic + new_aquatic,
    adult = surviving_adult + new_adult
  )

}

# simulate for a full timeseries, with optional multiple years of burnin
simulate_population <- function(
    hourly_conditions,
    lifehistory_functions,
    initial_adult = 100,
    initial_aquatic = 100,
    burnin_years = 1) {

  # add whole year of burnin
  n_times <- length(hourly_conditions$water_temperature)
  index <- rep(seq_len(n_times), burnin_years + 1)

  # pull out environment-dependent lifehistory parameter timeseries

  # mosquito (ie. aquatic) development rate
  mdr <- lifehistory_functions$mdr_temp(
    hourly_conditions$water_temperature[index]
  )

  # eggs per female per day
  efd <- lifehistory_functions$efd_temp(
    hourly_conditions$air_temperature[index]
  )

  # daily aquatic survival at zero density (ie. due to water temperature), which
  # we will modify it later by density, to avoid having to re-run
  # temperature-dependent survival calculations at every iteration
  das_zerodensity <- lifehistory_functions$das_temp(
    hourly_conditions$water_temperature[index]
  )

  # daily (adult) survival at this temperature and humidity
  ds <- lifehistory_functions$ds_temp_humid(
    temperature = hourly_conditions$air_temperature[index],
    humidity = hourly_conditions$humidity[index]
  )

  # relative amount of larval habitat available
  larval_habitat_area <- hourly_conditions$water_surface_area[index]

  # simulate the population
  n <- length(index)

  # vectors for tracking
  adult_states <- rep(NA, n)
  aquatic_states <- rep(NA, n)

  # current state
  state <- list(
    adult = initial_adult,
    aquatic = initial_aquatic
  )

  # pass in the daily aquatic survival density modification function, as it is
  # the only one that is dynamic (depends on the previous state, so introduces
  # density dependence)
  for (t in seq_len(n)) {
    # update the state
    state <- iterate_state(
      state = state,
      t = t,
      das_densmod_function = lifehistory_functions$das_densmod,
      larval_habitat_area = larval_habitat_area,
      mdr = mdr,
      efd = efd,
      das_zerodensity = das_zerodensity,
      ds = ds
    )
    # track the recent states
    aquatic_states[t] <- state$aquatic
    adult_states[t] <- state$adult
  }

  # keep only the final year (post burnin)
  keep_index <- tail(seq_along(index), n_times)

  # and return the vectors, discarding the burnin
  list(
    aquatic = aquatic_states[keep_index],
    adult = adult_states[keep_index]
  )

}
