# simulate mosquito population dynamics over time, from microclimate and water
# availability data

# functions largely copied over from
# https://github.com/idem-lab/anopheles_stephensi_expansion
# to centralise them in this package

# Build a simple two stage model (aquatic stages and adults), with effects of
# density dependence (daily aquatic survival; DAS), water temperature (DAS and
# aquatic development rate; MDR), air temperature (adult survival; DS, and egg
# laying; EFD), and humidity (DS). Construct as a dynamic matrix.

# construct the matrix appropriately, given the larval
# density in the previous timestep
create_matrix <- function(state,
                          das_densmod_function,
                          larval_habitat_area,
                          water_temperature,
                          mdr,
                          efd,
                          das_zerodensity,
                          ds,
                          timestep = 1 / 24) {

  # given the previous state, surface area, and survival at zero density (ie. as
  # a function of water temperature), compute daily aquatic survival (density- and
  # temperature-dependent)
  das <- das_densmod_function(surv_prob = das_zerodensity,
                              density = state[1] / larval_habitat_area)

  # convert all of these to the required timestep (survivals cumulative, rates
  # linear)
  das_step <- das ^ timestep
  ds_step <- ds ^ timestep
  mdr_step <- mdr * timestep
  efd_step <- efd * timestep

  # construct the matrix
  #     L                A
  # L   das * (1-mdr)    ds * efd
  # A   das * mdr        ds
  matrix(
    c(
      das_step * (1 - mdr_step), # top left
      das_step * (mdr_step), # bottom left
      ds_step * efd_step, # top right
      ds_step # bottom right
    ),
    nrow = 2,
    ncol = 2
  )

}

# iterate the state of the model
iterate_state <- function(state,
                          t,
                          das_densmod_function,
                          larval_habitat_area,
                          water_temperature,
                          mdr,
                          efd,
                          das_zerodensity,
                          ds) {
  mat <- create_matrix(state = state,
                       das_densmod_function = das_densmod_function,
                       larval_habitat_area = larval_habitat_area[t],
                       water_temperature = water_temperature[t],
                       mdr = mdr[t],
                       efd = efd[t],
                       das_zerodensity = das_zerodensity[t],
                       ds = ds[t])
  mat %*% state
}

# simulate for a full timeseries, with optional multiple years of burnin
simulate_population <- function(
    conditions,
    lifehistory_functions,
    larval_habitat_area = rep(pi, length(conditions$day)),
    initial_state = rep(100, 2),
    burnin_years = 1) {

  # add whole year of burnin
  n_times <- length(conditions$water_temperature)
  index <- rep(seq_len(n_times), burnin_years + 1)

  # pull out timeseries needed for simulating
  water_temperature <- conditions$water_temperature[index]
  mdr <- lifehistory_functions$mdr_temp(
    conditions$water_temperature[index])
  efd <- lifehistory_functions$efd_temp(
    conditions$air_temperature[index])
  # aquatic survival at zero density (based on temperature), so we can just
  # modify it later withoutnhaving to re-run temperature sruvival relationship
  # at every iteration
  das_zerodensity <- lifehistory_functions$das_temp(
    conditions$water_temperature[index])
  ds <- lifehistory_functions$ds_temp_humid(
    temperature = conditions$air_temperature[index],
    humidity = conditions$humidity[index])
  larval_habitat_area <- larval_habitat_area[index]

  # simulate the population
  n <- length(index)
  states <- matrix(0, n, 2)
  colnames(states) <- c("aquatic", "adult")
  state <- initial_state

  # pass int he daily aquatic survival function, as it is the only one that is
  # dynamic (depends on the previous state, so introduces density dependence)
  for (t in seq_len(n)) {
    state <- iterate_state(state,
                           t = t,
                           das_densmod_function = lifehistory_functions$das_densmod,
                           larval_habitat_area = larval_habitat_area,
                           water_temperature = water_temperature,
                           mdr = mdr,
                           efd = efd,
                           das_zerodensity = das_zerodensity,
                           ds = ds)
    states[t, ] <- state
  }

  # keep only the final year (post burnin)
  keep_index <- tail(seq_along(index), n_times)
  states[keep_index, ]
}
