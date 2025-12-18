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
                          das_function,
                          larval_habitat_area,
                          water_temperature,
                          mdr,
                          efd,
                          ds,
                          timestep = 1 / 24) {

  # given the previous state, surface area, and temperature, compute daily
  # aquatic survival (density- and temperature-dependent)
  das <- das_function(temperature = water_temperature,
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
                          das_function,
                          larval_habitat_area,
                          water_temperature,
                          mdr,
                          efd,
                          ds) {
  mat <- create_matrix(state = state,
                       das_function,
                       larval_habitat_area = larval_habitat_area[t],
                       water_temperature = water_temperature[t],
                       mdr = mdr[t],
                       efd = efd[t],
                       ds = ds[t])
  mat %*% state
}

# access a list of the lifehistory functions needed for the named species
get_lifehistory_functions <- function(
    species = c("An. stephensi", "An. gambiae"),
    storage_path = "data/life_history_params/dehydrated") {

  # enforce the species label
  species <- match.arg(species)

  # load the daily adult survival for either An. gambiae or An. stephensi
  ds_temp_humid = rehydrate_lifehistory_function(
    file.path(storage_path, "ds_temp_humid.RDS")
  )

  # subset to this species
  ds_function <- function(temperature, humidity) {
    ds_temp_humid(temperature, humidity, species = species)
  }

  # for the others, load the relevant RDS object
  species_suffix <- switch(species,
                           "An. stephensi" = "As",
                           "An. gambiae" = "Ag")

  # development rate of aquatic stages as a function of water temperature
  mdr_function = rehydrate_lifehistory_function(
    file.path(storage_path,
              sprintf("mdr_temp_%s.RDS", species_suffix))
  )

  # daily survival probability of aquatic stages as a function of water
  # temperature and density of aquatic stages
  das_function = rehydrate_lifehistory_function(
    file.path(storage_path,
              sprintf("das_temp_dens_%s.RDS", species_suffix))
  )

  # daily egg laying as a function of air temperature
  efd_function = rehydrate_lifehistory_function(
    file.path(storage_path,
              sprintf("efd_temp_%s.RDS", species_suffix))
  )

  # return as a named list of functions
  list(
    ds_function = ds_function,
    mdr_function = mdr_function,
    das_function = das_function,
    efd_function = efd_function
  )

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
  mdr <- lifehistory_functions$mdr_function(
    conditions$water_temperature[index])
  efd <- lifehistory_functions$efd_function(
    conditions$air_temperature[index])
  ds <- lifehistory_functions$ds_function(
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
                           das_function = lifehistory_functions$das_function,
                           larval_habitat_area = larval_habitat_area,
                           water_temperature = water_temperature,
                           mdr = mdr,
                           efd = efd,
                           ds = ds)
    states[t, ] <- state
  }

  # keep only the final year (post burnin)
  keep_index <- tail(seq_along(index), n_times)
  states[keep_index, ]
}
