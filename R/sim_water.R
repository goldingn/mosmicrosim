# simulate water availability over time, from microclimate data

# functions largely copied over from
# https://github.com/idem-lab/anopheles_stephensi_expansion
# to centralise them in this package

######
# larval habitat modelling

# assume each pool is a cone with 90 degree angle, and define functions to
# convert between measurements
cone_depth_to_radius <- function(depth) {
  depth * sin(pi / 2)
}

cone_depth_to_volume <- function(depth) {
  pi * cone_depth_to_radius(depth) ^ 2 * depth / 3
}

cone_volume_to_depth <- function(volume) {
  ((volume * 3) / (pi * sin(pi / 2) ^ 2)) ^ (1/3)
}

cone_volume_to_surface <- function(volume) {
  # depth <- cone_volume_to_depth(volume)
  # pi * cone_depth_to_radius(depth) ^ 2
  # with the right angle assumption, we have sin(pi / 2) = 1, which cancels
  # leaving this (which is a bit faster to evaluate in the larval habitat model)
  # pi * (3 / pi) ^ (2/3) * volume ^ (2/3)
  3.046474 * volume ^ 0.666667
}

# use the Buck equation to estimate the vapour pressure in kPa
# see: https://en.wikipedia.org/wiki/Vapour_pressure_of_water
# and: https://en.wikipedia.org/wiki/Arden_Buck_equation
saturation_vapour_pressure_kpa <- function(temperature_celsius) {
  0.61121 * exp((18.678 - temperature_celsius / 234.5) * (temperature_celsius / (257.14 + temperature_celsius)))
}

# air pressure as a function of altitude
air_pressure_kpa <- function(altitude_metres) {
  # At low altitudes above sea level, the pressure decreases by about 1.2 kPa
  # (12 hPa) for every 100  metres.
  p0 <- 101.325 # kPa at sea level
  p0 - 1.2 * altitude_metres / 100
}

# compute the amount of evaporation from a water body (kg per hour) as a function of
# temperature, humidity, altitude, windspeed, and surface area
evaporation_rate_kg_h <- function(temperature_c,
                                  relative_humidity,
                                  altitude_m,
                                  windspeed_m_s,
                                  surface_area_m2) {

  # https://www.engineeringtoolbox.com/evaporation-water-surface-d_690.html
  # rate of evaporation gh (kg/hour) from a water body is:
  # gh = theta * A * (xs - x)
  # theta = (25 + 19 v) = evaporation coefficient (kg/m2h)
  # v = velocity of air above the water surface (m/s)
  # A = surface_area (m^2)
  # xs = *maximum* humidity ratio of saturated air at the same temperature as
  #   the water surface (kg/kg)  (kg H2O in kg Dry Air)
  # x = humidity ratio of the air

  v <- windspeed_m_s
  theta <- (25 + 19 * v)
  A <- surface_area_m2

  # compute the maximum humidity ratio xs, given the saturation vapour pressure
  # p_ws and atmospheric pressure p_a
  # xs = 0.62198 p_ws / (p_a - p_ws)
  p_ws <- saturation_vapour_pressure_kpa(temperature_c)
  p_a <- air_pressure_kpa(altitude_m)
  xs <- 0.62198 * p_ws / (p_a - p_ws)

  # now compute the humidity ratio x for the actual air, from the air pressure
  # and actual vapour pressure

  # relative humidity is the ratio of the partial pressure of water vapour in
  # the air (p), to the saturation vapour pressure (p_s), so use it to convert
  # back to the partial pressure of water vapour p
  # RH = 100 * p_w / p_ws
  p_w <- p_ws * (relative_humidity / 100)
  x <- 0.62198 * p_w / (p_a - p_w)

  # put it all together
  evaporation_kg_hour <- theta * A * (xs - x)
  evaporation_kg_hour
}

# iterate the water volume in the cone, with accounting for evaporation and
# rainfall
iterate_cone_volume <- function(water_volume,
                                t,
                                rainfall_mm_h,
                                evaporation_rate_kg_h_m2,
                                # temperature_c,
                                # relative_humidity,
                                # altitude_m,
                                # windspeed_m_s,
                                max_cone_volume = pi/3,
                                catchment_area = pi) {

  # at each 1h timestep, compute evaporation before inflow
  surface_area_m2 <- cone_volume_to_surface(water_volume)
  loss_kg_h <- evaporation_rate_kg_h_m2[t, ] * surface_area_m2

  # this is already a 1h timestep, so just convert this to m^3
  # 1kg = 1,000cm3
  # 1m3 = 1,000,000 cm3, so
  # 1m3 = 1,000kg
  loss <- loss_kg_h / 1000

  # lose water, avoiding negative volume
  water_volume <- ensure_positive(water_volume - loss)
  # # if not running vectorised, this is faster:
  # water_volume <- max(water_volume - loss, 0)

  # 1000mm per metre, and catchment area is in m2, so convert to m3 (including
  # multiplier on inflow)
  gain_m_h <- rainfall_mm_h[t, ] / 1000
  gain <- catchment_area * gain_m_h

  # gain water, capping maximum volume
  water_volume <- enforce_max(water_volume + gain, max = max_cone_volume)
  # # if not running vectorised, this is faster:
  # water_volume <- min(water_volume + gain, max = max_cone_volume)

  water_volume

}

# given an hourly timeseries of conditions for a single location (including
# rainfall, windspeed and altitude), compute a timeseries of water surface area
# from a simple cone model. This is a thin single-location convenience wrapper
# around simulate_ephemeral_habitat_vectorised() (which the pipeline calls
# directly), so the two stay in lock-step. Burnin is on by default
# (burnin_years = 1): the first year of the timeseries is replayed burnin_years
# times before the period of interest, and that burnin discarded, so the
# returned timeseries starts at the water level implied by the climate rather
# than from the (arbitrary) initial_volume.
simulate_ephemeral_habitat <- function(hourly_climate,
                                       altitude,
                                       initial_volume = 0,
                                       burnin_years = 1,
                                       max_cone_depth = 1,
                                       inflow_multiplier = 1) {

  # reshape the single location's conditions into one-column matrices and defer
  # to the vectorised implementation
  as_column <- function(x) matrix(x, ncol = 1)

  water <- simulate_ephemeral_habitat_vectorised(
    rainfall_matrix = as_column(hourly_climate$rainfall),
    air_temperature_matrix = as_column(hourly_climate$air_temperature),
    humidity_matrix = as_column(hourly_climate$humidity),
    windspeed_matrix = as_column(hourly_climate$windspeed),
    altitude_vector = altitude,
    initial_volume = initial_volume,
    burnin_years = burnin_years,
    max_cone_depth = max_cone_depth,
    inflow_multiplier = inflow_multiplier
  )

  # return the surface-area timeseries for the single location as a vector
  water$water_surface_area[, 1]

}

# Vectorised (across multiple pixels) version of simulate_ephemeral_habitat().
# Given a time-by-pixel matrices of the hourly values of each of the
# microclimate conditions used to model water (rainfall, air temperature,
# humidity, and windspeed), a vector of the altitudes of the
# pixels (defaulting to sea-level), the initial water volume (defaulting to 0),
# a number of years of burnin (defaulting to 1), the maximum cone depth
# (defaulting to 1) and an inflow multiplier to modify run-on of water into each
# pool (defaulting to 1, no run-on), compute a timeseries of water surface area
# from a simple 90 degree angle cone model. With burnin_years >= 1 the first
# year of the timeseries is replayed that many times before the period of
# interest (and the burnin discarded), so the pools start at the water level
# implied by the climate.
simulate_ephemeral_habitat_vectorised <- function(
    rainfall_matrix,
    air_temperature_matrix,
    humidity_matrix,
    windspeed_matrix,
    altitude_vector = rep(0, ncol(rainfall_matrix)),
    initial_volume = 0,
    burnin_years = 1,
    max_cone_depth = 1,
    inflow_multiplier = 1) {

  if (burnin_years < 0) {
    stop("burnin_years must be non-negative, but was ", burnin_years,
         call. = FALSE)
  }

  # add whole years of burnin, repeating the first year for each. Guard against
  # a timeseries shorter than a year, which would otherwise index past the data
  n_times_obs <- nrow(rainfall_matrix)
  n_times_year <- min(365 * 24, n_times_obs)
  first_year_index <- seq_len(n_times_year)
  timeseries_index <- seq_len(n_times_obs)
  index <- c(rep(first_year_index, burnin_years),
             timeseries_index)

  # get other dimensions
  n_pixels <- ncol(rainfall_matrix)
  n_times_sim <- length(index)

  # pull out timeseries needed for simulating
  rainfall <- rainfall_matrix[index, , drop = FALSE]
  air_temperature <- air_temperature_matrix[index, , drop = FALSE]
  relative_humidity <- humidity_matrix[index, , drop = FALSE]
  windspeed <- windspeed_matrix[index, , drop = FALSE]

  # create a corresponding matrix of altitudes
  altitude <- matrix(data = 0,
                     nrow = n_times_sim,
                     ncol = n_pixels) |>
    sweep(MARGIN = 2,
          STATS = altitude_vector,
          FUN = "+")

  # set up the cone model

  # define the cone - assume it cone has a maximum depth of 1m, so a maximum
  # volume of 1.05m3, beyond which it cannot get more full
  max_cone_volume <- cone_depth_to_volume(max_cone_depth)

  # the catchment area is arbitrary (modelling only relative
  # abundance), but set it to the maximum cone surface area - pi!
  catchment_area <- inflow_multiplier * cone_volume_to_surface(max_cone_volume)

  # precalculate the evaporation rate per unit surface area
  loss_kg_h_m2 <- evaporation_rate_kg_h(
    temperature_c = air_temperature,
    relative_humidity = relative_humidity,
    altitude_m = altitude,
    windspeed_m_s = windspeed,
    surface_area_m2 = 1
  )

  # simulate the water volume, vectorised across pixels
  volumes <- matrix(data = NA,
                    nrow = n_times_sim,
                    ncol = n_pixels)

  volume <- rep(initial_volume, n_pixels)

  for (t in seq_len(n_times_sim)) {
    volume <- iterate_cone_volume(water_volume = volume,
                                  t = t,
                                  rainfall_mm_h = rainfall,
                                  evaporation_rate_kg_h_m2 = loss_kg_h_m2,
                                  max_cone_volume = max_cone_volume,
                                  catchment_area = catchment_area)
    volumes[t, ] <- volume
  }

  # keep only the simulation data (post burnin)
  keep_index <- tail(seq_along(index), n_times_obs)
  volumes <- volumes[keep_index, , drop = FALSE]

  # convert to (depth then) surface area
  surface_areas <- cone_volume_to_surface(volumes)

  # return as a named list
  list(
    water_surface_area = surface_areas
  )

}
