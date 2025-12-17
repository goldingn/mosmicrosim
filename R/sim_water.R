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
  # with the right angle asusmption, we have sin(pi / 2) = 1, which cancels
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
  loss_kg_h <- evaporation_rate_kg_h_m2[t] * surface_area_m2

  # this is already a 1h timestep, so just convert this to m^3
  # 1kg = 1,000cm3
  # 1m3 = 1,000,000 cm3, so
  # 1m3 = 1,000kg
  loss <- loss_kg_h / 1000

  # lose water, avoiding negative volume
  water_volume <- ensure_positive(water_volume - loss)

  # 1000mm per metre, and catchment area is in m2, so convert to m3 (including
  # multiplier on inflow)
  gain_m_h <- rainfall_mm_h[t] / 1000
  gain <- catchment_area * gain_m_h

  # gain water, capping maximum volume
  water_volume <- enforce_max(water_volume + gain, max = max_cone_volume)

  water_volume

}

# given an hourly timeseries of conditions (including rainfall, windspeed and
# altitude), compute a timeseries of water surface area from a simple cone model
simulate_ephemeral_habitat <- function(hourly_climate,
                                       altitude,
                                       initial_volume = 0,
                                       burnin_years = 0,
                                       max_cone_depth = 1,
                                       inflow_multiplier = 1) {

  # add whole year of burnin
  n_times <- length(hourly_climate$water_temperature)
  index <- rep(seq_len(n_times), burnin_years + 1)

  # pull out timeseries needed for simulating
  rainfall <- hourly_climate$rainfall[index]
  air_temperature <- hourly_climate$air_temperature[index]
  relative_humidity <- hourly_climate$humidity[index]
  windspeed <- hourly_climate$windspeed[index]

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

  # simulate the water volume
  n <- length(index)
  volumes <- numeric(length = length(index))
  volume <- 0

  for (t in seq_along(volumes)) {
    volume <- iterate_cone_volume(water_volume = volume,
                                  t = t,
                                  rainfall_mm_h = rainfall,
                                  evaporation_rate_kg_h_m2 = loss_kg_h_m2,
                                  # temperature_c = air_temperature,
                                  # relative_humidity = relative_humidity,
                                  # altitude_m = altitude,
                                  # windspeed_m_s = windspeed,
                                  max_cone_volume = max_cone_volume,
                                  catchment_area = catchment_area)
    volumes[t] <- volume
  }

  # keep only the final year (post burnin)
  keep_index <- tail(seq_along(index), n_times)
  volumes <- volumes[keep_index]

  # convert to (depth then) surface area
  surface_areas <- cone_volume_to_surface(volumes)

}
