# The water simulation used to start every pool empty (initial_volume = 0) with
# no burn-in, so the modelled water surface area - and hence the larval habitat
# area - was biased downwards at the start of the record. simulate_ephemeral_
# habitat_vectorised() (used by the pipeline) and its single-location wrapper
# simulate_ephemeral_habitat() now burn in the water volume by default, replaying
# the first year of climate before the period of interest so the pools start at
# the level implied by the climate. These tests pin that behaviour down.

# a synthetic seasonal climate whose record starts in the wet season (so the
# pool should already hold water at the start), with a hot dry season (so the
# pool dries out completely part way through), returned as time-by-pixel
# matrices with `n_pixels` identical columns
fake_climate_matrices <- function(n_days = 200, annual_rain_mm = 500,
                                  n_pixels = 1) {
  n_hours <- n_days * 24
  day <- (seq_len(n_hours) - 1) %/% 24 + 1
  hour <- (seq_len(n_hours) - 1) %% 24

  wet <- pmax(cos(2 * pi * (day - 1) / n_days), 0)

  set.seed(2024)
  rainfall <- rbinom(n_hours, 1, 0.05 * wet) * wet
  rainfall <- rainfall * annual_rain_mm / sum(rainfall)

  temperature <- 28 + 8 * cos(2 * pi * (day - n_days / 2) / n_days) +
    4 * sin(2 * pi * hour / 24)

  as_cols <- function(x) matrix(rep(x, n_pixels), ncol = n_pixels)
  list(
    rainfall = as_cols(rainfall),
    air_temperature = as_cols(temperature),
    humidity = as_cols(50 + 30 * wet),
    windspeed = as_cols(rep(1, n_hours))
  )
}

run_vectorised <- function(m, initial_volume = 0, burnin_years = 1,
                           altitude = 500) {
  simulate_ephemeral_habitat_vectorised(
    rainfall_matrix = m$rainfall,
    air_temperature_matrix = m$air_temperature,
    humidity_matrix = m$humidity,
    windspeed_matrix = m$windspeed,
    altitude_vector = rep(altitude, ncol(m$rainfall)),
    initial_volume = initial_volume,
    burnin_years = burnin_years
  )$water_surface_area
}

test_that("initial_volume changes the result when there is no burn-in", {
  m <- fake_climate_matrices()
  max_volume <- cone_depth_to_volume(1)

  empty <- run_vectorised(m, initial_volume = 0, burnin_years = 0)
  full <- run_vectorised(m, initial_volume = max_volume, burnin_years = 0)

  # a fuller starting pool means more water early on, and on average
  expect_gt(full[1, 1], empty[1, 1])
  expect_gt(mean(full), mean(empty))
})

test_that("burn-in makes the result independent of the initial volume", {
  m <- fake_climate_matrices()
  max_volume <- cone_depth_to_volume(1)

  # the cone dynamics contract, so after a burn-in the simulation has forgotten
  # where it started
  from_empty <- run_vectorised(m, initial_volume = 0, burnin_years = 1)
  from_full <- run_vectorised(m, initial_volume = max_volume, burnin_years = 1)

  expect_equal(from_empty, from_full, tolerance = 1e-8)
})

test_that("burn-in stops the record starting from a dry pool", {
  m <- fake_climate_matrices()

  # the record starts in the wet season, so the pool should already hold water.
  # starting empty with no burn-in biases the early record downwards
  burnt_in <- run_vectorised(m, burnin_years = 1)
  no_burnin <- run_vectorised(m, initial_volume = 0, burnin_years = 0)

  expect_identical(no_burnin[1, 1], 0)
  expect_gt(burnt_in[1, 1], 0)

  # and that downward bias shows up across the first month of surface areas
  first_month <- seq_len(24 * 30)
  expect_gt(mean(burnt_in[first_month, 1]), mean(no_burnin[first_month, 1]))
})

test_that("surface areas are finite, non-negative and within the cone bounds", {
  m <- fake_climate_matrices()
  max_area <- cone_volume_to_surface(cone_depth_to_volume(1))

  areas <- run_vectorised(m)

  expect_equal(dim(areas), c(nrow(m$rainfall), ncol(m$rainfall)))
  expect_false(anyNA(areas))
  expect_true(all(areas >= 0))
  expect_true(all(areas <= max_area + 1e-8))
  # this climate must actually dry the pool out, or the bounds test is hollow
  expect_true(any(areas == 0))
})

test_that("burnin_years must be non-negative", {
  m <- fake_climate_matrices()
  expect_error(
    run_vectorised(m, burnin_years = -1),
    "burnin_years must be non-negative"
  )
})

test_that("the scalar wrapper matches the vectorised implementation", {
  m <- fake_climate_matrices()

  # the single-location wrapper should defer exactly to the vectorised version
  vectorised <- run_vectorised(m, initial_volume = 0, burnin_years = 1)

  hourly_climate <- list(
    rainfall = m$rainfall[, 1],
    air_temperature = m$air_temperature[, 1],
    humidity = m$humidity[, 1],
    windspeed = m$windspeed[, 1],
    water_temperature = m$air_temperature[, 1]
  )
  scalar <- simulate_ephemeral_habitat(hourly_climate, altitude = 500,
                                       initial_volume = 0, burnin_years = 1)

  expect_equal(scalar, vectorised[, 1], tolerance = 1e-12)
})

test_that("burn-in is applied independently across pixels", {
  # two pixels sharing a climate, but the second only gets 30% of the rain
  m <- fake_climate_matrices(n_pixels = 2)
  m$rainfall[, 2] <- m$rainfall[, 2] * 0.3
  max_volume <- cone_depth_to_volume(1)

  from_empty <- run_vectorised(m, initial_volume = 0, burnin_years = 1)
  from_full <- run_vectorised(m, initial_volume = max_volume, burnin_years = 1)

  # each column burns in to the same trajectory regardless of the start
  expect_equal(from_empty, from_full, tolerance = 1e-8)
  # and the drier pixel ends up with less water
  expect_lt(mean(from_empty[, 2]), mean(from_empty[, 1]))
})
