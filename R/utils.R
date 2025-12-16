# a bunch of utility functions

# is this filepath writeable?
file_writeable <- function(file_path) {

  file_already_existed <- file.exists(file_path)

  # check if the file exists, and temporarily create it if not
  if (!file_already_existed) {
    # delete exit, if it didn't already exist
    on.exit(file.remove(file_path))
    file.create(file_path)
  }

  # then check it's writeable
  file.access(file_path, mode = 2) == 0

}

# return the expected path to the extdata directory of an installed package
package_extdata_dir <- function(package = "NicheMapR") {
  file.path(.libPaths()[1], package, "extdata")
}

# get a valid filename for the nichemapr data object
nichemapr_data_filename <- function() {

  # define a datafile in the package extdata directory
  file_name <- "nichemapr_global_data.RDS"
  file_path <- file.path(package_extdata_dir(), file_name)

  # the file needs to exist to be writeable...

  # if we can't write here, write to a tempfile and let the user know
  if (!file_writeable(file_path)) {
    file_path <- file.path(tempdir(), file_name)
    message("storing NicheMapR global climate data in a temporary file ",
            "valid only for this R session")
  }

  file_path

}

# convert a date vector to a numeric day of the year, indexing from 1
doy <- function(date) {
  as.integer(format(date, "%j"))
}

# check the dates are contiguous, increasing, and non-repeating
check_contiguous_dates <- function(dates) {

  contiguous_dates <- identical(dates,
                                dates[1] + seq_along(dates) - 1)

  if (!contiguous_dates) {
    stop("dates must be contiguous, increasing, and non-repeating ",
         "(a vector of consecutive days)",
         call. = FALSE)
  }

}

# compute mean vapour pressure from mean temperature and vapour pressure deficit (kPa)
vapour_pressure <- function(tmean, vpd_kpa) {

  # Vapour pressure (Pa)
  vapour_pressure <- NicheMapR::WETAIR(
    # drybulb temp (C)
    db = tmean,
    # max relative humidity
    rh = 100)$e

  # subtract vapour pressure deficit from hypothetical amount to get actual

  # convert units
  vpd_pa <- vpd_kpa * 1000
  vapour_pressure - vpd_pa

}

relative_humidity <- function(temperature, vapour_pressure) {
  # compute saturation vapour pressure (Pa)
  saturation_vapour_pressure <- NicheMapR::WETAIR(db = temperature,
                                                  rh = 100)$esat
  # compute and clamp relative humidity
  rh <- (vapour_pressure / saturation_vapour_pressure) * 100
  rh <- pmin(rh, 100)
  rh <- pmax(rh, 0.01)
  rh
}


# adjust the wind speed from the height at which it was measured to the height
# required for modelling, based on the terrain (default to level grass and
# NicheMapR terraclimate default heights).
adjust_wind_speed <- function(wind_speed,
                              height_required = 2,
                              height_measured = 10,
                              terrain = c("level grass",
                                          "open water",
                                          "crops/bushes",
                                          "trees/buildings/mountains")) {

  # From NicheMapR:
  # Correct for fact that wind is measured at 10m height
  # wind shear equation: v / vo = (h / ho)^a
  # where
  # v = the velocity at height h (m/s)
  # vo = the velocity at height ho (m/s)
  # a = the wind shear exponent
  # source http://www.engineeringtoolbox.com/wind-shear-d_1215.html

  # Terrain                       Wind Shear Exponent
  # -----------------------------|---------------------------
  # Open water                    0.1
  # Smooth, level, grass-covered  0.15 (or more commonly 1/7)
  # Row crops 	                  0.2
  # Low bushes with a few trees 	0.2
  # Heavy trees 	                0.25
  # Several buildings 	          0.25
  # Hilly, mountainous terrain 	  0.25

  terrain <- match.arg(terrain)
  exponent <- switch(terrain,
                     "open water" = 0.1,
                     "level grass" = 0.15,
                     "crops/bushes" = 0.2,
                     "trees/buildings/mountains" = 0.25)

  # wind shear equation: v / vo = (h / ho)^a
  wind_speed * (height_required / height_measured) ^ exponent

}

# Calculate the average percentage of the sky covered by cloud. solar_radiation
# is terraclimate's downward surface shortwave radiation (W/m^2).
# clear_sky_radiation is the solar radiation in a clear sky in the same units,
# from Nichemapr. multiplier is a manual muliplier on the output percentage. The
# results will then be clamped to between 0 and 100 percent.
cloud_cover <- function(solar_radiation, clear_sky_radiation, multiplier = 1) {
  cloud <- (1 - solar_radiation / clear_sky_radiation) * 100
  cloud <- cloud * multiplier
  cloud <- pmin(cloud, 100)
  cloud <- pmax(cloud, 0)
  cloud
}

# numerical stuff

# ensure all values of x are greater than or equal to min
# this is much faster than pmax(x, min)
ensure_gte <- function(x, min = .Machine$double.eps) {
  mask <- as.numeric(x >= min)
  x * mask + min * (1 - mask)
}

# this is much faster than pmax(0, x)
ensure_positive <- function(x) {
  x * as.numeric(x > 0)
}

# this is much faster than pmin(x, max)
enforce_max <- function(x, max) {
  too_big <- x > max
  x[too_big] <- max
  x
}

# calculate the mean of the variable only for the time when it has a positive
# value
mean_if_positive <- function(x) {
  mean(x[x > 0])
}
