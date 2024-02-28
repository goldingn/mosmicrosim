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
