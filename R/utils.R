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
