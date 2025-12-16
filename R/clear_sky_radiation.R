# compute clearsky solar radiation over 12 synoptic months for a given latitude
# and longitude, either using a precomputed lookup table, or with NicheMapR

# compute monthly mean clear sky solar radiation max/min for this location using
# NicheMapR microclimate model in solonly mode. the microclimate model gives
# hourly data for representative days. This is summarised across th
# representative days as the mean radiation *during daylight hours*, since it is
# not made explicit how these summary data are represented in terraclimate or
# the reanalysis model it is informed by. This assumption matches with cloud
# cover data for Perth from both Hulmes (core data in the NicheMapR package) and
# from the BOM's directly observed OKTA maps:
# http://www.bom.gov.au/climate/maps/averages/cloud/
#' @export
clear_sky_radiation_12mo <- function(latitude,
                                     longitude,
                                     which = c("lookup",
                                               "nichemapr")) {

  which <- match.arg(which)

  switch(which,
         "lookup" = clear_sky_radiation_12mo_lookup(
           latitude = latitude,
           longitude = longitude
         ),
         "nichemapr" = clear_sky_radiation_12mo_nichemapr(
           latitude = latitude,
           longitude = longitude
         ))

}

clear_sky_radiation_12mo_nichemapr <- function(latitude, longitude) {

  # first, mock up a full micro input to microclimate for 12 arbitrary days,
  # setting cloud cover to 0, and unused dummy vectors for all other inputs

  dummy_dates <- as.Date("2020-01-01") + 1:12
  micro <- create_micro(latitude,
                        longitude,
                        altitude_m = 0,
                        dates = dummy_dates,
                        daily_temp_max_c = rep(0, 12),
                        daily_temp_min_c = rep(0, 12),
                        daily_rh_max_perc = rep(0, 12),
                        daily_rh_min_perc = rep(0, 12),
                        daily_cloud_max_perc = rep(0, 12),
                        daily_cloud_min_perc = rep(0, 12),
                        daily_wind_max_ms = rep(0, 12),
                        daily_wind_min_ms = rep(0, 12),
                        daily_rainfall_mm = rep(0, 12))
  # edit to not run anything except solar radiation
  micro$microinput["solonly"] <- 1

  # don't run daily - use these 12 days of the year as the month midpoints
  micro$microinput["microdaily"] <- 0
  micro$doy <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)

  # run the microclimate model on this
  micro_res <- NicheMapR::microclimate(micro)

  # pull out the solar radiation for this clear sky scenario and aggregate to
  # monthly means
  means <- tapply(X = micro_res$metout[, "SOLR"],
                  INDEX = micro_res$metout[, "DOY"],
                  FUN = mean_if_positive)

  names(means) <- NULL
  means

}


# # load the tc-GADS hybrid lookup raster
tc_gads_file <- system.file("extdata",
                            "tc_gads.tif",
                            package = "mosmicrosim")
tc_gads <- terra::rast(
  tc_gads_file
)

# wrap it, so it can be serialised (e.g. passed to future) without losing the
# pointer
tc_gads_wrapped <- terra::wrap(tc_gads)

clear_sky_radiation_12mo_lookup <- function(latitude = latitude,
                                            longitude = longitude) {

  # lookup the tc-GADS hybrid cell number using the raster
  tc_gads <- terra::unwrap(tc_gads_wrapped)
  cell_ids <- terra::extract(tc_gads,
                             cbind(longitude, latitude),
                             cells = FALSE) |>
    dplyr::as_tibble()

  # pull out the clearsky monthly information for this cell
  cell_ids |>
    dplyr::left_join(
      clear_sky_lookup,
      by = "cell_id"
    ) |>
    dplyr::select(
      -cell_id
    ) |>
    # unpack list
    tidyr::unnest(
      clear_sky_SOLR
    ) |>
    dplyr::pull(
      clear_sky_SOLR
    ) |>
    # scrub names
    `names<-`(NULL)

}

# write a batched version of this

# # load a random location from the tc raster to test this
# tc_native_file <- system.file("extdata",
#                             "tc_native.tif",
#                             package = "mosmicrosim")
# tc_native <- terra::rast(
#   tc_native_file
# )
#
# i <- sample.int(terra::ncell(tc_native), 1)
# coords <- terra::xyFromCell(tc_native, i)
#
# lat <- coords[1, 2]
# lon <- coords[1, 1]
#
# bench::mark(
#   tmp_lookup <- mosmicrosim:::clear_sky_radiation_12mo_lookup(latitude = lat,
#                                                               longitude = lon),
#   tmp_nichemapr <- mosmicrosim:::clear_sky_radiation_12mo_nichemapr(latitude = lat,
#                                        longitude = lon),
#   check = TRUE
# )
# # lookup is ~4x faster
