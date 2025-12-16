# compute solar attenuation data for a given latitude and longitude,
# either using spatially-varying estimates from GADS, via a precomputed lookup
# table, or by running GADS code (attempting to use fortran version, and falling
# back on the R version if that fails), or using the non-spatial estimate of
# Elterman.
#' @export
solar_attenuation <- function(latitude,
                              longitude,
                              which = c("GADS_lookup",
                                        "GADS",
                                        "Elterman")) {

  which <- match.arg(which)

  solar_attenuation <- switch(which,
                              "GADS_lookup" = solar_attenuation_gads_lookup(
                                latitude = latitude,
                                longitude = longitude),
                              "GADS" = solar_attenuation_gads(
                                latitude = latitude,
                                longitude = longitude,
                                method = "fortran>R"),
                              "Elterman" = solar_attenuation_elterman)

  # return these
  data.frame(wavelengths_nm, solar_attenuation)

}

# NicheMapR wavelengths, in nanometres
wavelengths_nm <- c(290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360,
                    370, 380, 390, 400, 420, 440, 460, 480, 500, 520, 540,
                    560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760,
                    780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980,
                    1000, 1020, 1080, 1100, 1120, 1140, 1160, 1180, 1200,
                    1220, 1240, 1260, 1280, 1300, 1320, 1380, 1400, 1420,
                    1440, 1460, 1480, 1500, 1540, 1580, 1600, 1620, 1640,
                    1660, 1700, 1720, 1780, 1800, 1860, 1900, 1950, 2000,
                    2020, 2050, 2100, 2120, 2150, 2200, 2260, 2300, 2320,
                    2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600, 2700,
                    2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600,
                    3700, 3800, 3900, 4000)

# "The original profile from Elterman, L. 1970. Vertical-attenuation model with
# eight surface meteorological ranges 2 to 13 kilometers. U. S. Airforce
# Cambridge Research Laboratory, Bedford, Mass." From:
# https://github.com/mrke/NicheMapR/blob/f57a2000aaab9d14d1b6390eded823877010d856/R/micro_terra.R#L1475
# note there is another hardcoded non-GADS option floating around in the
# NicheMapR codebase which is described as 'appropriate for Melbourne/Adelaide",
# which we will ignore
solar_attenuation_elterman <- c(0.42, 0.415, 0.412, 0.408, 0.404, 0.4, 0.395,
                                0.388, 0.379, 0.379, 0.379, 0.375, 0.365, 0.345,
                                0.314, 0.3, 0.288, 0.28, 0.273, 0.264, 0.258,
                                0.253, 0.248, 0.243, 0.236, 0.232, 0.227, 0.223,
                                0.217, 0.213, 0.21, 0.208, 0.205, 0.202, 0.201,
                                0.198, 0.195, 0.193, 0.191, 0.19, 0.188, 0.186,
                                0.184, 0.183, 0.182, 0.181, 0.178, 0.177, 0.176,
                                0.175, 0.175, 0.174, 0.173, 0.172, 0.171, 0.17,
                                0.169, 0.168, 0.167, 0.164, 0.163, 0.163, 0.162,
                                0.161, 0.161, 0.16, 0.159, 0.157, 0.156, 0.156,
                                0.155, 0.154, 0.153, 0.152, 0.15, 0.149, 0.146,
                                0.145, 0.142, 0.14, 0.139, 0.137, 0.135, 0.135,
                                0.133, 0.132, 0.131, 0.13, 0.13, 0.129, 0.129,
                                0.128, 0.128, 0.128, 0.127, 0.127, 0.126, 0.125,
                                0.124, 0.123, 0.121, 0.118, 0.117, 0.115, 0.113,
                                0.11, 0.108, 0.107, 0.105, 0.103, 0.1)


# run GADS to compute
solar_attenuation_gads <- function(latitude,
                                   longitude,
                                   method = c("fortran>R", "R", "fortran")) {

  # switch between fortran (crashy) and R versions of GADS
  method <- match.arg(method)
  gads <- switch(method,
                 "fortran>R" = gads_try_fortran,
                 fortran = NicheMapR::rungads,
                 R = NicheMapR::gads.r)

  # Compute for both summer and winter, then average them. Just use 100% RH
  # because NicheMapr does (note this is for the full sky, not at ground level).
  summer <- gads(lat = latitude,
                    lon = longitude,
                    relhum = 1,
                    season = 0)
  winter <- gads(lat = latitude,
                 lon = longitude,
                 relhum = 1,
                 season = 1)

  optical_depth_data <- data.frame(lambda = summer[, 1],
                              opt_depth = (summer[, 2] + winter[, 2]) / 2)

  # fit a sixth-order polynomial to these to smooth them a bit :/
  poly_model <- lm(opt_depth ~ poly(lambda, 6, raw = TRUE),
                   data = optical_depth_data)

  # predict back to the NicheMapR wavelength and return
  predict(poly_model,
          data.frame(lambda = wavelengths_nm))

}

# attempt to run fortran GADS, and fall back on R GADS in the event of an error
gads_try_fortran <- function(lat, lon, relhum, season) {

  tryCatch(
    NicheMapR::rungads(lat = lat,
                      lon = lon,
                      relhum = relhum,
                      season = season),
    error = function(e) {
      NicheMapR::gads.r(lat = lat,
                         lon = lon,
                         relhum = relhum,
                         season = season)
    }
  )

}

# # load the GADS lookup raster (replace this filepath with a system.file() call
# # once the package is building)
ye_gads_file <- system.file("extdata",
                            "ye_gads.tif",
                            package = "mosmicrosim")
ye_gads <- terra::rast(
  ye_gads_file
)
# wrap it, so it can be serialised (e.g. passed to future) without losing the
# pointe
ye_gads_wrapped <- terra::wrap(ye_gads)

solar_attenuation_gads_lookup <- function(latitude = latitude,
                                          longitude = longitude) {

  # lookup the GADS cell number using the raster
  ye_gads <- terra::unwrap(ye_gads_wrapped)
  cell_ids <- terra::extract(ye_gads,
                             cbind(longitude, latitude),
                             cells = FALSE) |>
    dplyr::as_tibble()

  # pull out the solar attenuation information for this cell
  cell_ids |>
    dplyr::left_join(
      solar_attenuation_lookup,
      by = "cell_id"
    ) |>
    dplyr::select(
      -cell_id
    ) |>
    # unpack list
    unlist(
      solar_attenuation
    ) |>
    # scrub names
    `names<-`(1:111)

}

# to do next, make all of these function batched

# lat <- runif(1, -90, 90)
# lon <- runif(1, -180, 180)
#
# bench::mark(
#   tmp_lookup <- solar_attenuation_gads_lookup(latitude = lat,
#                                               longitude = lon),
#   tmp_either <- solar_attenuation_gads(latitude = lat,
#                                        longitude = lon,
#                                        method = "fortran>R"),
#   tmp <- solar_attenuation_gads(latitude = lat,
#                                 longitude = lon,
#                                 method = "fortran"),
#   tmp_R <- solar_attenuation_gads(latitude = lat,
#                                   longitude = lon,
#                                   method = "R"),
#   check = FALSE
# )
#
# # lookup is ~100x faster
#
# # lookup value identical to R, and minimal difference to fortran (different rounding/lookup in the fortran code?)
# max(abs(tmp_lookup - tmp_either))
# max(abs(tmp_lookup - tmp_R))
# max(abs(tmp_lookup - tmp))
