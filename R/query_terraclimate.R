# read in terraclimate data
library(ncdf4)
library(lubridate)
library(mgcv)

source("R/utils.R")

# note there is a google earth engine for terraclimate...
# https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE#bands

# get the URL to the terraclimate monthly summary data
terraclimate_url <- function(variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad")) {
  variable <- match.arg(variable)
  base_url <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_%s_1958_CurrentYear_GLOBE.nc#fillmismatch"
  sprintf(base_url, variable)
}

# open the connection to terraclimate monthly data
terraclimate_open <- function(variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad")) {
  url <- terraclimate_url(variable = variable)
  ncdf4::nc_open(url)
}

terraclimate_available_indices <- function() {
  # open a connection with the default cube (tmax), and close it again on
  # exiting this function
  con <- terraclimate_open()
  on.exit(ncdf4::nc_close(con))

  # find times and places for which there are data, and return
  list(
    times = ncvar_get(con, "time"),
    longitudes = ncvar_get(con, "lon"),
    latitudes = ncvar_get(con, "lat")
  )
}

# given a slice (index to initial time, lat and long, and number of elements in
# each of these dimensions), extract the array of the named variable from
# terraclimate
terraclimate_fetch <- function(
    slice,
    variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad")) {
  # open a connection with the required cube, and close it again on exiting
  # this function
  variable <- match.arg(variable)
  con <- terraclimate_open(variable)
  on.exit(ncdf4::nc_close(con))

  # fetch the required slice of data and return
  ncdf4::ncvar_get(con,
                   varid = variable,
                   start = slice$start,
                   count = slice$count)
}

# find the nearest lat and long in terra climate to the user's. Check the
# requested coordinate against the available set. If there isn't one within 1/48
# of a degree, get one slightly further (basing this on NicheMapR's logic)
nearest_coord <- function(value, available) {
  error <- abs(available - value)
  acceptable <- error < 1/48
  if (any(acceptable)) {
    index <- which(acceptable)
  } else {
    index <- which(error < 1/47.9)[1]
  }
  index
}

# all days since 1900 until the end of the year
days_since_1900 <- function() {
  start <- as.Date("1900-01-01")
  end <- Sys.Date()
  seq(start, end, "days")
}

# return a sequence of last days of the month, for all months encompassing the
# requested dates
monthly_summary_days <- function(dates, buffer_days = 183) {
  # start of month buffer_days before the first date
  min_date <- lubridate::floor_date(min(dates) - buffer_days, unit = "month")
  # start of month buffer_days after the last date
  max_date <- lubridate::floor_date(max(dates) + buffer_days, unit = "month")
  # sequence of all the first dates of the months in between
  seq(min_date, max_date, by = "1 month")
}

# error nicely if the user provides inaccessible dates
check_dates <- function(dates, terraclimate_available) {

  check_contiguous_dates(dates)

  # get the summary dates available and the time period they cover (the start of
  # the first date summarised to the end of the last)
  available_dates <- days_since_1900()[terraclimate_available$times+1]
  coverage_min <- min(available_dates)
  coverage_max <- lubridate::ceiling_date(max(available_dates),
                                        unit = "month") - 1

  # error if the request cannot be fulfilled
  if (min(dates) < coverage_min | max(dates) > coverage_max) {
    stop("the terraclimate archive covers ",
         coverage_min,
         " to ",
         coverage_max,
         "\n       you requested dates between ",
         min(dates),
         " and ",
         max(dates),
         call. = FALSE
         )
  }

}

terraclimate_build_slice <- function(longitude,
                                     latitude,
                                     dates) {

  # times and places for which there are data in terraclimate
  available <- terraclimate_available_indices()

  # check the user-provided dates are contiguous, increasing, and we can get
  # them from terraclimate
  check_dates(dates, available)

  # the timeseries is summarised by month and indexed by the first day of the
  # month, so find all summary dates relevant to our target dates, with some
  # buffering either side to make sure the splining is consistent for those dates
  required_summary_dates <- monthly_summary_days(dates)

  # the terraclimate temporal index is integer days since the first day of 1900
  # (ie. indexed from 0, unlike R's indexing from 1, so subtract 1) so find our
  # month-start dates as indices
  time_index_1900 <- match(required_summary_dates, days_since_1900()) - 1

  # Remove any of these months that are not available. The dates themselves are
  # all included (checked above), so this is just removing buffer dates we can't
  # access.
  indices_keep <- time_index_1900 %in% available$times
  required_summary_dates <- required_summary_dates[indices_keep]
  time_index_1900 <- time_index_1900[indices_keep]

  # now index the elements of the NCDF array
  time_index <- match(time_index_1900, available$times)

  # Find the nearest terraclimate lat longs to the ones we want
  longitude_index <- nearest_coord(longitude, available$longitudes)
  latitude_index <- nearest_coord(latitude, available$latitudes)

  # return a list with the initial time/place location and size of slice
  list(
    start = c(longitude_index, latitude_index, time_index[1]),
    count = c(1, 1, length(time_index)),
    dates = list(
      start = required_summary_dates,
      end = lubridate::ceiling_date(required_summary_dates, unit = "1 month") -1
    )
  )

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

# download the NicheMapR global climate zipfile
download_nichemapr_data <- function(zip_file_path) {

  # download the zipfile to the zip file path
  file_url <- "https://github.com/mrke/NicheMapR/releases/download/v2.0.0/global.climate.zip"
  download.file(file_url, zip_file_path, mode = "wb")

}

# unzip the NicheMapR global climate zipfile to an NCDF in the R session's
# temp directory
unzip_nichemapr_data <- function(zip_file_path, ncdf_file_path) {
  unzip(zip_file_path,
        files = basename(ncdf_file_path),
        exdir = dirname(ncdf_file_path))
}

# get the NicheMapR global climate data object, downloading it if needed
get_nichemapr_data <- function() {

  # work out where the zipfile should be (wither NicheMapR extdata folder or a
  # tempfile) and where the unzipped ncdf file should be (a tempfile)
  zip_file_path <- nichemapr_data_filename()
  ncdf_file_path <- file.path(tempdir(), "global_climate.nc")

  # download it if it doesn't exist
  if (!file.exists(zip_file_path)) {
    download_nichemapr_data(zip_file_path)
  }

  # unzip it if it doesn't exist
  if (!file.exists(ncdf_file_path)) {
    unzip_nichemapr_data(zip_file_path, ncdf_file_path)
  }

  # read it into memory
  global_climate <- terra::rast(ncdf_file_path)

  # give it names matching definition in NicheMapR micro_global, to query later,
  names(global_climate) <- c(
    "altitude",
    paste0("rainfall", 1:12),
    paste0("rainy_days", 1:12),
    paste0("windspeed", 1:12),
    paste0("temperature_min", 1:12),
    paste0("temperature_max", 1:12),
    paste0("rh_min", 1:12),
    paste0("rh_max", 1:12),
    paste0("cloud_cover", 1:12)
  )

  global_climate

}

# get and cache the NicheMapR global climate data object, using a tempfile
get_cloud_cover_raster <- function() {

  nichemapr_data <- get_nichemapr_data()

  # pull out the cloud cover layers
  cloud_cover_data <- nichemapr_data[[paste0("cloud_cover", 1:12)]]

  # divide by 10 to get percentages
  cloud_cover_data / 10

}

get_altitude_raster <- function() {

  nichemapr_data <- get_nichemapr_data()

  # pull out the altitude layer
  altitude_data <- nichemapr_data[["altitude"]]

  altitude_data

}

# lookup altitude in metres
altitude_m <- function(longitude, latitude) {
  altitude_raster <- get_altitude_raster()
  altitude_data <- terra::extract(altitude_raster,
                                  data.frame(longitude, latitude))
  altitude_data$altitude
}

# calculate the mean of the variable only for the time when it has a positive
# value
mean_if_positive <- function(x) {
  mean(x[x > 0])
}

# compute monthly mean clear sky solar radiation max/min for this location using
# NicheMapR microclimate model in solonly mode. the microclimate model gives
# hourly data for representative days. This is summarised across th
# representative days as the mean radiation *during daylight hours*, since it is
# not made explicit how these summary data are represented in terraclimate or
# the reanalysis model it is informed by. This assumption matches with cloud
# cover data for Perth from both Hulmes (core data in the NicheMapR package) and
# from the BOM's directly observed OKTA maps:
# http://www.bom.gov.au/climate/maps/averages/cloud/
clear_sky_radiation_12mo <- function(latitude, longitude) {

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
  micro_res <- microclimate(micro)

  # pull out the solar radiation for this clear sky scenario and aggregate to
  # monthly means
  means <- tapply(X = micro_res$metout[, "SOLR"],
                  INDEX = micro_res$metout[, "DOY"],
                  FUN = mean_if_positive)

  means

}

# compare this with running micro_global (slower) as called recursively in
# micro_global for microclima diffusity output
# clearsky_radiation_monthly <- clear_sky_radiation_12mo(latitude = latitude,
#                                                        longitude = longitude)
# solar_attenuation_data <- solar_attenuation(latitude, longitude)
# micro_clearsky <- micro_global(loc = c(longitude, latitude),
#                                clearsky = 1,
#                                run.gads = 0,
#                                TAI = solar_attenuation_data$solar_attenuation,
#                                timeinterval = 12,
#                                solonly = 1)
# micro_global_res <- tapply(X = micro_clearsky$metout[, "SOLR"],
#                            INDEX = micro_clearsky$metout[, "DOY"],
#                            FUN = mean_if_positive)
# plot(clearsky_radiation_monthly ~ micro_global_res)
# abline(0, 1)

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

# spline through dates, fitting a periodic spline to seasonality
spline_seasonal <- function(values, dates, dates_predict,
                            knots_per_year_trend = 2,
                            knots_per_year_season = 6) {

  min_date <- min(dates)

  df_fit <- data.frame(
    values = values,
    datenum = as.numeric(dates - min_date),
    doy = doy(dates)
  )

  df_predict <- data.frame(
    datenum = as.numeric(dates_predict - min_date),
    doy = doy(dates_predict)
  )

  # adjust knots based on length of timeseries
  duration_years <- max(df_fit$datenum) / 365

  # trend, based on number of years
  knots_trend <- ceiling(knots_per_year_trend * duration_years)

  # seasonal, based on length if less than one year
  fraction_of_year <- pmin(duration_years, 1)
  knots_season <- ceiling(knots_per_year_season * fraction_of_year)

  # but we need at least 2 knots on each
  knots_trend <- pmax(knots_trend, 2)
  knots_season <- pmax(knots_season, 2)

  m <- mgcv::gam(values ~ s(datenum, k = knots_trend) +
                   s(doy, bs = "cp", k = knots_season),
                 data = df_fit,
                 # cyclic on days of the year
                 knots = list(doy = c(0.5, 365.5)),
                 # overfit as much as possible
                 gamma = 0.1)

  predict(m, df_predict)

}


# Get weather for central Perth WA for 2020-2023
longitude <- 115.86
latitude <- -31.95
dates <- seq(as.Date("2010-03-01"),
             as.Date("2020-12-31"),
             by = "1 day")

# build the indexing slice
slice <- terraclimate_build_slice(longitude = longitude,
                                  latitude = latitude,
                                  dates = dates)

# get all variables
vars <- c("tmax", "tmin", "ppt", "ws", "vpd", "srad")
climate_monthly <- data.frame(
  start_date = slice$dates$start,
  end_date = slice$dates$end,
  mid_date = slice$dates$start + (slice$dates$end - slice$dates$start) / 2)
for(var in vars) {
  climate_monthly[, var] <- terraclimate_fetch(slice, var)
}

# compute relative humidity

# mean temperature and vapour pressure
climate_monthly$tmean <- (climate_monthly$tmax + climate_monthly$tmin) / 2
climate_monthly$vp <- vapour_pressure(climate_monthly$tmean, climate_monthly$vpd)

# RH at max/min temps
climate_monthly$rhmax <- relative_humidity(climate_monthly$tmax, climate_monthly$vp)
climate_monthly$rhmin <- relative_humidity(climate_monthly$tmin, climate_monthly$vp)

# adjust wind speed to measured height and split by max/min
climate_monthly$wsmax <- adjust_wind_speed(climate_monthly$ws)
climate_monthly$wsmin <- climate_monthly$wsmax * 0.1

# compute cloud cover from solar radiation (terraclimate) and expected clear sky
# radiation (microclimate code using coordinates and aerosol data) and manual
# multipliers to enforce daily variation, as used in NicheMapR

# monthly summaries of clear sky radiation (doesn't depend on weather)
clearsky_radiation_monthly <- clear_sky_radiation_12mo(latitude = latitude,
                                                       longitude = longitude)
monthly_index <- lubridate::month(climate_monthly$mid_date)

# convert to cloud cover percentage
climate_monthly$ccmax <- cloud_cover(climate_monthly$srad,
                                     clearsky_radiation_monthly[monthly_index],
                                     multiplier = 1.2)
climate_monthly$ccmin <- cloud_cover(climate_monthly$srad,
                                     clearsky_radiation_monthly[monthly_index],
                                     multiplier = 0.8)

# # compare with the Hulmes interpolated cloud cover raster at this site
# cloud_cover_data <- get_cloud_cover_raster()
# cloud_cover_hulmes <- terra::extract(cloud_cover_data,
#                                      data.frame(longitude, latitude))
# srad_sry <- tapply(climate_monthly$srad, monthly_index, FUN = mean)
#
# # our calculation from solar radiation is too high
#
# par(mfrow = c(1, 2))
#
# # from Hulmes:
# plot(t(cloud_cover_hulmes[, -1]), ylim = c(0, 100))
#
# # from our calculation:
# plot(cloud_cover(srad_sry, clearsky_radiation_monthly, multiplier = 1),
#      ylim = c(0, 100))

# spline these to the requested dates
climate_daily <- data.frame(
  date = dates
)
new_var <- c("tmax", "tmin",
             "rhmax", "rhmin",
             "wsmin", "wsmax",
             "ccmax", "ccmin")
for (var in new_var) {
  climate_daily[, var] <- spline_seasonal(values = climate_monthly[, var],
                                          dates = climate_monthly[, "mid_date"],
                                          dates_predict = climate_daily[, "date"])
}

# it would be nice if we could use the correct likelihood for disaggregation
# in these spline models, but this is cheap and probably not much different

par(mfrow = c(1, 1))
op <- par(mfrow = n2mfrow(length(new_var)),
          mar = c(3, 3, 2, 1))
for (var in new_var) {
  ylims <- range(c(climate_daily[, var], climate_monthly[, var]))
  xlims <- range(c(climate_daily$date, climate_monthly$mid_date))
  plot(climate_daily[, var] ~ climate_daily$date,
       type = "l",
       col = grey(0.4),
       lwd = 2,
       ylab = "",
       xlab = "",
       main = var,
       pch = 16,
       ylim = ylims,
       xlim = xlims)

  points(climate_monthly[, var] ~ climate_monthly$mid_date,
         lwd = 1.5,
         cex = 0.7)
}
par(op)

# need to also obtain:
#  altitude (DEM lookup, easy)

alt <- altitude_m(longitude = longitude, latitude = latitude)


# add a wind shear exponent adjustment to the microclimate simulation set up?

# then build a wrapper function that:

# 1. Downloads terraclimate data in tiles rather than pixels

# 2. Downloads, processes, and saves as monthly daily max and min (then spline
# later) summary data for required quantities

# 3. Caches these, and only downloads new pixels/times as needed (at least by
# month and pixel)



