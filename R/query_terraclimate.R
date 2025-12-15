# NOTE: this is just a script, dumped in a package, which breaks everything.
# this code needs to be wrapped up in functions and called in e.g. a vignette

# get an old install of NicheMapR that works for this code
# remotes::install_github("mrke/NicheMapR@v3.3.1")

# read in terraclimate data
library(ncdf4)
library(lubridate)
library(mgcv)
library(NicheMapR)
library(tidyverse)
library(terra)
library(abind)

source("R/utils.R")
source("R/create_micro.R")
source("R/solar_attenuation.R")

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

  # cache the available indices
  tc_avail_indices <- options()$.tc_avail_indices

  # download if required
  if (is.null(tc_avail_indices)) {

    # open a connection with the default cube (tmax), and close it again on
    # exiting this function
    con <- terraclimate_open()
    on.exit(ncdf4::nc_close(con))

    # find times and places for which there are data, and return
    tc_avail_indices <-list(
      times = ncvar_get(con, "time"),
      longitudes = ncvar_get(con, "lon"),
      latitudes = ncvar_get(con, "lat")
    )

    options(.tc_avail_indices = tc_avail_indices)
  }

  tc_avail_indices
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

# build a terraclimate slice for rectangular spatial extent, over the
# requested dates
terraclimate_build_slice_extent <- function(extent,
                                            dates) {

  # check the extent object
  if(!inherits(extent, "SpatExtent")) {
    stop("extent should be a terra SpatExtent object, e.g. created like this:
         extent <- ext(<xmin>, <xmax>, <ymin>, <ymax>)",
         call. = FALSE)
  }

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

  # Find the nearest terraclimate lat longs to the edges of the extent
  ext_list <- as.list(extent)
  xmin_index <- nearest_coord(ext_list$xmin, available$longitudes)
  xmax_index <- nearest_coord(ext_list$xmax, available$longitudes)
  # note: indices are ordered left -> right (like coordinates) and top -> bottom
  # (the reverse of coordinates), so we need to flip the latitudes
  ymin_index <- nearest_coord(ext_list$ymax, available$latitudes)
  ymax_index <- nearest_coord(ext_list$ymin, available$latitudes)

  xrange <- xmax_index - xmin_index + 1
  yrange <- ymax_index - ymin_index + 1

  # return a list with the initial time/place location and size of slice
  list(
    start = c(xmin_index, ymin_index, time_index[1]),
    count = c(xrange, yrange, length(time_index)),
    dates = list(
      start = required_summary_dates,
      end = lubridate::ceiling_date(required_summary_dates, unit = "1 month") -1
    ),
    longitudes = seq(ext_list$xmin, ext_list$xmax, length.out = xrange),
    # flip back the latitudes
    latitudes = seq(ext_list$ymax, ext_list$ymin, length.out = yrange)
  )

}

# build a terraclimate slice for a single point in space, over the requested
# dates
terraclimate_build_slice_point <- function(longitude,
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

# given a template raster object, build a terraclimate-harmonised 5km template.
# Ie. a raster covering all the non-NA cells in the template raster whose cells
# align with those in terraclimate
make_terraclimate_template <- function(template) {

  # make template binary
  template <- template * 0 + 1

  # trim, expand, then buffer these by one cell, to ensure we can interpolate as
  # much as possible when resampling back to this raster from the terraclimate one
  template <- template |>
    terra::trim() |>
    terra::extend(1) |>
    terra::focal(w = 3,
                 fun = max,
                 na.rm = TRUE)

  # now find the equivalent terraclimate raster cells

  # find the grid setup in terraclimate
  terraclimate_available <- terraclimate_available_indices()

  # create an empty raster that aligns with terraclimate and has extent no smaller
  # than our raster
  dims <- as.list(ext(template))
  x_valid <- terraclimate_available$longitudes >= dims$xmin &
    terraclimate_available$longitudes <= dims$xmax
  y_valid <- terraclimate_available$latitudes >= dims$ymin &
    terraclimate_available$latitudes <= dims$ymax
  tc_ext <- ext(c(range(terraclimate_available$longitudes[x_valid]),
                  range(terraclimate_available$latitudes[y_valid])))

  # find the non-NA elements in this raster by building a slice with one date and
  # extracting one value
  raster_slice <- terraclimate_build_slice_extent(tc_ext,
                                                  dates = as.Date("2020-01-01"))
  # this function is set up to extract extra dates to ensure we can spline to
  # these dates, but we only need one value for this so manually switch to the one
  # date at the middle of the time slice
  timesteps <- raster_slice$count[3]
  raster_slice$start[3] <- raster_slice$start[3] + ceiling(timesteps / 2)
  raster_slice$count[3] <- 1
  # get the values
  vals <- terraclimate_fetch(raster_slice)

  # put them in a raster
  tc_template <- rast(nrows = sum(y_valid),
                      ncols = sum(x_valid),
                      extent = tc_ext,
                      vals = t(vals))

  # resample the target raster to this one
  new_template_mask <- terra::resample(template,
                                       tc_template,
                                       method = "max")

  # and use it to mask the terraclimate values, set values to 0, and return
  tc_template_masked <- terra::mask(tc_template, new_template_mask)
  tc_template_masked <- tc_template_masked * 0

  # now trim it again, and return
  tc_template <- terra::trim(tc_template_masked)
  tc_template

}

# given a vector defining an extent object, and a terraclimate template raster,
# return a logical for whether there are any non-na values in that extent
check_tile <- function(vec, template) {
  this_tile <- ext(vec)
  this_raster <- terra::crop(template, this_tile)
  length(cells(this_raster)) > 0
}

# given a vector defining an extent object, and a terraclimate template raster,
# return a vector for a trimmed extent object to the smallest rectangle
# including only non-NA cells
trim_tile <- function(vec, template) {

  # make a raster with this extent, trim it, and convert back ot an extent
  tile <- terra::ext(vec)
  raster <- terra::crop(template, tile)
  raster_trimmed <- terra::trim(raster)
  tile_trimmed <- terra::ext(raster_trimmed)

  # convert to a spatvector, negative buffer by a quarter of a cell so each
  # cell centroid is in only one tile, and convert back ot an extent
  buffer_distance <- -1 * res(template) / 4
  spatvector_trimmed <- vect(tile_trimmed)
  spatvector_trimmed_shrunk <- terra::buffer(spatvector_trimmed,
                                             buffer_distance)
  tile_trimmed_shrunk <- terra::ext(spatvector_trimmed_shrunk)

  # return as a vector
  as.vector(tile_trimmed_shrunk)
}

# given a template raster denoting non-NA cells for processing, and a target
# number of tiles, return the dimensions for approximately that number of tiles
# as a matrix with each row giving the tile extents (xmin, xmax, ymin, ymax).
# Tiles will only be returned that contain non-NA cells, and will be trimmed to
# include only rows and columns that contain non-NA cells.
make_tiles <- function(template, target_n_tiles = 100) {

  # work out the fraction of non-missing cells in the raster, and adjust up the
  # number of tiles
  n_cells_all <- prod(dim(template)[c(1:2)])
  n_cells_valid <- length(cells(template))
  prop_cells_valid <- n_cells_valid / n_cells_all
  target_n_tiles_adj <- target_n_tiles / prop_cells_valid

  # aggregate a raster into approximately this number of square tiles
  agg_factor <- round(sqrt(n_cells_all / target_n_tiles_adj))
  tiles_raster <- terra:::aggregate(template, agg_factor)
  tiles <- terra::getTileExtents(template, tiles_raster)

  # weed out tiles with no values in them
  tiles_valid <- apply(tiles,
                       MARGIN = 1,
                       FUN = check_tile,
                       template = template)
  tiles <- tiles[tiles_valid, ]

  # trim them down to only the valid cells, and shrink them all by less than half
  # a cell to prevent excessive processing of cells
  tiles <- t(apply(tiles,
                   MARGIN = 1,
                   FUN = trim_tile,
                   template = template))

  tiles
}

# given an integer 'tile_number' for the tile to extract (indexing tiles), a
# matrix 'tiles' of tile extents (wach rwo giving the bounding box), a vector of
# 'dates' for which we need data, and a vector of 'variables' required, extract
# monthly terraclimate data on those variables, for all cells in that tile, for
# all available months, where possible extending beyond 'dates' to enable
# improved spline interpolation, and return as a tibble, identifying the cell by
# the latitude and longitude of its centroid.
terraclimate_extract_tile <- function(tile_number, tiles, dates, variables) {

  # make an extent object for this tile
  tile <- ext(tiles[tile_number, ])

  # build the slice for all pixels (including NAs) in this tile
  tile_slice <- terraclimate_build_slice_extent(
    extent = tile,
    dates = dates
  )

  # download each of the variables for this slice
  stack_list <- list()
  for(var in variables) {
    stack_list[[var]] <- terraclimate_fetch(tile_slice, var)
  }

  # combine into tidy long-format tibble, by making a 4D array with appropriate
  # dimension names, coercing to long form via a table, and then tidying up the
  # tibble
  stack <- do.call(abind,
                   c(stack_list, list(along = 4)))

  dimnames(stack) <- list(
    longitude = tile_slice$longitudes,
    latitude = tile_slice$latitudes,
    start = as.character(tile_slice$dates$start),
    variable = variables
  )

  stack |>
    as.data.frame.table(
      stringsAsFactors = FALSE
    ) |>
    as_tibble() |>
    rename(
      value = Freq
    ) |>
    filter(
      !is.na(value)
    ) |>
    mutate(
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude),
      start = as.Date(start)
    ) |>
    # add on the end dates
    left_join(
      as_tibble(tile_slice$dates),
      by = "start"
    ) |>
    relocate(
      end,
      .after = start
    )

}


# demo


# Get weather for central Perth WA for 2020-24
longitude <- 115.86
latitude <- -31.95
dates <- seq(as.Date("2020-01-01"),
             as.Date("2024-12-31"),
             by = "1 day")

# build the indexing slice at this location
slice <- terraclimate_build_slice_point(longitude = longitude,
                                        latitude = latitude,
                                        dates = dates)

# get all variables
vars <- c("tmax", "tmin",  # temperature
          "ppt",  # precipitation?
          "ws",  # wind speed
          "vpd",  # vapor pressure deficit (for rel humidity)
          "srad")  # solar radiation

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
climate_monthly$rhmax <- relative_humidity(climate_monthly$tmin, climate_monthly$vp)
climate_monthly$rhmin <- relative_humidity(climate_monthly$tmax, climate_monthly$vp)

# adjust wind speed to from measured height (10m) to match the other inputs (2m)
# and split by max/min
climate_monthly$wsmax <- adjust_wind_speed(climate_monthly$ws,
                                           height_required = 2,
                                           height_measured = 10)
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

# convert precipitation to log scale to spline
climate_monthly$rainfall_mm <- climate_monthly$ppt
climate_monthly$log_rainfall_mm <- log1p(climate_monthly$rainfall_mm)

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
             "ccmax", "ccmin",
             "log_rainfall_mm")
for (var in new_var) {
  climate_daily[, var] <- spline_seasonal(values = climate_monthly[, var],
                                          dates = climate_monthly[, "mid_date"],
                                          dates_predict = climate_daily[, "date"])
}

# convert log rainfall back
climate_daily$rainfall_mm <- expm1(climate_daily$log_rainfall_mm)

# it would be nice if we could use the correct likelihood for disaggregation
# in these spline models, but this is cheap and probably not much different

var_plot <- c("tmax", "tmin",
              "rhmax", "rhmin",
              "wsmin", "wsmax",
              "ccmax", "ccmin",
              "rainfall_mm")
par(mfrow = c(1, 1))
op <- par(mfrow = n2mfrow(length(var_plot)),
          mar = c(3, 3, 2, 1))
for (var in var_plot) {
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

altitude <- altitude_m(longitude = longitude, latitude = latitude)


# plug these into the microclimate simulation:

shade_proportion <- 0.95
adult_height <- 0.1

micro <- create_micro(latitude = latitude,
                      longitude = longitude,
                      altitude_m = altitude,
                      dates = dates,
                      daily_temp_max_c = climate_daily$tmax,
                      daily_temp_min_c = climate_daily$tmin,
                      daily_rh_max_perc = climate_daily$rhmax,
                      daily_rh_min_perc = climate_daily$rhmin,
                      daily_cloud_max_perc = climate_daily$ccmax,
                      daily_cloud_min_perc = climate_daily$ccmin,
                      daily_wind_max_ms = climate_daily$wsmax,
                      daily_wind_min_ms = climate_daily$wsmin,
                      daily_rainfall_mm = climate_daily$rainfall_mm,
                      weather_height_m = 2,
                      adult_height_m = adult_height,
                      even_rain = FALSE,
                      shade_prop = shade_proportion)

# profvis::profvis(
system.time(
  sim <- NicheMapR::microclimate(micro)
)
# )

# # it's frustrating that runshade=1 is required, contact Mike with a reprex and
# # to to ask why?
# micro2 <- micro
# micro2$microinput["runshade"] <- 0
# system.time(
#   sim2 <- NicheMapR::microclimate(micro2)
# )
# summary(sim2$shadmet[, "TALOC"])


# plot some examples of microclimate conditions, with the external conditions
# over the top

n_days <- 52 * 7
days <- 1:n_days
hours <- 0 * 24 + 1:(24 * n_days)
dates_plot <- dates[days] + lubridate::hours(1)
hours_plot <- dates[1] + lubridate::hours(hours)

par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1)

vals <- sim$shadmet[, "TALOC"][hours]
plot(vals ~ hours_plot,
     type = "l",
     ylab = "C",
     las = 1,
     lwd = 0.1,
     ylim = c(0, 40))
lines(climate_daily$tmax[days] ~ dates_plot,
      lty = 2)
lines(climate_daily$tmin[days] ~ dates_plot,
      lty = 2)
title(main = "Microclimate temperature")

vals <- sim$shadmet[, "RHLOC"][hours]
plot(vals ~ hours_plot,
     type = "l",
     ylab = "%",
     las = 1,
     lwd = 0.1,
     ylim = c(0, 100))
lines(climate_daily$rhmin[days] ~ dates_plot,
      lty = 2)
lines(climate_daily$rhmax[days] ~ dates_plot,
      lty = 2)
title(main = "Relative humidity")

vals <- sim$shadmet[, "VLOC"][hours]
plot(vals ~ hours_plot,
     type = "l",
     ylab = "m/s",
     las = 1,
     lwd = 0.1,
     ylim = c(0, 5))
lines(climate_daily$wsmax[days] ~ dates_plot,
      lty = 2)
lines(climate_daily$wsmin[days] ~ dates_plot,
      lty = 2)
title(main = "Wind speed")



# batch process in tiles


# load a template raster for Africa
template <- rast("~/Dropbox/github/ir_cube/data/clean/raster_mask.tif")

# create a corresponding terraclimate raster (aligned with terraclimate grid
# location, and excluding missing cells in terrclimate or more than 1 cell away
# from template cells with data)
tc_template <- make_terraclimate_template(template)

# now batch-process this by defining tiles covering the continent,
# and extracting whole slices for those tiles for the required times

# make some tiles
tiles <- make_tiles(tc_template, target_n_tiles = 100)

# # not too bad
# nrow(tiles)
#
# # check how they look
# par(mfrow = c(1, 1))
# plot(tc_template,
#      box = FALSE,
#      axes = FALSE)
# for(i in seq_len(n_tiles)) {
#   plot(ext(tiles[i, ]),
#        add = TRUE)
# }
# # beautiful.
#
# # zoom in to check the margins are non-overlapping and contain centoids, as
# # expected
# zoom <- ext(-7.5,
#             -7.2,
#             32.3,
#             32.6)
# plot(tc_template,
#      ext = zoom,
#      box = FALSE,
#      axes = FALSE)
# for(i in seq_len(n_tiles)) {
#   plot(ext(tiles[i, ]),
#        add = TRUE)
# }
# zoom_rast <- crop(tc_template, zoom)
# zoom_points <- xyFromCell(zoom_rast, cells(zoom_rast))
# points(zoom_points,
#        pch = 16)

# define all dates to extract per tile
dates <- seq(as.Date("2020-01-01"),
             as.Date("2024-12-31"),
             by = "1 day")

# extract all terracclimate data for this tile
tile_data <- terraclimate_extract_tile(tile_number = 1,
                                       tiles = tiles,
                                       dates = dates,
                                       variables = vars)

# plot some things
tile_data_plot <- tile_data |>
  # add a cell id for plotting
  group_by(
    latitude,
    longitude
  ) |>
  mutate(
    cell_id = cur_group_id(),
    .before = everything()
  ) |>
  ungroup()

tile_data_plot |>
  filter(
    cell_id %in% sample.int(n_distinct(tile_data_plot$cell_id), 5)
  ) |>
  mutate(
    date = start + (end - start) / 2,
    .after = end
  ) |>
  ggplot(
    aes(
      x = date,
      y = value,
      colour = factor(cell_id)
    )
  ) +
  geom_line() +
  facet_wrap(
    ~variable,
    scales = "free_y"
  ) +
  theme_minimal()


# To do:

# Build wrapper functions to:

# Create processing tiles - FUNCTION DONE

# Within each tile:

  # download and format all terraclimate data for the tile 2000-2025 - FUNCTION DONE

  # For each pixel in the tile:

    # run spline interpolation to get daily outdoor data 2000-2025

    # run NicheMapR to get hourly microclimate data 2000-2025

    # run the population dynamic models to get hourly population data 2000-2025

    # summarise the population dynamic outputs to monthly data 2000-2025

  # write the monthly summaries for this tile to disk as a CSV file

# load the template raster and CSVs to create (~300) monthly GeoTIFFs of monthly
# data, and 12 layers of synoptic values.




# also:

# need to add in the other microclimate parameters (stone substrate type, etc)
# to microclimate simulation

# add a wind shear exponent adjustment to the microclimate simulation set up?




# Note: we could use GPM IMERG remotely-sensed daily 12km precipitation data,
# rather than Terraclimate 5km (downscaled from CRU TS4.0 reanalysis of weather
# station data), as Tas did for her model. But those data have not been
# downscaled, and the stochastic nature of the rainfall might cause convergence
# issues with the population dynamics simulation (even on an hourly timestep).
# Given our aim is to capture broad-scale seasonality and spatial variation in
# climatic suitability, the monthly, but spatially downscaled, terraclimate
# data. This is also much easier to process, since there's no open THREDDS cube
# interface, and the daily layers are stored as single-day grids (on portals and
# MAP's GeoTIFF library)

# imerg_dir <- "/mnt/s3/mastergrids/Other_Global_Covariates/Rainfall/GPMM_IMerg_Daily/v07B_Total/12km"
# files <- list.files(imerg_dir, full.names = TRUE)
# file.size(files[1])
# # imerg <- rast(files)
#
# read_daily_rain <- function(year, day_of_year, res = c("5km", "12km")) {
#   res <- match.arg(res)
#   ddd <- sprintf("%03d", day_of_year)
#   rr  <- rast(
#     sprintf("/mnt/s3/mastergrids/Other_Global_Covariates/Rainfall/GPMM_IMerg_Daily/v07B_Total/%s/GPMM-IMerg-V07B-MM_Total.%d.%s.Data.%s.%s.tif",
#             res,
#             max(2001, year),
#             ddd,
#             res,
#             ifelse(res == "5km", "NN", "Data"))
#   )
#   rr
# }
#
# library(terra)
# rain_5 <- read_daily_rain(2024, 179)
# system.time(
#   # rain <- read_daily_rain(2024, 179)
#   rain <- read_daily_rain(2024, 179, "12km")
# )
# system.time(
#   res <- terra::extract(rain_12, cbind(longitude, latitude))
# )
