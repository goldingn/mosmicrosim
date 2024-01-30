# read in terraclimate data
library(ncdf4)
library(lubridate)
library(mgcv)

source("R/utils.R")

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

# spline these to the requested dates
climate_daily <- data.frame(
  date = dates
)
new_var <- c("tmax", "tmin", "rhmax", "rhmin", "ws", "srad")
for (var in new_var) {
  climate_daily[, var] <- spline_seasonal(values = climate_monthly[, var],
                                          dates = climate_monthly[, "mid_date"],
                                          dates_predict = climate_daily[, "date"])
}

# it would be nice if we could use the correct likelihood for disaggregation
# here, but this is cheap and probably not much different

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
#  cloud cover max/min (based on SRAD modification of existing cloudcover)

# from the helpfile for NicheMapR::get.global.climate() :

# Cloud cover comes from a bilinear interpolation of a lower resolution version
# of this dataset (New, M., M. Hulme, and P. D. Jones. 1999. Representing
# twentieth century space-time climate variability. Part 1: development of a
# 1961-90 mean monthly terrestrial climatology. Journal of Climate 12:829-856.).

# need to do some more transformations and wind min/max as per NicheMapR

# WNMAXX <- WIND
# WNMINN<-WNMAXX * 0.1 # impose diurnal cycle

# Correct for fact that wind is measured at 10m height
# wind shear equation: v / vo = (h / ho)^a
# where
# v = the velocity at height h (m/s)
# vo = the velocity at height ho (m/s)
# a = the wind shear exponent

# Terrain                       Wind Shear Exponent
# -----------------------------|---------------------------
# Open water                    0.1
# Smooth, level, grass-covered  0.15 (or more commonly 1/7)
# Row crops 	                  0.2
# Low bushes with a few trees 	0.2
# Heavy trees 	                0.25
# Several buildings 	          0.25
# Hilly, mountainous terrain 	  0.25

# source http://www.engineeringtoolbox.com/wind-shear-d_1215.html

# WNMINN <- WNMINN * (2 / 10) ^ 0.15
# WNMAXX <- WNMAXX * (2 / 10) ^ 0.15

# add this wind shear exponent to the microclimate simulation set up?


# then build a wrapper function that:

# 1. Downloads terraclimate data in tiles rather than pixels

# 2. Downloads, processes, and saves as monthly daily max and min (then spline
# later) summary data for required quantities

# 3. Caches these, and only downloads new pixels/times as needed (at least by
# month and pixel)



