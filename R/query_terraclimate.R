# read in terraclimate data
library(ncdf4)
library(lubridate)

# get the URL to the terraclimate data
terraclimate_url <- function(variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad", "soil")) {
  variable <- match.arg(variable)
  base_url <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_%s_1958_CurrentYear_GLOBE.nc#fillmismatch"
  sprintf(base_url, variable)
}

# open the connection to terraclimate
terraclimate_open <- function(variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad", "soil")) {
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
    variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad", "soil")) {
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

# find the nearest lat and long in terra climate the user's. Check the requested
# coordinate agains the available set. If there isn't one within 1/48 of a
# degree, get one slightly further (basing this on NicheMapR's logic)
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
monthly_summary_days <- function(dates) {
  # find the start of month immediately before the first date
  min_date <- lubridate::floor_date(min(dates), unit = "month")
  # find the start of month immediately after the last date
  max_date <- lubridate::ceiling_date(max(dates), unit = "month")
  # roll back to the ends of the previous months (how terraclimate summarises)
  seq(min_date, max_date, by = "1 month") - 1
}

# error nicely if the user provides inaccessible dates
check_dates <- function(dates, terraclimate_available) {

  check_contiguous_dates(dates)

  # get the summary dates available and the time period they cover (the start of
  # the first date summarised to the end of the last)
  available_dates <- days_since_1900()[terraclimate_available$times]
  coverage_min <- lubridate::floor_date(min(available_dates),
                                        unit = "month")
  coverage_max <- max(available_dates)

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

  # the timeseries is summarised by month and indexed by the last day of the
  # month, so find all summary dates relevant to our target dates
  required_summary_dates <- monthly_summary_days(dates)

  # the terraclimate temporal index is integer days since the start of 1900 so
  # find our month-end dates as indices
  time_index_1900 <- match(required_summary_dates, days_since_1900())

  # I'm insecure
  if (!all(time_index_1900 %in% available$times)) {
    stop("Some date error. It's Nick's fault, please tell him.")
  }

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
      start = lubridate::floor_date(required_summary_dates, unit = "1 month"),
      end = required_summary_dates
    )
  )

}

# Get weather for central Perth WA for 2020-2023
longitude <- 115.86
latitude <- -31.95
dates <- seq(as.Date("2020-01-13"),
             as.Date("2021-07-20"),
             by = "1 day")

# build the indexing slice
slice <- terraclimate_build_slice(longitude = longitude,
                                  latitude = latitude,
                                  dates = dates)

# get all variables
vars <- c("tmax", "tmin", "ppt", "ws", "vpd", "srad", "soil")
results <- list()
for(var in vars) {
  results[[var]] <- terraclimate_fetch(slice, var)
}

climate_data <- cbind(
  data.frame(
    start_date = slice$dates$start,
    end_date = slice$dates$end
  ),
  as.data.frame(results)
)

climate_data$date <- with(climate_data,
                          start_date + round((end_date - start_date) / 2))
plot(tmax ~ date,
     data = climate_data,
     ylab = "temperature",
     ylim = c(0, max(tmax)),
     pch = 16,
     col = "red")
points(tmin ~ date,
       data = climate_data,
       pch = 16,
       col = "blue")

# come back and build a version that pulls these out for blocks of pixels at a
# time, and stores them locally


