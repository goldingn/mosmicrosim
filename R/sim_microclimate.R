# functions for modelling microclimates at a sub-month resolution, to inform
# population dynamics

# given a tibble of monthly climate data for NicheMapR (as created by
# process_terraclimate_tile_vars()) return a
# tibble of daily climate data interpolated for all the days in those dates those dates
daily_from_monthly_climate <- function(monthly_climate) {

  # all dates to predict to
  dates_predict <- seq(from = min(monthly_climate$start),
                       to = max(monthly_climate$end),
                       by = 1)

  monthly_climate |>

    # maybe transform daily rainfall for splining, then convert back later
    dplyr::mutate(
      rainfall_daily = log1p(rainfall_daily)
    ) |>

    tidyr::pivot_longer(
      cols = !any_of(c("mid_date", "start", "end")),
      names_to = "variable",
      values_to = "value"
    ) |>
    # for each variable, run spline_seasonal to interpolate to daily data
    dplyr::group_by(
      variable
    ) |>
    dplyr::summarise(
      value = list(spline_seasonal(values = value,
                                   dates = mid_date,
                                   dates_predict = dates_predict)),
      date = list(dates_predict),
      .groups = "drop"
    ) |>
    tidyr::unnest(
      c(date, value)
    ) |>
    tidyr::pivot_wider(
      names_from = variable,
      values_from = value
    ) |>
    # transform back the daily rainfall, and force to be non-negative
    dplyr::mutate(
      rainfall_daily = pmax(0, expm1(rainfall_daily))
    )

}


# given a tibble of (spline-interpolated) daily climate data for a given
# location, the location's altitude, and some microclimate parameters, run
# NicheMapR's microclimate model to model hourly microclimate conditions at that
# location
hourly_from_daily_climate <- function(latitude,
                                            longitude,
                                            altitude,
                                            daily_climate,
                                            microclimate) {

  micro <- create_micro(
    latitude = latitude,
    longitude = longitude,
    altitude = altitude,
    dates = daily_climate$date,
    daily_temp_max_c = daily_climate$tmax,
    daily_temp_min_c = daily_climate$tmin,
    daily_rh_max_perc = daily_climate$rhmax,
    daily_rh_min_perc = daily_climate$rhmin,
    daily_cloud_max_perc = daily_climate$ccmax,
    daily_cloud_min_perc = daily_climate$ccmin,
    daily_wind_max_ms = daily_climate$wsmax,
    daily_wind_min_ms = daily_climate$wsmin,
    daily_rainfall_mm = daily_climate$rainfall_daily,
    adult_height_m = microclimate$adult_height,
    shade_prop = microclimate$shade_proportion)

  # run NicheMapR microclimate model on it
  sim <- NicheMapR::microclimate(micro)

  # get hourly rainfall (assume even, so just daily divided by 24)
  hourly_rainfall <- dplyr::tibble(
    date = daily_climate$date,
    rainfall = daily_climate$rainfall_daily / 24
  )

  # get the hourly micro (shaded) habitat results, and pull out useful variables
  hourly_climate <- dplyr::tibble(
    date = rep(daily_climate$date, each = 24),
    hour = rep(0:23, length(daily_climate$date)),
    air_temperature = sim$shadmet[, "TALOC"],
    humidity = sim$shadmet[, "RHLOC"],
    water_temperature = sim$shadsoil[, "D0cm"],
    windspeed = sim$shadmet[, "VLOC"]
  ) |>
    # append rainfall
    dplyr::left_join(
      hourly_rainfall,
      by = "date"
    )

  hourly_climate

}



# Hourly interpolation of air temperatures from vector on daily maxima
# (daily_max) and minima (daily_min). This is based on the nichemapper SINEC
# routine, employing a sine curve for temperature from sunrise (when the
# temperature is at the daily minimum) to sunset each day, an exponential decay
# from sunset to midnight, and a linear trend from midnight to the next sunrise.
# The timing of sunrise (time_sunrise, default 6am), sunset (time_sunset,
# default 6pm), and the timing of the daily max temperature (time_daily_max,
# default 1pm) can be specified as units of hours (6 is 6am, 13
# is 1pm, 24 is midnight, etc.).
interpolate_daily_air_temp <- function(daily_max, daily_min,
                                       time_sunrise = 6,
                                       time_sunset = 18,
                                       time_daily_max = 13) {

  # rename inputs to scheme
  temp_max_day <- daily_max
  temp_min_day <- daily_min

  # multiply times by 100 to match nichemapper units and parameter choices
  time_sunrise <- time_sunrise * 100
  time_sunset <- time_sunset * 100
  time_daily_max <- time_daily_max * 100

  # set up times for interpolation
  n_days <- length(temp_max_day)
  times <- rep(100 * seq_len(24), n_days)

  # compute a reference time, half way from sunrise to max temperature
  time_reference <- (time_daily_max - time_sunrise) / 2 + time_sunrise

  # parameter of overnight exponential decay
  tau <- 3 / (2400 - time_sunset + time_sunrise)

  # Z parameter for overnight exponential decay
  Z_night <- sin((360 * (time_sunset - time_reference) / (2 * (time_daily_max - time_sunrise))) / 57.29577)

  # tomorrow's minimum (repeating the last one)
  temp_min_tomorrow_day <- c(temp_min_day[-1], temp_min_day[n_days])

  # compute the A parameter (half the temperature range, NOT the average temperature
  # as stated in Mike's doc)
  A_day <- (temp_max_day - temp_min_day) / 2

  # compute each day's sunset temperature
  temp_sunset_day <- A_day * Z_night + temp_min_day + A_day

  # get yesterday's temp_sunset (repeating the first one)
  temp_sunset_yesterday_day <- c(temp_sunset_day[1], temp_sunset_day[-n_days])

  # expand vectors to hours, for vectorised lookup
  temp_max_hour <- rep(temp_max_day, each = 24)
  temp_min_hour <- rep(temp_min_day, each = 24)
  A_hour <- rep(A_day, each = 24)
  temp_sunset_hour <- rep(temp_sunset_day, each = 24)
  temp_min_tomorrow_hour <- rep(temp_min_tomorrow_day, each = 24)
  temp_sunset_yesterday_hour <- rep(temp_sunset_yesterday_day, each = 24)

  # for each hour, compute the *next minimum temperature*, which is the same day's
  # minimum (before or at sunrise) or the next day's minimum (after sunrise)
  after_sunrise <- as.numeric(times > time_sunrise)
  temp_min_next_hour <- temp_min_hour * (1 - after_sunrise) +
    temp_min_tomorrow_hour * after_sunrise

  # for each hour, compute the *last sunset temperature*, which is the same day's
  # sunset temperature (after or at sunset) or the previous day's sunset
  # temperature (before sunset)
  before_sunset <- as.numeric(times < time_sunset)
  temp_sunset_last_hour <- temp_sunset_hour * (1 - before_sunset) +
    temp_sunset_yesterday_hour * before_sunset

  # mask for daylight vs nighttime hours
  is_light <- as.numeric(times >= time_sunrise & times <= time_sunset)

  # now we can start

  # make a blank vector to fill
  temp_hour <- rep(NA, n_days * 24)

  # for daylight hours, use the sine curve
  Z_hour_light <- sin((360 * (times - time_reference) / (2 * (time_daily_max - time_sunrise))) / 57.29577)
  temp_hour_light <- A_hour * Z_hour_light + temp_min_hour + A_hour

  # for nighttime, use exponential decay from the last sunset temperature to the
  # next sunrise
  time_since_ss <- times - time_sunset + 2400 * (times < 1800)
  E <- tau * time_since_ss
  temp_hour_night = (temp_sunset_last_hour - temp_min_next_hour) / exp(E) + temp_min_next_hour

  # combine these
  temp_hour <- temp_hour_light * is_light + temp_hour_night * (1 - is_light)

  # now apply linear interpolation from each midnight to sunrise, to remove
  # pre-dawn dips

  # get the previous day's midnight temperature and apply to hours
  temp_midnight_day <- temp_hour[times == 2400]
  temp_midnight_yesterday_day <- c(temp_midnight_day[1], temp_midnight_day[-n_days])
  temp_midnight_yesterday_hour <- rep(temp_midnight_yesterday_day, each = 24)

  temp_hour_predawn <- temp_min_hour + (temp_midnight_yesterday_hour - temp_min_hour) * (1 - times / time_sunrise)
  is_predawn <- times < time_sunrise
  temp_hour[is_predawn] <- temp_hour_predawn[is_predawn]

  # return the temperature and times
  data.frame(
    time = times / 100,
    temperature = temp_hour
  )

}

# # demo
# n_days <- 12
# max_temps <- runif(n_days, 25, 32)
# min_temps <- runif(n_days, 11, 15)
#
# res <- interpolate_daily_air_temp(
#   daily_max = max_temps,
#   daily_min = min_temps
# )
#
# # plot interpolation and maxima/minima
# plot(res$temperature,
#      type = "l")
# points(max_temps ~ which(res$time == 13),
#        pch = 16, col = "red")
# points(min_temps ~ which(res$time == 6),
#        pch = 16, col = "blue")




# Hourly interpolation of relative humidity from vectors of daily maxima
# (daily_max) and minima (daily_min). This is based on the nichemapper VSINE
# routine, with linear interpolation from midnight (when humidity is halfway
# between the previous day's daily maximum and this day's minimum) to sunrise
# (when humidity is at its daily maximum) to the timing of daily minimum. back
# to midnight. The timing of sunrise (time_sunrise, default 6am) and of the
# daily minimum humidity (time_daily_min, default 1pm) can be specified in units
# of  hours (6 is 6am, 13 is 1pm, 24 is midnight, etc.).
interpolate_daily_humidity <- function(daily_max, daily_min,
                                       time_sunrise = 6,
                                       time_daily_min = 13) {

  result <- vsine_interpolate_daily(daily_max = daily_max,
                                    daily_min = daily_min,
                                    time_sunrise = time_sunrise,
                                    time_inflection = time_daily_min,
                                    sunrise_value = "max")
  names(result) <- c("time", "humidity")
  result

}

# # demo
# n_days <- 6
# max_humids <- runif(n_days, 85, 92)
# min_humids <- runif(n_days, 51, 55)
#
# res <- interpolate_daily_humidity(
#   daily_max = max_humids,
#   daily_min = min_humids
# )
#
# # plot interpolation and maxima/minima
# plot(res$humidity,
#      type = "l")
# points(max_humids ~ which(res$time == 6),
#        pch = 16, col = "red")
# points(min_humids ~ which(res$time == 13),
#        pch = 16, col = "blue")

# Hourly interpolation of wind speed from vectors of daily maxima (daily_max)
# and minima (daily_min). This is based on the nichemapper VSINE routine, with
# linear interpolation from midnight (when wind speed is halfway between the
# previous day's daily maximum and this day's daily minimum) to sunrise (when
# wind speed is at its daily minimum) to the timing of daily maximum, back to
# midnight. The timing of sunrise (time_sunrise, default 6am) and of the daily
# maximum windspeed (time_daily_max, default 1pm) can be specified in units of
# hours (6 is 6am, 13 is 1pm, 24 is midnight, etc.).
interpolate_daily_windspeed <- function(daily_max, daily_min,
                                       time_sunrise = 6,
                                       time_daily_max = 13) {

  result <- vsine_interpolate_daily(daily_max = daily_max,
                                    daily_min = daily_min,
                                    time_sunrise = time_sunrise,
                                    time_inflection = time_daily_max,
                                    sunrise_value = "min")
  names(result) <- c("time", "windspeed")
  result

}

# # demo
# n_days <- 6
# max_wind <- runif(n_days, 3, 5)
# min_wind <- runif(n_days, 1, 1.5)
#
# res <- interpolate_daily_windspeed(
#   daily_max = max_wind,
#   daily_min = min_wind
# )
#
# # plot interpolation and maxima/minima
# plot(res$windspeed,
#      type = "l")
# points(max_wind ~ which(res$time == 13),
#        pch = 16, col = "red")
# points(min_wind ~ which(res$time == 6),
#        pch = 16, col = "blue")

# Hourly interpolation of cloud cover from vectors of daily maxima (daily_max)
# and minima (daily_min). This is based on the nichemapper VSINE routine, with
# linear interpolation from midnight (when cloud cover is halfway between the
# previous day's daily minimum and this day's daily maximum) to sunrise (when
# cloud cover is at its daily maximum) to the timing of daily minimum, back to
# midnight. The timing of sunrise (time_sunrise, default 6am) and of the daily
# minimum cloud cover (time_daily_max, default 1pm) can be specified in units of
# hours (6 is 6am, 13 is 1pm, 24 is midnight, etc.).
interpolate_daily_cloudcover <- function(daily_max, daily_min,
                                         time_sunrise = 6,
                                         time_daily_min = 13) {

  result <- vsine_interpolate_daily(daily_max = daily_max,
                                    daily_min = daily_min,
                                    time_sunrise = time_sunrise,
                                    time_inflection = time_daily_min,
                                    sunrise_value = "max")
  names(result) <- c("time", "cloudcover")
  result

}

# # demo
# n_days <- 6
# max_cloud <- runif(n_days, 80, 95)
# min_cloud <- runif(n_days, 10, 15)
#
# res <- interpolate_daily_cloudcover(
#   daily_max = max_cloud,
#   daily_min = min_cloud
# )
#
# # plot interpolation and maxima/minima
# plot(res$cloudcover,
#      type = "l")
# points(max_cloud ~ which(res$time == 6),
#        pch = 16, col = "red")
# points(min_cloud ~ which(res$time == 13),
#        pch = 16, col = "blue")

# given a vector of times, the time of sunrise (max humidity, min wind), and the
# other inflection point (min humidity, max wind), all in units of hours,
# return a corresponding vector of weights for the VSINE routine
vsine_weights <- function(times, time_sunrise = 6, time_inflection = 13) {

  # make an unscaled curve from 0 to 1, hitting 0.5 at midnight.

  # midnight to sunrise, 0.5 to 1
  mask_mn_to_sr <- as.numeric(times <= time_sunrise)
  rel_mn_to_sr <- (time_sunrise - times) / time_sunrise
  weight_mn_to_sr <- (1 - 0.5 * pmax(0, rel_mn_to_sr)) * mask_mn_to_sr

  # sunrise to inflection, 1 to 0
  mask_sr_to_inf <- as.numeric(times > time_sunrise & times <= time_inflection)
  rel_sr_to_inf <- (times - time_sunrise) / (time_inflection - time_sunrise)
  weight_sr_to_inf <- (1 - pmin(1, pmax(0, rel_sr_to_inf))) * mask_sr_to_inf

  # min to midnight, 0 to 0.5
  mask_inf_to_mn <- as.numeric(times > time_inflection)
  rel_inf_to_mn <- (times - time_inflection) / (24 - time_inflection)
  weight_inf_to_mn <- 0.5 * pmin(1, pmax(0, rel_inf_to_mn)) * mask_inf_to_mn

  # combined weights
  weight <- weight_mn_to_sr + weight_sr_to_inf + weight_inf_to_mn
  weight

}

# Hourly interpolation of climate variables from vectors of daily maxima
# (daily_max) and minima (daily_min). This is based on the nichemapper VSINE
# routine, with linear interpolation between: midnight, sunrise, and an
# inflection point around the time air temperature peaks. This function can
# either model the maximum at sunrise (sunrise_value = "max") and the minimum at the inflection point,
# or vice versa (sunrise_value = "min"). In either case, at midnight, the variable is halfway between
# the previous day's inflection point value and this day's sunrise value. The
# timing of sunrise (time_sunrise, default 6am) and of the inflection point
# (time_inflection, default 1pm) can be specified in units of hours (6
# is 6am, 13 is 1pm, 24 is midnight, etc.).
vsine_interpolate_daily <- function(daily_max, daily_min,
                                    time_sunrise = 6,
                                    time_inflection = 13,
                                    sunrise_value = c("max", "min")) {

  # get the sunrise value to decide the direction
  sunrise_value <- match.arg(sunrise_value)


  # some duplicate (actually, inverted) code here from the humidity function
  # that could be abstracted into a helper function

  # rename inputs to match naming scheme
  max_day <- daily_max
  min_day <- daily_min

  # set up times for interpolation
  n_days <- length(max_day)
  times <- rep(seq_len(24), n_days)

  # expand min and max vectors to hours, for vectorised lookup
  max_hour <- rep(max_day, each = 24)
  min_hour <- rep(min_day, each = 24)

  # make an unscaled curve of weights from 1 (time_daily_max) to 0
  # (time_sunrise), passing through 0.5 at midnight
  weight <- vsine_weights(times = times,
                          time_sunrise = time_sunrise,
                          time_inflection = time_inflection)

  # if the variable is at its minimum at sunrise, negate the weights to have
  # value '1' then and value '0' at time_inflection
  if (sunrise_value == "min") {
    weight <- 1 - weight
  }

  # now adjust the definition of max and min across times to ensure smooth
  # continuity in linear interpolations. We need to switch around the min and
  # max values depending on which peaks at sunrise

  # if the variable is at its minimum at sunrise
  if (sunrise_value == "min") {
    sunval_day <- min_day
    infval_day <- max_day
  } else {
    sunval_day <- max_day
    infval_day <- min_day
  }

  # for each day, if we are before sunrise use yesterday's inflection value,
  # otherwise use today's

  # repeat today's value at inflection across hours
  infval_hour <- rep(infval_day, each = 24)

  # get yesterday's value at inflection
  infval_yesterday_day <- c(infval_day[1], infval_day[-n_days])
  infval_yesterday_hour <- rep(infval_yesterday_day, each = 24)

  # find out if we are before sunrise
  before_sunrise <- as.numeric(times < time_sunrise)

  # reassign these
  infval_last_hour <- infval_hour * (1 - before_sunrise) +
    infval_yesterday_hour * before_sunrise

  # now the same for the value at sunrise

  # for each day, if we are before time_daily_max (peak) use today's min,
  # otherwise use tomorrow's min

  # repeat today's value at sunrise across hours
  sunval_hour <- rep(sunval_day, each = 24)

  # get tomorrow's sunrise values
  sunval_tomorrow_day <- c(sunval_day[-1], sunval_day[n_days])
  sunval_tomorrow_hour <- rep(sunval_tomorrow_day, each = 24)

  # find out if we are after inflection
  after_daily_max <- as.numeric(times > time_daily_max)

  # reassign these
  sunval_next_hour <- sunval_hour * (1 - after_daily_max) +
    sunval_tomorrow_hour * after_daily_max

  # now scale these according to the weights curve
  if (sunrise_value == "min") {
    diff_hour <- infval_last_hour - sunval_next_hour
    value_hour <- sunval_next_hour + weight * diff_hour
  } else {
    diff_hour <- sunval_next_hour - infval_last_hour
    value_hour <- infval_last_hour + weight * diff_hour
  }

  # return the values and times
  data.frame(
    time = times,
    value = value_hour
  )

}

# # demo
# set.seed(4)
# n_days <- 3
# max <- runif(n_days, 80, 95)
# min <- runif(n_days, 10, 15)
#
# res <- vsine_interpolate_daily(
#   daily_max = max,
#   daily_min = min,
#   sunrise_value = "max"
# )
#
# # plot interpolation and maxima/minima
# plot(res$value,
#      type = "l")
# points(max ~ which(res$time == 6),
#        pch = 16, col = "red")
# points(min ~ which(res$time == 13),
#        pch = 16, col = "blue")


