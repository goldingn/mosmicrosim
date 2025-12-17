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
# default 1pm) can be specified as integers in military time (600 is 6am, 1300
# is 1pm, 2400 is midnight, etc.).
interpolate_daily_air_temp <- function(daily_max, daily_min,
                                       time_sunrise = 600,
                                       time_sunset = 1800,
                                       time_daily_max = 1300) {

  # rename inputs to scheme
  temp_max_day <- daily_max
  temp_min_day <- daily_min

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
    time = times,
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
# points(max_temps ~ which(res$time == 1300),
#        pch = 16, col = "red")
# points(min_temps ~ which(res$time == 600),
#        pch = 16, col = "blue")

# Hourly interpolation of relative humidity from vectors of daily maxima
# (daily_max) and minima (daily_min). This is based on the nichemapper VSINE
# routine, with linear interpolation from midnight (when humidity is halfway
# between the previous day's daily maximum and min) to sunrise (when humidity is
# at its daily maximum) to the timing of daily minimum. back to midnight. The
# timing of sunrise (time_sunrise, default 6am), sunset (time_sunset, default
# 6pm), and the timing of the daily minimum humidity (time_daily_min, default
# 1pm) can be specified in hours (600 is 6am, 1300 is 1pm, 2400 is midnight,
# etc.).
interpolate_daily_humidity <- function(daily_max, daily_min,
                                       time_sunrise = 600,
                                       time_sunset = 1800,
                                       time_daily_min = 1300) {

  # rename inputs to scheme
  humid_max_day <- daily_max
  humid_min_day <- daily_min

  # set up times for interpolation
  n_days <- length(humid_max_day)
  times <- rep(100 * seq_len(24), n_days)

  # expand min and diff vectors to hours, for vectorised lookup
  humid_min_hour <- rep(humid_min_day, each = 24)
  humid_max_hour <- rep(humid_max_day, each = 24)

  # make an unscaled curve from 0 to 1, hitting 0.5 at midnight. Export this to
  # another function

  # midnight to sunrise, 0.5 to 1
  mask_mn_to_sr <- as.numeric(times <= time_sunrise)
  rel_mn_to_sr <- (time_sunrise - times) / time_sunrise
  weight_mn_to_sr <- (1 - 0.5 * pmax(0, rel_mn_to_sr)) * mask_mn_to_sr

  # sunrise to min, 1 to 0
  mask_sr_to_min <- as.numeric(times > time_sunrise & times <= time_daily_min)
  rel_sr_to_min <- (times - time_sunrise) / (time_daily_min - time_sunrise)
  weight_sr_to_min <- (1 - pmin(1, pmax(0, rel_sr_to_min))) * mask_sr_to_min

  # min to midnight, 0 to 0.5
  mask_min_to_mn <- as.numeric(times > time_daily_min)
  rel_min_to_mn <- (times - time_daily_min) / (2400 - time_daily_min)
  weight_min_to_mn <- 0.5 * pmin(1, pmax(0, rel_min_to_mn)) * mask_min_to_mn

  # combined weights
  weight <- weight_mn_to_sr + weight_sr_to_min + weight_min_to_mn


  # there are jumps between successive days' midnight values, so blend them
  # between the two by shifting the max and min values on hours

  # for each midnight, the minimum should be the previous day's and the maximum
  # the next day's


  # midnight to sunrise (peak)
  #   should be previous day's min, this day's max

  # sunrise to daily min (trough)
  #   should be that day's min and max

  # daily min to midnight
  #   should be this day's min, next day's max


  # this day's max goes from this day's sunrise to next day's sunrise

  # for each day, if we are before sunrise (peak) use yesterday's min, otherwise
  # use today's

  # get yesterday's min
  humid_min_yesterday_day <- c(humid_min_day[1], humid_min_day[-n_days])
  humid_min_yesterday_hour <- rep(humid_min_yesterday_day, each = 24)

  # find out if we are before sunrise
  before_sunrise <- as.numeric(times < time_sunrise)

  # reassign these
  humid_min_last_hour <- humid_min_hour * (1 - before_sunrise) +
    humid_min_yesterday_hour * before_sunrise


  # this day's max goes from the last day's min time to this day's min time

  # for each day, if we are before time_daily_min (trough) use today's max,
  # otherwise use tomorrow's max

  # get tomorrow's max
  humid_max_tomorrow_day <- c(humid_max_day[-1], humid_max_day[n_days])
  humid_max_tomorrow_hour <- rep(humid_max_tomorrow_day, each = 24)

  # find out if we are before sunrise
  after_daily_min <- as.numeric(times > time_daily_min)

  # reassign these
  humid_max_next_hour <- humid_max_hour * (1 - after_daily_min) +
    humid_max_tomorrow_hour * after_daily_min

  humid_diff_hour <- humid_max_next_hour - humid_min_last_hour

  # scale these according to to max and min
  humid_hour <- humid_min_last_hour + weight * humid_diff_hour

  # return the humidity and times
  data.frame(
    time = times,
    humidity = humid_hour
  )

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
# points(max_humids ~ which(res$time == 600),
#        pch = 16, col = "red")
# points(min_humids ~ which(res$time == 1300),
#        pch = 16, col = "blue")
