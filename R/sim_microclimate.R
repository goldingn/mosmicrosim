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
