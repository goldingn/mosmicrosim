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
