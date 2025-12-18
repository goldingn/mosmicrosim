# given climate data 'values' at dates 'dates', interpolate values at
# 'dates_predict'. Where dates_predict are within the timeseries, interpolate
# linearly. Where they fall outside, Fit a GAM with a long-term spline and a
# periodic spline (for within-year seasonality) to predict those. Blend the two
# estimates over the edges, to ensure no jumps
spline_seasonal <- function(values, dates, dates_predict,
                            knots_per_year_trend = 2,
                            knots_per_year_season = 6) {

  min_date <- min(dates)

  # training and prediction data
  df_fit <- data.frame(
    values = values,
    datenum = as.numeric(dates - min_date),
    doy = doy(dates)
  )

  df_predict <- data.frame(
    datenum = as.numeric(dates_predict - min_date),
    doy = doy(dates_predict)
  )

  # do linear interpolation to all dates first
  values_pred <- approx(x = df_fit$datenum,
                        y = df_fit$values,
                        xout = df_predict$datenum)$y

  # is any extrapolation needed?
  extrapolation <- dates_predict > max(dates) |
    dates_predict < min(dates)

  # if so, fit a spline model
  if (any(extrapolation)) {

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

    # we need to smooth between the linear interpolation and the spline
    # interpolation, so work out where the spline interpolation is being used,
    # and only predict to those elements (gam prediction is costly).

    # Compute weights for linear vs spline interpolation, to smooth between
    # these two estimates. Set this so that predictions between two observed
    # dates neither of which are the first and last observed dates get a weight
    # of 1 (linear interpolation) predictions beyond the first and last observed
    # dates get a weight of 0 (spline interpolation), and predictions in between
    # get weights from 0 to 1, depending on their distance between the
    # penultimate and ultimate observed date.

    # find the key dates
    first_observed <- min(dates)
    last_observed <- max(dates)
    second_observed <- min(dates[dates > first_observed])
    penultimate_observed <- max(dates[dates < last_observed])

    # find the durations of the shoulder periods
    time_first_to_second <- as.numeric(second_observed - first_observed)
    time_penultimate_to_last <- as.numeric(last_observed - penultimate_observed)

    # compute time since first observed and time before last observed
    time_since_first <- as.numeric(dates_predict - first_observed)
    time_before_last <- as.numeric(last_observed - dates_predict)

    # compute these *relative* to the shoulders
    rel_time_since_first <- time_since_first / time_first_to_second
    rel_time_before_last <- time_before_last / time_penultimate_to_last

    # clamp this to be 0 if before first (needs interpolating) and 1
    # if after second
    forward_weight <- pmax(0, pmin(1, rel_time_since_first))

    # clamp this to be 0 if after last (needs interpolating) and 1
    # if before penultimate
    backward_weight <- pmax(0, pmin(1, rel_time_before_last))

    # combine into one weight
    weight <- forward_weight * backward_weight

    # for linear, set any extrapolations to 0 (instead of NA)
    values_pred_linear <- values_pred
    values_pred_linear[extrapolation] <- 0

    # predict the gam only where the weight on the spline is nonzero
    values_pred_spline <- values_pred_linear * 0
    to_insert <- (1 - weight) > 0
    values_pred_spline[to_insert] <- as.numeric(predict(m, df_predict[to_insert, ]))

    # combine the two methods to get the new values
    values_pred <- values_pred_linear * weight +
      values_pred_spline * (1 - weight)

  }

  # return the interpolated values
  values_pred

}
