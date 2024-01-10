# create a 'micro' file for input to NicheMapR::microclimate(), given vectors of
# daily macroclimate data

#' @param latitude,longitude,altitude_m scalars for the coordinates and altitude
#'   at these coordinates, to compute solar radiation and air pressure
#' @param adult_height_m scalar height above ground experienced by adult
#'   mosquitoes in metres
#' @param shade_prop proportion of full shade experienced by mosquitoes
#' @param even_rain whether to spread daily rainfall out evenly over 24. If
#'   FALSE, have it all fall at midnight.
#' @param dates vector of dates corresponding to data vectors - these must be
#'   the dates in Greenwich Mean Time (so that longitude can be used to compute
#'   solar noon)
#' @param daily_temp_max_c vector of daily maximum temperatures in centigrade
#' @param daily_temp_min_c vector of daily minimum temperatures in centigrade
#' @param daily_rainfall_mm vector of daily rainfall in millimetres
create_micro <- function(
    latitude,
    longitude,
    altitude_m,
    adult_height_m = 1,
    shade_prop = 0,
    even_rain = TRUE,
    dates,
    daily_temp_max_c,
    daily_temp_min_c,
    daily_rainfall_mm) {

  # check the dates are contiguous and integer days
  contiguous_dates <- identical(dates,
                                seq(dates[1],
                                    dates[length(dates)],
                                    by = 1))
  if (!contiguous_dates) {
    stop("dates must be contiguous, increasing, and non-repeating ",
         "(a vector of consecutive days)")
  }

  n_days <- length(dates)

  # check the data vectors are of the same length as the dates
  good_vectors <- length(daily_temp_max_c) == n_days &
    length(daily_temp_min_c) == n_days &
    length(daily_rainfall_mm) == n_days

  if (!good_vectors) {
    stop("all temperature and rainfall vectors must have the same number of ",
         "elements as dates")
  }

  # set up the micro objects

  # microinput is a hideous vector of doubles controlling behaviour of the
  # fortran code, create it from defaults and these args
  microinput <- create_microinput(
    n_days = n_days,
    latitude = latitude,
    longitude = longitude,
    altitude_m = altitude_m,
    shade_prop = shade_prop,
    adult_height_m = adult_height_m,
    mean_temp_c = mean(c(daily_temp_max_c, daily_temp_min_c)),
    even_rain = even_rain)

  # set the shade to a percentage, and only run once
  minshade <- maxshade <- 100 * shade_prop

  stop("not yet implemented")
  micro <- list(microinput = microinput,
                ...)

  micro
}

# create the vector of control arguments for the fortran code
create_microinput <- function(
    # number of days we are running for
    n_days,
    # coordinates
    latitude,
    longitude,
    # Elevation, if to be user specified (m)
    altitude_m,
    # shade proportion
    shade_prop,
    # Local height (m) at which air temperature, wind speed and humidity are to
    # be computed for organism of interest
    adult_height_m,
    # average temperature over the timeseries, to compute the deep soil
    # temperature
    mean_temp_c,
    even_rain = TRUE) {

  # find the hemisphere (1 is north)
  hemisphere <- ifelse(latitude < 0, 2, 1)

  # absolute degrees and minutes in this hemisphere
  lat_degree <- abs(trunc(latitude))
  lat_minute <- 60 * (abs(latitude) - lat_degree)
  lon_degree <- abs(trunc(longitude))
  lon_minute <- 60 * (abs(longitude) - lon_degree)

  # define horizon angles (assume no buildings/hills, so 0) and convert to
  # radians to calc view factors
  hori <- rep(0, 24)
  view_factor <- 1 - mean(sin(hori * pi / 180))

  # return a horrifying vector - order matters, names do not :(
  microinput <- c(ndays = n_days,
                  RUF = microinput_defaults$RUF,
                  ERR = microinput_defaults$ERR,
                  Usrhyt = adult_height_m,
                  Refhyt = microinput_defaults$Refhyt,
                  Numtyps = microinput_defaults$Numtyps,
                  Z01 = microinput_defaults$Z01,
                  Z02 = microinput_defaults$Z02,
                  ZH1 = microinput_defaults$ZH1,
                  ZH2 = microinput_defaults$ZH2,
                  # index for start day and number of days to run for
                  idayst = 1,
                  ida = n_days,
                  HEMIS = hemisphere,
                  # decimal degree lat/long in degrees and minutes
                  ALAT = lat_degree,
                  AMINUT = lat_minute,
                  ALONG = lon_degree,
                  ALMINT = lon_minute,
                  # set the solar hour angle to the longitude (dates in GMT, so
                  # 0 at noon on meridian)
                  ALREF = lon_degree,
                  slope = microinput_defaults$slope,
                  azmuth = microinput_defaults$aspect,
                  ALTT = altitude_m,
                  CMH2O = microinput_defaults$CMH2O,
                  microdaily = microinput_defaults$microdaily,
                  tannul = mean_temp_c,
                  EC = microinput_defaults$EC,
                  VIEWF = view_factor,
                  snowtemp = microinput_defaults$snowtemp,
                  snowdens = microinput_defaults$snowdens,
                  snowmelt = microinput_defaults$snowmelt,
                  undercatch = microinput_defaults$undercatch,
                  rainmult = microinput_defaults$rainmult,
                  # don't run multiple times with different shades, just run it
                  # at the minimum shade level
                  runshade = 0,
                  # don't run the soil moisture model
                  runmoist = 0,
                  maxpool = microinput_defaults$maxpool,
                  evenrain = ifelse(even_rain, 1, 0),
                  # don't run the snow model
                  snowmodel = 0,
                  rainmelt = microinput_defaults$rainmelt,
                  writecsv = microinput_defaults$writecsv,
                  densfun = microinput_defaults$densfun,
                  hourly = microinput_defaults$hourly,
                  rainhourly = microinput_defaults$rainhourly,
                  lamb = microinput_defaults$lamb,
                  IUV = microinput_defaults$IUV,
                  RW = microinput_defaults$RW,
                  PC = microinput_defaults$PC,
                  RL = microinput_defaults$RL,
                  SP = microinput_defaults$SP,
                  R1 = microinput_defaults$R1,
                  IM = microinput_defaults$IM,
                  MAXCOUNT = microinput_defaults$MAXCOUNT,
                  IR = microinput_defaults$IR,
                  message = microinput_defaults$message,
                  fail = microinput_defaults$fail,
                  snowcond = microinput_defaults$snowcond,
                  # snow interception fraction for when there's shade (0-1)
                  intercept = shade_prop * 0.3,
                  grasshade = microinput_defaults$grasshade,
                  solonly = microinput_defaults$solonly,
                  ZH = microinput_defaults$ZH,
                  D0 = microinput_defaults$D0,
                  TIMAXS = microinput_defaults$TIMAXS,
                  TIMINS = microinput_defaults$TIMINS,
                  spinup = microinput_defaults$spinup,
                  # these next two undefined and hardcoded in R
                  dewrain = 0,
                  timestep = 360,
                  maxsurf = microinput_defaults$maxsurf)

  microinput

}


# nichemapr microclimate model defaults:
microinput_defaults <- list(
  # Integrator error tolerance for soil temperature calculations
  ERR = 1.5,
  # Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass
  # may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001
  # (snow) - 0.02 m.
  RUF = 0.004,
  # Reference height (m), reference height at which air temperature, wind speed
  # and relative humidity input data are measured
  Refhyt = 2,
  # number of substrate types (depths) for soil modelling
  Numtyps = 2,
  # Top (1st) segment roughness height(m) - IF NO EXPERIMENTAL WIND PROFILE DATA
  # SET THIS TO ZERO! (then RUF and Refhyt used)
  Z01 = 0,
  # 2nd segment roughness height(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET
  # THIS TO ZERO! (then RUF and Refhyt used).
  Z02 = 0,
  # Top of (1st) segment, height above surface(m) - IF NO EXPERIMENTAL WIND
  # PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used).
  ZH1 = 0,
  # 2nd segment, height above surface(m) - IF NO EXPERIMENTAL WIND PROFILE DATA
  # SET THIS TO ZERO! (then RUF and Refhyt used).
  ZH2 = 0,
  # Eccentricity of the earth's orbit (current value 0.0167238, ranges between
  # 0.0034 to 0.058)
  EC = 0.0167238,
  # Slope and aspect in degrees (aspect 0 = north); could compute from DEM, but
  # ignore
  slope = 0,
  aspect = 0,
  # Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air
  # conditions; 2.0 = humid, tropical conditions (note this is for the whole
  # atmospheric profile, not just near the ground)
  CMH2O = 1,
  # run microclimate model where one iteration of each day occurs and last day
  # gives initial conditions for present day with an initial 3 day burn in
  microdaily = 1,
  # Temperature (°C) at which precipitation falls as snow
  snowtemp = 1.5,
  # snow density (mg/m3), overridden by densfun
  snowdens = 0.375,
  # proportion of calculated snowmelt that doesn't refreeze
  snowmelt = 1,
  # undercatch multipier for converting rainfall to snow
  undercatch = 1,
  # Rain multiplier for surface soil moisture (-), used to induce runon
  rainmult = 1,
  # Max depth for water pooling on the surface (mm), to account for runoff
  maxpool = 10000,
  # parameter in equation that melts snow with rainfall as a function of air
  # temp
  rainmelt = 0.0125,
  # Make Fortran code write output as csv files? 1=yes, 0=no
  writecsv = 0,
  # slope and intercept of model of snow density as a linear function of
  # snowpack age if first two values are nonzero, and following the exponential
  # function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394 if all values are
  # non-zero; if it is c(0,0,0,0) then fixed density used
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  # unclear what these do, set to 0 mostly
  hourly = 0,
  rainhourly = 0,
  # Return wavelength-specific solar radiation output?
  lamb = 0,
  # Use gamma function for scattered solar radiation? (computationally
  # intensive))
  IUV = 0,
  # root radius, m
  R1 = 0.001,
  # resistance per unit length of root, m3 kg-1 s-1
  RW = 2.5e+10,
  # resistance per unit length of leaf, m3 kg-1 s-1
  RL = 2e+6,
  # critical leaf water potential for stomatal closure, J kg-1
  PC = -1500,
  # stability parameter for stomatal closure equation, -
  SP = 10,
  # maximum allowable mass balance error, kg
  IM = 1e-06,
  # maximum iterations for mass balance, -
  MAXCOUNT = 500,
  # undefined:
  IR = 0,
  # allow the Fortran integrator to output warnings? (1) or not (0)
  message = 0,
  # how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)
  fail = 24 * 365,
  # effective snow thermal conductivity W/mC (if zero, uses inbuilt function of
  # density)
  snowcond = 0,
  # if 1, means shade is removed when snow is present, because shade is cast by
  # grass/low shrubs
  grasshade = 0,
  # Only run SOLRAD to get solar radiation? 1=yes, 0=no
  solonly = 0,
  # heat transfer roughness height (m) for Campbell and Norman air
  # temperature/wind speed profile (invoked if greater than 0, 0.02 * canopy
  # height in m if unknown)
  ZH = 0,
  # zero plane displacement correction factor (m) for Campbell and Norman air
  # temperature/wind speed profile (0.6 * canopy height in m if unknown)
  D0 = 0,
  # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to
  # solar noon, humidity and cloud cover max's relative to sunrise
  TIMAXS = c(1, 1, 0, 0),
  # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to
  # sunrise, humidity and cloud cover min's relative to solar noon
  TIMINS = c(0, 0, 1, 1),
  # always set to this when initial soil temperatures are not given (uses air
  # temperature to initialise)
  spinup = 1,
  # undefined
  maxsurf = 95)


