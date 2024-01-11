#' Create a microclimate config file
#'
#' Create a 'micro' file for input to NicheMapR::microclimate(), given vectors of
#' daily macroclimate data
#'
#' @param latitude,longitude,altitude_m scalars for the coordinates and altitude
#'   at these coordinates, to compute solar radiation and air pressure
#' @param dates vector of dates corresponding to data vectors - these must be
#'   the dates in Greenwich Mean Time (so that longitude can be used to compute
#'   solar noon)
#' @param daily_temp_max_c vector of daily maximum temperature in centigrade
#' @param daily_temp_min_c vector of daily minimum temperature in centigrade
#' @param daily_rh_max_perc vector of daily maximum relative humidity in percent
#' @param daily_rh_min_perc vector of daily minimum relative humidity in percent
#' @param daily_cloud_max_perc vector of daily maximum cloud cover in percent
#' @param daily_cloud_min_perc vector of daily minimum cloud cover in percent
#' @param daily_wind_max_ms vector of daily maximum wind speed in m/s
#' @param daily_wind_min_ms vector of daily minimum wind speed in m/s
#' @param daily_rainfall_mm vector of daily rainfall in millimetres
#' @param weather_height_m scalar height in metres at which the weather data
#'   were recorded
#' @param adult_height_m scalar height above ground experienced by adult
#'   mosquitoes in metres
#' @param shade_prop proportion of full shade experienced by mosquitoes
#' @param even_rain whether to spread daily rainfall out evenly over 24. If
#'   FALSE, have it all fall at midnight.
#' @param soil_reflectance soil solar reflectance, decimal %
#' @param solar_attenuation_data dataframe giving on solar attenuation at the
#'   111 wavelengths for the given location, as created using
#'   `solar_attentuation()`
#' @export
create_micro <- function(
    latitude,
    longitude,
    altitude_m,
    dates,
    daily_temp_max_c,
    daily_temp_min_c,
    daily_rh_max_perc,
    daily_rh_min_perc,
    daily_cloud_max_perc,
    daily_cloud_min_perc,
    daily_wind_max_ms,
    daily_wind_min_ms,
    daily_rainfall_mm,
    weather_height_m = 2,
    adult_height_m = 1,
    shade_prop = 0,
    even_rain = TRUE,
    soil_reflectance = 0.15,
    solar_attenuation_data = solar_attenuation(latitude, longitude)) {

  n_days <- length(dates)

  # check the dates are contiguous, increasing, and integer days
  contiguous_dates <- identical(dates,
                                seq(dates[1],
                                    dates[n_days],
                                    by = 1))
  if (!contiguous_dates) {
    stop("dates must be contiguous, increasing, and non-repeating ",
         "(a vector of consecutive days)")
  }

  # check the data vectors are of the same length as the dates
  vectors <- list(
    daily_temp_max_c,
    daily_temp_min_c,
    daily_rh_max_perc,
    daily_rh_min_perc,
    daily_cloud_max_perc,
    daily_cloud_min_perc,
    daily_wind_max_ms,
    daily_wind_min_ms,
    daily_rainfall_mm
  )

  good_vectors <- vapply(vectors,
                         function(x) {length(x) == n_days},
                         FUN.VALUE = logical(1))

  if (!all(good_vectors)) {
    stop("all weather data vectors must be the same length as the vector of ",
         "dates")
  }

  # set up the micro objects

  # mean air temperature, used for deep and initial soil temperatures
  mean_temp_c <- mean(c(daily_temp_max_c, daily_temp_min_c))

  # horizon angles (assume not blocked by hills or buildings)
  hori <- rep(0, 24)

  # microinput is a hideous vector of doubles controlling behaviour of the
  # fortran code. Create it from these args, and a large set of defaults
  microinput <- create_microinput(
    n_days = n_days,
    latitude = latitude,
    longitude = longitude,
    altitude_m = altitude_m,
    shade_prop = shade_prop,
    adult_height_m = adult_height_m,
    mean_temp_c = mean_temp_c,
    hori = hori,
    even_rain = even_rain,
    weather_height_m = weather_height_m)

  # empty matrix of tides since we don't use them
  tides <- matrix(data = 0,
                  nrow = 24 * n_days,
                  ncol = 3)

  # These vectors for hourly data are all blanks, since we tell the microclimate
  # model to impute them with sine curves from daily max and mins. I'm not sure
  # whether the fortran code is writing over their memory slots or defining new
  # ones, so I'm creating separate objects (with separate memory allocation)
  # instead of copying a blank vector around, to avoid potential memory issues
  # later.
  TAIRhr <- rep(0, 24 * n_days)
  RHhr <- rep(0, 24 * n_days)
  WNhr <- rep(0, 24 * n_days)
  CLDhr <- rep(0, 24 * n_days)
  SOLRhr <- rep(0, 24 * n_days)
  RAINhr <- rep(0, 24 * n_days)
  ZENhr <- rep(-1, 24 * n_days)
  IRDhr <- rep(-1, 24 * n_days)

  # soil data
  soil_objects <- make_soil_data(n_days = n_days,
                                 mean_temp_c = mean_temp_c)

  micro <- list(

    # blank tides data, since we aren't using that
    tides = tides,

    # the microinput control vector created above
    microinput = microinput,

    # day of the year
    doy = doy(dates),

    # substrate longwave IR emissivity (decimal %), typically close to 1
    SLES = rep(0.95, n_days),

    # soil stuff
    DEP = soil_objects$soil_depths_cm,
    Nodes = soil_objects$soil_nodes,

    # set the shade to a percentage and repeat to a vector (both the same as we
    # will only run once)
    MAXSHADES = rep(100 * shade_prop, n_days),
    MINSHADES = rep(100 * shade_prop, n_days),

    # pass in the weather data vectors as maxima and minima
    TMAXX = daily_temp_max_c,
    TMINN = daily_temp_min_c,
    RHMAXX = daily_rh_max_perc,
    RHMINN = daily_rh_min_perc,
    CCMAXX = daily_cloud_max_perc,
    CCMINN = daily_cloud_min_perc,
    WNMAXX = daily_wind_max_ms,
    WNMINN = daily_wind_min_ms,

    # pass in blank vectors for hourly weather data (since the hourly mode is
    # turned off)
    TAIRhr = TAIRhr,
    RHhr = RHhr,
    WNhr = WNhr,
    CLDhr = CLDhr,
    SOLRhr = SOLRhr,
    RAINhr = RAINhr,
    ZENhr = ZENhr,
    IRDhr = IRDhr,

    # soil solar reflectance
    REFLS = rep(soil_reflectance, n_days),

    # daily percentage of the ground surface area that is covered in water -
    # overridden by soil moisture model, so set to 0
    PCTWET = rep(0, n_days),

    # set the initial soil temperatures to the mean of the air temperature over
    # the simulation
    soilinit = soil_objects$soil_temp_c_init,

    # horizon angles, to compute solar gain over time
    hori = hori,

    # solar attentuation at a given latitude and longitude (from GADS)
    TAI = solar_attenuation_data$solar_attenuation,

    # properties and initial moisture of the soil
    soilprops = soil_objects$soil_properties,
    moists = soil_objects$soil_moisture,

    # daily rainfall vector
    RAINFALL = daily_rainfall_mm,

    # temperature of the soil at depth
    tannulrun = soil_objects$deep_soil_temp_c,

    # more soil properties:
    PE = soil_objects$PE,
    KS = soil_objects$KS,
    BB = soil_objects$BB,
    BD = soil_objects$BD,
    DD = soil_objects$DD,
    L = soil_objects$L,

    # leaf area index, used to partition traspiration/evaporation from PET
    LAI = rep(0.1, n_days)
  )

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
    # angle (degrees) up to the horizon in 15 degree rotational angles, starting
    # at 0 (north) - horixzon angle of 0 implies no buildings/hills
    hori = rep(0, 24),
    even_rain = TRUE,
    # reference height in metres at which air temperature, wind speed and
    # relative humidity input data are measured
    weather_height_m = 2) {

  # find the hemisphere (1 is north)
  hemisphere <- ifelse(latitude < 0, 2, 1)

  # absolute degrees and minutes in this hemisphere
  lat_degree <- abs(trunc(latitude))
  lat_minute <- 60 * (abs(latitude) - lat_degree)
  lon_degree <- abs(trunc(longitude))
  lon_minute <- 60 * (abs(longitude) - lon_degree)

  # convert horizon angles to radians to calc view factors
  view_factor <- 1 - mean(sin(hori * pi / 180))

  # return a horrifying vector - order matters, names do not :(
  microinput <- c(ndays = n_days,
                  RUF = microinput_defaults$RUF,
                  ERR = microinput_defaults$ERR,
                  Usrhyt = adult_height_m,
                  # reference height in metres at which air temperature, wind
                  # speed and relative humidity input data are measured
                  Refhyt = weather_height_m,
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
                  # run the soil moisture model
                  runmoist = 1,
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
  # and relative humidity input data are measured. Note all input vectors
  # assumed to be at this height
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


#' Make all the objects related to the soil structure being simulated
#' @param n_days the number of days in the simulation
#' @param mean_temp_c the average air temperature over the simulation (used for
#'   initial values and deep soil temperatures)
#' @param organic_cap is there an organic cap present on the soil surface? (cap
#'   has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)
#' @param bulk_density Soil bulk density (Mg/m3), single value or vector of 10
#'   specific to each depth
#' @param density Soil minerals density, single value or vector of 10 specific
#'   to each depth (Mg/m3)
#' @param thermal_cond Soil minerals thermal conductivity, single value or vector of
#'   10 specific to each depth (W/mK)
#' @param spec_heat Soil minerals specific heat, single value or vector of 10 specific
#'   to each depth (J/kg-K)
#' @noRd
make_soil_data <- function(n_days,
                           mean_temp_c,
                           organic_cap = TRUE,
                           bulk_density = 1.3,
                           density = 2.56,
                           thermal_cond = 2.5,
                           spec_heat = 870) {

  # depth in cm: "must be 10 values starting from 0, and more closely
  # spaced near the surface". fine.
  soil_depths_cm <- c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200)

  # density of roots in the soil at these different depths (m/m3),
  root_density_vector <- 10000 * c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0,
                            1.8, 0.9, 0.6, 0.8, 0.4, 0.4, 0, 0)

  # define the depth (integer giving the deepest node) for each substrate type.
  # We must pass a matrix of dimension 10 (number of substrate types to model)
  # into microclimate, but we are only using 2 (optional organic cap and then
  # main soil type; setting Numtyps = 2 in the microinput), so leave rest blank
  soil_nodes <- matrix(0, nrow = 10, ncol = n_days)
  soil_nodes[1, ] <- 3
  soil_nodes[2, ] <- 9

  # populate a matrix of soil properties (5 columns), again for up to 10 types
  # (rows) but only using 2 of them
  soil_properties <- matrix(0, nrow = 10, ncol = 5)
  # soil bulk density
  soil_properties[1:2, 1] <- bulk_density
  # saturated water content
  soil_properties[1:2, 2] <- min(0.26, 1 - bulk_density / density)
  # thermal conductivity
  soil_properties[1:2, 3] <- thermal_cond
  # heat capacity
  soil_properties[1:2, 4] <- spec_heat
  # mineral density
  soil_properties[1:2, 5] <- density

  # if there is an 'organic cap' on the soil, overwrite the thermal conductivity
  # and heat capacity for profile 1
  if (organic_cap) {
    soil_properties[1, 3] <- 0.2
    soil_properties[1, 4] <- 1920
  }

  # initial soil water content at each soil node, m3/m3
  soil_moisture_vector <- c(0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.3, 0.3, 0.3, 0.3)
  # make empty matrix for soil moisture values through time
  soil_moisture <- matrix(0, nrow = 10, ncol = n_days)
  # fill with the initial soil moisture values
  soil_moisture[1:10, ] <- soil_moisture_vector
  # enforce a maximum (density = bulk_density + water?)
  soil_moisture <- pmin(soil_moisture, 1 - bulk_density / density)

  # output a list of all these for setting up the microclimate model
  list(

    soil_depths_cm = soil_depths_cm,
    soil_nodes = soil_nodes,
    soil_properties = soil_properties,
    soil_moisture = soil_moisture,

    # use the mean air temperature to set the initial soil temperatures and the
    # deep soil temperature
    soil_temp_c_init = rep(mean_temp_c, 20),
    deep_soil_temp_c = rep(mean_temp_c, n_days),

    # various other soil properties (vector of 19 values descending through soil
    # for specified soil nodes in parameter DEP and points half way between)

    # Air entry potential (J/kg)
    PE = rep(1.1, 19),

    # Saturated conductivity (kg s/m3)
    KS = rep(0.0037, 19),

    # Campbell's soil 'b' parameter (unitless)
    BB = rep(4.5, 19),

    # Soil bulk density (Mg/m3)
    BD = rep(bulk_density, 19),

    # Soil density (Mg/m3)
    DD = rep(density, 19),

    # Root density (m/m3),
    L = root_density_vector

  )

  # NOTE: the following can be queried from soilgrids: PE, BB, BD, KS,
  # BulkDensity, based on the soil types at the 10 different depths. This then
  # requires Numtyps = 10 in microinput, and setting soil_nodes values to 1:10

}

