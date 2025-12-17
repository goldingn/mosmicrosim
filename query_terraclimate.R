# script to mess around with mapping microclimates using terraclimate data
# - needs to be turned into a vignette eventually

# get an old install of NicheMapR that works for this code
# remotes::install_github("mrke/NicheMapR@v3.3.1")

# load all internal package functions, for hacking purposes
pkgload::load_all()

# load major dependencies
library(NicheMapR)
library(terra)
library(ggplot2)

# demo: pointwise analysis

# # Get weather for central Perth WA for 2020-24
# longitude <- 115.86
# latitude <- -31.95
# dates <- seq(as.Date("2020-01-01"),
#              as.Date("2024-12-31"),
#              by = "1 day")
#
# # build the indexing slice at this location
# slice <- terraclimate_build_slice_point(longitude = longitude,
#                                         latitude = latitude,
#                                         dates = dates)
#
# # get all variables
# vars <- c("tmax", "tmin",  # temperature
#           "ppt",  # precipitation
#           "ws",  # wind speed
#           "vpd",  # vapor pressure deficit (for rel humidity)
#           "srad")  # solar radiation
#
# climate_monthly <- data.frame(
#   start_date = slice$dates$start,
#   end_date = slice$dates$end,
#   mid_date = slice$dates$start + (slice$dates$end - slice$dates$start) / 2)
#
# for(var in vars) {
#   climate_monthly[, var] <- terraclimate_fetch(slice, var)
# }
#
# # compute relative humidity
#
# # mean temperature and vapour pressure
# climate_monthly$tmean <- (climate_monthly$tmax + climate_monthly$tmin) / 2
# climate_monthly$vp <- vapour_pressure(climate_monthly$tmean, climate_monthly$vpd)
#
# # RH at max/min temps
# climate_monthly$rhmax <- relative_humidity(climate_monthly$tmin, climate_monthly$vp)
# climate_monthly$rhmin <- relative_humidity(climate_monthly$tmax, climate_monthly$vp)
#
# # adjust wind speed to from measured height (10m) to match the other inputs (2m)
# # and split by max/min
# climate_monthly$wsmax <- adjust_wind_speed(climate_monthly$ws,
#                                            height_required = 2,
#                                            height_measured = 10)
# climate_monthly$wsmin <- climate_monthly$wsmax * 0.1
#
# # compute cloud cover from solar radiation (terraclimate) and expected clear sky
# # radiation (microclimate code using coordinates and aerosol data) and manual
# # multipliers to enforce daily variation, as used in NicheMapR
#
# # monthly summaries of clear sky radiation (doesn't depend on weather)
# clearsky_radiation_monthly <- clear_sky_radiation_12mo(latitude = latitude,
#                                                        longitude = longitude)
# monthly_index <- lubridate::month(climate_monthly$mid_date)
#
# # convert to cloud cover percentage
# climate_monthly$ccmax <- cloud_cover(climate_monthly$srad,
#                                      clearsky_radiation_monthly[monthly_index],
#                                      multiplier = 1.2)
# climate_monthly$ccmin <- cloud_cover(climate_monthly$srad,
#                                      clearsky_radiation_monthly[monthly_index],
#                                      multiplier = 0.8)
#
# # convert precipitation to log scale to spline
# climate_monthly$rainfall_mm <- climate_monthly$ppt
# climate_monthly$log_rainfall_mm <- log1p(climate_monthly$rainfall_mm)
#
# # # compare with the Hulmes interpolated cloud cover raster at this site
# # cloud_cover_data <- get_cloud_cover_raster()
# # cloud_cover_hulmes <- terra::extract(cloud_cover_data,
# #                                      data.frame(longitude, latitude))
# # srad_sry <- tapply(climate_monthly$srad, monthly_index, FUN = mean)
# #
# # # our calculation from solar radiation is too high
# #
# # par(mfrow = c(1, 2))
# #
# # # from Hulmes:
# # plot(t(cloud_cover_hulmes[, -1]), ylim = c(0, 100))
# #
# # # from our calculation:
# # plot(cloud_cover(srad_sry, clearsky_radiation_monthly, multiplier = 1),
# #      ylim = c(0, 100))
#
# # spline these to the requested dates
# climate_daily <- data.frame(
#   date = dates
# )
# new_var <- c("tmax", "tmin",
#              "rhmax", "rhmin",
#              "wsmin", "wsmax",
#              "ccmax", "ccmin",
#              "log_rainfall_mm")
# for (var in new_var) {
#   climate_daily[, var] <- spline_seasonal(values = climate_monthly[, var],
#                                           dates = climate_monthly[, "mid_date"],
#                                           dates_predict = climate_daily[, "date"])
# }
#
# # convert log rainfall back
# climate_daily$rainfall_mm <- expm1(climate_daily$log_rainfall_mm)
#
# # it would be nice if we could use the correct likelihood for disaggregation
# # in these spline models, but this is cheap and probably not much different
#
# var_plot <- c("tmax", "tmin",
#               "rhmax", "rhmin",
#               "wsmin", "wsmax",
#               "ccmax", "ccmin",
#               "rainfall_mm")
# par(mfrow = c(1, 1))
# op <- par(mfrow = n2mfrow(length(var_plot)),
#           mar = c(3, 3, 2, 1))
# for (var in var_plot) {
#   ylims <- range(c(climate_daily[, var], climate_monthly[, var]))
#   xlims <- range(c(climate_daily$date, climate_monthly$mid_date))
#   plot(climate_daily[, var] ~ climate_daily$date,
#        type = "l",
#        col = grey(0.4),
#        lwd = 2,
#        ylab = "",
#        xlab = "",
#        main = var,
#        pch = 16,
#        ylim = ylims,
#        xlim = xlims)
#
#   points(climate_monthly[, var] ~ climate_monthly$mid_date,
#          lwd = 1.5,
#          cex = 0.7)
# }
# par(op)
#
# altitude <- altitude_m(longitude = longitude, latitude = latitude)
#
#
# # plug these into the microclimate simulation:
#
# shade_proportion <- 0.95
# adult_height <- 0.1
#
# micro <- create_micro(latitude = latitude,
#                       longitude = longitude,
#                       altitude_m = altitude,
#                       dates = dates,
#                       daily_temp_max_c = climate_daily$tmax,
#                       daily_temp_min_c = climate_daily$tmin,
#                       daily_rh_max_perc = climate_daily$rhmax,
#                       daily_rh_min_perc = climate_daily$rhmin,
#                       daily_cloud_max_perc = climate_daily$ccmax,
#                       daily_cloud_min_perc = climate_daily$ccmin,
#                       daily_wind_max_ms = climate_daily$wsmax,
#                       daily_wind_min_ms = climate_daily$wsmin,
#                       daily_rainfall_mm = climate_daily$rainfall_mm,
#                       weather_height_m = 2,
#                       adult_height_m = adult_height,
#                       even_rain = FALSE,
#                       shade_prop = shade_proportion)
#
# # profvis::profvis(
# system.time(
#   sim <- NicheMapR::microclimate(micro)
# )
# # )
#
# # # it's frustrating that runshade=1 is required, contact Mike with a reprex and
# # # to to ask why?
# # micro2 <- micro
# # micro2$microinput["runshade"] <- 0
# # system.time(
# #   sim2 <- NicheMapR::microclimate(micro2)
# # )
# # summary(sim2$shadmet[, "TALOC"])
#
#
# # plot some examples of microclimate conditions, with the external conditions
# # over the top
#
# n_days <- 52 * 7
# days <- 1:n_days
# hours <- 0 * 24 + 1:(24 * n_days)
# dates_plot <- dates[days] + lubridate::hours(1)
# hours_plot <- dates[1] + lubridate::hours(hours)
#
# par(mfrow = c(2, 2), mar = c(5, 4, 4, 2) + 0.1)
#
# vals <- sim$shadmet[, "TALOC"][hours]
# plot(vals ~ hours_plot,
#      type = "l",
#      ylab = "C",
#      las = 1,
#      lwd = 0.1,
#      ylim = c(0, 40))
# lines(climate_daily$tmax[days] ~ dates_plot,
#       lty = 2)
# lines(climate_daily$tmin[days] ~ dates_plot,
#       lty = 2)
# title(main = "Microclimate temperature")
#
# vals <- sim$shadmet[, "RHLOC"][hours]
# plot(vals ~ hours_plot,
#      type = "l",
#      ylab = "%",
#      las = 1,
#      lwd = 0.1,
#      ylim = c(0, 100))
# lines(climate_daily$rhmin[days] ~ dates_plot,
#       lty = 2)
# lines(climate_daily$rhmax[days] ~ dates_plot,
#       lty = 2)
# title(main = "Relative humidity")
#
# vals <- sim$shadmet[, "VLOC"][hours]
# plot(vals ~ hours_plot,
#      type = "l",
#      ylab = "m/s",
#      las = 1,
#      lwd = 0.1,
#      ylim = c(0, 5))
# lines(climate_daily$wsmax[days] ~ dates_plot,
#       lty = 2)
# lines(climate_daily$wsmin[days] ~ dates_plot,
#       lty = 2)
# title(main = "Wind speed")
#


# demo: batch process analysis in tiles

# load a template raster for Africa
template <- rast("~/Dropbox/github/ir_cube/data/clean/raster_mask.tif")

# create a corresponding terraclimate raster (aligned with terraclimate grid
# location, and excluding missing cells in terraclimate or more than 1 cell away
# from template cells with data)
tc_template <- make_terraclimate_template(template)

# make an altitude raster on the right pixels by disaggregating and resampling
# the lower-res NicheMapR one, and interpolating coastal regions
altitude_nichemapr <- get_altitude_raster()
agg <- round(mean(res(altitude_nichemapr) / res(tc_template)))
altitude <- altitude_nichemapr |>
  terra::disagg(agg) |>
  terra::crop(tc_template) |>
  terra::resample(tc_template,
                  method = "bilinear",
                  na.rm = TRUE) |>
  terra::mask(tc_template)

# extend out this layer to match
altitude_blur <- terra::focal(altitude,
                              w = 5,
                              fun = "mean",
                              na.policy = "only")
altitude_blur[!is.na(altitude)] <- altitude
altitude <- terra::mask(altitude_blur, tc_template)

# why does the point fall outside the terraclimate template?


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

# define all dates to extract per tile (earlier ones will be extracted to
# enable burnin)
dates <- seq(as.Date("2000-01-01"),
             as.Date("2024-12-31"),
             by = "1 day")

# set up microclimate parameters
microclimate_params <- list(
  shade_proportion = 1,
  adult_height = 1
)

# variables required for running NicheMapR sims

# extract all terraclimate data for this tile
tile_data <- terraclimate_extract_tile(tile_number = 1,
                                       tiles = tiles,
                                       dates = dates) |>
  # subset to only those cells inside the template raster
  clean_tile_data(tc_template)


# subset this for now for testing
tile_data_sub <- tile_data |>
  dplyr::filter(
    longitude %in% longitude[1:2],
    latitude %in% latitude[1:2]
  )


system.time(
# process these variables for input to NicheMapR
sims <- tile_data_sub |>

  # convert all terraclimate variables to the required NicheMapR input types,
  # stored in a per-location tibble in the 'monthly_climate' column
  process_terraclimate_tile_vars() |>

  # append altitude of each site based on the coordinates
  dplyr::mutate(
    altitude = altitude_m(
      longitude = longitude,
      latitude = latitude,
      altitude_raster = altitude
    ),
    .before = monthly_climate
  ) |>

  # this sub-monthly modelling section to be replaced with a single wrapper
  # function taking monthly climate data and converting it into monthly vector
  # abundance data

  # interpolate daily climate data from monthly climate data
  dplyr::group_by(
    longitude,
    latitude,
    altitude
  ) |>
  dplyr::summarise(
    daily_climate = list(
      daily_from_monthly_climate(
        monthly_climate = monthly_climate[[1]]
      )
    ),
    .groups = "drop"
  ) |>

  # simulate hourly (micro-)climate data from daily climate data
  dplyr::group_by(
    longitude,
    latitude,
    altitude,
  ) |>
  dplyr::summarise(
    hourly_climate = list(
      hourly_from_daily_climate(
        latitude = latitude,
        longitude = longitude,
        altitude = altitude,
        daily_climate = daily_climate[[1]],
        microclimate = microclimate_params
      )
    ),
    .groups = "drop"
  ) |>
  dplyr::group_by(
    longitude,
    latitude,
    altitude,
    hourly_climate
  ) |>
  dplyr::summarise(
    water_surface_area = list(
      simulate_ephemeral_habitat(
        hourly_climate = hourly_climate[[1]],
        altitude = altitude,
        initial_volume = 0,
        burnin_years = 0,
        max_cone_depth = 1,
        inflow_multiplier = 1
      )
    ),
    .groups = "drop"
  )
)

# ~5 mins per pixel runtime

# ncell(tc_template) / 64
# # 45000 cells per CPU
# 5 * 45000 / 60
# 3750 hours
# 3750/ 24
#
# 156 days

# we need a faster alternative to nichemapr:

# hourly interpolation of air temperature, humidity, and windspeed.

# use this to solve for water volume

# use air temperature and humidity to drive adult vector model (shaded ambient
# conditions)

# work out how to approximate water temperature:

# lose heat from evaporation

# at each timestep, given the mass of water evaporated E (in kg) and the Latent
# Heat of Vaporization λ (in kJ/kg, for water approximately 2260 kJ/kg at
# standard pressure) representing the energy to change liquid to gas. Then can
# calculate Heat Loss Q_{evap} as Q_{evap} = E λ (in kJ).

# gain radiant heat from air (given water and air temperatures): Q_{rad}
# following Stefan-Boltzmann Law and current temperatures of water and air

# Q_{rad}= sigma epsilon A (T_{water}^{4} - T_{air}^{4})
# sigma: Stefan-Boltzmann constant (5.67x10^-8 W / m^2 K ^4)
# epsilon: Emissivity (water is ~0.97-0.99, air gases are complex)

# For a fully shaded pool (no direct solar radiation)
# get Heat exchange: Q = Q_{rad} - Q_{evap}

# For an unshaded pool (solar radiation, accounting for cloud cover - c)
# get Heat exchange: Q = Q_{rad} - Q_{solar} - Q_{evap}

# we can model either a shaded (radiation from air) or an unshaded puddle
# (radiation plus solar gain, from clearsky and terraclimate-derived
# ccmin/ccmax) https://codes.ecmwf.int/grib/param-db/169, assuming low or no
# albedo)

# Calculate heat gain from solar radiation Q_{solar}:
# Incoming solar radiation ISR = clearsky_SRAD (in m^2) * cloud_cover_fraction
# Absorbed solar radiation ASR = ISR * (1 - albedo) (assume very low albedo for water)
# Surface area (A)
# Q_{solar} = ASR × A (in Watts or Joules/second)

# then compute temperature change as:
#   delta_{T} = Q / (m c)
# where m (in kg) is current mass of water and c is specific heat capacity of
# water https://en.wikipedia.org/wiki/Specific_heat_capacity c = 4184 J/kg/K, so
# convert into kJ and solve for change in temperature

# the just iterate!



df <- sims$hourly_climate[[1]]

df |>
  dplyr::filter(
    date >= as.Date("2024-06-14"),
    date < as.Date("2024-06-16"),
  ) |>
  dplyr::mutate(
    day = as.numeric(date - min(date)),
    hour = hour + day * 24
  ) |>
  dplyr::select(
    -date,
    -day
  ) |>
  tidyr::pivot_longer(
    cols = !any_of("hour"),
    names_to = "variable",
    values_to = "value"
  ) |>
  ggplot(
    aes(x = hour,
        y = value)
  ) +
  geom_line() +
  facet_wrap(
    ~variable,
    scales = "free_y"
  ) +
  theme_minimal()


tile_data_sub |>
  process_terraclimate_tile_vars() |>
  dplyr::pull(monthly_climate) |>
  unlist(recursive = FALSE) |>
  dplyr::as_tibble() |>
  dplyr::filter(
    start == as.Date("2024-06-01")
  )

# microclimate air temperature is lower, humidity is higher

# what if we skipped the nichemapr stuff?

# spline interpolation of day/night temps, cloud cover, etc.

# see how nichemapr puts the peaks and troughs for hourly data and align






par(mfrow = c(2, 2))
tile_data_sub |>
  dplyr::filter(
    start > min(dates),
    variable == "ppt"
  ) |>
  dplyr::pull(value) |>
  plot(type = "l",
       main = "monthly rainfall")

sims$hourly_climate[[1]] |>
  dplyr::filter(date > min(dates)) |>
  dplyr::pull(rainfall) |>
  plot(type = "l",
       main = "hourly rainfall")
sims$water_surface_area[[1]][sims$hourly_climate[[1]]$date > min(dates)] |>
  plot(type = "l",
       main = "water surface area")
sims$hourly_climate[[1]] |>
  dplyr::filter(date > min(dates)) |>
  dplyr::pull(air_temperature) |>
  plot(type = "l",
       main = "hourly temperature")

# need to copy over and apply the an gambiae population dynamics

# wrapper function:

# # need to add microclimate information (like in An stephensi version)


# input: tile_data_for_nichemapr

# group by locations

# then in the major parallel wrapper, for each row:
# - take the monthly data tibble and compute a tibble of the daily data (using
#    spline seasonal)
# - feed this and the other row information into create_micro()
# - feed the output of create_micro() into NicheMapR::microclimate()
# - feed the output of NicheMapR::microclimate() into the waterbody simulation
#    (to bring over from Anopheles stephensi work)
# - feed the output of NicheMapR::microclimate() and the waterbody simulation
#    into the population dynamics simulation (to bring over from Anopheles
#    stephensi work)

# need wrapper function to do the following for monthly data at each location:
# 1. spline interpolation of all variables (to daily)
# 2. microclimate simulation (to hourly)
# 3. water body simulation (hourly)
# 4. population dynamics (hourly)
# 5. summarise population dynamics back to monthly

# then execute this in parallel across all cells in a tile

# To do:

# Build wrapper functions to:

# Create processing tiles - DONE

# create clearsky_rad raster for terraclimate - DONE

# Within each tile:

# download all terraclimate data for the tile 2000-2025 - DONE

  # For each pixel in the tile:

    # format terraclimate to get NicheMapR inputs - DONE

    # run spline interpolation to get daily outdoor data 2000-2025 - DONE

    # run NicheMapR to get hourly microclimate data 2000-2025 - DONE

    # run cone model to get hourly water surface area 2000-2025 - DONE

    # run the population dynamic models to get hourly population data 2000-2025

    # summarise the population dynamic and water surface area outputs to monthly
    # data 2000-2025

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
