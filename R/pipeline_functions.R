# Pipeline functions for modelling vector metrics (population sizes and
# lifehistory parameters) from terraclimate data and microclimate models


# given a template raster aligned with the terraclimate grid, make
# (approximately target_n_tiles) processing tiles, and return a tibble with a
# column of tile numbers and a list-column of extents
#   tiles <- make_tiles(tc_template, target_n_tiles = 100)
make_tiles <- function(tc_template, target_n_tiles = 100) {

  tc_template |>
    tile_template_raster(
      target_n_tiles = 100
    ) |>
    dplyr::as_tibble() |>
    dplyr::rowwise() |>
    dplyr::summarise(
      extent = list(
        terra::ext(c(xmin, xmax, ymin, ymax))
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      tile = dplyr::row_number(),
      .before = everything()
    )
}

# then loop across the tiles in this, probably not in a pipe, but a
# for-loop saving results, or split across machines

# for each extent in the tibble of tiles, and specified vector of dates run:
#   pixel_terraclimate_data <- extract_terraclimate_tile(extent, dates)
# to return a tibble of pixels, with latitude, longitude, and a
# list-column of the monthly terraclimate variables extracted
extract_terraclimate_tile <- function(extent, dates, tc_template) {

  # build the slice for all pixels (including NAs) in this tile
  tile_slice <- terraclimate_build_slice_extent(
    extent = extent,
    dates = dates
  )

  variables <- c("tmax", "tmin", "ppt", "ws", "vpd", "srad")

  # download each of the variables for this slice
  stack_list <- list()
  for(var in variables) {
    stack_list[[var]] <- terraclimate_fetch(tile_slice, var)
  }

  # combine into tidy long-format tibble, by making a 4D array with appropriate
  # dimension names, coercing to long form via a table, and then tidying up the
  # tibble
  stack <- do.call(abind::abind,
                   c(stack_list, list(along = 4)))

  dimnames(stack) <- list(
    longitude = tile_slice$longitudes,
    latitude = tile_slice$latitudes,
    start = as.character(tile_slice$dates$start),
    variable = variables
  )

  stack |>
    as.data.frame.table(
      stringsAsFactors = FALSE
    ) |>
    dplyr::as_tibble() |>
    dplyr::rename(
      value = Freq
    ) |>
    dplyr::filter(
      !is.na(value)
    ) |>
    dplyr::mutate(
      longitude = as.numeric(longitude),
      latitude = as.numeric(latitude),
      start = as.Date(start)
    ) |>
    # add on the end dates
    dplyr::left_join(
      dplyr::as_tibble(tile_slice$dates),
      by = "start"
    ) |>
    dplyr::relocate(
      end,
      .after = start
    ) |>
    # remove pixels where the template has no data
    clean_tile_data(tc_template) |>
    # convert to a tibble for each pixel
    tidyr::pivot_wider(
      names_from = variable,
      values_from = value
    ) |>
    dplyr::group_by(
      latitude,
      longitude
    ) |>
    dplyr::summarise(
      monthly_terraclimate = list(
        dplyr::tibble(
          start, end, tmax, tmin, ppt, ws, vpd, srad
        )
      ),
      .groups = "drop"
    )

}

# then, run the following on the whole tibble returned by
# extract_terraclimate_tile:
#   pixel_monthly_climate <- terraclimate_to_monthly_climate(
#     pixel_terraclimate_data
#   )
# to return a per-pixel tibble with: latitude, longitude, *altitude*, and
# list-column of tibbles with the monthly climate data
terraclimate_to_monthly_climate <- function(pixel_terraclimate_data) {

  # get the rasters of terraclimate-aligned elevation, and terraclimate-GADS
  # hybrid cell numbers (for clear sky radiation)
  tc_elev <- terra::unwrap(tc_elev_wrapped)
  tc_gads <- terra::unwrap(tc_gads_wrapped)

  # find unique coordinates in this tile
  unique_coords <- pixel_terraclimate_data |>
    dplyr::distinct(
      longitude,
      latitude
    )

  # create a lookup from the coordinates to their clearsky values
  clearsky_lookup <- unique_coords |>
    dplyr::mutate(
      cell_id = terra::extract(tc_gads,
                               unique_coords,
                               ID = FALSE)[, 1]
    ) |>
    # and the clearsky monthly information for these cells
    dplyr::left_join(
      clear_sky_lookup,
      by = "cell_id"
    ) |>
    dplyr::select(
      -cell_id
    ) |>
    # add on a month ID
    dplyr::mutate(
      month_id = list(1:12),
      .before = clear_sky_SOLR
    ) |>
    # unnest lists
    tidyr::unnest(
      c(month_id, clear_sky_SOLR)
    ) |>
    # get rid of column-ness
    dplyr::mutate(
      clear_sky_SOLR = as.numeric(clear_sky_SOLR)
    )

  # process the data
  pixel_terraclimate_data |>
    # look up the altitude for all pixels
    dplyr::mutate(
      altitude = altitude_m(
        longitude = longitude,
        latitude = latitude,
        altitude_raster = tc_elev
      ),
      .before = monthly_terraclimate
    ) |>
    # for each pixel, convert monthly terraclimate to monthly climate - the
    # variables needed for modelling microclimates
    dplyr::group_by(
      latitude,
      longitude,
      altitude
    ) |>
    dplyr::summarise(
      monthly_climate = list(
        process_terraclimate_pixel_vars(
          monthly_terraclimate[[1]],
          latitude,
          longitude,
          clearsky_lookup
        )
      ),
      .groups = "drop"
    )

}

# then run the following on the whole tibble:
#   pixel_daily_climate <- interpolate_daily_climate(pixel_monthly_climate)
# to return a per-pixel tibble with: latitude, longitude, altitude, and
# list-column of tibbles with the daily climate data
interpolate_daily_climate <- function(pixel_monthly_climate) {

  pixel_monthly_climate |>
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
    )

}

# then run the following on the whole tibble:
#   pixel_hourly_microclimate <- simulate_hourly_microclimate(
#     pixel_daily_climate
#   )
# to return a per-pixel tibble with: latitude, longitude, altitude, and
# list-column of tibbles with the hourly microclimate data (not including water
# temperature)
simulate_hourly_microclimate <- function(pixel_daily_climate) {

  pixel_daily_climate |>
    dplyr::group_by(
      longitude,
      latitude,
      altitude,
    ) |>
    dplyr::summarise(
      hourly_microclimate = list(
        hourly_from_daily_climate_ambient(
          daily_climate = daily_climate[[1]]
        )
      ),
      .groups = "drop"
    )

}


# then run the following on the whole tibble:
#   pixel_hourly_conditions <- simulate_hourly_conditions(
#     pixel_hourly_microclimate,
#     model_water_temperature = FALSE,
#     water_shade_proportion = 1
#   )
# to return a per-pixel tibble with: latitude, longitude, (*not altitude*), and
# list-column of tibbles with the hourly condition data, including microclimate
# (air), water surface area, and water temperature (altitude is needed to model
# evaporation). At first, fix model_water_temperature = FALSE and just copy the
# air temperature, then when we implement a basic model of water temperature, at
# first fix water_shade_proportion = 1 and ignore solar gain in water bodies,
# then later implement the water temperature model with solar gain, and with
# clearsky radiation multiplied by (1 - water_shade_proportion) * (1 -
# cloud_cover)
simulate_hourly_conditions <- function(
    pixel_hourly_microclimate,
    model_water_temperature = FALSE,
    water_shade_proportion = 1
) {

  # warn the user if they tried to set a model parameter for water shade
  # proportion, but forgot to run the model
  unused_setting <- water_shade_proportion != 1 & !model_water_temperature
  if (unused_setting) {
    warning("water temperature is not being modelled, ",
            "so water_shade_proportion argument is ignored",
            call. = FALSE)
  }

  if (model_water_temperature) {
    stop("model of water temperature not yet implemented, ",
         "model_water_temperature must be set to FALSE",
         call. = FALSE)
  }

  pixel_hourly_microclimate |>
    dplyr::group_by(
      longitude,
      latitude,
    ) |>
    dplyr::summarise(
      hourly_conditions = list(
        dplyr::tibble(
          # drop some conditions that don't directly affect vectors, but keep
          # the rest
          dplyr::select(
            hourly_microclimate[[1]],
            -cloudcover,
            -windspeed,
            -rainfall
          ),
          # model and add the water surface area
          water_surface_area = simulate_ephemeral_habitat(
            hourly_climate = hourly_microclimate[[1]],
            altitude = altitude,
            initial_volume = 0,
            burnin_years = 0,
            max_cone_depth = 1,
            inflow_multiplier = 1
          )
        )
      ),
      .groups = "drop"
    )
}

# then run the following on the whole tibble:
  # pixel_hourly_lifehistory <- simulate_hourly_lifehistory(
  #   pixel_hourly_conditions,
  #   species = "An. gambiae"
  # )
# to return a per-pixel tibble with: latitude, longitude, and list-column of
# tibbles with the hourly lifehistory and water surface area simulations, and
# the species name
simulate_hourly_lifehistory <- function(
  pixel_hourly_conditions,
  species = c("An. gambiae", "An. stephensi")
) {
  species <- match.arg(species)

  lifehistory <- switch(
    species,
    "An. gambiae" = lifehistory_functions$An_gambiae,
    "An. stephensi" = lifehistory_functions$An_stephensi
  )

  pixel_hourly_conditions |>
    dplyr::mutate(
      species = species
    ) |>
    dplyr::group_by(
      longitude,
      latitude,
      species
    ) |>
    dplyr::summarise(
      # compute the pre-computed (environmentally-determined) lifehistory
      # variables
      hourly_lifehistory = list(
        simulate_vector_lifehistory(
          hourly_conditions[[1]],
          lifehistory = lifehistory
        )
      ),
      # also attach the daily aquatic survival density modification functions,
      # for use in dynamically solving density dependent aquatic-stage survival
      # effects
      das_densmod_function = list(
        lifehistory$das_densmod
      ),
      .groups = "drop"
    )

}

# then run the following on the whole tibble:
#   pixel_hourly_vectors <- simulate_hourly_vectors(pixel_hourly_lifehistory)
# to return a per-pixel tibble with: latitude, longitude, and list-column of
# tibbles with the hourly vector information (number of adults, number of
# aquatics, and the lifehistory parameters, water surface area simulations and
# the species name)
simulate_hourly_vectors <- function(pixel_hourly_lifehistory) {

  # run the population simulation code, for each pixel separately
  pixel_hourly_lifehistory |>
    dplyr::group_by(
      longitude,
      latitude,
      species
    ) |>
    dplyr::summarise(
      # return a combined tibble of vector populations and environment-dependent
      # lifehistory parameters
      hourly_vector = list(
        dplyr::bind_cols(
          # put date and hour first
          dplyr::select(
            hourly_lifehistory[[1]],
            date,
            hour
          ),
          # simulated population sizes
          simulate_population(
            # precomputed lifehistory parameters
            hourly_lifehistory = hourly_lifehistory[[1]],
            # and the corresponding density dependence modifier on aqautic-stage
            # survival
            das_densmod_function = das_densmod_function[[1]]
          ),
          # and environment-dependent lifehistory parameters
          dplyr::select(
            hourly_lifehistory[[1]],
            -date,
            -hour
          )
        )
      ),
      .groups = "drop"
    )

}

# then run the following on the whole tibble:
#   pixel_vectors <- summarise_vectors(pixel_hourly_vectors)
# to return a per-pixel tibble with: latitude, longitude, and list-column of
# tibbles with the monthly (or daily, if aggregate_by = "day") vector
# information (average number of adults, average number of aquatics, average
# lifehistory parameters and water surface areas and the species name)
summarise_vectors <- function(pixel_hourly_vectors,
                              aggregate_by = c("month", "day")) {

  aggregate_by <- match.arg(aggregate_by)

  # loop across the pixels, for each tibble of hourly vector data, group by
  # either the date or the month and summarise by the mean for the useful
  # quantities

  pixel_hourly_vectors |>
    dplyr::group_by(
      longitude,
      latitude,
      species
    ) |>
    dplyr::summarise(
      pixel_vectors = list(
        aggregate_vectors(
          hourly_vector[[1]],
          aggregate_by = aggregate_by
        )
      ),
      .groups = "drop"
    )
}

# aggregate the vector data to a given date or a given month
aggregate_vectors <- function(vector_tibble,
                              aggregate_by = c("day", "month")) {

  aggregate_by <- match.arg(aggregate_by)

  vector_variables <- c(
    "adult",
    "aquatic",
    "mdr",
    "efd",
    "das_zerodensity",
    "ds"
  )

  if (aggregate_by == "day") {
    aggregated_tibble <- vector_tibble |>
      dplyr::group_by(
        date
      ) |>
      dplyr::summarise(
        across(
          all_of(vector_variables),
          mean
        ),
        .groups = "drop"
      )
  } else if (aggregate_by == "month") {
    aggregated_tibble <- vector_tibble |>
      dplyr::mutate(
        start = lubridate::floor_date(date,
                                      unit = "month"),
        end = lubridate::ceiling_date(date,
                                      unit = "month") - 1,
      ) |>
      dplyr::group_by(
        start,
        end
      ) |>
      dplyr::summarise(
        across(
          all_of(vector_variables),
          mean
        ),
        .groups = "drop"
      )

  }
  # grouped_tibble <- switch(
  #   aggregate_by
  # )

  aggregated_tibble
}

