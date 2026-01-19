# Pipeline functions for modelling vector metrics (population sizes and
# lifehistory parameters) from terraclimate data and microclimate models

# given a named list of tidy format tibbles of variables, each with the columns:
# 'time_index', 'pixel_index', 'value', and with the list element names giving
# the name of the variable, efficiently combine them into a single tibble with
# each variable's value under its name and side by side. This is around 1e4 times
# faster than using dplyr::pivot_wider (or the dtplyr equivalent), like this:
# variable_list |>
#   dplyr::bind_rows(
#     .id = "variable"
#   ) |>
# tidyr::pivot_wider(
#   names_from = variable,
#   values_from = value
# )
recombine_variables <- function(variable_list) {
  # pull out a copy of the indices
  indices <- variable_list[[1]] |>
    dplyr::select(
      time_index,
      pixel_index
    )
  # get a list of the values for each variables, as a vector
  values_list <- lapply(variable_list,
                        function(x) {x$value})
  # name them
  names(values_list) <- names(variable_list)
  # combine them
  dplyr::bind_cols(
    indices, values_list
  )
}

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

  result <- stack |>
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
    clean_tile_data(
      tc_template
    ) |>
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

  # if we have extracted the final available month in the terraclimate data,
  # then pad the end with an additional synthetic month of data (to enable
  # spline interpolation to the back half of the last month) computed as the
  # average of the last three of those months. Ie. if we end on a December,
  # impute an January from the next year based on up to the last three January
  # datapoints in the extracted data
  last_month_start <- max(result$monthly_terraclimate[[1]]$start)
  end_of_the_line <- last_month_start == last_available_month_start()

  if (end_of_the_line) {
    result <- result |>
      dplyr::group_by(
        latitude,
        longitude
      ) |>
      dplyr::mutate(
        monthly_terraclimate = list(
          pad_last_month(
            monthly_terraclimate[[1]],
            n_previous = 3
          )
        )
      ) |>
      dplyr::ungroup()
  }

  result

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

  # the climate variables we care about for modelling water, or for modelling
  # mosquito lifehistory parameters later
  climate_vars <- c(
    "rainfall",
    "air_temperature",
    "water_temperature",
    "humidity",
    "windspeed"
  )

  # dimensions of things
  n_pixels <- nrow(pixel_hourly_microclimate)
  n_times <- nrow(pixel_hourly_microclimate$hourly_microclimate[[1]])

  # add a pixel index to the hourly microclimate data
  data_with_pixel <- pixel_hourly_microclimate |>
    dplyr::mutate(
      pixel_index = dplyr::row_number(),
      .before = everything()
    )

  # get a lookup from this to the pixel info, so we can add the info back later
  pixel_info <- data_with_pixel |>
    dplyr::select(
      -hourly_microclimate
    )

  # create a lookup between a time index and the dates and htimes
  times_info <- data_with_pixel$hourly_microclimate[[1]] |>
    dplyr::select(
      date,
      hour
    ) |>
    dplyr::mutate(
      time_index = dplyr::row_number(),
      .before = everything()
    )

  # subset to wide format with only the pixel/time indices
  variables_sub <- data_with_pixel |>
    # drop unneeded pixel info (to prevent duplication and memory clog)
    dplyr::select(
      pixel_index,
      hourly_microclimate
    ) |>
    # unnest the microclimate data and add on the time index
    tidyr::unnest(
      hourly_microclimate
    ) |>
    dplyr::left_join(
      times_info,
      by = c("date", "hour")
    ) |>
    # drop unneeded time and climate info (to prevent duplication and memory clog)
    dplyr::select(
      all_of(
        c("pixel_index",
          "time_index",
          climate_vars)
      )
    )

  # create a named list of time-by-pixel matrices of the values of each variable
  variable_list <- lapply(climate_vars,
                          function(var_name) {
                            variables_sub |>
                              dplyr::select(
                                pixel_index,
                                time_index,
                                !!var_name
                              ) |>
                              # force arrangement in the right order
                              dplyr::arrange(
                                pixel_index,
                                time_index
                              ) |>
                              dplyr::pull(
                                !!var_name
                              ) |>
                              matrix(
                                nrow = n_times,
                                ncol = n_pixels
                              )
                          })
  names(variable_list) <- climate_vars

  # NOTE: this is inelegant. The code previously did a pivot longer, grouped,
  # did a group_split, and then pivoted wider, but that was much slower

  # now run vectorised ephemeral habitat simulation to return a named list with
  # time-by-pixel matrices of water_surface_area and water_temperature
  water_list <- simulate_ephemeral_habitat_vectorised(
    rainfall_matrix = variable_list$rainfall,
    air_temperature_matrix = variable_list$air_temperature,
    humidity_matrix = variable_list$humidity,
    windspeed_matrix = variable_list$windspeed,
    altitude_vector = pixel_info$altitude,
    initial_volume = 0,
    burnin_years = 0,
    max_cone_depth = 1,
    inflow_multiplier = 1
  )

  # now convert this list of pixel-by-time condition information back into a
  # tibble with a list column of water conditions per pixel, with just the pixel
  # and time indices
  water_sub <- water_list |>
    lapply(
      function(matrix) {
        # convert matrix to tidy format tibble with indices
        n_pixels <- ncol(matrix)
        matrix |>
          `colnames<-`(
            seq_len(n_pixels)
          ) |>
          dplyr::as_tibble() |>
          dplyr::mutate(
            time_index = dplyr::row_number(),
            .before = everything()
          ) |>
          tidyr::pivot_longer(
            cols = !any_of("time_index"),
            names_to = "pixel_index",
            values_to = "value"
          ) |>
          dplyr::mutate(
            pixel_index = as.numeric(pixel_index)
          )
      }
    ) |>
    # stack these side by side (with lapply - much faster than bind_rows and
    # pivot_wider)
    recombine_variables()

  # drop index information to recombine this with the microclimate information
  # (faster than left_join)
  water_sub_noindex <- water_sub |>
    dplyr::arrange(
      pixel_index,
      time_index
    ) |>
    dplyr::select(
      -time_index,
      -pixel_index
    )

  # recombine with the microclimate variables, add time information, nest by
  # pixel, and add back on the pixel information
  variables_sub |>
    dplyr::arrange(
      pixel_index,
      time_index
    ) |>
    dplyr::bind_cols(
      water_sub_noindex
    ) |>
    # # faster than this
    # dplyr::left_join(
    #   water_sub,
    #   by = c("time_index", "pixel_index")
    # ) |>
    # add on the time information and drop the index
    dplyr::left_join(
      times_info,
      by = "time_index"
    ) |>
    dplyr::select(
      -time_index
    ) |>
    # turn the useful conditions information (not rainfall or windspeed) into a
    # list column
    dplyr::group_by(
      pixel_index
    ) |>
    dplyr::summarise(
      hourly_conditions = list(
        dplyr::tibble(
          date,
          hour,
          water_surface_area,
          air_temperature,
          humidity,
          water_temperature
        )
      ),
      .groups = "drop"
    ) |>
    # now add back on the pixel information
    dplyr::left_join(
      pixel_info,
      by = "pixel_index"
    ) |>
    dplyr::relocate(
      hourly_conditions,
      .after = everything()
    ) |>
    dplyr::select(
      -pixel_index
    )

#   # combine this list with the microclimate information that is needed for
#   # modelling mosquito lifehistory parameters
#   final_list <- c(variable_list,
#                   water_list)
#
#   # now convert this list of pixel-by-time condition information back into a
#   # tibble with a list column of conditions per pixel (with time info added back
#   # on), and add back on the pixel info
#   final_list |>
#     lapply(
#       function(matrix) {
#         # convert matrix to tidy format tibble with indices
#         n_pixels <- ncol(matrix)
#         matrix |>
#           `colnames<-`(
#             seq_len(n_pixels)
#           ) |>
#           dplyr::as_tibble() |>
#           dplyr::mutate(
#             time_index = dplyr::row_number(),
#             .before = everything()
#           ) |>
#           tidyr::pivot_longer(
#             cols = !any_of("time_index"),
#             names_to = "pixel_index",
#             values_to = "value"
#           ) |>
#           dplyr::mutate(
#             pixel_index = as.numeric(pixel_index)
#           )
#       }
#     ) |>
#     # stack these side by side (with lapply - much faster than bind_rows and
#     # pivot_wider)
#     recombine_variables() |>
# #     dplyr::bind_rows(
# #       .id = "variable"
# #     ) |>
# #     tidyr::pivot_wider(
# #       names_from = variable,
# #       values_from = value
# #     ) |>
#     # add on the time information and drop the index
#     dplyr::left_join(
#       times_info,
#       by = "time_index"
#     ) |>
#     dplyr::select(
#       -time_index
#     ) |>
#     # turn the useful conditions information (not rainfall or windspeed) into a
#     # list column
#     dplyr::group_by(
#       pixel_index
#     ) |>
#     dplyr::summarise(
#       hourly_conditions = list(
#         dplyr::tibble(
#           date,
#           hour,
#           water_surface_area,
#           air_temperature,
#           humidity,
#           water_temperature
#         )
#       ),
#       .groups = "drop"
#     ) |>
#     # now add back on the pixel information
#     dplyr::left_join(
#       pixel_info,
#       by = "pixel_index"
#     ) |>
#     dplyr::relocate(
#       hourly_conditions,
#       .after = everything()
#     ) |>
#     dplyr::select(
#       -pixel_index
#     )

}

# then run the following on the whole tibble:
  # pixel_hourly_lifehistory <- simulate_hourly_lifehistory(
  #   pixel_hourly_conditions,
  #   species = "An. gambiae"
  # )
# to return a per-pixel tibble with: latitude, longitude, and list-column of
# tibbles with the hourly lifehistory and water surface area simulations, and
# the species name. If fast_ds_temp_humid = TRUE (the default) then the
# ds_temp_humid() function for that species (daily adult survival as a fucntion
# of temperature and humidity) is replaced with a much faster but still very
# accurate bilinear interpolation function
simulate_hourly_lifehistory <- function(
  pixel_hourly_conditions,
  species = c("An. gambiae", "An. stephensi"),
  fast_ds_temp_humid = TRUE
) {
  species <- match.arg(species)

  lifehistory <- switch(
    species,
    "An. gambiae" = lifehistory_functions$An_gambiae,
    "An. stephensi" = lifehistory_functions$An_stephensi
  )

  # optionally replace the ds_temp_humid function with a much faster emulator
  # (bilinear interpolation)
  if (fast_ds_temp_humid) {
    lifehistory$ds_temp_humid <- make_temp_humid_interpolator(
      lifehistory$ds_temp_humid
    )
  }


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

  # add a pixel index to the hourly microclimate data
  data_with_pixel <- pixel_hourly_lifehistory |>
    dplyr::mutate(
      pixel_index = dplyr::row_number(),
      .before = everything()
    )

  # get a lookup from this to the pixel info, so we can add the info back later
  pixel_info <- data_with_pixel |>
    dplyr::select(
      -hourly_lifehistory
    )

  # check they are all the same species, so we can use the same density
  # dependence modification function across all pixels
  all_same_species <- all(pixel_info$species == pixel_info$species[1])
  if (!all_same_species) {
    stop(
      "All pixel calculations must be for the same species when running
      simulate_hourly_vectors(). Consider doing `group_by(species)`.",
    call. = FALSE)
  }

  das_densmod_function <- pixel_info$das_densmod_function[[1]]

  # create a lookup between a time index and the dates and htimes
  times_info <- data_with_pixel$hourly_lifehistory[[1]] |>
    dplyr::select(
      date,
      hour
    ) |>
    dplyr::mutate(
      time_index = dplyr::row_number(),
      .before = everything()
    )

  parameter_names <- c(
    "mdr",
    "efd",
    "das_zerodensity",
    "ds",
    "larval_habitat_area"
  )

  # dimensions of things
  n_pixels <- nrow(pixel_hourly_lifehistory)
  n_times <- nrow(pixel_hourly_lifehistory$hourly_lifehistory[[1]])

  # subset to wide format with only the pixel/time indices
  parameters_sub <- data_with_pixel |>
    # drop unneeded pixel info (to prevent duplication and memory clog)
    dplyr::select(
      pixel_index,
      hourly_lifehistory
    ) |>
    # unnest the lifehistory data and add on the time index
    tidyr::unnest(
      hourly_lifehistory
    ) |>
    dplyr::left_join(
      times_info,
      by = c("date", "hour")
    ) |>
    # drop unneeded time and climate info (to prevent duplication and memory clog)
    dplyr::select(
      all_of(
        c("pixel_index",
          "time_index",
          parameter_names)
      )
    )

  # create a named list of time-by-pixel matrices of the values of each
  # parameter
  parameter_list <- lapply(parameter_names,
                          function(param_name) {
                            parameters_sub |>
                              dplyr::select(
                                pixel_index,
                                time_index,
                                !!param_name
                              ) |>
                              # force arrangement in the right order
                              dplyr::arrange(
                                pixel_index,
                                time_index
                              ) |>
                              dplyr::pull(
                                !!param_name
                              ) |>
                              matrix(
                                nrow = n_times,
                                ncol = n_pixels
                              )
                          })
  names(parameter_list) <- parameter_names

  # NOTE: this is inelegant. The code previously did a pivot longer, grouped,
  # did a group_split, and then pivoted wider, but that was much slower

  # now run vectorised population simulation to return a named list with
  # time-by-pixel matrices of adult and aquatic population sizes
  population_list <- simulate_population_vectorised(
    # precomputed lifehistory parameters in time-by-pixel matrices
    larval_habitat_area_matrix = parameter_list$larval_habitat_area,
    mdr_matrix = parameter_list$mdr,
    efd_matrix = parameter_list$efd,
    das_zerodensity_matrix = parameter_list$das_zerodensity,
    ds_matrix = parameter_list$ds,
    # and the density dependence modifier on aquatic-stage
    # survival
    das_densmod_function = das_densmod_function
  )

  # now convert this list of pixel-by-time population information back into a
  # tibble with avalues by pixel and time index
  population_sub <- population_list |>
    lapply(
      function(matrix) {
        # convert matrix to tidy format tibble with indices
        n_pixels <- ncol(matrix)
        matrix |>
          `colnames<-`(
            seq_len(n_pixels)
          ) |>
          dplyr::as_tibble() |>
          dplyr::mutate(
            time_index = dplyr::row_number(),
            .before = everything()
          ) |>
          tidyr::pivot_longer(
            cols = !any_of("time_index"),
            names_to = "pixel_index",
            values_to = "value"
          ) |>
          dplyr::mutate(
            pixel_index = as.numeric(pixel_index)
          )
      }
    ) |>
    # stack these side by side (with lapply - much faster than bind_rows and
    # pivot_wider)
    recombine_variables()


  # # drop index information to recombine this with the parameter information
  # # (faster than left_join)
  population_sub_noindex <- population_sub |>
    dplyr::arrange(
      pixel_index,
      time_index
    ) |>
    dplyr::select(
      -time_index,
      -pixel_index
    )

  # recombine with the microclimate variables, add time information, nest by
  # pixel, and add back on the pixel information
  parameters_sub |>
    dplyr::arrange(
      pixel_index,
      time_index
    ) |>
    dplyr::bind_cols(
      population_sub_noindex
    ) |>
    # # faster than this
    # dplyr::left_join(
    #   population_sub,
    #   by = c("time_index", "pixel_index")
    # ) |>
    # add on the time information and drop the index
    dplyr::left_join(
      times_info,
      by = "time_index"
    ) |>
    dplyr::select(
      -time_index
    ) |>
    # turn the useful conditions information (not rainfall or windspeed) into a
    # list column
    dplyr::group_by(
      pixel_index
    ) |>
    dplyr::summarise(
      hourly_vector = list(
        dplyr::tibble(
          date,
          hour,
          adult,
          aquatic,
          larval_habitat_area,
          ds,
          mdr,
          efd,
          das_zerodensity
        )
      ),
      .groups = "drop"
    ) |>
    # now add back on the pixel information
    dplyr::left_join(
      pixel_info,
      by = "pixel_index"
    ) |>
    dplyr::relocate(
      hourly_vector,
      .after = everything()
    ) |>
    dplyr::select(
      -pixel_index
    )


  # # combine this list with the microclimate information that is needed for
  # # modelling mosquito lifehistory parameters
  # final_list <- c(parameter_list,
  #                 population_list)

  # # now convert this list of pixel-by-time condition information back into a
  # # tibble with a list column of conditions per pixel (with time info added back
  # # on), and add back on the pixel info
  # final_list |>
  #   lapply(
  #     function(matrix) {
  #       # convert matrix to tidy format tibble with indices
  #       n_pixels <- ncol(matrix)
  #       matrix |>
  #         `colnames<-`(
  #           seq_len(n_pixels)
  #         ) |>
  #         dplyr::as_tibble() |>
  #         dplyr::mutate(
  #           time_index = dplyr::row_number(),
  #           .before = everything()
  #         ) |>
  #         tidyr::pivot_longer(
  #           cols = !any_of("time_index"),
  #           names_to = "pixel_index",
  #           values_to = "value"
  #         ) |>
  #         dplyr::mutate(
  #           pixel_index = as.numeric(pixel_index)
  #         )
  #     }
  #   ) |>
  #   # stack these side by side (with lapply - much faster than bind_rows and
  #   # pivot_wider)
  #   recombine_variables() |>
  #   # add on the time information and drop the index
  #   dplyr::left_join(
  #     times_info,
  #     by = "time_index"
  #   ) |>
  #   dplyr::select(
  #     -time_index
  #   ) |>
  #   # turn the useful conditions information (not rainfall or windspeed) into a
  #   # list column
  #   dplyr::group_by(
  #     pixel_index
  #   ) |>
  #   dplyr::summarise(
  #     hourly_vector = list(
  #       dplyr::tibble(
  #         date,
  #         hour,
  #         adult,
  #         aquatic,
  #         larval_habitat_area,
  #         ds,
  #         mdr,
  #         efd,
  #         das_zerodensity
  #       )
  #     ),
  #     .groups = "drop"
  #   ) |>
  #   # now add back on the pixel information
  #   dplyr::left_join(
  #     pixel_info,
  #     by = "pixel_index"
  #   ) |>
  #   dplyr::relocate(
  #     hourly_vector,
  #     .after = everything()
  #   ) |>
  #   dplyr::select(
  #     -pixel_index
  #   )

}

# then run the following on the whole tibble:
#   pixel_vectors <- summarise_vectors(pixel_hourly_vectors)
# to return a per-pixel tibble with: latitude, longitude, and list-column of
# tibbles with the monthly (or daily, if aggregate_by = "day") vector
# information (average number of adults, average number of aquatics, average
# lifehistory parameters and water surface areas and the species name).
# included_dates can be specified to enable subsetting of dates before
# summarisation
summarise_vectors <- function(pixel_hourly_vectors,
                              aggregate_by = c("month", "day"),
                              included_dates = NULL) {

  aggregate_by <- match.arg(aggregate_by)

  # if the included dates are null, return all dates
  if (is.null(included_dates)) {
    included_dates <- unique(pixel_hourly_vectors$hourly_vector[[1]]$date)
  }

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
          aggregate_by = aggregate_by,
          included_dates = included_dates
        )
      ),
      .groups = "drop"
    )
}

# aggregate the vector data to a given date or a given month. included_dates can
# be specified to enable subsetting of dates before summarisation
aggregate_vectors <- function(vector_tibble,
                              aggregate_by = c("day", "month"),
                              included_dates = NULL) {

  aggregate_by <- match.arg(aggregate_by)

  # if the included dates are null, return all dates
  if (is.null(included_dates)) {
    included_dates <- unique(vector_tibble$date)
  }

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
      dplyr::filter(
        date %in% included_dates
      ) |>
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
      dplyr::filter(
        date %in% included_dates
      ) |>
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

  aggregated_tibble
}


# overall function to run the entire pipeline to monthly vector data. Input the
# pixel-level terracliamate data, as produced by extract_terraclimate_tile() and
# stored in saved tile info, the species name, and the dates to be included in
# the final monthly summaries. Return monthly-aggregated pixel-level vector
# timeseries information.
tc_to_vectors <- function(pixel_terraclimate_data,
                          species = c("An. gambiae",
                                      "An. stephensi"),
                          included_dates = NULL) {

  species <- match.arg(species)

  # get terraclimate data for this batch of pixels
  pixel_terraclimate_data |>
    # convert terraclimate data into monthly variables needed for microclimate
    # modelling
    terraclimate_to_monthly_climate() |>
    # interpolate these to daily max/min data
    interpolate_daily_climate() |>
    # interpolate microclimates on an hourly timestep
    simulate_hourly_microclimate() |>
    # model conditions experienced by vectors in the microclimate (microclimate,
    # plus water surface area and water temperature) on an hourly timestep
    simulate_hourly_conditions() |>
    # model vector lifehistory parameters
    simulate_hourly_lifehistory(
      species = species
    ) |>
    # model vector populations and transmission-relevant parameters
    simulate_hourly_vectors() |>
    # aggregate by month
    summarise_vectors(
      aggregate_by = "month",
      included_dates = included_dates
    )

}


# Create vector rasters from processed tile data. Takes a template raster to
# which to align the outputs, the species name and the vector variable name
# required, the directory 'vector_tile_dir' containing the RDS files of monthly
# vector summary information for different tiles (as produced by tc_to_vectors()
# ), the directory 'vector_raster_dir' for the output rasters to be saved into
# (this will be created if it does not exist). Separate single-band GeoTiff
# rasters will be saved into 'vector_raster_dir', with the filename formatted
# approximately as: `<species>.<year>.<month>.<variable>.tif`
create_vector_rasters <- function(
    template = template,
    species = c("An. gambiae", "An. stephensi"),
    variable = "adult",
    vector_tile_dir = "processing/vectors",
    vector_raster_dir = "processing/vector_rasters"
) {

  # make the terraclimate raster
  terraclimate_template <- make_terraclimate_template(template)

  # list the slices (years and months)
  slices <- get_slices(vector_tile_dir)

  n_slices <- nrow(slices)

  # loop through the slices (years and months)
  for (slice_idx in seq_len(n_slices)) {

    print(slice_idx)

    slice <- slices[slice_idx, ]

    # call create_vector_raster_slice(), and then wrap in furrr::future_map()
    create_vector_raster_slice(
      year = slice$year,
      month = slice$month,
      template = template,
      terraclimate_template = terraclimate_template,
      species = species,
      variable = variable,
      vector_tile_dir = vector_tile_dir,
      vector_raster_dir = vector_raster_dir
    )

  }

}

# Create vector rasters from processed tile data for single temporal slice (year
# and month). Takes a year, month, template raster to which to align the
# outputs, corresponding terraclimate-aligned template raster, the species name
# and the vector variable name required, the directory 'vector_tile_dir'
# containing the RDS files of monthly vector summary information for different
# tiles (as produced by tc_to_vectors() ), the directory 'vector_raster_dir' for
# the output rasters to be saved into (this will be created if it does not
# exist). Separate single-band GeoTiff rasters will be saved into
# 'vector_raster_dir', with the filename formatted approximately as:
# `<species>.<year>.<month>.<variable>.tif`. Note that the template and
# terraclimate template must be passed 'wrapped' (using terra::wrap), so that
# they can be passed to future for parallel processing
create_vector_raster_slice <- function(
    year,
    month,
    template_wrapped,
    terraclimate_template_wrapped,
    species = c("An. gambiae", "An. stephensi"),
    variable = "adult",
    vector_tile_dir = "processing/vectors",
    vector_raster_dir = "processing/vector_rasters"
) {

  # unwrap the rasters, so we can pass them to future
  template <- terra::unwrap(
    template_wrapped
  )

  terraclimate_template <- terra::unwrap(
    terraclimate_template_wrapped
  )

  # enforce the species name spelling
  species <- match.arg(species)

  # simplify the species name for filenames
  species_simple <- tolower(gsub("An. ", "an_", species))

  # create the output raster path, if needed
  if (!dir.exists(vector_raster_dir)) {
    dir.create(vector_raster_dir,
               recursive = TRUE)
  }

  # create a filepath for the final output raster
  slice_fname <- sprintf("%s.%i.%02i.%s.tif",
                         species_simple,
                         year,
                         month,
                         variable)

  slice_fpath <- file.path(vector_raster_dir, slice_fname)

  # list the vector tile files in the given repository
  vector_tile_fpaths <- list.files(vector_tile_dir,
                                   full.names = TRUE)

  # create a temporary terraclimate-aligned raster for this slice. This will be
  # repeatedly called and updated, but discarded after each slice is processed,
  # once the final version has been resampled to the original template raster
  tc_slice_fpath <- tempfile(
    pattern = sprintf("tc_slice_%i_%i_",
                      year,
                      month),
    fileext = ".tif"
  )

  # put the empty terraclimate template raster in this file
  terra::writeRaster(terraclimate_template,
                     tc_slice_fpath,
                     overwrite = TRUE)

  # loop through the tiles, pulling out the information for this slice and
  # putting it in the slice raster
  for (tile_fpath in vector_tile_fpaths) {

    # load the tile data and subset to the species of interest
    tile_data <- tile_fpath |>
      readRDS() |>
      dplyr::filter(
        species == species
      )

    # from the first pixel, find an index to the correct times
    first_pixel <- tile_data |>
      dplyr::pull(
        pixel_vectors
      ) |>
      dplyr::first()

    row_index <- first_pixel |>
      dplyr::mutate(
        matches_year = lubridate::year(start) == year,
        matches_month = lubridate::month(start) == month,
        matches_slice = matches_year & matches_month
      ) |>
      dplyr::pull(
        matches_slice
        ) |>
      which()

    # extract the required variable values, in the slice year and month, for
    # this tile
    values <- tile_data |>
      dplyr::select(
        pixel_vectors
      ) |>
      dplyr::mutate(
        pixel_index = dplyr::row_number(),
        pixel_vectors = lapply(
          pixel_vectors,
          function(pixel_vectors) {
            # subset to required row and column
            pixel_vectors[row_index, variable]
          }
        )
      ) |>
      tidyr::unnest(
        pixel_vectors
      ) |>
      # check they are in the original order
      dplyr::arrange(
        pixel_index
      ) |>
      dplyr::pull(
        !!variable
      )

    # read in the temporary slice raster
    tc_slice_rast <- terra::rast(tc_slice_fpath)

    # match coordinates to cells in template
    coords_mat <- tile_data |>
      dplyr::select(
        longitude,
        latitude
      ) |>
      as.matrix()

    cell_idx <- terra::cellFromXY(
      tc_slice_rast,
      coords_mat
    )

    # put these values in the raster
    tc_slice_rast[cell_idx] <- values

    # save it to disk again
    terra::writeRaster(tc_slice_rast,
                       tc_slice_fpath,
                       overwrite = TRUE)
  }

  # now we are done inserting all the tile date, read in the temporary slice
  # raster again
  tc_slice_rast <- terra::rast(tc_slice_fpath)

  # resample it to the original template
  slice_rast <- terra::resample(tc_slice_rast,
                                template,
                                method = "bilinear")

  # and save this to disk
  terra::writeRaster(slice_rast,
                     slice_fpath,
                     overwrite = TRUE)

  # delete the old one to free up disk space
  file.remove(tc_slice_fpath)

  # return something so we can all it with furrr:future_map() without returning
  # a raster
  invisible(TRUE)

}

# given the path to the directory where the tiled vector summaries are, pull out
# the different slices (years and months) from the first tile, and return as a
# tibble
get_slices <- function(vector_tile_dir = "processing/vectors") {

  vector_tile_files <- list.files(vector_tile_dir,
                                  full.names = TRUE)

  # get the dates from the first tile and format as distinct slices
  vector_tile_files[1] |>
    readRDS() |>
    dplyr::pull(pixel_vectors) |>
    dplyr::first() |>
    dplyr::select(start) |>
    dplyr::mutate(
      year = lubridate::year(start),
      month = lubridate::month(start),
      .before = everything()
    ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      slice = dplyr::row_number(),
      .before = everything()
    )

}

