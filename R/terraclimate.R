# internal functions to access and process terraclimate data, either at a single
# point or a tile to batch-process covariates

# note there is a google earth engine for terraclimate...
# https://developers.google.com/earth-engine/datasets/catalog/IDAHO_EPSCOR_TERRACLIMATE#bands

# get the URL to the terraclimate monthly summary data
terraclimate_url <- function(variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad")) {
  variable <- match.arg(variable)
  base_url <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_%s_1958_CurrentYear_GLOBE.nc#fillmismatch"
  sprintf(base_url, variable)
}

# open the connection to terraclimate monthly data
terraclimate_open <- function(variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad")) {
  url <- terraclimate_url(variable = variable)
  ncdf4::nc_open(url)
}

terraclimate_available_indices <- function() {

  # cache the available indices
  tc_avail_indices <- options()$.tc_avail_indices

  # download if required
  if (is.null(tc_avail_indices)) {

    # open a connection with the default cube (tmax), and close it again on
    # exiting this function
    con <- terraclimate_open()
    on.exit(ncdf4::nc_close(con))

    # find times and places for which there are data, and return
    tc_avail_indices <-list(
      times = ncdf4::ncvar_get(con, "time"),
      longitudes = ncdf4::ncvar_get(con, "lon"),
      latitudes = ncdf4::ncvar_get(con, "lat")
    )

    options(.tc_avail_indices = tc_avail_indices)
  }

  tc_avail_indices
}

# given a slice (index to initial time, lat and long, and number of elements in
# each of these dimensions), extract the array of the named variable from
# terraclimate
terraclimate_fetch <- function(
    slice,
    variable = c("tmax", "tmin", "ppt", "ws", "vpd", "srad")) {
  # open a connection with the required cube, and close it again on exiting
  # this function
  variable <- match.arg(variable)
  con <- terraclimate_open(variable)
  on.exit(ncdf4::nc_close(con))

  # fetch the required slice of data and return
  ncdf4::ncvar_get(con,
                   varid = variable,
                   start = slice$start,
                   count = slice$count)
}

# find the nearest lat and long in terra climate to the user's. Check the
# requested coordinate against the available set. If there isn't one within 1/48
# of a degree, get one slightly further (basing this on NicheMapR's logic)
nearest_coord <- function(value, available) {
  error <- abs(available - value)
  acceptable <- error < 1/48
  if (any(acceptable)) {
    index <- which(acceptable)
  } else {
    index <- which(error < 1/47.9)[1]
  }
  index
}

# all days since 1900 until the end of the year
days_since_1900 <- function() {
  start <- as.Date("1900-01-01")
  end <- Sys.Date()
  seq(start, end, "days")
}

# return a sequence of last days of the month, for all months encompassing the
# requested dates
monthly_summary_days <- function(dates, buffer_days = 183) {
  # start of month buffer_days before the first date
  min_date <- lubridate::floor_date(min(dates) - buffer_days, unit = "month")
  # start of month buffer_days after the last date
  max_date <- lubridate::floor_date(max(dates) + buffer_days, unit = "month")
  # sequence of all the first dates of the months in between
  seq(min_date, max_date, by = "1 month")
}

# error nicely if the user provides inaccessible dates
check_dates <- function(dates, terraclimate_available) {

  check_contiguous_dates(dates)

  # get the summary dates available and the time period they cover (the start of
  # the first date summarised to the end of the last)
  available_dates <- days_since_1900()[terraclimate_available$times+1]
  coverage_min <- min(available_dates)
  coverage_max <- lubridate::ceiling_date(max(available_dates),
                                          unit = "month") - 1

  # error if the request cannot be fulfilled
  if (min(dates) < coverage_min | max(dates) > coverage_max) {
    stop("the terraclimate archive covers ",
         coverage_min,
         " to ",
         coverage_max,
         "\n       you requested dates between ",
         min(dates),
         " and ",
         max(dates),
         call. = FALSE
    )
  }

}

# build a terraclimate slice for rectangular spatial extent, over the
# requested dates
terraclimate_build_slice_extent <- function(extent,
                                            dates) {

  # check the extent object
  if(!inherits(extent, "SpatExtent")) {
    stop("extent should be a terra SpatExtent object, e.g. created like this:
         extent <- ext(<xmin>, <xmax>, <ymin>, <ymax>)",
         call. = FALSE)
  }

  # times and places for which there are data in terraclimate
  available <- terraclimate_available_indices()

  # check the user-provided dates are contiguous, increasing, and we can get
  # them from terraclimate
  check_dates(dates, available)

  # the timeseries is summarised by month and indexed by the first day of the
  # month, so find all summary dates relevant to our target dates, with some
  # buffering either side to make sure the splining is consistent for those dates
  required_summary_dates <- monthly_summary_days(dates)

  # the terraclimate temporal index is integer days since the first day of 1900
  # (ie. indexed from 0, unlike R's indexing from 1, so subtract 1) so find our
  # month-start dates as indices
  time_index_1900 <- match(required_summary_dates, days_since_1900()) - 1

  # Remove any of these months that are not available. The dates themselves are
  # all included (checked above), so this is just removing buffer dates we can't
  # access.
  indices_keep <- time_index_1900 %in% available$times
  required_summary_dates <- required_summary_dates[indices_keep]
  time_index_1900 <- time_index_1900[indices_keep]

  # now index the elements of the NCDF array
  time_index <- match(time_index_1900, available$times)

  # Find the nearest terraclimate lat longs to the edges of the extent
  ext_list <- as.list(extent)
  xmin_index <- nearest_coord(ext_list$xmin, available$longitudes)
  xmax_index <- nearest_coord(ext_list$xmax, available$longitudes)
  # note: indices are ordered left -> right (like coordinates) and top -> bottom
  # (the reverse of coordinates), so we need to flip the latitudes
  ymin_index <- nearest_coord(ext_list$ymax, available$latitudes)
  ymax_index <- nearest_coord(ext_list$ymin, available$latitudes)

  xrange <- xmax_index - xmin_index + 1
  yrange <- ymax_index - ymin_index + 1

  # return a list with the initial time/place location and size of slice
  list(
    start = c(xmin_index, ymin_index, time_index[1]),
    count = c(xrange, yrange, length(time_index)),
    dates = list(
      start = required_summary_dates,
      end = lubridate::ceiling_date(required_summary_dates, unit = "1 month") -1
    ),
    longitudes = seq(ext_list$xmin, ext_list$xmax, length.out = xrange),
    # flip back the latitudes
    latitudes = seq(ext_list$ymax, ext_list$ymin, length.out = yrange)
  )

}

# build a terraclimate slice for a single point in space, over the requested
# dates
terraclimate_build_slice_point <- function(longitude,
                                           latitude,
                                           dates) {

  # times and places for which there are data in terraclimate
  available <- terraclimate_available_indices()

  # check the user-provided dates are contiguous, increasing, and we can get
  # them from terraclimate
  check_dates(dates, available)

  # the timeseries is summarised by month and indexed by the first day of the
  # month, so find all summary dates relevant to our target dates, with some
  # buffering either side to make sure the splining is consistent for those dates
  required_summary_dates <- monthly_summary_days(dates)

  # the terraclimate temporal index is integer days since the first day of 1900
  # (ie. indexed from 0, unlike R's indexing from 1, so subtract 1) so find our
  # month-start dates as indices
  time_index_1900 <- match(required_summary_dates, days_since_1900()) - 1

  # Remove any of these months that are not available. The dates themselves are
  # all included (checked above), so this is just removing buffer dates we can't
  # access.
  indices_keep <- time_index_1900 %in% available$times
  required_summary_dates <- required_summary_dates[indices_keep]
  time_index_1900 <- time_index_1900[indices_keep]

  # now index the elements of the NCDF array
  time_index <- match(time_index_1900, available$times)

  # Find the nearest terraclimate lat longs to the ones we want
  longitude_index <- nearest_coord(longitude, available$longitudes)
  latitude_index <- nearest_coord(latitude, available$latitudes)

  # return a list with the initial time/place location and size of slice
  list(
    start = c(longitude_index, latitude_index, time_index[1]),
    count = c(1, 1, length(time_index)),
    dates = list(
      start = required_summary_dates,
      end = lubridate::ceiling_date(required_summary_dates, unit = "1 month") -1
    )
  )

}

# given a template raster object, build a terraclimate-harmonised 5km template.
# Ie. a raster covering all the non-NA cells in the template raster whose cells
# align with those in terraclimate
make_terraclimate_template <- function(template) {

  # make template binary
  template <- template * 0 + 1

  # trim, expand, then buffer these by one cell, to ensure we can interpolate as
  # much as possible when resampling back to this raster from the terraclimate one
  template <- template |>
    terra::trim() |>
    terra::extend(1) |>
    terra::focal(w = 3,
                 fun = max,
                 na.rm = TRUE)

  # now find the equivalent terraclimate raster cells

  # find the grid setup in terraclimate
  terraclimate_available <- terraclimate_available_indices()

  # create an empty raster that aligns with terraclimate and has extent no smaller
  # than our raster
  dims <- as.list(ext(template))
  x_valid <- terraclimate_available$longitudes >= dims$xmin &
    terraclimate_available$longitudes <= dims$xmax
  y_valid <- terraclimate_available$latitudes >= dims$ymin &
    terraclimate_available$latitudes <= dims$ymax
  tc_ext <- ext(c(range(terraclimate_available$longitudes[x_valid]),
                  range(terraclimate_available$latitudes[y_valid])))

  # find the non-NA elements in this raster by building a slice with one date and
  # extracting one value
  raster_slice <- terraclimate_build_slice_extent(tc_ext,
                                                  dates = as.Date("2020-01-01"))
  # this function is set up to extract extra dates to ensure we can spline to
  # these dates, but we only need one value for this so manually switch to the one
  # date at the middle of the time slice
  timesteps <- raster_slice$count[3]
  raster_slice$start[3] <- raster_slice$start[3] + ceiling(timesteps / 2)
  raster_slice$count[3] <- 1
  # get the values
  vals <- terraclimate_fetch(raster_slice)

  # put them in a raster
  tc_template <- rast(nrows = sum(y_valid),
                      ncols = sum(x_valid),
                      extent = tc_ext,
                      vals = t(vals))

  # resample the target raster to this one
  new_template_mask <- terra::resample(template,
                                       tc_template,
                                       method = "max")

  # and use it to mask the terraclimate values, set values to 0, and return
  tc_template_masked <- terra::mask(tc_template, new_template_mask)
  tc_template_masked <- tc_template_masked * 0

  # now trim it again, and return
  tc_template <- terra::trim(tc_template_masked)
  tc_template

}

# given a vector defining an extent object, and a terraclimate template raster,
# return a logical for whether there are any non-na values in that extent
check_tile <- function(vec, template) {
  this_tile <- ext(vec)
  this_raster <- terra::crop(template, this_tile)
  length(cells(this_raster)) > 0
}

# given a vector defining an extent object, and a terraclimate template raster,
# return a vector for a trimmed extent object to the smallest rectangle
# including only non-NA cells
trim_tile <- function(vec, template) {

  # make a raster with this extent, trim it, and convert back ot an extent
  tile <- terra::ext(vec)
  raster <- terra::crop(template, tile)
  raster_trimmed <- terra::trim(raster)
  tile_trimmed <- terra::ext(raster_trimmed)

  # convert to a spatvector, negative buffer by a quarter of a cell so each
  # cell centroid is in only one tile, and convert back ot an extent
  buffer_distance <- -1 * res(template) / 4
  spatvector_trimmed <- vect(tile_trimmed)
  spatvector_trimmed_shrunk <- terra::buffer(spatvector_trimmed,
                                             buffer_distance)
  tile_trimmed_shrunk <- terra::ext(spatvector_trimmed_shrunk)

  # return as a vector
  as.vector(tile_trimmed_shrunk)
}

# given a template raster denoting non-NA cells for processing, and a target
# number of tiles, return the dimensions for approximately that number of tiles
# as a matrix with each row giving the tile extents (xmin, xmax, ymin, ymax).
# Tiles will only be returned that contain non-NA cells, and will be trimmed to
# include only rows and columns that contain non-NA cells.
make_tiles <- function(template, target_n_tiles = 100) {

  # work out the fraction of non-missing cells in the raster, and adjust up the
  # number of tiles
  n_cells_all <- prod(dim(template)[c(1:2)])
  n_cells_valid <- length(cells(template))
  prop_cells_valid <- n_cells_valid / n_cells_all
  target_n_tiles_adj <- target_n_tiles / prop_cells_valid

  # aggregate a raster into approximately this number of square tiles
  agg_factor <- round(sqrt(n_cells_all / target_n_tiles_adj))
  tiles_raster <- terra:::aggregate(template, agg_factor)
  tiles <- terra::getTileExtents(template, tiles_raster)

  # weed out tiles with no values in them
  tiles_valid <- apply(tiles,
                       MARGIN = 1,
                       FUN = check_tile,
                       template = template)
  tiles <- tiles[tiles_valid, ]

  # trim them down to only the valid cells, and shrink them all by less than half
  # a cell to prevent excessive processing of cells
  tiles <- t(apply(tiles,
                   MARGIN = 1,
                   FUN = trim_tile,
                   template = template))

  tiles
}

# given an integer 'tile_number' for the tile to extract (indexing tiles), a
# matrix 'tiles' of tile extents (each row giving the bounding box), a vector of
# 'dates' for which we need data, and a vector of 'variables' required (by
# default those we can convert o NicheMapR inputs), extract monthly terraclimate
# data on those variables, for all cells in that tile, for all available months,
# where possible extending beyond 'dates' to enable improved spline
# interpolation, and return as a tibble, identifying the cell by the latitude
# and longitude of its centroid.
terraclimate_extract_tile <- function(tile_number,
                                      tiles,
                                      dates,
                                      variables = c("tmax",
                                                    "tmin",
                                                    "ppt",
                                                    "ws",
                                                    "vpd",
                                                    "srad")) {

  # make an extent object for this tile
  tile <- ext(tiles[tile_number, ])

  # build the slice for all pixels (including NAs) in this tile
  tile_slice <- terraclimate_build_slice_extent(
    extent = tile,
    dates = dates
  )

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
    )

}

# clean extracted terraclimate tile data, by removing all cells that are outside
# the raster template
clean_tile_data <- function(tile_data, tc_template) {

  # get all coordinates in tile data
  all_coords <- tile_data |>
    dplyr::distinct(
      longitude, latitude
    )

  # find the valid ones
  valid_coords <- all_coords |>
    dplyr::mutate(
      value = terra::extract(
        tc_template,
        all_coords,
        ID = FALSE
      )[, 1]
    ) |>
    dplyr::filter(
      !is.na(value)
    ) |>
    dplyr::select(
      -value
    )

  # subset to these
  valid_coords |>
    dplyr::left_join(
      tile_data
    )

}


# Given a tibble of extracted monthly terraclimate data in wide format
# per-variable (i.e. rows for observations but separate columns for tmax, tmin,
# etc.), compute additional variables used by NicheMapR: min/max relative
# humidity, min/max windspeed, min/max cloud cover, and rename precipitation as
# rainfall (assumed no snow in areas we consider for mosquito population
# dynamics). This involves execution of clear_sky_radiation_12mo(), which
# involves running NicheMapR in solonly mode for each pixel, and can be slow. the user m
process_terraclimate_tile_vars <- function(tile_data) {

  # batch append the clearsky monthly synoptic data from the in-package
  # information

  # get the raster of tc-GADS hybrid cell numbers
  tc_gads <- terra::unwrap(tc_gads_wrapped)

  # find unique coordinates in this tile
  unique_coords <- tile_data |>
    dplyr::distinct(
      longitude,
      latitude
    )

  # add on the cell ID
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

  # process terraclimate data into format required for NicheMapR micro function
  tile_data |>
    # add a midpoint date, for looking up solar radiation max/min, and later
    # for splining, plus a number of days for computing daily rainfall
    dplyr::mutate(
      mid_date = start + (end - start) / 2,
      n_days = as.numeric(end - start),
      .before = start
    ) |>
    # pivot wider to compute things
    tidyr::pivot_wider(
      names_from = variable,
      values_from = value
    ) |>
    # compute relative humidity max/min, by first computing vapour pressure
    # from mean temp and vp deficit then computing RH at max/min temps
    dplyr::mutate(
      tmean = (tmax + tmin) / 2,
      vp = vapour_pressure(tmean, vpd),
      rhmax = relative_humidity(tmin, vp),
      rhmin = relative_humidity(tmax, vp)
    ) |>
    # adjust wind speed from measured height (10m) to match the other inputs
    # (2m) and split by max/min (using NicheMapR's ad-hoc adjusment for min)
    dplyr::mutate(
      wsmax = adjust_wind_speed(ws, height_required = 2, height_measured = 10),
      wsmin = wsmax * 0.1
    ) |>
    # append the expected clear sky radiation (microclimate code using
    # coordinates and aerosol data) and use this to compute cloud cover
    # percentage from solar radiation (terraclimate) and and manual multipliers
    # to enforce daily variation, as used in NicheMapR
    dplyr::mutate(
      month_id = lubridate::month(mid_date)
    ) |>
    dplyr::left_join(
      clearsky_lookup,
      by = c("latitude", "longitude", "month_id")
    ) |>
    # compute cloud cover
    dplyr::mutate(
      ccmax = cloud_cover(srad, clear_sky_SOLR, multiplier = 1.2),
      ccmin = cloud_cover(srad, clear_sky_SOLR, multiplier = 0.8)
    ) |>
    # convert monthly to daily precipitation (assumed rainfall in non-snowfall areas
    # we care about)
    dplyr::mutate(
      rainfall_daily = ppt / n_days
    ) |>
    # group by lat and long and store the monthly data we need as a list column
    dplyr::group_by(
      longitude,
      latitude
    ) |>
    dplyr::summarise(
      monthly_climate = list(
        tibble::tibble(
          start,
          end,
          mid_date,
          tmin, tmax,
          rhmin, rhmax,
          wsmin, wsmax,
          ccmin, ccmax,
          rainfall_daily
        )
      ),
      .groups = "drop"
    )

}
