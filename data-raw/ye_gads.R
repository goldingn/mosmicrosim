# A lot of overhead in running NicheMapR comes from querying GADS aerosol data
# to estimate solar radiation, and then using it to compute clear sky radiation
# for given pixels. NicheMapR can query GADS in its native 5 degree resolution
# to obtain solar attenuation curves, either using Fortran (microclimate
# function in solonly mode) or a slower R approximation (gads.r). NicheMapR can
# then use microclimate to compute clear sky solar radiation (5 degree
# horizontal resolution, user-provided latitude) from this. However, since the
# underlying GADS data for this are on such a course grid, pre-computing and
# then looking up solar attenuation and clear sky radiation values is faster,
# especially when we'll be running this for a lot of cells. We do this here, and
# store the precomputed values in this package for later use.

# For solar attenuation we create a raster to index the GADS cells, and a lookup
# table to return the solar attenuation vector data at each of those cells.
# Functions in solar_attenuation.R query this lookup whenever solar attenuation
# data are needed.

# For clear sky radiation, we precompute monthly synoptic values of total
# daytime clear sky solar radiation at locations that match the native
# resolution of the terraclimate data. We do this by lowering the horizontal
# resolution of these data to match GADS to speed up computation and reduce file
# storage.

library(NicheMapR)
library(tidyverse)
library(terra)
library(furrr)
library(mosmicrosim)

# source("R/create_micro.R")
# source("R/solar_attenuation.R")
# source("R/clear_sky_radiation.R")
# source("R/utils.R")

# create an empty raster to index GADS data locations

# GADS data are defined at these nodes, as stated at the top ofs NicheMapR's R
# implementation of GADS
gads_res <- 5
gads_longitudes <- seq(-180, 180, by = gads_res)
gads_latitudes <- seq(-90, 90, by = gads_res)

# NicheMapR finds the nearest node to a given set of coordinates to get the
# corresponding solar attenuation information. So make our raster a grid with
# these as the centroids
gads_pad <- gads_res / 2

ye_gads <- rast(
  nrows = length(gads_latitudes),
  ncols = length(gads_longitudes),
  xmin = min(gads_longitudes) - gads_pad,
  xmax = max(gads_longitudes) + gads_pad,
  ymin = min(gads_latitudes) - gads_pad,
  ymax = max(gads_latitudes) + gads_pad
)

# put the cell number in each cell, so we can use extract to do a lookup
ye_gads[] <- seq_len(ncell(ye_gads))

names(ye_gads) <- "cell_id"

# create a lookup from these cell indices to the coordinates at the cell
# centroids. These centroids will be used to access the fortran solar
# attenuation data.

gads_lookup_coords <- terra::xyFromCell(
  ye_gads,
  cells(ye_gads)
) |>
  dplyr::as_tibble() |>
  dplyr::rename(
    longitude = x,
    latitude = y,
  ) |>
  # add on cell IDs (the slow way to make sure they are in the right order)
  dplyr::mutate(
    cell_id = terra::extract(
      ye_gads,
      cells(ye_gads)
    )[, 1],
    .before = everything()
  )

# now compute solar attenuation at each cell in this raster, computed in
# parallel using future (takes a few minutes)
plan(multisession, workers = 8)

# compute solar attenuation for each cell, and store as a list column
solar_attenuation_lookup <- gads_lookup_coords |>
  # compute the solar attenuation curves
  dplyr::mutate(
    solar_attenuation = future_map2(.x = latitude,
                                    .y = longitude,
                                    .f = mosmicrosim:::solar_attenuation_gads)
  ) |>
  # drop the coordinates, now we have the lookup
  dplyr::select(
    -latitude,
    -longitude
  )

# now make a raster for terraclimate's vertical (latitudinal) resolution, and
# approximately (ie. integer multiples of terraclimate resolution) GADS'
# horizontal (longitudinal) resolution

# first make the terraclimate native resolution raster
terraclimate_available <- terraclimate_available_indices()
tc_longitudes <- terraclimate_available$longitudes
tc_latitudes <- terraclimate_available$latitudes

tc_horiz_res <- abs(mean(diff(tc_longitudes)))
tc_vert_res <- abs(mean(diff(tc_latitudes)))

tc_horiz_pad <- tc_horiz_res / 2
tc_vert_pad <- tc_vert_res / 2

tc_native <- rast(
  nrows = length(tc_latitudes),
  ncols = length(tc_longitudes),
  xmin = min(tc_longitudes) - tc_horiz_pad,
  xmax = max(tc_longitudes) + tc_horiz_pad,
  ymin = min(tc_latitudes) - tc_vert_pad,
  ymax = max(tc_latitudes) + tc_vert_pad
)

tc_native[] <- seq_len(ncell(tc_native))
names(tc_native) <- "cell_id"

# now make a version with GADs' horizontal resolution
tc_gads <- rast(
  nrows = length(tc_latitudes),
  ncols = length(gads_longitudes),
  xmin = min(gads_longitudes) - gads_pad,
  xmax = max(gads_longitudes) + gads_pad,
  ymin = min(tc_latitudes) - tc_vert_pad,
  ymax = max(tc_latitudes) + tc_vert_pad
)

tc_gads[] <- seq_len(ncell(tc_gads))
names(tc_gads) <- "cell_id"

# now compute synoptic monthly clear sky solar radiation for each cell in
# tc_gads

tc_gads_lookup_coords <- terra::xyFromCell(
  tc_gads,
  cells(tc_gads)
) |>
  dplyr::as_tibble() |>
  dplyr::rename(
    longitude = x,
    latitude = y,
  ) |>
  # add on cell IDs (the slow way to make sure they are in the right order)
  dplyr::mutate(
    cell_id = terra::extract(
      tc_gads,
      cells(tc_gads)
    )[, 1],
    .before = everything()
  )

plan(multisession, workers = 8)

# compute monthly synoptic celar sky radiation values for the world, at hybrid
# GADS and terraclimate resolution. Takes about 30 mins in parallel on my MBP
clear_sky_lookup <- tc_gads_lookup_coords |>
  dplyr::mutate(
    clear_sky_SOLR = future_map2(.x = latitude,
                                 .y = longitude,
                                 .f = mosmicrosim:::clear_sky_radiation_12mo_nichemapr,
                                 # NicheMapR seems to be using a random seed
                                 # somewhere, and future does not like it. But
                                 # these calculations are deterministic
                                 .options = furrr_options(seed = NULL))
  ) |>
  # drop the coordinates, now we have the lookup
  dplyr::select(
    -latitude,
    -longitude
  )

# NOTE: the above code is circular and non-reproducible. Solar attenuation
# lookup needs to be saved and stored in the package before clearsky can run.
# However they need to be run together to go into the same internal sysdata.rda
# file. to work around this, just comment out the clearsky_lookup save on the
# first run through, rebuild the package, and on the second run through either
# rebuild both, or skip running solar_attenuation_lookup the second time and
# reload with:
#   solar_attenuation_lookup <- mosmicrosim:::solar_attenuation_lookup
# before re-saving the files below.

# save the solar attenuation and clear sky lookups in sysdata.R for internal
# package use in lookup functions
usethis::use_data(solar_attenuation_lookup,
                  clear_sky_lookup,
                  internal = TRUE,
                  overwrite = TRUE,
                  compress = "xz")

# save the GADS raster, the tc native raster, and the tc-GADS hybrid raster in
# inst/extdata to access later
terra::writeRaster(ye_gads,
                   "inst/extdata/ye_gads.tif",
                   overwrite = TRUE)

terra::writeRaster(tc_native,
                   "inst/extdata/tc_native.tif",
                   overwrite = TRUE)

terra::writeRaster(tc_gads,
                   "inst/extdata/tc_gads.tif",
                   overwrite = TRUE)
