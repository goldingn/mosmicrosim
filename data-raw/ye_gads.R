# a lot of overhead in running NicheMapR comes from querying GADS aerosol data
# to estimate solar radiation. NicheMapR can do this using Fortran or a slower R
# approximation. However, the underlying data for this are on a course grid, so
# pre-computing a looking up values is faster, especially when we'll be running
# this for a lot of cells. Here we create a raster to index those cells, and a
# lookup table to return the solar attenuation data at those lookups. Functions
# in solar_attenuation.R query this lookup whenever solar attenuation data are
# needed.

library(NicheMapR)
library(tidyverse)
library(terra)
library(furrr)

source("R/solar_attenuation.R")

# create an empty raster to look up GADS data

# GADS data are defined at these nodes, as stated at the top ofs NicheMapR's R
# implementation of GADS
res <- 5
longitudes <- seq(-180, 175, by = res)
latitudes <- seq(-90, 90, by = res)

# NicheMapR finds the nearest node to a given set of coordinates to get the
# corresponding solar attenuation information. So make our raster a grid with
# these as the centroids
pad <- res / 2

ye_gads <- rast(
  nrows = length(latitudes),
  ncols = length(longitudes),
  xmin = min(longitudes) - pad,
  xmax = max(longitudes) + pad,
  ymin = min(latitudes) - pad,
  ymax = max(latitudes) + pad
)

# put the cell number in each cell, so we can use extract to do a lookup
ye_gads[] <- seq_len(ncell(ye_gads))

names(ye_gads) <- "cell_id"

# now compute solar attenuation at each cell in this raster, computed in
# parallel using future (takes a minute or so)
plan(multisession, workers = 8)

# compute solar attenuation for each cell, and store as a list column
solar_attenuation_lookup <- terra::xyFromCell(
  ye_gads,
  cells(ye_gads)
) |>
  as_tibble() |>
  rename(
    longitude = x,
    latitude = y,
  ) |>
  # add on cell IDs (the slow way to make sure they are in the right order)
  mutate(
    cell_id = terra::extract(
      ye_gads,
      cells(ye_gads)
    )[, 1],
    .before = everything()
  ) |>
  # compute the solar attenuation curves
  mutate(
    solar_attenuation = future_map2(.x = latitude,
                                    .y = longitude,
                                    .f = solar_attenuation_gads)
  ) |>
  # drop the coordinates, now we have the lookup
  select(
    -latitude,
    -longitude
  )

# save the lookup in sysdata.R (for internal package use)
usethis::use_data(solar_attenuation_lookup,
                  internal = TRUE,
                  overwrite = TRUE,
                  compress = "xz")

# save the gads raster in inst/extdata
terra::writeRaster(ye_gads,
                   "inst/extdata/ye_gads.tif",
                   overwrite = TRUE)
