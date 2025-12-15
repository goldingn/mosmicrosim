# a lot of overhead in running NicheMapR comes from querying GADS aerosol data
# to estimate solar radiation. NicheMapR can do this using Fortran or a slower R
# approximation. However, the underlying data for this are on a course grid, so
# pre-computing a looking up values is faster, especially when we'll be running
# this for a lot of cells. Here we create a raster to index those cells, and a
# lookup table to return the solar attenuation data at those lookups. Functions
# in solar_attenuation.R query this lookup whenever solar attenuation data are
# needed.

library(tidyverse)
library(terra)
library(furrr)

# create an empty raster to look up GADS data

# the dimensions of the underlying GADS raster are stated in NicheMapR's R
# implementation of GADS
ye_gads <- rast(
  nrows = 36,
  ncols = 71,
  xmin = -180,
  xmax = 175,
  ymin = -90,
  ymax = 90
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
  # add on cell IDs
  mutate(
    cell_id = terra::extract(
      ye_gads,
      cells(ye_gads)
      # ID = FALSE
    )[, 1],
    .before = everything()
  ) |>
  mutate(
    solar_attenuation = future_map2(.x = latitude,
                                    .y = longitude,
                                    .f = solar_attenuation_gads)
  ) |>
  rowwise() |>
  # add on wavelengths, since map2 is apparently disposing of them
  mutate(
    wavelength_nm = list(wavelengths_nm),
    .before = solar_attenuation
  ) |>
  unnest(
    c(wavelength_nm, solar_attenuation)
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
