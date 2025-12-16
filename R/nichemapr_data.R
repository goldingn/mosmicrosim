# internal functions to access NicheMapR's internal global climate data

# download the NicheMapR global climate zipfile
download_nichemapr_data <- function(zip_file_path) {

  # download the zipfile to the zip file path
  file_url <- "https://github.com/mrke/NicheMapR/releases/download/v2.0.0/global.climate.zip"
  download.file(file_url, zip_file_path, mode = "wb")

}

# unzip the NicheMapR global climate zipfile to an NCDF in the R session's
# temp directory
unzip_nichemapr_data <- function(zip_file_path, ncdf_file_path) {
  unzip(zip_file_path,
        files = basename(ncdf_file_path),
        exdir = dirname(ncdf_file_path))
}

# get the NicheMapR global climate data object, downloading it if needed
get_nichemapr_data <- function() {

  # work out where the zipfile should be (wither NicheMapR extdata folder or a
  # tempfile) and where the unzipped ncdf file should be (a tempfile)
  zip_file_path <- nichemapr_data_filename()
  ncdf_file_path <- file.path(tempdir(), "global_climate.nc")

  # download it if it doesn't exist
  if (!file.exists(zip_file_path)) {
    download_nichemapr_data(zip_file_path)
  }

  # unzip it if it doesn't exist
  if (!file.exists(ncdf_file_path)) {
    unzip_nichemapr_data(zip_file_path, ncdf_file_path)
  }

  # read it into memory
  global_climate <- terra::rast(ncdf_file_path)

  # give it names matching definition in NicheMapR micro_global, to query later,
  names(global_climate) <- c(
    "altitude",
    paste0("rainfall", 1:12),
    paste0("rainy_days", 1:12),
    paste0("windspeed", 1:12),
    paste0("temperature_min", 1:12),
    paste0("temperature_max", 1:12),
    paste0("rh_min", 1:12),
    paste0("rh_max", 1:12),
    paste0("cloud_cover", 1:12)
  )

  global_climate

}

# load the NicheMapR global climate raster
global_climate_file <- system.file("extdata",
                                   "nichemapr_global_climate.tif",
                                   package = "mosmicrosim")
global_climate <- terra::rast(
  global_climate_file
)
# wrap it, so it can be serialised (e.g. passed to future) without losing the
# pointer
global_climate_wrapped <- terra::wrap(global_climate)

# get and cache the NicheMapR global climate data object, using a tempfile
get_cloud_cover_raster <- function() {

  nichemapr_data <- terra::unwrap(global_climate_wrapped)

  # pull out the cloud cover layers
  cloud_cover_data <- nichemapr_data[[paste0("cloud_cover", 1:12)]]

  # divide by 10 to get percentages
  cloud_cover_data / 10

}

get_altitude_raster <- function() {

  nichemapr_data <- terra::unwrap(global_climate_wrapped)

  # pull out the altitude layer
  altitude_data <- nichemapr_data[["altitude"]]

  altitude_data

}

# lookup altitude in metres, defaulting to the nichemapr version
altitude_m <- function(longitude, latitude,
                       altitude_raster = get_altitude_raster()) {
  altitude_data <- terra::extract(altitude_raster,
                                  data.frame(longitude, latitude),
                                  ID = FALSE)
  altitude_data[, 1]
}
