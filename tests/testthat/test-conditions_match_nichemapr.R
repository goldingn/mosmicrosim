test_that("manual inputs match nichemap outputs", {

  # 1. use NicheMapR to get macroclimatic inputs and run the model

  # 2. compare with the same model applied to the reconsituted 'micro' config
  # file to check we are able to reload all the data

  # 3. pass tthe daily climate data into our new shimming function and run
  # microclimate to check it passes all arguments correctly

  library(NicheMapR)

  # Get weather for central Perth WA for 2021-2022
  loc <- c(115.86, -31.95)
  year_start <- 2021
  year_end <- 2022

  # 1. using NicheMapR, extract macroclimatic data from terraclimate, write the
  # microclimate config file data to CSVs, and run NicheMapR microclimate model
  # on it to define results to compare against

  # temporarily move us to a temporary directory, so that nichemapr writes these
  # files there

  # where we started
  owd <- getwd()
  # when we leave this function, return us there
  on.exit(setwd(owd))
  # put us in a new, temporary working directory
  nwd <- tempdir()
  setwd(nwd)

  # note: NicheMapR must be loaded in the namespace to find the table
  # CampNormTbl9_1, since it appears not to be visible within the package.
  res_nmr <- NicheMapR::micro_terra(loc = loc,
                                    ystart = year_start,
                                    yfinish = year_end,
                                    # daily data (highest resolution)
                                    timeinterval = 365,
                                    # rain falls evenly over the month
                                    rainfrac = 0,
                                    # rain falls evenly over the day
                                    evenrain = 1,
                                    # don't run the snow model
                                    snowmodel = 0,
                                    # run the moisture model
                                    runmoist = 1,
                                    # 50% shade (max has to be higher, even
                                    # though it isn't used)
                                    minshade = 0,
                                    maxshade = 100,
                                    # force the intercept to match minshade, not
                                    # maxshade
                                    intercept = 0.5 * 0.3,
                                    # don't run it twice for different shade
                                    # levels
                                    runshade = 0,
                                    # organism at 1m
                                    Usrhyt = 1,
                                    # use R version of GADS bc of stochastically
                                    # crashing fortran version
                                    run.gads = 2,
                                    # write extracted macroclimate data to CSV
                                    # files
                                    write_input = 1)

  # now reload macroclimatic conditions and push them into our function
  files <- list.files("micro csv input",
                      full.names = TRUE)

  load_micro <- function(filename) {
    table <- read.csv(filename)
    table_sub <- table[, -1]
    if (inherits(table_sub, "data.frame")) {
      table_sub <- as.matrix(table_sub)
    }
    table_sub
  }
  micro_data <- lapply(files, load_micro)
  names(micro_data) <- gsub(".csv$", "", basename(files))

  # manually reenter and tidy these as the names and order change
  micro <- with(micro_data,
                list(tides = tides,
                     microinput = microinput,
                     doy = doy,
                     SLES = SLES,
                     DEP = DEP,
                     Nodes = Nodes,
                     MAXSHADES = Maxshades,
                     MINSHADES = Minshades,
                     TMAXX = TMAXX,
                     TMINN = TMINN,
                     RHMAXX = RHMAXX,
                     RHMINN = RHMINN,
                     CCMAXX = CCMAXX,
                     CCMINN = CCMINN,
                     WNMAXX = WNMAXX,
                     WNMINN = WNMINN,
                     TAIRhr = TAIRhr,
                     RHhr = RHhr,
                     WNhr = WNhr,
                     CLDhr = CLDhr,
                     SOLRhr = SOLRhr,
                     RAINhr = RAINhr,
                     ZENhr = ZENhr,
                     IRDhr = IRDhr,
                     REFLS = REFLS,
                     PCTWET = PCTWET,
                     soilinit = soilinit,
                     hori = hori,
                     TAI = TAI,
                     soilprops = soilprop,
                     moists = moists,
                     RAINFALL = rain,
                     tannulrun = tannulrun,
                     PE = PE,
                     KS = KS,
                     BB = BB,
                     BD = BD,
                     DD = DD,
                     L = L,
                     LAI = LAI)
  )
  # put us back in the right working directory
  setwd(owd)

  # 2. rerun directly with this rebuilt config file, to check we have access to
  # all the right objects and arguments
  res <- NicheMapR::microclimate(micro)

  # 3. now pass in the main ones to our new shimming function with a simpler
  # interface to see if we can successfully recreate the outputs

  # recreate the dates
  dates <- as.Date(sprintf("%s-01-01", year_start)) + seq_along(micro$doy) - 1

  # extract the altitude (this is the right element)
  altitude <- micro$microinput[21]

  # now rerun via our shimming function:
  micro_new <- create_micro(
    latitude = loc[2],
    longitude = loc[1],
    altitude_m = altitude,
    dates = dates,
    daily_temp_max_c = micro$TMAXX,
    daily_temp_min_c = micro$TMINN,
    daily_rh_max_perc = micro$RHMAXX,
    daily_rh_min_perc = micro$RHMINN,
    daily_cloud_max_perc = micro$CCMAXX,
    daily_cloud_min_perc = micro$CCMINN,
    daily_wind_max_ms = micro$WNMAXX,
    daily_wind_min_ms = micro$WNMINN,
    daily_rainfall_mm = micro$RAINFALL)

  res_new <- NicheMapR::microclimate(micro_new)

  # check they yield approximately the same microclimate forcing parameters
  # (assuming a loss of precision due to writing and reading the parameters) we
  # use metout variables: TAREF (air temp in C) RH (relative humidity in %) VLOC
  # (wind speed in m/s)
  cols <- c("TAREF", "RH", "VLOC")

  # compute the absolute difference as a percentage of the value
  perc_diff <- function(a, b) {
    100 * abs(a - b) / abs(a)
  }

  res_diff <- perc_diff(res_nmr$metout[, cols], res$metout[, cols])
  new_diff <- perc_diff(res_nmr$metout[, cols], res_new$metout[, cols])

  # ensure it is very small (less than 0.001%)

  # check we can rebuild the original microclimate config file
  testthat::expect_lt(max(res_diff), 1e-3)

  # check we can pass it all in correctly via our shimming function and defaults
  testthat::expect_lt(max(new_diff), 1e-3)

})
