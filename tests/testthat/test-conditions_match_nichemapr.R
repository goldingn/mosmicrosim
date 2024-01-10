test_that("manual inputs match nichemap outputs", {

  # use a nichemapr function to get macroclimatic inputs
  # compare them with the results of the equivalent inbuilt nichemapr function
  # use default
  library(NicheMapR)

  # Central Perth WA for 2021-2022
  loc <- c(115.86, -31.95)
  year_start <- 2021
  year_end <- 2022

  # extract macrocliamtic data from terraclimate (writing to a CSV) and run
  # nichemapr microclimate model on it

  # temporarily move us to a temporary directory, so that nichemapr writes these
  # files their

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
                                    # spread out rainfall over months and days
                                    evenrain = 1,
                                    rainfrac = 0,
                                    snowmodel = 0,
                                    runmoist = 0,
                                    runshade = 0,
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

  res <- NicheMapR::microclimate(micro)

  # n_years <- length(year_start:year_end)
  # # most elements of the input give daily values
  # length(micro$TMAXX) / n_years
  # # each row of the output gives hourly values
  # nrow(res$metout) / (365 * n_years)

  # check they yield approximately the same microclimate forcing parameters
  # (assuming a loss of precision due to writing and reading the parameters) we
  # use metout variables: TAREF (air temp in C) RH (relative humidity in %) VLOC
  # (wind speed in m/s)
  cols <- c("TAREF", "RH", "VLOC")

  # compute the absolute difference as a percentage of the value
  perc_diff <- function(a, b) {
    100 * abs(a - b) / abs(a)
  }

  metout_diff <- perc_diff(res_nmr$metout[, cols], res$metout[, cols])

  # ensure it is very small (less than 0.001%)
  testthat::expect_lt(max(metout_diff), 1e-3)

})
