# convert a date vector to a numeric day of the year, indexing from 1
doy <- function(date) {
  as.integer(format(date, "%j"))
}

# check the dates are contiguous, increasing, and non-repeating
check_contiguous_dates <- function(dates) {

  contiguous_dates <- identical(dates,
                                dates[1] + seq_along(dates) - 1)

  if (!contiguous_dates) {
    stop("dates must be contiguous, increasing, and non-repeating ",
         "(a vector of consecutive days)",
         call. = FALSE)
  }

}
