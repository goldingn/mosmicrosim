# convert a date vector to a numeric day of the year, indexing from 1
doy <- function(date) {
  as.numeric(format(dates, "%j"))
}
