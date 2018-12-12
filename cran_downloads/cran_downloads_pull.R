# http://cran-logs.rstudio.com/

setwd("~/git/snfa-scripts/cran_downloads")

rm(list = ls())

start = as.Date('2018-12-01')
today = as.Date('2018-12-11')

all.days = as.character(seq(start, today, by = 'day'))
existing.days = gsub(".csv.gz$", "", grep(".csv.gz$", list.files("logs"), value = TRUE))
missing.days = setdiff(all.days, existing.days)

for (day in missing.days){
  year = as.POSIXlt(day)$year + 1900
  url = paste0('http://cran-logs.rstudio.com/', year, '/', day, '.csv.gz')
  out.file = paste0("logs/", day, ".csv.gz")
  download.file(url, out.file)
}

