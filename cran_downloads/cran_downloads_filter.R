library(dplyr)

# http://cran-logs.rstudio.com/

setwd("~/git/snfa-scripts/cran_downloads")

rm(list = ls())

package.name = "snfa"

log.files = list.files("logs")

df = read.csv(gzfile(paste0("logs/", log.files[1])))
df = df %>% filter(package == package.name)

for (log.file in log.files[-1]){
  temp.df = read.csv(gzfile(paste0("logs/", log.file)))
  temp.df = temp.df %>% filter(package == package.name)
  df = rbind(df, temp.df)
  rm(temp.df)
}

write.csv(df, "log_filter.csv", row.names = FALSE)
