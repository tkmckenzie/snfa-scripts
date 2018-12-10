library(dplyr)
library(ggplot2)

# http://cran-logs.rstudio.com/

setwd("~/git/snfa-scripts/cran_downloads")

rm(list = ls())

df = read.csv("log_filter.csv")

country.df = df %>%
  group_by(country) %>%
  summarize(count = n())
ggplot(country.df, aes(country, count)) +
  geom_bar(stat = "identity") +
  theme(axis.text = element_text(angle = 90, hjust = 0, vjust = 0.5))

day.df = df %>%
  group_by(date) %>%
  summarize(count = n()) %>%
  mutate(date = as.Date(date))
ggplot(day.df, aes(date, count)) +
  geom_line()