library(dplyr)
library(lubridate)

setwd("~/git/snfa/macro_data")

rm(list = ls())

df.Y = read.csv("Y.csv", stringsAsFactors = FALSE)
df.K = read.csv("K.csv", stringsAsFactors = FALSE)
df.K.price = read.csv("K_price.csv", stringsAsFactors = FALSE)
df.L = read.csv("L.csv", stringsAsFactors = FALSE)
df.L.price = read.csv("L_price.csv", stringsAsFactors = FALSE)
df.gdpdef = read.csv("GDPDEF.csv", stringsAsFactors = FALSE)

names(df.Y) = c("DATE", "Y")
names(df.K) = c("DATE", "K")
names(df.K.price) = c("DATE", "K.price")
names(df.L) = c("DATE", "L")
names(df.L.price) = c("DATE", "L.price")
names(df.gdpdef) = c("DATE", "gdpdef")

df.Y$Year = year(df.Y$DATE)
df.K$Year = year(df.K$DATE)
df.K.price$Year = year(df.K.price$DATE)
df.L$Year = year(df.L$DATE)
df.L.price$Year = year(df.L.price$DATE)
df.gdpdef$Year = year(df.gdpdef$DATE)

USMacro = df.Y %>%
  full_join(df.K, by = c("DATE", "Year")) %>%
  full_join(df.K.price, by = c("DATE", "Year")) %>%
  full_join(df.L, by = c("DATE", "Year")) %>%
  full_join(df.L.price, by = c("DATE", "Year")) %>%
  full_join(df.gdpdef, by = c("DATE", "Year"))

#Adjust to 2010 dollars
gdpdef.1983 = unlist(USMacro %>% filter(Year == 1983) %>% select(gdpdef))
gdpdef.2009 = unlist(USMacro %>% filter(Year == 2009) %>% select(gdpdef))
gdpdef.2010 = unlist(USMacro %>% filter(Year == 2010) %>% select(gdpdef))
gdpdef.2011 = unlist(USMacro %>% filter(Year == 2011) %>% select(gdpdef))

USMacro$Y = USMacro$Y * (gdpdef.2010 / gdpdef.2009)
USMacro$K = USMacro$K * (gdpdef.2010 / gdpdef.2011)
USMacro$L.price = USMacro$L.price * (gdpdef.2010 / gdpdef.1983)

#Convert to billions of dollars annually
USMacro$K = USMacro$K / 1000

#Convert labor price to be price for 1000 workers annually
USMacro$L.price = USMacro$L.price * 52 * 1000

#Convert capital price to be price of $1 billion capital annually
USMacro$K.price = 1e9 * (USMacro$K.price + 100) / 100

#Select variables
USMacro = USMacro %>%
  select(Year, Y, K, K.price, L, L.price)

#Save
USMacro = as.data.frame(USMacro)
save(USMacro, file = "../snfa/data/USMacro.RData")
