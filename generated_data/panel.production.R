setwd("~/git/snfa/generated_data")

rm(list = ls())

set.seed(100)

num.firms <- 20
num.inputs <- 3
num.years <- 10

beta <- runif(num.inputs, 0, 1)
TFP.trend = 0.25
TFP <- cumsum(rnorm(num.years)) + TFP.trend * (1:num.years)

sd.measurement <- 0.05
sd.inefficiency <- 0.01

f <- function(X){
  return(TFP + X %*% beta)
}
gen.firm.data <- function(i){
  X = matrix(runif(num.years * num.inputs, 1, 10), ncol = num.inputs)
  y = f(X) +
    rnorm(num.years, sd = sd.measurement) -
    abs(rnorm(num.years, sd = sd.inefficiency))
  firm.df <- data.frame(Firm = i,
                        Year = 1:num.years,
                        X = exp(X),
                        y = exp(y))
}

panel.production = Reduce(rbind, lapply(1:num.firms, gen.firm.data))
panel.production$Firm = as.factor(panel.production$Firm)

save(panel.production, file = "../snfa/data/panel.production.RData")
