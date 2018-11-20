setwd("~/git/snfa/generated_data")

set.seed(100)

N <- 50
x <- runif(N, 10, 100)
y <- sapply(x, function(x) 500 * x^0.25 - dnorm(x, mean = 70, sd = 10) * 8000) - abs(rnorm(N, sd = 20))
y <- y - min(y) + 10
df <- data.frame(x, y)

univariate <- df
save(univariate, file = "../snfa/data/univariate.RData")
