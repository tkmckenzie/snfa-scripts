#' Randomly generated panel of production data
#' 
#' A dataset for illustrating technical and efficiency
#' changes using smooth non-parametric frontiers.
#' 
#' @format A data frame with 200 observations of six variables.
#' \describe{
#'   \item{Firm}{Firm identifier}
#'   \item{Year}{Year of observation}
#'   \item{X.1}{Input 1}
#'   \item{X.2}{Input 2}
#'   \item{X.3}{Input 3}
#'   \item{y}{Output}
#' }
#' @details 
#' Generated with the following code:
#' \preformatted{
#' set.seed(100)
#' 
#' num.firms <- 20
#' num.inputs <- 3
#' num.years <- 10
#' 
#' beta <- runif(num.inputs, 0, 1)
#' TFP.trend = 0.25
#' TFP <- cumsum(rnorm(num.years)) + TFP.trend * (1:num.years)
#' 
#' sd.measurement <- 0.05
#' sd.inefficiency <- 0.01
#' 
#' f <- function(X){
#'   return(TFP + X %*% beta)
#' }
#' gen.firm.data <- function(i){
#'   X = matrix(runif(num.years * num.inputs, 1, 10), ncol = num.inputs)
#'   y = f(X) +
#'     rnorm(num.years, sd = sd.measurement) -
#'     abs(rnorm(num.years, sd = sd.inefficiency))
#'   firm.df <- data.frame(Firm = i,
#'                         Year = 1:num.years,
#'                         X = exp(X),
#'                         y = exp(y))
#' }
#' 
#' panel.production = Reduce(rbind, lapply(1:num.firms, gen.firm.data))
#' panel.production$Firm = as.factor(panel.production$Firm)
#' }
#' 
"panel.production"

#' Randomly generated univariate data
#' 
#' A dataset for illustrating univariate non-parametric boundary
#' regressions and various constraints.
#' 
#' @format A data frame with 50 observations of two variables.
#' \describe{
#'   \item{x}{Input}
#'   \item{y}{Output}
#' }
#' @details 
#' Generated with the following code:
#' \preformatted{
#' set.seed(100)
#'
#' N <- 50
#' x <- runif(N, 10, 100)
#' y <- sapply(x, function(x) 500 * x^0.25 - dnorm(x, mean = 70, sd = 10) * 8000) - abs(rnorm(N, sd = 20))
#' y <- y - min(y) + 10
#' df <- data.frame(x, y)
#' }
#' 
"univariate"

#' US Macroeconomic Data
#' 
#' A dataset of real output, labor force, capital stock,
#' wages, and interest rates for the U.S. between 1929 and 2014, 
#' as available. All nominal values converted to 2010 U.S. dollars
#' using GDP price deflator.
#' 
#' @format A data frame with 89 observations of four variables.
#' \describe{
#'   \item{Year}{Year}
#'   \item{Y}{Real GDP, in billions of dollars}
#'   \item{K}{Capital stock, in billions of dollars}
#'   \item{K.price}{Annual cost of $1 billion of capital, using 10-year treasury}
#'   \item{L}{Labor force, in thousands of people}
#'   \item{L.price}{Annual wage for one thousand people}
#' }
#' @source \url{https://fred.stlouisfed.org/}
#' 
"USMacro"
