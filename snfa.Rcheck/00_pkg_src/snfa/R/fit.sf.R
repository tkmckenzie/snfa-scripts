#' Non-parametric stochastic frontier
#' 
#' Fits stochastic frontier of data with kernel smoothing, imposing monotonicity and/or concavity constraints.
#' 
#' @param X Matrix of inputs
#' @param y Vector of outputs
#' @param X.constrained Matrix of inputs where constraints apply
#' @param H.inv Inverse of the smoothing matrix (must be positive definite); defaults to rule of thumb
#' @param H.mult Scaling factor for rule of thumb smoothing matrix
#' @param method Constraints to apply; "u" for unconstrained, "m" for monotonically increasing, and "mc" for monotonically increasing and concave
#' @param scale.constraints Boolean, whether to scale constraints by their average value, can help with convergence
#' 
#' @return Returns a list with the following elements
#' \item{y.fit}{Estimated value of the frontier at X.fit}
#' \item{gradient.fit}{Estimated gradient of the frontier at X.fit}
#' \item{mean.efficiency}{Average efficiency for X, y as a whole}
#' \item{mode.efficiency}{Modal efficiencies for each observation in X, y}
#' \item{X.eval}{Matrix of inputs used for fitting}
#' \item{X.constrained}{Matrix of inputs where constraints apply}
#' \item{X.fit}{Matrix of inputs where curve is fit}
#' \item{H.inv}{Inverse smoothing matrix used in fitting}
#' \item{method}{Method used to fit frontier}
#' \item{scaling.factor}{Factor by which constraints are multiplied before quadratic programming}
#' 
#' @details
#' This method fits non-parametric stochastic frontier models. The data-generating process
#' is assumed to be of the form
#' 
#' \deqn{\ln y_i = \ln f(x_i) + v_i - u_i,}
#' 
#' where \eqn{y_i} is the \eqn{i}th observation of output, \eqn{f} is a continuous
#' function, \eqn{x_i} is the \eqn{i}th observation of input, \eqn{v_i} is a
#' normally-distributed error term (\eqn{v_i\sim N(0, \sigma_v^2)}), and \eqn{u_i} is a
#' normally-distributed error term truncated below at zero (\eqn{u_i\sim N^+(0, \sigma_u)}). 
#' Aigner et al. developed methods to decompose
#' \eqn{\varepsilon_i = v_i - u_i} into its basic components.
#' 
#' This procedure first fits the mean of the data using \code{fit.mean},
#' producing estimates of output \eqn{\hat{y}}. Log-proportional errors are calculated as
#' 
#' \deqn{\varepsilon_i = \ln(y_i / \hat{y}_i).}
#' 
#' Following Aigner et al. (1977), parameters of one- and two-sided error distributions
#' are estimated via maximum likelihood. First,
#' 
#' \deqn{\hat{\sigma}^2 = \frac1N \sum_{i=1}^N \varepsilon_i^2.}
#' 
#' Then, \eqn{\hat{\lambda}} is estimated by solving
#' 
#' \deqn{\frac1{\hat{\sigma}^2} \sum_{i=1}^N \varepsilon_i\hat{y}_i + \frac{\hat{\lambda}}{\hat{\sigma}} \sum_{i=1}^N \frac{f_i^*}{1 - F_i^*}y_i = 0,}
#' 
#' where \eqn{f_i^*} and \eqn{F_i^*} are standard normal density and distribution function,
#' respectively, evaluated at \eqn{\varepsilon_i\hat{\lambda}\hat{\sigma}^{-1}}. Parameters of
#' the one- and two-sided distributions are found by solving the identities
#' 
#' \deqn{\sigma^2 = \sigma_u^2 + \sigma_v^2}
#' \deqn{\lambda = \frac{\sigma_u}{\sigma_v}.}
#' 
#' Mean efficiency over the sample is given by 
#' 
#' \deqn{\exp\left(-\frac{\sqrt{2}}{\sqrt{\pi}}\right)\sigma_u,}
#' 
#' and modal efficiency for each observation is given by
#' 
#' \deqn{-\varepsilon(\sigma_u^2/\sigma^2).}
#' 
#' @examples
#' data(USMacro)
#' 
#' USMacro <- USMacro[complete.cases(USMacro),]
#' 
#' #Extract data
#' X <- as.matrix(USMacro[,c("K", "L")])
#' y <- USMacro$Y
#' 
#' #Fit frontier
#' fit.sf <- fit.sf(X, y,
#'                  X.constrained = X,
#'                  method = "mc")
#'
#' print(fit.sf$mean.efficiency)
#' # [1] 0.9772484
#' 
#' #Plot efficiency over time
#' library(ggplot2)
#' 
#' plot.df <- data.frame(Year = USMacro$Year,
#'                       Efficiency = fit.sf$mode.efficiency)
#' 
#' ggplot(plot.df, aes(Year, Efficiency)) +
#'   geom_line()
#'   
#' @references
#' \insertRef{AignerLovellSchmidt}{snfa}\cr\cr
#' \insertRef{ParmeterRacine}{snfa}
#' 
#' @export

fit.sf <-
  function(X, y, X.constrained = NA, H.inv = NA, H.mult = 1, method = "u", scale.constraints = TRUE){
    # require(quadprog)
    # require(abind)
    
    #X.eval, y.eval used for curve fitting
    #Constraints only applied at X.constrained, y.constrained
    #Curve only fit at X.fit
    
    #H.mult widens bandwidth, only used when H.inv not specified
    
    #Strategy: Use fit.mean to fit conditional mean, take difference from observed data, fit to one-sided error model
    
    if (nrow(X) != length(y)) stop("X must have same number of observations as y.")
    
    N <- nrow(X)
    
    #Reflect data for fitting
    reflected.data <- reflect.data(X, y)
    
    X.eval <- reflected.data$X.reflected
    y.eval <- reflected.data$y.reflected
    
    X.fit <- X
    
    if (any(is.na(X.constrained))){
      X.constrained <- X.fit
    }
    
    if (any(is.na(H.inv))){
      H.inv <- H.inv.select(X, H.mult = H.mult)
    }
    
    mean.estimate <- fit.mean(X.eval, y.eval, X.constrained, X.fit, H.inv = H.inv, method = method, scale.constraints = scale.constraints)
    if (!mean.estimate$solution) stop("Curve fitting failed.")
    
    errors <- log(y / mean.estimate$y.fit)
    
    #Estimate sigma.sq and lambda from system of equations in Aigner, Lovell, Schmidt (1977)
    sigma.sq.est <- sum(errors^2) / N
    sigma.est <- sqrt(sigma.sq.est)
    
    lambda.obj <- function(lambda.log){
      lambda <- exp(lambda.log)
      eval.point <- errors * lambda / sigma.est
      normal.eval <- exp(stats::dnorm(eval.point, log = TRUE) - stats::pnorm(eval.point, log = TRUE, lower.tail = FALSE))
      return(sum(errors * mean.estimate$y.fit) / sigma.sq.est + sum(normal.eval * y) * lambda / sigma.est)
    }
    lambda.est <- exp(rootSolve::multiroot(lambda.obj, start = 0)$root)
    
    #Back out sd.u and sd.v
    var.v <- sigma.sq.est / (1 + lambda.est^2)
    sd.v <- sqrt(var.v)
    sd.u <- lambda.est * sd.v
    
    mean.efficiency <- exp(-sd.u * sqrt(2) / sqrt(pi))
    mode.efficiency <- ifelse(errors < 0, 1, exp(-errors * (sd.u^2 / (sd.u^2 + sd.v^2))))
    
    return(list(y.fit = mean.estimate$y.fit, #Estimates of y at X.fit
                gradient.fit = mean.estimate$gradient.fit, #Estimates of gradient at X.fit
                mean.efficiency = mean.efficiency,
                mode.efficiency = mode.efficiency,
                X.eval = X.eval, X.constrained = X.constrained, X.fit = X.fit,
                H.inv = H.inv, method = method, scaling.factor = mean.estimate$scaling.factor))
  }
