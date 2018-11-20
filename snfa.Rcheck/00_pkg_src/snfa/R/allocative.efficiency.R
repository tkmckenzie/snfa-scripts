#' Allocative efficiency estimation
#' 
#' Fits frontier to data and estimates technical and allocative efficiency
#' 
#' @param X Matrix of inputs
#' @param y Vector of outputs
#' @param X.price Matrix of input prices
#' @param y.price Vector of output prices
#' @param X.constrained Matrix of inputs where constraints apply
#' @param H.inv Inverse of the smoothing matrix (must be positive definite); defaults to rule of thumb
#' @param H.mult Scaling factor for rule of thumb smoothing matrix
#' @param model Type of frontier to use; "br" for boundary regression, "sf" for stochastic frontier
#' @param method Constraints to apply; "u" for unconstrained, "m" for monotonically increasing, and "mc" for monotonically increasing and concave
#' @param scale.constraints Boolean, whether to scale constraints by their average value, can help with convergence
#' 
#' @return Returns a list with the following elements
#' \item{y.fit}{Estimated value of the frontier at X.fit}
#' \item{gradient.fit}{Estimated gradient of the frontier at X.fit}
#' \item{technical.efficiency}{Estimated technical efficiency}
#' \item{log.overallocation}{Estimated log-overallocation of each input for each observation}
#' \item{X.eval}{Matrix of inputs used for fitting}
#' \item{X.constrained}{Matrix of inputs where constraints apply}
#' \item{H.inv}{Inverse smoothing matrix used in fitting}
#' \item{method}{Method used to fit frontier}
#' \item{scaling.factor}{Factor by which constraints are multiplied before quadratic programming}
#' 
#' @details 
#' This function estimates allocative inefficiency using the methodology in McKenzie
#' (2018). The estimation process is a non-parametric analogue of Schmidt and Lovell
#' (1979). First, the frontier is fit using either a boundary regression or stochastic
#' frontier as in Racine et al. (2009), from which technical efficiency is estimated.
#' Then, gradients and price ratios are computed for each observation and compared to 
#' determine the extent of misallocation. Specifically, log-overallocation is computed as
#' 
#' \deqn{\log\left(\frac{w_i^j}{p_i}\right) - \log\left(\phi_i\frac{\partial f(x_i)}{\partial x^j}\right),}
#' 
#' where \eqn{\phi_i} is the efficiency of observation \eqn{i},
#' \eqn{\partial f(x_i) / \partial x^j} is the marginal productivity of input \eqn{j}
#' at observation \eqn{i}, \eqn{w_i^j} is the cost of input \eqn{j} for observation
#' \eqn{i}, and \eqn{p_i} is the price of output for observation \eqn{i}.
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
#' X.price <- as.matrix(USMacro[,c("K.price", "L.price")])
#' y.price <- rep(1e9, nrow(USMacro)) #Price of $1 billion of output is $1 billion
#' 
#' #Run model
#' efficiency.model <- allocative.efficiency(X, y,
#'                                          X.price, y.price,
#'                                          X.constrained = X,
#'                                          model = "br",
#'                                          method = "mc")
#' 
#' #Plot technical/allocative efficiency over time
#' library(ggplot2)
#' 
#' technical.df <- data.frame(Year = USMacro$Year,
#'                           Efficiency = efficiency.model$technical.efficiency)
#' 
#' ggplot(technical.df, aes(Year, Efficiency)) +
#'   geom_line()
#'
#' allocative.df <- data.frame(Year = rep(USMacro$Year, times = 2),
#'                             log.overallocation = c(efficiency.model$log.overallocation[,1],
#'                                                    efficiency.model$log.overallocation[,2]),
#'                             Variable = rep(c("K", "L"), each = nrow(USMacro)))
#' 
#' ggplot(allocative.df, aes(Year, log.overallocation)) +
#'   geom_line(aes(color = Variable))
#'
#' #Estimate average overallocation across sample period
#' lm.model <- lm(log.overallocation ~ 0 + Variable, allocative.df)
#' summary(lm.model)
#'   
#' @references
#' \insertRef{AignerLovellSchmidt}{snfa}\cr\cr
#' \insertRef{McKenzie}{snfa}\cr\cr
#' \insertRef{ParmeterRacine}{snfa}\cr\cr
#' \insertRef{SchmidtLovell}{snfa}
#' 
#' @export

allocative.efficiency <-
  function(X, y, X.price, y.price, X.constrained = NA, H.inv = NA, H.mult = 1, model = "br", method = "u", scale.constraints = TRUE){
    # require(quadprog)
    # require(abind)
    
    #X.eval, y.eval used for curve fitting
    #Constraints only applied at X.constrained, y.constrained
    #Curve only fit at X.fit
    
    #H.mult widens bandwidth, only used when H.inv not specified
    
    if (!(model %in% c("br", "sf"))) stop("model must be \"br\" or \"sf\".")
    
    if ((nrow(X) != length(y)) | (nrow(X) != nrow(X.price)) | (nrow(X) != length(y.price))) stop("X, y, X.price, and y.price must have same number of observations.")
    if (ncol(X) != ncol(X.price)) stop("X must have same number of columns as X.price.")
    
    N <- nrow(X)
    
    #Reflect data for fitting
    reflected.data <- reflect.data(X, y)
    
    X.eval <- reflected.data$X.reflected
    y.eval <- reflected.data$y.reflected
    
    if (any(is.na(X.constrained))){
      X.constrained <- X
    }
    
    if (any(is.na(H.inv))){
      H.inv <- H.inv.select(X, H.mult = H.mult)
    }
    
    #Fit frontier
    if (model == "br"){
      m <- fit.boundary(X.eval, y.eval, X, y, X.constrained, X, y, H.inv, H.mult, method, scale.constraints)
      
      technical.efficiency <- m$efficiency
      gradient <- m$gradient.fit
    } else{
      m <- fit.sf(X, y, X.constrained, H.inv, H.mult, method, scale.constraints)
      
      technical.efficiency <- m$mode.efficiency
      gradient <- m$gradient.fit
    }
    
    #Compute price ratios
    price.ratio <- apply(X.price, 2, function(col) col / y.price)
    
    #Compute mraginal productivities
    marginal.productivities <- apply(gradient, 2, function(col) col * technical.efficiency)
    
    log.overallocation <- log(price.ratio) - log(marginal.productivities)
    
    
    return(list(y.fit = m$y.fit, #Estimates of y at X.fit
                gradient.fit = gradient, #Estimates of gradient at X.fit
                technical.efficiency = technical.efficiency,
                log.overallocation = log.overallocation,
                X.eval = X.eval, X.constrained = X.constrained,
                H.inv = H.inv, method = method, scaling.factor = m$scaling.factor))
  }
