#' Kernel smoothing with additional constraints
#' 
#' Fits conditional mean of data with kernel smoothing, imposing monotonicity and/or concavity constraints.
#' 
#' @param X.eval Matrix of inputs used for fitting
#' @param y.eval Vector of outputs used for fitting
#' @param X.constrained Matrix of inputs where constraints apply
#' @param X.fit Matrix of inputs where curve is fit; defaults to X.constrained
#' @param H.inv Inverse of the smoothing matrix (must be positive definite); defaults to rule of thumb
#' @param H.mult Scaling factor for rule of thumb smoothing matrix
#' @param method Constraints to apply; "u" for unconstrained, "m" for monotonically increasing, and "mc" for monotonically increasing and concave
#' @param scale.constraints Boolean, whether to scale constraints by their average value, can help with convergence
#' 
#' @return Returns a list with the following elements
#' \item{y.fit}{Estimated value of the frontier at X.fit}
#' \item{gradient.fit}{Estimated gradient of the frontier at X.fit}
#' \item{solution}{Boolean; TRUE if frontier successfully estimated}
#' \item{X.eval}{Matrix of inputs used for fitting}
#' \item{X.constrained}{Matrix of inputs where constraints apply}
#' \item{X.fit}{Matrix of inputs where curve is fit}
#' \item{H.inv}{Inverse smoothing matrix used in fitting}
#' \item{method}{Method used to fit frontier}
#' \item{scaling.factor}{Factor by which constraints are multiplied before quadratic programming}
#' 
#' @details
#' This method uses kernel smoothing to fit the mean of the data
#' while imposing specified monotonicity and concavity constraints. The procedure is
#' derived from Racine et al. (2009), which develops kernel smoothing methods with
#' bounding, monotonicity and concavity constraints. Specifically, the smoothing procedure
#' involves finding optimal weights for a Nadaraya-Watson estimator of the form 
#' 
#' \deqn{\hat{y} = m(x) = \sum_{i=1}^N p_i A(x, x_i) y_i,}
#' 
#' where \eqn{x} are inputs, \eqn{y} are outputs, \eqn{p} are weights, subscripts
#' index observations, and 
#' 
#' \deqn{A(x, x_i) = \frac{K(x, x_i)}{\sum_{h=1}^N K(x, x_h)}}
#' 
#' for a kernel \eqn{K}. This method uses a multivariate normal kernel of the form
#' 
#' \deqn{K(x, x_h) = \exp\left(-\frac12 (x - x_h)'H^{-1}(x - x_h)\right),}
#' 
#' where \eqn{H} is a bandwidth matrix. Bandwidth selection is performed via Silverman's
#' (1986) rule-of-thumb, in the function \code{H.inv.select}.
#' 
#' Optimal weights \eqn{\hat{p}} are selected by solving the quadratic programming problem
#' 
#' \deqn{\min_p \mbox{\ \ }-\mathbf{1}'p + \frac12 p'p.}
#' 
#' Monotonicity constraints of the following form can be imposed at 
#' specified points:
#' 
#' \deqn{\frac{\partial m(x)}{\partial x^j} = \sum_{h=1}^N p_h \frac{\partial A(x, x_h)}{\partial x^j} y_h \geq 0 \mbox{\ \ \ \ }\forall x, j,}
#' 
#' where superscripts index inputs. Finally concavity constraints of the following form can also be imposed using Afriat's
#' (1967) conditions:
#' 
#' \deqn{m(x) - m(z) \leq \nabla_x m(z) \cdot (x - z) \mbox{\ \ \ \ }\forall x, z.}
#' 
#' The gradient of the estimated curve at a point \eqn{x} is given by 
#' 
#' \deqn{\nabla_x m(x) = \sum_{i=1}^N \hat{p}_i \nabla_x A(x, x_i) y_i,}
#' 
#' where \eqn{\hat{p}_i} are estimated weights.
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
#' #Reflect data for fitting
#' reflected.data <- reflect.data(X, y)
#' X.eval <- reflected.data$X
#' y.eval <- reflected.data$y
#' 
#' #Fit frontier
#' fit.mc <- fit.mean(X.eval, y.eval, 
#'                    X.constrained = X,
#'                    X.fit = X,
#'                    method = "mc")
#'
#' #Plot input productivities over time
#' library(ggplot2)
#' plot.df <- data.frame(Year = rep(USMacro$Year, times = 2),
#'                       Elasticity = c(fit.mc$gradient.fit[,1] * X[,1] / y,
#'                                      fit.mc$gradient.fit[,2] * X[,2] / y),
#'                       Variable = rep(c("Capital", "Labor"), each = nrow(USMacro)))
#' 
#' ggplot(plot.df, aes(Year, Elasticity)) +
#'   geom_line() +
#'   facet_grid(Variable ~ ., scales = "free_y")
#'   
#' @references
#' \insertRef{ParmeterRacine}{snfa}
#' 
#' @export

fit.mean <-
  function(X.eval, y.eval, X.constrained = NA, X.fit = NA, H.inv = NA, H.mult = 1, method = "u", scale.constraints = TRUE){
    # require(quadprog)
    # require(abind)
    
    #X.eval, y.eval used for curve fitting
    #Constraints only applied at X.constrained, y.constrained
    #Curve only fit at X.fit
    
    #H.mult widens bandwidth, only used when H.inv not specified
    #H.inv is created using X.fit; if curve is being fit at points that are not data, user should specify H.inv
    
    stopifnot(method %in% c("u", "m", "mc")) #(u)nconstrained, (m)onotonically increasing, (c)oncave
    
    if (any(is.na(X.constrained))){
      X.constrained <- X.eval
    }
    
    if (any(is.na(X.fit))){
      X.fit <- X.constrained
    }
    
    N.eval <- nrow(X.eval)
    N.constrained <- nrow(X.constrained)
    N.fit <- nrow(X.fit)
    if ((ncol(X.eval) == ncol(X.constrained)) & (ncol(X.constrained) == ncol(X.fit))){
      k <- ncol(X.eval)
    } else{
      stop("X.eval, X.constrained, and X.fit must have same number of columns.")
    }
    
    #Rule of thumb for multivariate:
    if (any(is.na(H.inv))){
      H.inv <- H.inv.select(X.constrained, H.mult = H.mult)
    }
    
    #K[i,j] = K(X[i,], X.fit[j,])
    K <- apply(X.fit, 1, function(X.0) apply(X.eval, 1, function(X.i) kernel.func(X.0, X.i, H.inv)))
    
    #A[i,j] = A(X[i,], X[j,])
    A <- t(K) / colSums(K)
    
    A.y <- t(t(A) * y.eval)
    
    #K.grad[i,j,k] = dK(X[i,], X[j,]) / dX_k
    K.grad <- aperm(abind::abind(lapply(1:N.eval, function(index.0) kernel.grad.mult(index.0, H.inv, X.eval, X.fit, K)), along = 3), c(3, 2, 1))
    
    #A.grad[i,j,k] = dA(X[i,], X[j,]) / dX_k
    K.sum <- colSums(K)
    K.grad.sum <- apply(K.grad, c(2, 3), sum)
    A.grad <- aperm(abind::abind(lapply(1:N.eval, function(j) t((K.grad[j,,] * K.sum - K.grad.sum * K[j,]) / (K.sum)^2)), along = 3), c(2, 3, 1))
    
    if (method == "u"){
      scaling.factor <- 1
      p.hat <- rep(1, N.eval)
    } else{
      # monotonicity.constraints <- t(t(apply(A.grad, 1:2, sum)) * y)
      monotonicity.constraints <- Reduce(rbind, lapply(1:k, function(j) t(t(A.grad[,,j]) * y.eval)))
      
      if (scale.constraints){
        log.scaling.factor <- floor(-mean(log(abs(monotonicity.constraints), base = 10)))
        if (log.scaling.factor < 0){
          log.scaling.factor <- 0
        }
        scaling.factor <- 10^log.scaling.factor
      } else{
        scaling.factor <- 1
      }
      
      if (method == "m"){
        constraints <- monotonicity.constraints
        
        p.hat <- try(quadprog::solve.QP(Dmat = diag(N.eval),
                                        dvec = rep(1, N.eval),
                                        Amat = scaling.factor * t(constraints),
                                        bvec = scaling.factor * c(rep(0, N.fit * k)),
                                        meq = 0)$solution, 
                     silent = TRUE)
      } else{
        concavity.constraints <- t(Reduce(cbind, lapply(1:N.fit, function(i) concavity.matrix(i, X.fit, y.eval, A, A.grad))))
        trivial.constraints <- seq(1, N.fit^2, by = N.fit + 1)
        if (length(trivial.constraints) > 0) concavity.constraints = concavity.constraints[-trivial.constraints,] #Remove 0 * p >= 0 constraints
        
        constraints <- rbind(monotonicity.constraints, concavity.constraints)
        
        if (scale.constraints){
          log.scaling.factor <- floor(-mean(log(abs(constraints), base = 10)))
          if (log.scaling.factor < 0){
            log.scaling.factor <- 0
          }
          scaling.factor <- 10^log.scaling.factor
        } else{
          scaling.factor <- 1
        }
        
        p.hat <- try(quadprog::solve.QP(Dmat = diag(N.eval),
                                        dvec = rep(1, N.eval),
                                        Amat = scaling.factor * t(constraints),
                                        bvec = scaling.factor * c(rep(0, N.fit * (N.fit + k - 1))),
                                        meq = 0)$solution, 
                     silent = TRUE)
      }
    }
    
    if (mode(p.hat) == "numeric"){
      # print(mode(p.hat))
      y.fit <- A.y %*% p.hat
      
      gradient.fit <- sapply(1:k, function(d) t(t(A.grad[,,d]) * y.eval) %*% p.hat)
      
      return(list(y.fit = y.fit, #Estimates of y at X.fit
                  gradient.fit = gradient.fit, #Estimates of gradient at X.fit
                  # efficiency = y / y.hat[1:N],
                  solution = TRUE, #Indicator for convergence
                  X.eval = X.eval, X.constrained = X.constrained, X.fit = X.fit,
                  H.inv = H.inv, method = method, scaling.factor = scaling.factor))
    } else{
      return(list(y.fit = NA, #Estimates of y at X.fit
                  gradient.fit = NA, #Estimates of gradient at X.fit
                  # efficiency = NA,
                  solution = FALSE,
                  X.eval = X.eval, X.constrained = X.constrained, X.fit = X.fit,
                  H.inv = H.inv, method = method, scaling.factor = scaling.factor))
    }
  }
