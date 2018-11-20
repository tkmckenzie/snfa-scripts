#' Multivariate smooth boundary fitting with additional constraints
#' 
#' Fits boundary of data with kernel smoothing, imposing monotonicity and/or concavity constraints.
#' 
#' @param X.eval Matrix of inputs used for fitting
#' @param y.eval Vector of outputs used for fitting
#' @param X.bounded Matrix of inputs where bounding constraints apply
#' @param y.bounded Vector of outputs where bounding constraints apply
#' @param X.constrained Matrix of inputs where monotonicity/concavity constraints apply
#' @param X.fit Matrix of inputs where curve is fit; defaults to X.constrained
#' @param y.fit.observed Vector of outputs corresponding to observations in X.fit; used for efficiency calculation
#' @param H.inv Inverse of the smoothing matrix (must be positive definite); defaults to rule of thumb
#' @param H.mult Scaling factor for rule of thumb smoothing matrix
#' @param method Constraints to apply; "u" for unconstrained, "m" for monotonically increasing, and "mc" for monotonically increasing and concave
#' @param scale.constraints Boolean, whether to scale constraints by their average value, can help with convergence
#' 
#' @return Returns a list with the following elements
#' \item{y.fit}{Estimated value of the frontier at X.fit}
#' \item{gradient.fit}{Estimated gradient of the frontier at X.fit}
#' \item{efficiency}{Estimated efficiencies of y.fit.observed}
#' \item{solution}{Boolean; TRUE if frontier successfully estimated}
#' \item{X.eval}{Matrix of inputs used for fitting}
#' \item{X.constrained}{Matrix of inputs where monotonicity/concavity constraints apply}
#' \item{X.fit}{Matrix of inputs where curve is fit}
#' \item{H.inv}{Inverse smoothing matrix used in fitting}
#' \item{method}{Method used to fit frontier}
#' \item{scaling.factor}{Factor by which constraints are multiplied before quadratic programming}
#' 
#' @details
#' This method fits a smooth boundary of the data (with all data points below the boundary)
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
#' This method always imposes bounding constraints as specified points, given by
#' 
#' \deqn{m(x_i) - y_i = \sum_{h=1}^N p_h A(x_i, x_h) y_h - y_i \geq 0 \mbox{\ \ \ \ }\forall i.}
#' 
#' Additionally, monotonicity constraints of the following form can be imposed at 
#' specified points:
#' 
#' \deqn{\frac{\partial m(x)}{\partial x^j} = \sum_{h=1}^N p_h \frac{\partial A(x, x_h)}{\partial x^j} y_h \geq 0 \mbox{\ \ \ \ }\forall x, j,}
#' 
#' where superscripts index inputs. Finally concavity constraints of the following form can also be imposed using Afriat's
#' (1967) conditions:
#' 
#' \deqn{m(x) - m(z) \leq \nabla_x m(z) \cdot (x - z) \mbox{\ \ \ \ }\forall x, z.}
#' 
#' The gradient of the frontier at a point \eqn{x} is given by 
#' 
#' \deqn{\nabla_x m(x) = \sum_{i=1}^N \hat{p}_i \nabla_x A(x, x_i) y_i,}
#' 
#' where \eqn{\hat{p}_i} are estimated weights.
#' 
#' @examples
#' data(univariate)
#' 
#' #Set up data for fitting
#' 
#' X <- as.matrix(univariate$x)
#' y <- univariate$y
#' 
#' N.fit <- 100
#' X.fit <- as.matrix(seq(min(X), max(X), length.out = N.fit))
#' 
#' #Reflect data for fitting
#' reflected.data <- reflect.data(X, y)
#' X.eval <- reflected.data$X
#' y.eval <- reflected.data$y
#' 
#' #Fit frontiers
#' frontier.u <- fit.boundary(X.eval, y.eval, 
#'                            X.bounded = X, y.bounded = y,
#'                            X.constrained = X.fit,
#'                            X.fit = X.fit,
#'                            method = "u")
#'                           
#' frontier.m <- fit.boundary(X.eval, y.eval, 
#'                            X.bounded = X, y.bounded = y,
#'                            X.constrained = X.fit,
#'                            X.fit = X.fit,
#'                            method = "m")
#'                           
#' frontier.mc <- fit.boundary(X.eval, y.eval, 
#'                             X.bounded = X, y.bounded = y,
#'                             X.constrained = X.fit,
#'                             X.fit = X.fit,
#'                             method = "mc")
#'
#' #Plot frontier
#' library(ggplot2)
#' 
#' frontier.df <- data.frame(X = rep(X.fit, times = 3),
#'                           y = c(frontier.u$y.fit, frontier.m$y.fit, frontier.mc$y.fit),
#'                           model = rep(c("u", "m", "mc"), each = N.fit))
#' 
#' ggplot(univariate, aes(X, y)) +
#'   geom_point() +
#'   geom_line(data = frontier.df, aes(color = model))
#'
#' #Plot slopes
#' slope.df <- data.frame(X = rep(X.fit, times = 3),
#'                        slope = c(frontier.u$gradient.fit,
#'                                  frontier.m$gradient.fit,
#'                                  frontier.mc$gradient.fit),
#'                        model = rep(c("u", "m", "mc"), each = N.fit))
#'
#' ggplot(slope.df, aes(X, slope)) +
#'   geom_line(aes(color = model))
#'   
#' @references
#' \insertRef{ParmeterRacine}{snfa}
#' 
#' @export

fit.boundary <-
  function(X.eval, y.eval, X.bounded, y.bounded, X.constrained = NA, X.fit = NA, y.fit.observed = NA, H.inv = NA, H.mult = 1, method = "u", scale.constraints = TRUE){
    # require(quadprog)
    # require(abind)
    # require(prodlim)
    
    #X.eval, y.eval used for curve fitting
    #Bounding constraints only applied at X.bounded, y.bounded
    #Monotonicity/concavity constraints applied at X.constrained and X.fit (no y.constrained needed)
    #Curve only fit at X.fit
    
    #H.mult widens bandwidth, only used when H.inv not specified
    
    stopifnot(method %in% c("u", "m", "mc")) #(u)nconstrained, (m)onotonically increasing, (c)oncave
    
    if (any(is.na(X.constrained))){
      X.constrained <- X.bounded
    }
    
    if (any(is.na(X.fit))){
      X.fit <- X.constrained
    }
    
    if (!any(is.na(y.fit.observed)) & (length(y.fit.observed) != nrow(X.fit))){
      stop("y.fit.observed must have same number of observations as X.fit.")
    }
    
    N.eval <- nrow(X.eval)
    N.bounded <- nrow(X.bounded)
    N.constrained <- nrow(X.constrained)
    N.fit <- nrow(X.fit)
    if ((ncol(X.eval) == ncol(X.constrained)) & (ncol(X.constrained) == ncol(X.fit)) & (ncol(X.fit) == ncol(X.bounded))){
      k <- ncol(X.eval)
    } else{
      stop("X.eval, X.constrained, and X.fit must have same number of columns.")
    }
    
    #Combine X.fit and X.constrained
    #First eliminate observations in X.constrained that appear in X.fit
    duplicated.rows <- prodlim::row.match(as.data.frame(X.fit), as.data.frame(X.constrained))
    duplicated.rows <- duplicated.rows[!is.na(duplicated.rows)]
    
    if (length(duplicated.rows) == 0){
      X.fit <- rbind(X.fit, X.constrained)
    } else{
      X.fit <- rbind(X.fit, X.constrained[-duplicated.rows,,drop = FALSE])
    }
    
    N.fit.new <- nrow(X.fit)
    
    #Rule of thumb for multivariate:
    #Fit to X.constrained since that is typically raw data
    if (any(is.na(H.inv))){
      H.inv <- H.inv.select(X.constrained, H.mult = H.mult)
    }
    
    #Bounding constraints:
    # K <- apply(X.eval, 1, function(X.0) apply(X.constrained, 1, function(X.i) kernel.func(X.0, X.i, H.inv)))
    K <- apply(X.eval, 1, function(X.0) apply(X.bounded, 1, function(X.i) kernel.func(X.0, X.i, H.inv)))
    A <- K / rowSums(K)
    bounding.constraints <- t(t(A) * y.eval)
    
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
      if (scale.constraints){
        log.scaling.factor <- floor(-mean(log(abs(bounding.constraints), base = 10)))
        if (log.scaling.factor < 0){
          log.scaling.factor <- 0
        }
        scaling.factor <- 10^log.scaling.factor
      } else{
        scaling.factor <- 1
      }
      
      p.hat <- try(quadprog::solve.QP(Dmat = diag(N.eval),
                                      dvec = rep(1, N.eval),
                                      Amat = scaling.factor * t(bounding.constraints),
                                      bvec = scaling.factor * y.bounded,
                                      meq = 0)$solution,
                   silent = TRUE)
    } else{
      # monotonicity.constraints <- t(t(apply(A.grad, 1:2, sum)) * y)
      monotonicity.constraints <- Reduce(rbind, lapply(1:k, function(j) t(t(A.grad[,,j]) * y.eval)))
      
      if (method == "m"){
        constraints <- rbind(bounding.constraints, monotonicity.constraints)
        
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
                                        bvec = scaling.factor * c(y.bounded, rep(0, (N.fit.new) * k)),
                                        meq = 0)$solution, 
                     silent = TRUE)
      } else{
        concavity.constraints <- t(Reduce(cbind, lapply(1:(N.fit.new), function(i) concavity.matrix(i, X.fit, y.eval, A, A.grad))))
        trivial.constraints <- seq(1, (N.fit.new)^2, by = (N.fit.new) + 1)
        if (length(trivial.constraints) > 0) concavity.constraints = concavity.constraints[-trivial.constraints,] #Remove 0 * p >= 0 constraints
        
        constraints <- rbind(bounding.constraints, monotonicity.constraints, concavity.constraints)
        
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
                                        bvec = scaling.factor * c(y.bounded, rep(0, (N.fit.new) * (N.fit.new + k - 1))),
                                        meq = 0)$solution,
                     silent = TRUE)
      }
    }
    
    if (mode(p.hat) == "numeric"){
      # print(mode(p.hat))
      y.fit <- A.y %*% p.hat
      
      gradient.fit <- sapply(1:k, function(d) t(t(A.grad[,,d]) * y.eval) %*% p.hat)
      
      return(list(y.fit = y.fit[1:N.fit], #Estimates of y at X.fit
                  gradient.fit = gradient.fit[1:N.fit,], #Estimates of gradient at X.fit
                  efficiency = y.fit.observed / y.fit[1:N.fit],
                  solution = TRUE, #Indicator for convergence
                  X.eval = X.eval, X.constrained = X.constrained, X.fit = X.fit,
                  H.inv = H.inv, method = method, scaling.factor = scaling.factor))
    } else{
      return(list(y.fit = NA, #Estimates of y at X.fit
                  gradient.fit = NA, #Estimates of gradient at X.fit
                  efficiency = NA,
                  solution = FALSE,
                  X.eval = X.eval, X.constrained = X.constrained, X.fit = X.fit,
                  H.inv = H.inv, method = method, scaling.factor = scaling.factor))
    }
  }
