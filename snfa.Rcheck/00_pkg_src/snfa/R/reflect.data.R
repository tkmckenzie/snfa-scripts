#' Data reflection for kernel smoothing
#' 
#' This function reflects data below minimum and above maximum for use in reducing endpoint bias in kernel smoothing.
#' 
#' @param X Matrix of inputs
#' @param y Vector of outputs
#' 
#' @return Returns a list with the following elements
#' \item{X.reflected}{Reflected values of X}
#' \item{y.reflected}{Reflected values of y}
#' 
#' @examples
#' data(univariate)
#' 
#' #Extract data
#' X <- as.matrix(univariate$x)
#' y <- univariate$y
#' 
#' #Reflect data
#' reflected.data <- reflect.data(X, y)
#' 
#' X.reflected <- reflected.data$X
#' y.reflected <- reflected.data$y
#' 
#' #Plot
#' library(ggplot2)
#' 
#' plot.df <- data.frame(X = X.reflected,
#'                       y = y.reflected,
#'                       data = rep(c("reflected", "actual", "reflected"), each = nrow(X)))
#' 
#' ggplot(plot.df, aes(X, y)) +
#'   geom_point(aes(color = data))
#' 
#' @export

reflect.data <-
  function(X, y){
    #Rotates data below and above max values for boundary correction
    m <- cbind(X, y)
    k <- ncol(m)
    
    rotation.origin.lower <- apply(m, 2, min)
    rotation.origin.upper <- apply(m, 2, max)
    
    m.translated.lower <- sapply(1:ncol(m), function(i) m[,i] - rotation.origin.lower[i])
    m.translated.upper <- sapply(1:ncol(m), function(i) m[,i] - rotation.origin.upper[i])
    
    m.translated.rotated.lower <- -m.translated.lower
    m.translated.rotated.upper <- -m.translated.upper
    
    m.rotated.lower <- sapply(1:ncol(m), function(i) m.translated.rotated.lower[,i] + rotation.origin.lower[i])
    m.rotated.upper <- sapply(1:ncol(m), function(i) m.translated.rotated.upper[,i] + rotation.origin.upper[i])
    
    X.reflected <- rbind(m.rotated.lower[,-k,drop = FALSE], X, m.rotated.upper[,-k,drop = FALSE])
    y.reflected <- c(m.rotated.lower[,k], y, m.rotated.upper[,k])
    
    return(list(X.reflected = X.reflected, y.reflected = y.reflected))
  }
