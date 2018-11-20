#' Bandwidth matrix selection
#' 
#' Computes inverse of bandwidth matrix using rule-of-thumb from Silverman (1986).
#' 
#' @param X Matrix of inputs
#' @param H.mult Scaling factor for rule-of-thumb smoothing matrix
#' 
#' @return Returns inverse bandwidth matrix
#' 
#' @examples 
#' data(USMacro)
#' 
#' USMacro <- USMacro[complete.cases(USMacro),]
#' 
#' #Extract data
#' X <- as.matrix(USMacro[,c("K", "L")])
#' 
#' #Generate bandwidth matrix
#' print(H.inv.select(X))
#' #              [,1]         [,2]
#' # [1,] 3.642704e-08 0.000000e+00
#' # [2,] 0.000000e+00 1.215789e-08
#' 
#' @details 
#' This method performs selection of (inverse) multivariate bandwidth matrices using
#' Silverman's (1986) rule-of-thumb. Specifically, Silverman recommends setting the bandwidth
#' matrix to
#' 
#' \deqn{H_{jj}^{1/2} = \left(\frac{4}{M + 2}\right)^{1 / (M + 4)}\times N^{-1 / (M + 4)}\times \mbox{sd}(x^j) \mbox{\ \ \ \ for }j=1,...,M}
#' \deqn{H_{ab} = 0\mbox{\ \ \ \ for }a\neq b}
#' 
#' where \eqn{M} is the number of inputs, \eqn{N} is the number of observations, and
#' \eqn{\mbox{sd}(x^j)} is the sample standard deviation of input \eqn{j}.
#' 
#' @references
#' \insertRef{Silverman}{snfa}
#' 
#' @export

H.inv.select <-
  function(X, H.mult = 1){
    N <- nrow(X)
    k <- ncol(X)
    
    col.mult <- (4 / (k + 2))^(1 / (k + 4)) * N^(-1 / (k + 4)) #Silverman's
    # col.mult <- N.bounded^(-1 / (k + 4)) #Scott's
    H.diag <- apply(X, 2, function(col) H.mult * (col.mult * stats::sd(col))^2)
    if (ncol(X) > 1){
      H <- diag(H.diag)
    } else{
      H <- diag(1) * H.diag
    }
    H.inv <- solve(H)
    
    return(H.inv)
  }
