kernel.func <-
  function(X.0, X.i, H.inv){
    d <- X.0 - X.i
    return(c(exp(-0.5 * (t(d) %*% H.inv %*% d))))
  }
