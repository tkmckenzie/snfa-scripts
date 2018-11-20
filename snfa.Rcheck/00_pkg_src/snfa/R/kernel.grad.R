kernel.grad <-
  function(X.0, X.i, H.inv){
    return(kernel.func(X.0, X.i, H.inv) * (-H.inv %*% (X.0 - X.i)))
  }
