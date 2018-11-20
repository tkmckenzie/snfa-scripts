kernel.grad.mult <-
  function(index.0, H.inv, X, X.fit, K){
    return(t(K[index.0,] * t(-H.inv %*% (t(X.fit) - X[index.0,]))))
  }
