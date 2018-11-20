concavity.matrix <-
  function(index.X, X.fit, y, A, A.grad){
    return((A.grad[index.X,,] %*% (t(X.fit) - X.fit[index.X,]) - t(A) + A[index.X,]) * y)
  }
