LogLikSizeShape = function(X, mu, sigma, p, M, A) {
     ## Computes the size and shape likelihood of a matched matrix X
     ## from a mean mu inside a volume of a cube A.
     m = dim(X)[1]
     Q = (p-1) * m
     RDs = X[, 1:(p-1)] %*% t(mu[, 1:(p-1)]) / sigma^2
     ss = SvdNew(RDs)
     ss$d = c(1, 1, sign(det(ss$v %*% t(ss$u)))) * ss$d
     lik = RotSimLaplaceC(ss$d)[1]
     
     res = (p - M) * log(A) - (Q / 2) * log(2 * pi * sigma^2) - 
          (Norm(X[, 1:(p-1)])^2 + Norm(mu[, 1:(p-1)])^2) / (2 * sigma^2) + 
          lik
     return(res)
}


