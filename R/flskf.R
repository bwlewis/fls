`flskf` <-
function (A, b, mu=1, ncap=length(b))
{
  m <- nrow (A)
  n <- ncol (A)
  X <- array (0,c(n,ncap))
  R <- array (0,c(n,n))
  diag (R) <- 1/mu
  a <- t(A[1,,drop=FALSE])
  phi <- c((1/mu) * crossprod(a) + 1)
  w <- (1/phi) * R %*% a
  rho <- b[1]
  X[,1] <- rho * w
  for (j in 2:ncap) {
    R = R - phi * tcrossprod(w)
    diag (R) <- diag (R) + 1/mu
    a <- t(A[j,,drop=FALSE])
    phi <- c(t(a) %*% R %*% a)
    phi <- phi + 1
    rho <- c(b[j] - crossprod (a, X[,j-1]))
    w <- (1/phi) * R %*% a
    X[,j] <- X[,j-1] + rho*w
  }
  X
}

