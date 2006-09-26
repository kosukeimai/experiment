###
### Calculate the CACE with covariates and optional clustering using
### two-stage least squares
###

CACEcov <- function(Y, D, Z, X, grp = NULL, data = parent.frame(),
                    robust = FALSE){ 

  call <- match.call()
  Y <- matrix(eval(call$Y, data), ncol = 1)
  N <- nrow(Y)
  D <- matrix(eval(call$D, data), ncol = 1)
  Z <- matrix(eval(call$Z, data), nrow = N)
  X <- cbind(model.matrix(X, data = data), D)
  if (!is.null(grp)) {
    sgrp <- sort(grp, index.return = TRUE)
    grp <- grp[sgrp$ix]
    X <- X[sgrp$ix,]
    Z <- Z[sgrp$ix,]
    D <- D[sgrp$ix,]
    Y <- Y[sgrp$ix,]
  }
  
  Pz <- Z %*% solve(t(Z) %*% Z) %*% t(Z) 
  beta <- solve(t(X) %*% Pz %*% X) %*% t(X) %*% Pz %*% Y
  epsilon <- Y - X %*% beta
  est <- beta[length(beta)]
  
  if (is.na(grp)) {
    if (robust) { 
      tmp <- solve(t(X) %*% Pz %*% X)
      var <- tmp %*% (t(X) %*% Z %*% solve(t(Z) %*% Z)
                      %*% (t(Z) %*% diag(epsilon^2) %*% Z)
                      %*% solve(t(Z) %*% Z) %*% t(Z) %*% X) %*% tmp
    } else {
      sig2 <- t(epsilon) %*% epsilon / N
      var <- sig2 * solve(t(X) %*% Pz %*% X)
    }
  } else {
      tmp <- solve(t(X) %*% Pz %*% X)
      Omega <- matrix(0, ncol = N, nrow = N)
      n.grp <- length(unique(grp))
      counter <- 1
      for (i in 1:n.grp) {
        n.grp.obs <- sum(grp == unique(grp)[i])
        Omega[counter:(counter+n.grp.obs-1),counter:(counter+n.grp.obs-1)] <-
          epsilon[grp == unique(grp)[i]] %*% t(epsilon[grp == unique(grp)[i]])
        counter <- counter + n.grp.obs
      }
      var <- tmp %*% (t(X) %*% Z %*% solve(t(Z) %*% Z)
                      %*% (t(Z) %*% Omega %*% Z)
                      %*% solve(t(Z) %*% Z) %*% t(Z) %*% X) %*% tmp
    }

  return(list(est = est, var = var[nrow(var),ncol(var)]))
}
