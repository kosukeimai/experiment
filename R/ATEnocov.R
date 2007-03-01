###
### Calculate the ATE from cluster randomized experiments
###

ATEcluster <- function(Y, Z, data = parent.frame(), grp = NULL, 
                       match = NULL, size = NULL, unbiased = TRUE){

  ## getting the data
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  size <- eval(call$size, envir = data)
  match <- eval(call$match, envir = data)

  ## Organize the data in terms of group sums
  if (is.null(grp)) { # aggregate-mean data
    if (is.null(size)) {
      stop("`size' should be specified when `grp' is NULL")
    } else {
      M <- length(Y)
      Ysum <- Y*size
    }
  } else { # individual-level data
    ugrp <- unique(grp)
    M <- length(ugrp)
    if (is.null(size)) {
      tmp <- Z
      Ysum <- Z <- size <- rep(NA, M)
      for (i in 1:M) {
        size[i] <- sum(grp == ugrp[i])
        Ysum[i] <- sum(Y[grp == ugrp[i]])
        if (length(unique(tmp[grp == ugrp[i]])) != 1)
          stop("all units in the same cluster should have the same value of `Z'")
        else
          Z[i] <- unique(tmp[grp == ugrp[i]]) 
      }
      if (!is.null(match)) {
        tmp <- match
        match <- rep(NA, M)
        umatch <- unique(match)
        if (length(umatch) != M/2)
          stop("invalid input for `match'")
        match <- rep(NA, M)
        for (i in 1:M) {
          if (length(unique(tmp[grp == ugrp[i]])) != 1)
            stop("all units in the same cluster should have the same value of `match'")
          else
            match[i] <- unique(tmp[grp == ugrp[i]])
        }
      }
    } else {
      stop("`size' should be NULL when `grp' is specified")
    }
  }

  N <- sum(size)
  m1 <- sum(Z)
  m0 <- M-m1
  if (is.null(match)) { # without matching
     ATE.est <- M*(sum(Ysum[Z==1])/m1 - sum(Ysum[Z==0])/m0)/N
     ATE.var <- (M^2)*(m0*var(Ysum[Z==1])+m1*var(Ysum[Z==0]))/(m1*m0*(N^2))
   } else { # with matching
    ## check match
    if (length(unique(table(match))) > 1)
      stop("invalid input for `match'")
    if (unique(table(match)) != 2)
      stop("invalid input for `match'")
    ATE.est <- 2*(sum(Ysum[Z==1]) - sum(Ysum[Z==0]))/N
    umatch <- unique(match)
    diff <- rep(NA, M/2)
    for (i in 1:(M/2))
      diff[i] <- Ysum[(Z == 1) & (match == umatch[i])] -
        Ysum[(Z == 0) & (match == umatch[i])] 
    ATE.var <- 2*M*var(diff)/(N^2)
  }
  
  return(list(est = ATE.est, var = ATE.var))
}

