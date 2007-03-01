###
### Calculate the ATE from cluster randomized experiments
###
### unbiased = TRUE: my own estimator 
### unbiased = FALSE: standard estimator using the ``average method''
###            For matched-pair designs, see Donner and Klar
###            (Statistics in Medicine, 1987; Journal of Clinical Epidemiology, 1993)
###            A simple weight in Donner and Klar (1993) is used.
###

ATEnocov <- function(Y, Z, data = parent.frame(), grp = NULL,
                     match = NULL, size = NULL, unbiased = TRUE){

  ## an internal function that checks match and returns diff
  match.check <- function(Y, Z, match) { 
    n <- length(Y)
    if ((n %/% 2) != 0)
      stop("pair randomization requires the even number of observations")
    if (length(unique(table(match))) > 1)
      stop("invalid input for `match'")
    if (unique(table(match)) != 2)
      stop("invalid input for `match'")
    umatch <- unique(match)
    diff <- rep(NA, n/2)
    for (i in 1:(n/2))
      diff[i] <- Y[(Z == 1) & (match == umatch[i])] -
        Y[(Z == 0) & (match == umatch[i])] 
    return(diff)
  }
  
  ## getting the data
  call <- match.call()
  Y <- eval(call$Y, envir = data)
  Z <- eval(call$Z, envir = data)
  grp <- eval(call$grp, envir = data)
  size <- eval(call$size, envir = data)
  match <- eval(call$match, envir = data)

  ## checking Z
  if (sum(sort(unique(Z)) != c(0,1)) != 2)
    stop("`Z' should be binary taking the value of 0 or 1")

  ## ATE for unit randomization 
  if (is.null(grp) && is.null(size)) { 
    ATE.est <- mean(Y[Z==1])-mean(Y[Z==0])
    if (is.null(match)) { # without matching
      ATE.var <- var(Y[Z==1])/sum(Z==1)+var(Y[Z==0])/sum(Z==0)
    } else { # with matching
      diff <- match.check(Y, Z, match)
      ATE.var <- var(diff)/(length(Y)/2)
    }
    return(list(est = ATE.est, var = ATE.var, Y = Y, Z = Z, match = match))
  }

  ## ATE for cluster randomizaion
  if (is.null(grp)) { # aggregate data input
    if (!unbiased && !is.null(match)) { 
      stop("the input should be individual-level data for this estimator")
    } else {
      M <- length(Y)
      Ysum <- Y*size
    } 
  } else if (unbiased || !is.null(match)) { # individual data input
    ugrp <- unique(grp)
    M <- length(ugrp)
    if (is.null(size)) {
      tmp <- Z
      Ysum <- Z <- size <- rep(NA, M)
      for (i in 1:M) {
        size[i] <- sum(grp == ugrp[i])
        Ysum[i] <- sum(Y[grp == ugrp[i]])
        Zvalue <- unique(tmp[grp == ugrp[i]]) 
        if (length(Zvalue) != 1)
          stop("all units in the same cluster should have the same value of `Z'")
        else
          Z[i] <- Zvalue
      }
      if (!is.null(match)) {
        tmp <- match
        umatch <- unique(match)
        if (length(umatch) != M/2)
          stop("invalid input for `match'")
        match <- rep(NA, M)
        for (i in 1:M) {
          Mvalue <- unique(tmp[grp == ugrp[i]])
          if (length(Mvalue) != 1)
            stop("all units in the same cluster should have the same value of `match'")
          else
            match[i] <- Mvalue
        }
      }
    } else {
      stop("`size' should be NULL when `grp' is specified")
    }
  }
  
  if (unbiased) { ## my method
    N <- sum(size)
    m1 <- sum(Z)
    m0 <- M-m1
    if (is.null(match)) { # without matching
      ATE.est <- M*(sum(Ysum[Z==1])/m1 - sum(Ysum[Z==0])/m0)/N
      ATE.var <- (M^2)*(m0*var(Ysum[Z==1])+m1*var(Ysum[Z==0]))/(m1*m0*(N^2))
    } else { # with matching
      ATE.est <- 2*(sum(Ysum[Z==1]) - sum(Ysum[Z==0]))/N
      diff <- match.check(Ysum, Z, match)
      ATE.var <- 2*M*var(diff)/(N^2)
    }
    return(list(est = ATE.est, var = ATE.var, Ysum = Ysum, Z = Z,
                size = size, match = match, unbiased = unbiased))  
  } else { ## Donner and Klar
    if (is.null(match)) { # average method 
      ATE.est <- mean(Y[Z==1]) - mean(Y[Z==0])
      ATE.var <- varCluster(Y[Z==1], grp[Z==1])$var + varCluster(Y[Z==0], grp[Z==0])$var
      return(list(est = ATE.est, var = ATE.var, Y = Y, Z = Z, grp = grp))         
    } else { # with matching
      Y <- Ysum/size
      Y1 <- Y[Z==1]
      Y0 <- Y[Z==0]
      w1 <- size[Z==1]
      w0 <- size[Z==0]
      ATE.est <- weighted.mean(Y1, w1) - weighted.mean(Y0, w0)
      ind0 <- sort(match[Z==0], index.return = TRUE)
      ind1 <- sort(match[Z==1], index.return = TRUE)
      if (sum(ind0$x == ind1$x) != sum(Z==1))
        stop("invalid input for `match'.")
      D <- Y1[ind1$ix] - Y0[ind0$ix]
      w <- w1[ind1$ix] * w0[ind0$ix]/(w1[ind1$ix] + w0[ind0$ix])
      ATE.var <- weighted.var(D, w)*sum(w^2)/(sum(w)^2)
      return(list(est = ATE.est, var = ATE.var, Y = Y, Z = Z, grp = grp,
                  match = match))
    }  
  }
}

