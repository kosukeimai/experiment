bNegBinMixed <- function(Y, X, Z, grp, beta.start, gamma.start,
                         sig2.start, Phi.start, beta0, A0, a0,
                         b0, tau0, T0, varb, vars, varg, sims) {
  
  ## this code assumes the equal number of obs within each group
  counter <- rep(0, 2)
  counterg <- rep(0, length(table(grp))*2)
  res <- .C("R2bnegbinMixedMCMC", as.double(Y), as.double(X),
            as.double(Z), as.integer(grp), as.double(beta.start),
            as.double(gamma.start), as.double(sig2.start),
            as.double(solve(Psi.start)),
            as.integer(nrow(X)), as.integer(ncol(X)), as.integer(ncol(Z)),
            as.integer(length(table(grp))), 
            as.integer(max(table(grp))), as.double(beta0), as.double(A0),
            as.double(a0), as.double(b0), as.integer(tau0),
            as.double(T0), as.double(varb), as.double(vars),
            as.double(varg), as.integer(sims),
            counter = integer(2), counterg = integer(length(table(grp))),
            betaStore = double(sims*ncol(X)),
            gammaStore = double(sims*ncol(Z)*ngrp),
            sig2Store = double(sims),
            PsiStore = double(sims*ncol(Z)*(ncol(Z)+1)/2),
            PACKAGE = "experiment")

  return(list(beta = matrix(res$betaStore, byrow = TRUE, ncol = ncol(X)),
              sig2 = res$sig2Store,
              gamma = array(res$gammaStore, dim = c(ncol(Z), ngrp, sims)),
              Psi = matrix(res$PsiStore, byrow = TRUE, nrow = sims),
              accept = res$counter/sims,
              acceptg = res$counterg/sims))
}
