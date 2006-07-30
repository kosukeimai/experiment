bNormalMixed <- function(Y, X, Z, grp, beta.start, gamma.start,
                         Phi.start, beta0, A0, nu0, s0, tau0, T0,
                         imp = TRUE, sims) {
  
  ## this code assumes the equal number of obs within each group
  res <- .C("R2bNormalMixedGibbs", as.double(Y), as.double(X),
            as.double(Z), as.integer(grp), as.double(beta.start),
            as.double(gamma.start), as.double(solve(Psi.start)),
            as.integer(nrow(X)), as.integer(ncol(X)), as.integer(ncol(Z)),
            as.integer(ngrp), as.integer(table(ngrp)),
            as.integer(max(table(grp))), as.double(beta0), as.double(A0),
            as.integer(imp), as.integer(nu0), as.integer(s0),
            as.integer(tau0), as.double(T0), as.integer(sims),
            betaStore = double(sims*ncol(X)),
            gammaStore = double(sims*ncol(Z)*ngrp),
            sig2Store = double(sims),
            PsiStore = double(sims*ncol(Z)*(ncol(Z)+1)/2),
            PACKAGE = "are")

  return(list(beta = matrix(res$betaStore, byrow = TRUE, ncol = ncol(X)),
              gamma = array(res$gammaStore, dim = c(ncol(Z), ngrp, sims)),
              Psi = matrix(res$PsiStore, byrow = TRUE, nrow = sims)))
}
