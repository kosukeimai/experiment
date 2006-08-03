bprobitMixed <- function(Y, X, Z, grp, beta.start, 
                         Phi.start, beta0, A0, df, T0, mda =
                         FALSE, sims) {
  
  ## this code assumes the equal number of obs within each group
  res <- .C("R2bprobitMixedGibbs", as.integer(Y), as.double(X),
            as.double(Z), as.integer(grp), as.double(beta.start),
            as.double(Psi.start),
            as.integer(nrow(X)), as.integer(ncol(X)), as.integer(ncol(Z)),
            as.integer(length(table(grp))), as.integer(table(grp)),
            as.integer(max(table(grp))), as.double(beta0), as.double(A0),
            as.integer(df), as.double(T0), as.integer(mda),
            as.integer(sims), betaStore = double(sims*ncol(X)),
            gammaStore = double(sims*ncol(Z)*ngrp),
            PsiStore = double(sims*ncol(Z)*(ncol(Z)+1)/2),
            PACKAGE = "are")

  return(list(beta = matrix(res$betaStore, byrow = TRUE, ncol = ncol(X)),
              gamma = array(res$gammaStore, dim = c(ncol(Z), ngrp, sims)),
              Psi = matrix(res$PsiStore, byrow = TRUE, nrow = sims)))
}
