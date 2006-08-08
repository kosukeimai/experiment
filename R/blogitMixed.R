blogitMixed <- function(Y, X, Z, grp, beta.start,
                        Phi.start, beta0, A0, df, T0, sims, tune.fixed,
                        tune.random) {
  
  ## this code assumes the equal number of obs within each group
  res <- .C("R2logitMixedMetro", as.integer(Y), as.double(X),
            as.double(Z), as.integer(grp), as.double(beta.start),
            as.double(Psi.start), as.integer(nrow(X)),
            as.integer(max(Y)),
            as.integer(ncol(X)), as.integer(ncol(Z)),
            as.integer(length(table(grp))),
            as.integer(max(table(grp))),
            as.double(beta0), as.double(A0),
            as.integer(df), as.double(T0),
            as.double(tune.fixed), as.double(tune.random),
            as.integer(sims),
            acc_fixed = integer(ncol(X)*max(Y)),
            acc_random = integer(length(table(grp))*max(Y)),
            betaStore = double(sims*ncol(X)*max(Y)),
            PsiStore = double(sims*max(Y)*ncol(Z)*(ncol(Z)+1)/2),
            PACKAGE = "experiment")

  return(list(beta = matrix(res$betaStore, byrow = TRUE,
                ncol = ncol(X)*max(Y)),
              Psi = matrix(res$PsiStore, byrow = TRUE,
                ncol = length(table(grp))*ncol(Z)*(ncol(Z)+1)/2),
              beta.ratio = res$acc_fixed/sims,
              gamma.ratio = res$acc_random/sims))
}
