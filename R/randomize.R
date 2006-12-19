randomize <- function(data, group = c("Treat", "Control"), ratio =
                      NULL, indx = NULL, block = NULL, n.block = NULL,
                      match = NULL){  

  ## getting index
  n <- nrow(data)
  if (!is.null(indx))
    indx <- eval(indx, data)
  else if (is.null(rownames(data)))
    indx <- 1:n
  else
    indx <- rownames(data)
  
  ## groups
  m <- length(group)
  if (is.null(ratio))
    ratio <- 1/m
  ratio <- ratio/sum(ratio)
  if (sum(ratio <= 0) > 0)
    stop("invalid input for `ratio'.")
  
  ## blocking and matching variable
  if (is.null(block) & is.null(match)) {
    block <- rep(1, n)
    ttt <- sample(group, n, replace = TRUE, prob = ratio)
  } else if (!is.null(block) && !is.null(match)) {
    stop("invalid inputs for `block' and `match'.")
  } else if (is.null(match)) { ## blocking
    if (is.formula(block)) {
      tm <- terms(block)
      attr(tm, "intercept") <- 0
      data <- model.frame(tm, data = data, na.action = na.fail)
      X <- model.matrix(tm, data = data)
      block <- mahalanobis(X, apply(X, 2, mean), var(X))
    } else {
      block <- eval(block, data)
    }
    block.id <- rep(0, n)
    if (is.null(n.block)) {
      tmp <- unique(block)
      n.block <- length(n.block)
      for (i in 1:n.block)
        block.id[block == tmp[i]] <- i 
    } else {
      tmp <- quantile(block, (0:(n.block-1))/n.block)
      for (i in 1:n.block)
        block.id[block >= tmp[i]] <- block.id[block >= tmp[i]] + 1
    }
    ttt <- rep(NA, n)
    names(ttt) <- names(block.id) <- indx
    for (i in 1:n.block)
      ttt[block.id == i] <- sample(group, sum(block.id == i),
            replace = TRUE, prob = ratio) 
  } else { ## matching
    if (length(group) != 2)
      stop("2 groups are required for matching.")
    if (is.formula(match)) {
      tm <- terms(match)
      attr(tm, "intercept") <- 0
      data <- model.frame(tm, data = data, na.action = na.fail)
      X <- model.matrix(tm, data = data)
      match <- mahalanobis(X, apply(X, 2, mean), var(X))
    } else {
      match <- eval(match, data)
    }
    match.id <- ttt <- rep(NA, n)
    names(match.id) <- names(ttt) <- indx
    counter <- 1
    while (sum(is.na(match.id)) > 1) {
      unit <- sample(indx[is.na(match.id)], 1, replace = FALSE)
      diff <- abs(match[is.na(match.id)]-match[unit])
      mindiff <- names(sort(diff[diff>0]))[1]
      match.id[unit] <- match.id[mindiff] <- counter
      tmp <- sample(group, 2, replace = FALSE)
      ttt[unit] <- tmp[1]
      ttt[mindiff] <- tmp[2]
      counter <- counter + 1
    }
  }

  ## return the results
  return(list(treatment = ttt, data = data, block = block, match =
              match, block.id = block.id, match.id = match.id))
}
