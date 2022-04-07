#' Randomization of the Treatment Assignment for Conducting Experiments
#' 
#' This function can be used to randomize the treatment assignment for
#' randomized experiments. In addition to the complete randomization, it
#' implements randomized-block and matched-pair designs.
#' 
#' Randomized-block designs refer to the complete randomization of the
#' treatment within the pre-specified blocks which contain multiple
#' observations. Matched-pair designs refer to the randomization of the binary
#' treatment variable within the pre-specified pair of observations.
#' 
#' @aliases randomize Randomize
#' @param data A data frame containing the observations to which the treatments
#' are randomly assigned.
#' @param group A numerical or character vector indicating the
#' treatment/control groups. The length of the vector equals the total number
#' of such groups. The default specifies two groups called \dQuote{Treat} and
#' \dQuote{Control}.
#' @param ratio An optional numerical vector which specifies the proportion of
#' the treatment/control groups within the sample. The length of the vector
#' should equal the number of groups. The default is the equal allocation.
#' @param indx An optional variable name in the data frame to be used as the
#' names of the observations. If not specified, the row names of the data frame
#' will be used so long as they are available. If the row names are not
#' available, the integer sequence starting from 1 will be used.
#' @param block An optional variable name in the data frame or a formula to be
#' used as the blocking variables for randomized-block designs. If a variable
#' name is specified, then the unique values of that variable will form blocks
#' unless \code{n.block} is specified (see below). If a formula is specified,
#' it will be evaluated using \code{data} and then blocking will be based on
#' the \code{mahalanobis} distance of the resulting model matrix. In this case,
#' users may want to specify \code{n.block} to avoid creating blocks that have
#' too few observations.
#' @param n.block An optional scalar specifying the number of blocks to be
#' created for randomized block designs. If unspecified, the unique values of
#' the blocking variable will define blocks. If specified, the blocks of
#' roughly equal size will be created based on the \code{quantile} of the
#' blocking variable.
#' @param match An optional variable name in the data frame or a formula to be
#' used as the matching variables for matched-pair designs. This input is
#' applicable only to the case where there are two groups. Pairs of
#' observations will be formed based on the similar values of the matching
#' variable. If a formula is specified, the \code{mahalanobis} distance of the
#' resulting model matrix will be used.
#' @param complete logical. If it equals \code{TRUE} (default), then complete
#' randomization will be performed (within each block if randomized block
#' designs are used). Otherwise, simple randomization will be implemented. For
#' matched-pair designs, \code{complete} has to equal \code{TRUE}.
#' @return A list of class \code{randomize} which contains the following items:
#' \item{call}{ the matched call.  } \item{treatment}{ The vector of randomized
#' treatments.  } \item{data}{ The data frame that was used to conduct the
#' randomization.  } \item{block}{ The blocking variable that was used to
#' implement randomized-block designs.  } \item{match}{ The matching variable
#' that was used to implement matched-pair designs.  } \item{block.id}{ The
#' variable indicating which observations belong to which blocks in
#' randomized-block designs.  } \item{match.id}{ The variable indicating which
#' observations belong to which pairs in matched-pair designs.  }
#' @author Kosuke Imai, Department of Government and Department of Statistics, Harvard University
#' \email{imai@@Harvard.Edu}, \url{https://imai.fas.harvard.edu};
#' @keywords design
#' @export randomize
randomize <- function(data, group = c("Treat", "Control"), ratio =
                      NULL, indx = NULL, block = NULL, n.block = NULL,
                      match = NULL, complete = TRUE){  

  ## call
  call <- match.call()
  ## data 
  m <- length(group)
  if ((!is.null(call$block)) && (!is.null(call$match))) {
    stop("invalid inputs for `block' and `match'.")
  } else if (!is.null(call$block)) { ## blocking
    if ("formula" %in% class(block)) {
      tm <- terms(block)
      attr(tm, "intercept") <- 0
      data <- model.frame(tm, data = data, na.action = na.fail)
      X <- model.matrix(tm, data = data)
      block <- mahalanobis(X, apply(X, 2, mean), var(X))
    } else {
      block <- eval(call$block, envir = data)
    }
  } else { ## matching
    if (m != 2)
      stop("2 groups are required for matching.")
    if (!is.null(ratio))
      warning("`ratio' will be ignored.")
    if (complete) {
      if ("formula" %in% class(match)) {
        tm <- terms(match)
        attr(tm, "intercept") <- 0
        data <- model.frame(tm, data = data, na.action = na.fail)
        X <- model.matrix(tm, data = data)
        match <- mahalanobis(X, apply(X, 2, mean), var(X))
      } else {
        match <- eval(call$match, data)
      }
    } else {
      stop("`complete' should be TRUE for matching.")
    }
  }
    
  ## getting index
  n <- nrow(data)
  if (!is.null(indx))
    indx <- eval(indx, data)
  else if (is.null(rownames(data)))
    indx <- 1:n
  else
    indx <- rownames(data)
  
  ## groups
  if (is.null(ratio))
    ratio <- rep(1/m, m)
  ratio <- ratio/sum(ratio)
  if (sum(ratio < 0) > 1)
    stop("invalid input for `size'.")

  ## output
  res <- list(call = call, ratio = ratio)

  ## blocking and matching variable
  if (is.null(block) && is.null(match)) { 
    if (complete) { # complete randomization
      tmp <- ratio2size(n, ratio, group)
      ttt <- sample(tmp$vector, n, replace = FALSE)
    } else { # simple randomization
      ttt <- sample(group, n, replace = TRUE, prob = ratio)
    }
    names(ttt) <- indx
  } else if (is.null(match)) { ## blocking
    block.id <- rep(NA, n)
    if (is.null(n.block)) {
      tmp <- unique(block)
      n.block <- length(tmp)
      for (i in 1:n.block)
        block.id[block == tmp[i]] <- i 
    } else {
        tmp <- quantile(block, (0:(n.block-1))/n.block)
        block.id <- rep(0, n)
        for (i in 1:n.block) {
            block.id[block >= tmp[i]] <- block.id[block >= tmp[i]] + 1
        }
      if (sum(table(block.id) < m) > 0)
        stop("some blocks have too few observations.")
    }
    ttt <- rep(NA, n)
    names(ttt) <- names(block.id) <- indx
    for (i in 1:n.block) {
      howmany <- sum(block.id == i)
      if (complete) { # comlete randomization
        tmp <- ratio2size(howmany, ratio, group)
        ttt[block.id == i] <- sample(tmp$vector, howmany, replace = FALSE)
      } else {
        ttt[block.id == i] <- sample(group, sum(block.id == i),
            replace = TRUE, prob = ratio)
      }
    }
    res$block <- block
    res$block.id <- block.id
  } else { ## matching
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
    res$match <- match
    res$match.id <- match.id
  }

  ## return the results
  res$treatment <- ttt
  res$data <- data
  class(res) <- "randomize"
  
  return(res)
}

###
### This converts ratio into size while randomly allocating remainders 
###

ratio2size <- function(n, ratio, group) {
  m <- length(ratio)

  size <- round(ratio * n, digits = 0)
  if (sum(size) > n) {
    tmp <- sample(1:length(size), sum(size)-n, replace = FALSE)
    size[tmp] <- size[tmp] - 1
  }
  if (sum(size) < n) {
    tmp <- sample(1:length(size), n-sum(size), replace = FALSE)
    size[tmp] <- size[tmp] + 1
  }
  allgroup <- NULL
  for (i in 1:m)
    allgroup <- c(allgroup, rep(group[i], size[i]))
  return(list(size = size, vector = allgroup))
  
}
