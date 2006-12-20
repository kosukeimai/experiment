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
