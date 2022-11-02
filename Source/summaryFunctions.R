col_mean <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- mean(x[[i]])
  }
  output
}
col_sd <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- sd(x[[i]])
  }
  output
}
col_cv <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- sd(x[[i]]/mean(x[[i]]))
  }
  output
}
col_lci <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- qs(x[[i]],0.025)
  }
  output
}
col_uci <- function(x,y) {
  output <- vector("double", length(x))
  for (i in seq_along(y)) {
    output[[i]] <- qs(x[[i]],0.975)
  }
  output
}
qs <- function(x,y){as.numeric(quantile(x,y))}
