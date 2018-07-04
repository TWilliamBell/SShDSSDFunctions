library(mise)
library(matrixStats)

suppressWarnings(origpar <- par())

suppressWarnings(origpar$cin <- origpar$cra <- origpar$csi <- origpar$cxy <- origpar$din <- origpar$page <- NULL)

clear <- function(vars=FALSE,...) {
  mise(vars = vars,...)
}

euclidean <- function(vector) { # I know there's a norm function but this one only has one argument,
  # so it is simpler.  The other one doesn't default to the Euclidean norm.
  (sum(vector^2))^0.5
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

clear()