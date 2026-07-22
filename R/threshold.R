#' Threshold
#'
#' Gives the coefficient magnitude at which the spike overtakes the slab
#'
#' @param ss
#' @param sparsity
#'
#' @returns
#' @export
#'
#' @examples
# threshold <- function(ss, sparsity){
#   if(length(ss)!=2) stop("ss must be a vector of length 2")
#   if(ss[1]>ss[2]){
#     ss<- sort(ss)
#     print("ss warning, scale values supplied out of order. Spike scale taken as ss[1] and slab scale taken as ss[2].")
#   }
#
#   ss0 <- ss[1]
#   ss1 <- ss[2]
#
#   lambda0 <- 1/ss0
#   lambda1 <- 1/ss1
#
#   delta = (1/(lambda0-lambda1))*log(lambda0/lambda1 * (1-sparsity)/sparsity)
#   delta
# }


intersection_point <- function(lambda0, lambda1, theta){
  1/(lambda0-lambda1)*log(lambda0/lambda1 * (1-theta)/theta)
}

lambda1_range <- function(n, p){
  min_lambda1 = sqrt(n)/p
  max_lambda1 = 4*sqrt(n*log(p))
  ret <- c(min_lambda1, max_lambda1)
  ret
}

# lambda0_range <- function(n, ehat, theta, lambda1){
#   lambda1*(exp(1/(2*n)*(ehat-lambda1)^2)-1)*theta/(1-theta)
# }

lambda0_lower <- function(n, lambda1){
  sqrt(n)/2 + lambda1
}

threshold_upper <- function(n, lambda0, lambda1, theta){
  sqrt(2*n*log(1 + lambda0/lambda1 *theta/(1-theta)))+ lambda1
}

threshold_lower <- function(n, lambda0, lambda1, theta, d){
  sqrt(2*n*log(1 + lambda0/lambda1 *theta/(1-theta))-d)+lambda1
}

pstar_zero <- function(lambda0, lambda1, theta){
 1/( 1 + lambda0/lambda1 * (1-theta)/theta)
}

d_bound <- function(lambda0, lambda1, n){
  2*n -(1/(lambda0-lambda1) - sqrt(2*n))^2
}

c_plus <- function(n, lambda0, lambda1){
  0.5*(1+sqrt(1-4*n/(lambda0-lambda1)^2))
}


pstar <- function(x, lambda0, lambda1, theta){
  denom <- 1 + (1-theta)/theta *
           (lambda0/lambda1)*exp(-abs(x)*(lambda0-lambda1))
  1/denom
}

solve_for_x <- function(n, lambda0, lambda1, theta) {
  # 1. Calculate the target value C from c_plus
  C <- c_plus(n, lambda0, lambda1)

  # Check if C is NaN (happens if the term under the sqrt is negative)
  if (is.nan(C)) {
    warning("c_plus yielded NaN. Check if 4*n <= (lambda0 - lambda1)^2")
    return(c(NA, NA))
  }

  # 2. Define intermediate variables A and B
  A <- ((1 - theta) / theta) * (lambda0 / lambda1)
  B <- lambda0 - lambda1

  # 3. Calculate the inner value of the logarithm
  inner_val <- (1 - C) / (A * C)

  # The value inside the log must be strictly positive for a real solution
  if (inner_val <= 0) {
    warning("No real solution exists for x given these parameters.")
    return(c(NA, NA))
  }

  # 4. Solve for |x|
  abs_x <- -log(inner_val) / B

  # Return both the positive and negative roots
  return(c(abs_x, -abs_x))
}


s1_range <- function(n, p){
  max_s1 = p/sqrt(n)
  min_s1 = 1/(4*sqrt(n*log(p)))
  ret <- c(min_s1, max_s1)
  ret
}

s0_ladder <- function(s1, p, step){
  lambda1 <- 1/s1
  top_of_range <- p^2
  bottom_of_range <- lambda1 + step

  if (bottom_of_range > top_of_range) {
    return(numeric(0))
  }

  1/seq(from = bottom_of_range, to = top_of_range, by = step)
}

s0_ladder(10, 100, 100)
