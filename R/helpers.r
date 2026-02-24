#Laplace (double exponential) distribution functions in base R

# Density: f(x) = 1/(2*b) * exp(-|x - mu| / b)
dlaplace <- function(x, mu = 0, b = 1, log = FALSE) {
  if (any(b <= 0)) stop("scale 'b' must be > 0")
  z <- abs(x - mu) / b
  logd <- -log(2 * b) - z
  if (log) logd else exp(logd)
}

# CDF: for x < mu: 0.5 * exp((x - mu) / b)
#      for x >= mu: 1 - 0.5 * exp(-(x - mu) / b)
plaplace <- function(q, mu = 0, b = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(b <= 0)) stop("scale 'b' must be > 0")
  # ensure vector recycling like base R
  q <- as.numeric(q)
  mu <- as.numeric(mu)
  b <- as.numeric(b)
  # recycle
  n <- max(length(q), length(mu), length(b))
  q <- rep(q, length.out = n)
  mu <- rep(mu, length.out = n)
  b <- rep(b, length.out = n)

  p <- numeric(n)
  left <- q < mu
  # left side
  p[left] <- 0.5 * exp((q[left] - mu[left]) / b[left])
  # right side (including equality)
  p[!left] <- 1 - 0.5 * exp(-(q[!left] - mu[!left]) / b[!left])

  if (!lower.tail) p <- 1 - p
  if (log.p) log(p) else p
}

# Quantile function (inverse CDF)
qlaplace <- function(p, mu = 0, b = 1, lower.tail = TRUE, log.p = FALSE) {
  if (any(b <= 0)) stop("scale 'b' must be > 0")
  p <- as.numeric(p)
  mu <- as.numeric(mu)
  b <- as.numeric(b)
  n <- max(length(p), length(mu), length(b))
  p <- rep(p, length.out = n)
  mu <- rep(mu, length.out = n)
  b <- rep(b, length.out = n)

  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p

  # validate p
  if (any(p < 0 | p > 1, na.rm = TRUE)) stop("p must be in [0,1]")

  q <- numeric(n)
  # handle extremes
  q[p == 0] <- -Inf
  q[p == 1] <- Inf

  mid <- (p > 0) & (p < 1)
  if (any(mid)) {
    pm <- p[mid]
    mub <- mu[mid]
    bb <- b[mid]
    left <- pm < 0.5
    # for p < 0.5: mu + b * log(2p)
    q[mid][left] <- mub[left] + bb[left] * log(2 * pm[left])
    # for p >= 0.5: mu - b * log(2*(1-p))
    q[mid][!left] <- mub[!left] - bb[!left] * log(2 * (1 - pm[!left]))
  }

  q
}

# Random generation via inverse transform
rlaplace <- function(n, mu = 0, b = 1) {
  if (length(n) != 1 || n < 0) stop("'n' must be a non-negative integer scalar")
  if (any(b <= 0)) stop("scale 'b' must be > 0")
  n <- as.integer(n)
  u <- runif(n)
  qlaplace(u, mu = mu, b = b)
}
