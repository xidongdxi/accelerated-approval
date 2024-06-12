# Function to calculate the one-sided Type I error rate
type1_fun <- function(alpha, r, rho, p_S, p_F) {
  p_S + pmvnorm(lower = c(qnorm(1 - p_F), qnorm(1 - alpha)),
               upper = c(qnorm(1 - p_S), Inf),
               mean = rep(0, 2),
               corr = matrix(c(1, rho * sqrt(r), rho * sqrt(r), 1), nrow = 2),
               algorithm = GenzBretz(maxpts = 25000, abseps = 1e-6, releps = 0)
  )
}

# Function to calculate the difference between the Type I error and alpha
type1_fun_alpha <- function(alpha, r, rho, p_S, p_F) {
  type1_fun(alpha, r, rho, p_S, p_F) - alpha
}

# Function to search for the futility boundaryp_F to control the Type I error
# at level alpha for all possible values of the correlation between ORR and OS
type1_search <- function(alpha, r, p_S, p_F) {
  optimize(type1_fun, c(-1, 1), maximum = T, tol = .Machine$double.eps^0.5,
           alpha = alpha, r = r, p_S = p_S, p_F = p_F)$objective - alpha
}

# Function to calculate the ORR difference
orr_fun <- function(x, orr_c, nper, alpha) {
  (x - orr_c) / sqrt(x * (1 - x) / nper + orr_c * (1 - orr_c) / nper) - qnorm(1 - alpha)
}
