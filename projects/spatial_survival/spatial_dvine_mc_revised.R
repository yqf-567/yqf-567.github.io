# ============================================================================
# Spatial trivariate sequential-survival D-vine Monte Carlo simulation
# Revised implementation matching Sections 2-3 of spatial_survival.pdf.
#
# Dependencies: install.packages("mvtnorm")
# Uses cross-platform PSOCK parallel workers; no numDeriv dependency.
#
# QUICK START:
#   Rscript spatial_dvine_mc_revised.R
#   Change RUN_MODE <- "full" below for all 12 scenarios.
#
# IMPORTANT: The working stage-1 likelihood for later gaps is intentionally
# retained as in Equation (5). It can be biased under informative entry.
# No claim is made that the output reproduces the PDF's old tables exactly.
# ============================================================================

if (!requireNamespace("mvtnorm", quietly = TRUE)) {
  stop("Package 'mvtnorm' is required. Run install.packages('mvtnorm').")
}

# ========================= USER-EDITABLE OPTIONS ============================
RUN_MODE <- "full"  # "quick" (small diagnostic run), or "full" (PDF design)

CFG <- list(
  sample_sizes  = c(120L, 240L, 1200L),
  spatial_tau   = c(0.5, 0.9),
  censor_rates  = c(0.20, 0.40),
  replications  = 500L,
  cores         = 2L,          # PSOCK works on Windows, Linux and macOS
  master_seed   = 20261102L,
  calibration_n = 100000L,
  calculate_MT  = TRUE,       # TRUE: calculate cluster Murphy-Topel SEs
  score_step    = 1e-4,       # central-difference score step
  hessian_step  = 5e-3,       # larger step stabilizes second differences
  tvpack_eps    = 1e-8,       # deterministic multivariate Gaussian CDF
  maxit_stage1  = 160L,
  maxit_stage2  = 110L,
  output_prefix = "spatial_dvine_mc",
  alpha         = c(0, 0, 0),
  beta          = c(0.45, -0.35, 0.65),
  sigma         = c(0.85, 0.85, 0.85),
  rho           = c(0.55, 0.45, -0.35), # (rho12, rho23, rho13|2)
  positions     = rbind(c(0, 0), c(0.31, 0), c(0.10, 0.35)),
  spatial_range = 0.45
)

if (identical(RUN_MODE, "quick")) {
  CFG$sample_sizes <- 120L
  CFG$spatial_tau  <- 0.5
  CFG$censor_rates <- 0.2
  CFG$replications <- 5L
  CFG$output_prefix <- "spatial_dvine_mc_quick"
} else if (!identical(RUN_MODE, "full")) {
  stop("RUN_MODE must be 'quick' or 'full'.")
}
# ============================================================================

# Transform unconstrained conditional correlations into a valid 3x3 matrix.
# For a Gaussian D-vine: rho13 = rho12*rho23 + rho13|2*sqrt(...).
vine_cor <- function(z) {
  r <- tanh(z)
  r13 <- r[1L] * r[2L] + r[3L] *
    sqrt((1 - r[1L]^2) * (1 - r[2L]^2))
  matrix(c(1, r[1L], r13,
           r[1L], 1, r[2L],
           r13, r[2L], 1), nrow = 3L)
}

# Observed-data log likelihood of a jointly normal vector.
# 'exact' are log-times with observed events, 'cens' are lower limits
# (right-censored on the original gap-time scale). Unentered gaps have
# exact=FALSE and cens=FALSE, contributing no density or survival term.
# The normal integration dimension is <= 3 because a subject can have
# at most one censored gap across the three consecutive events.
loglik_censored <- function(mu, V, y, lo, exact, cens, cdf_eps) {
  e <- which(exact)
  cidx <- which(cens)
  if (length(e) == 0L && length(cidx) == 0L) return(0.0)

  ans <- 0.0
  if (length(e) > 0L) {
    W <- V[e, e, drop = FALSE]
    delta <- y[e] - mu[e]
    U <- chol(W)
    q <- backsolve(U, delta, transpose = TRUE)
    ans <- -0.5 * (length(e) * log(2 * pi) +
                   2 * sum(log(diag(U))) + sum(q^2))
  }
  if (length(cidx) == 0L) return(ans)

  if (length(e) > 0L) {
    # M = Cov(c,e) %*% Var(e)^(-1), evaluated by solving a linear system.
    M <- t(solve(W, V[e, cidx, drop = FALSE]))
    cond_mu <- as.vector(mu[cidx] + M %*% (y[e] - mu[e]))
    cond_V <- V[cidx, cidx, drop = FALSE] -
      M %*% V[e, cidx, drop = FALSE]
  } else {
    cond_mu <- mu[cidx]
    cond_V <- V[cidx, cidx, drop = FALSE]
  }
  cond_V <- 0.5 * (cond_V + t(cond_V))

  if (length(cidx) == 1L) {
    log_pr <- stats::pnorm(lo[cidx], mean = cond_mu,
                           sd = sqrt(cond_V[1L, 1L]),
                           lower.tail = FALSE, log.p = TRUE)
  } else if (length(cidx) <= 3L) {
    # TVPACK is deterministic for the present 2D/3D semi-infinite integral.
    pr <- as.numeric(mvtnorm::pmvnorm(
      lower = lo[cidx], upper = rep(Inf, length(cidx)),
      mean = cond_mu, sigma = cond_V,
      algorithm = mvtnorm::TVPACK(abseps = cdf_eps)
    ))
    if (!is.finite(pr)) stop("Nonfinite multivariate normal probability")
    log_pr <- log(pmax(pmin(pr, 1), 1e-300))
  } else {
    stop("At most 3 censored coordinates per block are expected.")
  }
  ans + log_pr
}

# Make one simulated sample. Crucially, subject rows are ordered as
# block 1 site 1, block 1 site 2, block 1 site 3, block 2 site 1, ...
# (as.vector(t(matrix))) rather than column-major site grouping.
simulate_blocks <- function(n, tau, censor_upper, cfg) {
  m <- nrow(cfg$positions)
  B <- n %/% m
  K <- cfg$K
  X <- matrix(stats::rnorm(n), nrow = B, ncol = m)
  xvec <- as.vector(t(X))

  # Each row is a subject; correlation between gaps is a Gaussian D-vine.
  residuals <- mvtnorm::rmvnorm(n, sigma = vine_cor(atanh(cfg$rho)))
  residuals <- sweep(residuals, 2L, cfg$sigma, `*`)

  log_gap <- matrix(NA_real_, nrow = n, ncol = 3L)
  for (j in seq_len(3L)) {
    # Across sites within the same block, gap-specific frailties have tau^2*K.
    Wj <- mvtnorm::rmvnorm(B, sigma = tau^2 * K)
    wvec <- as.vector(t(Wj))
    log_gap[, j] <- cfg$alpha[j] + cfg$beta[j] * xvec +
      residuals[, j] + wvec
  }

  G <- exp(log_gap)
  remaining <- stats::runif(n, 0, censor_upper) # one end time per subject
  observed_y <- matrix(NA_real_, n, 3L)
  lower_y <- matrix(NA_real_, n, 3L)
  exact <- matrix(FALSE, n, 3L)
  cens <- matrix(FALSE, n, 3L)
  entered <- rep(TRUE, n)

  for (j in seq_len(3L)) {
    event <- entered & (G[, j] <= remaining)
    censored <- entered & !event
    observed_y[event, j] <- log_gap[event, j]
    lower_y[censored, j] <- log(pmax(remaining[censored],
                                    .Machine$double.xmin))
    exact[event, j] <- TRUE
    cens[censored, j] <- TRUE
    remaining <- pmax(remaining - G[, j], 0)
    entered <- event
  }

  dat <- lapply(seq_len(B), function(b) {
    ix <- ((b - 1L) * m + 1L):(b * m)
    list(x = xvec[ix],
         y = observed_y[ix, , drop = FALSE],
         lo = lower_y[ix, , drop = FALSE],
         exact = exact[ix, , drop = FALSE],
         cens = cens[ix, , drop = FALSE])
  })

  list(blocks = dat,
       realized_censor = mean(!exact[, 3L]),
       entry_gap2 = mean(exact[, 1L]),
       entry_gap3 = mean(exact[, 2L]),
       event_gap3 = mean(exact[, 3L]))
}

# Marginal block log likelihood for a single gap (working stage 1).
block_marginal_loglik <- function(p, j, d, cfg) {
  V <- exp(2 * p[3L]) * diag(nrow(cfg$K)) +
    exp(2 * p[4L]) * cfg$K
  mu <- p[1L] + p[2L] * d$x
  loglik_censored(mu, V, d$y[, j], d$lo[, j],
                  d$exact[, j], d$cens[, j], cfg$tvpack_eps)
}

marginal_nll <- function(p, j, dat, cfg) {
  V <- exp(2 * p[3L]) * diag(nrow(cfg$K)) +
    exp(2 * p[4L]) * cfg$K
  -sum(vapply(dat, function(d) {
    mu <- p[1L] + p[2L] * d$x
    loglik_censored(mu, V, d$y[, j], d$lo[, j],
                    d$exact[, j], d$cens[, j], cfg$tvpack_eps)
  }, numeric(1L)))
}

# Site-major ordering: (site1 gap1, gap2, gap3, site2 gap1, ...).
# Independent gap-specific spatial frailties and a within-subject D-vine.
joint_covariance <- function(th, z, cfg) {
  P <- matrix(th, nrow = 4L, ncol = 3L)
  sig <- exp(P[3L, ])
  tau <- exp(P[4L, ])
  m <- nrow(cfg$K)
  kronecker(cfg$K, diag(tau^2, 3L)) +
    kronecker(diag(m), diag(sig, 3L) %*% vine_cor(z) %*%
                diag(sig, 3L))
}

# Full observed-history block log likelihood, Eq. (6)-(8) in the PDF.
block_joint_loglik <- function(th, z, d, cfg) {
  P <- matrix(th, nrow = 4L, ncol = 3L)
  mu <- as.vector(t(sweep(outer(d$x, P[2L, ]), 2L, P[1L, ], `+`)))
  loglik_censored(mu, joint_covariance(th, z, cfg),
                  as.vector(t(d$y)), as.vector(t(d$lo)),
                  as.vector(t(d$exact)), as.vector(t(d$cens)),
                  cfg$tvpack_eps)
}

joint_nll <- function(z, th, dat, cfg) {
  P <- matrix(th, nrow = 4L, ncol = 3L)
  V <- joint_covariance(th, z, cfg)  # same covariance for every block
  -sum(vapply(dat, function(d) {
    mu <- as.vector(t(sweep(outer(d$x, P[2L, ]), 2L, P[1L, ], `+`)))
    loglik_censored(mu, V, as.vector(t(d$y)),
                    as.vector(t(d$lo)), as.vector(t(d$exact)),
                    as.vector(t(d$cens)), cfg$tvpack_eps)
  }, numeric(1L)))
}

# In the PDF L2 = full observed-data likelihood / product(stage-1 marginal
# likelihoods). The denominator does not depend on z (the vine parameters),
# so optimizing joint_nll is IDENTICAL to optimizing Eq. (8) in z. Mixed
# z/theta derivatives and z scores used for Murphy-Topel are also identical.

# Numerical derivatives based on deterministic likelihoods. Central
# differences avoid the much higher number of calls of Richardson defaults.
fd_gradient <- function(f, x, hscale) {
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    h <- hscale * max(1, abs(x[i]))
    xp <- xm <- x
    xp[i] <- xp[i] + h
    xm[i] <- xm[i] - h
    out[i] <- (f(xp) - f(xm)) / (2 * h)
  }
  out
}

fd_hessian <- function(f, x, hscale) {
  p <- length(x)
  H <- matrix(0, p, p)
  f0 <- f(x)
  h <- hscale * pmax(1, abs(x))
  for (i in seq_len(p)) {
    xp <- xm <- x
    xp[i] <- xp[i] + h[i]
    xm[i] <- xm[i] - h[i]
    H[i, i] <- (f(xp) - 2 * f0 + f(xm)) / h[i]^2
    if (i < p) for (j in (i + 1L):p) {
      pp <- pm <- mp <- mm <- x
      pp[c(i, j)] <- x[c(i, j)] + c(h[i], h[j])
      pm[c(i, j)] <- x[c(i, j)] + c(h[i], -h[j])
      mp[c(i, j)] <- x[c(i, j)] + c(-h[i], h[j])
      mm[c(i, j)] <- x[c(i, j)] - c(h[i], h[j])
      H[i, j] <- (f(pp) - f(pm) - f(mp) + f(mm)) /
        (4 * h[i] * h[j])
      H[j, i] <- H[i, j]
    }
  }
  H
}

fd_jacobian <- function(f, x, hscale) {
  f0 <- f(x)
  J <- matrix(NA_real_, nrow = length(f0), ncol = length(x))
  for (j in seq_along(x)) {
    h <- hscale * max(1, abs(x[j]))
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    J[, j] <- (f(xp) - f(xm)) / (2 * h)
  }
  J
}

# Convert 12 unconstrained marginal parameters and three vine parameters
# to the 15 natural parameters shown in the paper.
natural_parameters <- function(th, z) {
  P <- matrix(th, 4L, 3L)
  c(P[1L, ], P[2L, ], exp(P[3L, ]), exp(P[4L, ]), tanh(z))
}

# Cluster-stacked Murphy-Topel covariance; spatial blocks are independent.
# Score vectors have 12 stage-1 components + 3 stage-2 components.
mt_standard_errors <- function(th, z, dat, cfg) {
  B <- length(dat)
  S <- matrix(NA_real_, nrow = B, ncol = 15L)
  H11 <- matrix(0, 12L, 12L)

  for (j in seq_len(3L)) {
    idx <- (4L * (j - 1L) + 1L):(4L * j)
    p <- th[idx]
    H11[idx, idx] <- fd_hessian(
      function(v) marginal_nll(v, j, dat, cfg),
      p, cfg$hessian_step)
  }

  H22 <- fd_hessian(
    function(v) joint_nll(v, th, dat, cfg),
    z, cfg$hessian_step)
  H21 <- fd_jacobian(function(a) {
    fd_gradient(function(v) joint_nll(v, a, dat, cfg),
                z, cfg$score_step)
  }, th, cfg$hessian_step)

  for (b in seq_len(B)) {
    d <- dat[[b]]
    for (j in seq_len(3L)) {
      idx <- (4L * (j - 1L) + 1L):(4L * j)
      S[b, idx] <- fd_gradient(
        function(v) block_marginal_loglik(v, j, d, cfg),
        th[idx], cfg$score_step)
    }
    S[b, 13L:15L] <- fd_gradient(
      function(v) block_joint_loglik(th, v, d, cfg),
      z, cfg$score_step)
  }

  # A = negative derivative of summed stacked scores.
  A <- rbind(cbind(H11, matrix(0, 12L, 3L)),
             cbind(H21, H22))
  if (any(!is.finite(A)) || any(!is.finite(S))) {
    stop("Nonfinite entries in stacked scores / derivative matrix")
  }
  Ainv <- solve(A)
  V <- Ainv %*% crossprod(S) %*% t(Ainv)
  V <- 0.5 * (V + t(V))

  # Delta-method Jacobian: (alpha1:3, beta1:3, sigma1:3,
  # tau1:3, rho12, rho23, rho13|2).
  J <- matrix(0, 15L, 15L)
  P <- matrix(th, 4L, 3L)
  for (j in seq_len(3L)) {
    k <- 4L * (j - 1L)
    J[j, k + 1L] <- 1
    J[j + 3L, k + 2L] <- 1
    J[j + 6L, k + 3L] <- exp(P[3L, j])
    J[j + 9L, k + 4L] <- exp(P[4L, j])
    J[j + 12L, j + 12L] <- 1 - tanh(z[j])^2
  }
  vars <- diag(J %*% V %*% t(J))
  if (any(!is.finite(vars)) || any(vars < -1e-8)) {
    stop("Murphy-Topel covariance has invalid diagonal entries")
  }
  sqrt(pmax(vars, 0))
}

# Independent calibration: P(T_3 > C), C~Uniform(0,upper), should equal
# target. T_3 is the sum of all three uncensored gap times.
# Marginally at a site, the three spatial frailties are independent N(0,tau^2).
calibrate_censoring <- function(tau, target, cfg, calibration_seed) {
  set.seed(calibration_seed)
  N <- cfg$calibration_n
  x <- stats::rnorm(N)
  eps <- mvtnorm::rmvnorm(N, sigma = vine_cor(atanh(cfg$rho)))
  eps <- sweep(eps, 2L, cfg$sigma, `*`)
  W <- matrix(stats::rnorm(3L * N, sd = tau), N, 3L)
  means <- sweep(outer(x, cfg$beta), 2L, cfg$alpha, `+`)
  total <- rowSums(exp(means + eps + W))
  fraction <- function(upper) mean(pmin(total / upper, 1))
  lower <- 1e-5
  upper <- max(10, stats::quantile(total, 0.99) * 10)
  while (fraction(upper) > target) upper <- upper * 2
  stats::uniroot(function(u) fraction(u) - target,
                 interval = c(lower, upper), tol = 1e-8)$root
}

# All Monte Carlo output is kept even if optim reports a nonzero code,
# as long as it supplies finite parameter values. MT errors do not delete
# parameter estimates; they are reported separately.
run_replication <- function(task, cfg, censor_upper) {
  assign(".Random.seed", task$seed, envir = .GlobalEnv)
  sim <- simulate_blocks(task$n, task$tau, censor_upper, cfg)
  dat <- sim$blocks
  th <- rep(NA_real_, 12L)
  z <- rep(NA_real_, 3L)
  conv1 <- rep(NA_integer_, 3L)
  conv2 <- NA_integer_
  issues <- character(0L)

  for (j in seq_len(3L)) {
    idx <- (4L * (j - 1L) + 1L):(4L * j)
    start <- c(cfg$alpha[j], cfg$beta[j],
               log(cfg$sigma[j]), log(task$tau))
    fit <- tryCatch(stats::optim(
      par = start,
      fn = function(v) marginal_nll(v, j, dat, cfg),
      method = "L-BFGS-B",
      lower = c(-5, -4, -2.3, -2.3),
      upper = c(5, 4, 1.7, 1.7),
      control = list(maxit = cfg$maxit_stage1)
    ), error = function(e) e)
    if (inherits(fit, "error")) {
      issues <- c(issues, paste0("stage1_gap", j, ": ",
                                conditionMessage(fit)))
    } else {
      th[idx] <- fit$par
      conv1[j] <- fit$convergence
      if (fit$convergence != 0L) {
        issues <- c(issues, paste0("stage1_gap", j,
                                  "_optim_code=", fit$convergence))
      }
    }
  }

  if (all(is.finite(th))) {
    fit <- tryCatch(stats::optim(
      par = atanh(0.7 * cfg$rho),
      fn = function(v) joint_nll(v, th, dat, cfg),
      method = "L-BFGS-B",
      lower = rep(-4.5, 3L), upper = rep(4.5, 3L),
      control = list(maxit = cfg$maxit_stage2)
    ), error = function(e) e)
    if (inherits(fit, "error")) {
      issues <- c(issues, paste0("stage2: ", conditionMessage(fit)))
    } else {
      z <- fit$par
      conv2 <- fit$convergence
      if (fit$convergence != 0L) {
        issues <- c(issues, paste0("stage2_optim_code=", fit$convergence))
      }
    }
  } else {
    issues <- c(issues, "stage2_not_attempted_due_to_stage1_error")
  }

  est <- natural_parameters(th, z)
  se <- rep(NA_real_, 15L)
  mt_status <- if (cfg$calculate_MT) "not_available" else "disabled"
  if (cfg$calculate_MT && all(is.finite(c(th, z)))) {
    se_try <- tryCatch(mt_standard_errors(th, z, dat, cfg),
                       error = function(e) e)
    if (inherits(se_try, "error")) {
      mt_status <- paste0("failed: ", conditionMessage(se_try))
    } else {
      se <- se_try
      mt_status <- "ok"
    }
  }

  list(replication = task$replication,
       est = est, se = se,
       stage1_codes = conv1, stage2_code = conv2,
       mt_status = mt_status,
       errors = paste(issues, collapse = " | "),
       realized_censor = sim$realized_censor,
       entry_gap2 = sim$entry_gap2,
       entry_gap3 = sim$entry_gap3,
       event_gap3 = sim$event_gap3)
}

# A top-level function, exported to PSOCK workers with all its dependencies.
safe_replication <- function(task, cfg, censor_upper) {
  tryCatch(run_replication(task, cfg, censor_upper),
           error = function(e) {
             list(replication = task$replication,
                  est = rep(NA_real_, 15L), se = rep(NA_real_, 15L),
                  stage1_codes = rep(NA_integer_, 3L),
                  stage2_code = NA_integer_, mt_status = "not_available",
                  errors = paste0("replication_error: ", conditionMessage(e)),
                  realized_censor = NA_real_, entry_gap2 = NA_real_,
                  entry_gap3 = NA_real_, event_gap3 = NA_real_)
           })
}

validate_config <- function(cfg) {
  m <- nrow(cfg$positions)
  if (m != 3L || ncol(cfg$positions) != 2L)
    stop("This trivariate design requires three 2D spatial site coordinates.")
  if (any(cfg$sample_sizes %% m != 0L) || any(cfg$sample_sizes < 6L))
    stop("Each sample size must be >= 6 and divisible by three sites.")
  if (any(cfg$censor_rates <= 0 | cfg$censor_rates >= 1))
    stop("Censoring rates must lie strictly between 0 and 1.")
  if (any(cfg$spatial_tau <= 0)) stop("spatial_tau must be positive.")
  if (length(cfg$rho) != 3L || any(abs(cfg$rho) >= 1))
    stop("The three partial vine correlations must lie in (-1,1).")
  if (length(cfg$beta) != 3L || length(cfg$alpha) != 3L ||
      length(cfg$sigma) != 3L || any(cfg$sigma <= 0))
    stop("alpha, beta, sigma must be length three; sigma must be positive.")
  if (!is.finite(cfg$cores) || cfg$cores < 1L)
    stop("cores must be a positive integer.")
  if (cfg$replications < 1L) stop("replications must be >= 1.")
  if (cfg$calibration_n < 1000L) stop("calibration_n must be >= 1000.")
  invisible(TRUE)
}

make_truth <- function(tau, cfg) {
  setNames(c(cfg$alpha, cfg$beta, cfg$sigma, rep(tau, 3L), cfg$rho),
           c(paste0("alpha", 1:3), paste0("beta", 1:3),
             paste0("sigma", 1:3), paste0("tau", 1:3),
             "rho12", "rho23", "rho13_given_2"))
}

make_outputs <- function(out, n, tau, censor_target, censor_upper, cfg) {
  truth <- make_truth(tau, cfg)
  pars <- names(truth)
  E <- t(vapply(out, `[[`, numeric(15L), "est"))
  SE <- t(vapply(out, `[[`, numeric(15L), "se"))
  colnames(E) <- colnames(SE) <- pars

  parameter_rows <- lapply(seq_along(truth), function(k) {
    e <- E[, k]
    se <- SE[, k]
    usable <- is.finite(e)
    with_se <- usable & is.finite(se) & se > 0
    data.frame(
      n = n, tau = tau, censor_target = censor_target,
      censor_upper = censor_upper,
      attempted = length(out), parameter = pars[k], true = truth[k],
      estimate_count = sum(usable), mt_se_count = sum(with_se),
      bias = if (any(usable)) mean(e[usable] - truth[k]) else NA_real_,
      empirical_sd = if (sum(usable) >= 2L) stats::sd(e[usable]) else NA_real_,
      mean_MT_se = if (any(with_se)) mean(se[with_se]) else NA_real_,
      coverage_95 = if (any(with_se)) {
        mean(abs(e[with_se] - truth[k]) <= 1.96 * se[with_se])
      } else NA_real_,
      row.names = NULL
    )
  })

  diagnostic_rows <- lapply(out, function(r) {
    data.frame(n = n, tau = tau, censor_target = censor_target,
               replication = r$replication,
               realized_censor = r$realized_censor,
               entry_gap2 = r$entry_gap2,
               entry_gap3 = r$entry_gap3,
               event_gap3 = r$event_gap3,
               stage1_gap1_code = r$stage1_codes[1L],
               stage1_gap2_code = r$stage1_codes[2L],
               stage1_gap3_code = r$stage1_codes[3L],
               stage2_code = r$stage2_code,
               mt_status = r$mt_status,
               errors = r$errors,
               row.names = NULL)
  })

  replication_rows <- lapply(seq_along(out), function(i) {
    data.frame(n = n, tau = tau, censor_target = censor_target,
               replication = out[[i]]$replication,
               parameter = pars, true = unname(truth),
               estimate = out[[i]]$est, MT_se = out[[i]]$se,
               row.names = NULL)
  })

  list(summary = do.call(rbind, parameter_rows),
       diagnostics = do.call(rbind, diagnostic_rows),
       replication_estimates = do.call(rbind, replication_rows))
}

# Entire simulation, reusing one PSOCK cluster across all scenarios.
run_study <- function(cfg) {
  validate_config(cfg)
  cfg$K <- exp(-as.matrix(stats::dist(cfg$positions)) / cfg$spatial_range)
  grid <- expand.grid(n = cfg$sample_sizes,
                      tau = cfg$spatial_tau,
                      censor_target = cfg$censor_rates,
                      KEEP.OUT.ATTRS = FALSE)
  ncore <- min(as.integer(cfg$cores), as.integer(cfg$replications))

  # User-chosen number of workers; do not assume detectCores() is reliable.
  cl <- parallel::makeCluster(ncore, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    if (!requireNamespace("mvtnorm", quietly = TRUE))
      stop("mvtnorm must be installed for every parallel worker")
    NULL
  })
  exported <- c("vine_cor", "loglik_censored", "simulate_blocks",
                "block_marginal_loglik", "marginal_nll", "joint_covariance",
                "block_joint_loglik", "joint_nll", "fd_gradient",
                "fd_hessian", "fd_jacobian", "natural_parameters",
                "mt_standard_errors", "run_replication", "safe_replication")
  parallel::clusterExport(cl, varlist = exported,
                          envir = environment(safe_replication))

  # L'Ecuyer-CMRG independent per-replication RNG streams mean that results
  # do not change with the number of PSOCK workers or load balancing.
  RNGkind("L'Ecuyer-CMRG")
  sums <- diags <- repsout <- vector("list", nrow(grid))
  for (s in seq_len(nrow(grid))) {
    n <- as.integer(grid$n[s])
    tau <- grid$tau[s]
    target <- grid$censor_target[s]
    cat(sprintf("Scenario %d/%d: n=%d, tau=%.2f, censor=%.0f%%\n",
                s, nrow(grid), n, tau, 100 * target))
    seed_cal <- as.integer((as.double(cfg$master_seed) + 1009 * s) %%
                             .Machine$integer.max)
    upper <- calibrate_censoring(tau, target, cfg, seed_cal)
    cat(sprintf("  Calibrated Uniform(0, %.5f) administrative end time\n",
                upper))

    seed_rep <- as.integer((as.double(cfg$master_seed) + 9973 * s) %%
                             .Machine$integer.max)
    set.seed(seed_rep)
    seed <- .Random.seed
    tasks <- vector("list", cfg$replications)
    for (i in seq_len(cfg$replications)) {
      tasks[[i]] <- list(replication = i, n = n, tau = tau, seed = seed)
      seed <- parallel::nextRNGStream(seed)
    }
    out <- parallel::parLapplyLB(cl, tasks, safe_replication,
                                 cfg = cfg, censor_upper = upper)
    outputs <- make_outputs(out, n, tau, target, upper, cfg)
    sums[[s]] <- outputs$summary
    diags[[s]] <- outputs$diagnostics
    repsout[[s]] <- outputs$replication_estimates

    # Save after every scenario so completed scenarios survive interruptions.
    utils::write.csv(do.call(rbind, sums[seq_len(s)]),
                     paste0(cfg$output_prefix, "_summary.csv"),
                     row.names = FALSE)
    utils::write.csv(do.call(rbind, diags[seq_len(s)]),
                     paste0(cfg$output_prefix, "_diagnostics.csv"),
                     row.names = FALSE)
    utils::write.csv(do.call(rbind, repsout[seq_len(s)]),
                     paste0(cfg$output_prefix, "_replications.csv"),
                     row.names = FALSE)
    usable <- sum(vapply(out, function(a) all(is.finite(a$est)), logical(1)))
    mt_ok <- sum(vapply(out, function(a) identical(a$mt_status, "ok"),
                        logical(1)))
    cat(sprintf("  Finished: %d/%d full estimates, %d MT SE sets.\n",
                usable, length(out), mt_ok))
  }
  invisible(do.call(rbind, sums))
}

# Optional simple deterministic indexing check independent of estimation.
# Test whether flattening a B-by-3 block/site matrix preserves site order.
stopifnot(identical(as.vector(t(matrix(1:6, nrow = 2L, byrow = TRUE))),
                    1:6))

if (sys.nframe() == 0L) {
  run_study(CFG)
}
